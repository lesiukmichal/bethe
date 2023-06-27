#include <cmath>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <sstream>
#include <limits>
#include <cassert>
#include <chrono>

#include "boost/math/quadrature/gauss.hpp"
#include "boost/math/constants/constants.hpp"

#include "../headers/main.hpp"
#include "../headers/ints.hpp"
#include "../headers/scf.hpp"
#include "../headers/cphf.hpp"
#include "../headers/bethe.hpp"
#include "../headers/powell.hpp"
#include "omp.h"

using namespace std;

real_ driver::OrbitalsGradientTarget() {
  vector <int> n_prim_orbital;
  vector <int> n_principal_orbital;

  for(int l=0; l<settings::global_lmax; l++) {
    if( settings::orbital_momentum[l] != 0 ) {
      n_prim_orbital.push_back( settings::orbital_momentum[l] );
      n_principal_orbital.push_back( l + 1 );
    };
  };

  vector <vector <real_>> A( n_prim_orbital.size(), vector <real_>( settings::n_orbital_legendre ) );
  for(int l=0; l<settings::global_lmax; l++)
    if( settings::orbital_momentum[l] != 0 )
      for(int k=0; k<settings::n_orbital_legendre; k++)
        A[l][k] = settings::orbital_basis_params[l][k];

  basisSTO <real_> basis( A, n_prim_orbital, n_principal_orbital );

  for(uint i=0; i<settings::extra_cphf_n.size(); i++) {
    basis.max_l = max( basis.max_l, settings::extra_cphf_l[i] );
    basis.basis.push_back( orbitalSTO <real_> (
      settings::extra_cphf_n[i],
      settings::extra_cphf_l[i],
      settings::extra_cphf_a[i] ) );
  };

  std_mtx C, S;
  C = scf::final_C;

  int n_extra_fun, n_extra_shl, n_basis_hf;
  n_basis_hf  = C.rows();
  n_extra_shl = settings::extra_cphf_n.size();
  n_extra_fun = 0;
  for(int l=0; l<n_extra_shl; l++)
    n_extra_fun += 2 * settings::extra_cphf_l[l] + 1;

  basis.n_shl += n_extra_shl;
  basis.n_fun += n_extra_fun;

  C.conservativeResize( basis.n_fun, basis.n_fun );
  ints::EngineUnified <real_> ( basis, basis, S, ints::IntegralType::Overlap );

  if( n_extra_fun > 0 ) {
    C.bottomRightCorner( n_extra_fun, n_extra_fun ) = std_mtx::Identity( n_extra_fun, n_extra_fun );
    C.topRightCorner( n_basis_hf, n_extra_fun ) = -( ( S * C ).bottomLeftCorner( n_extra_fun, n_basis_hf )
      * C.transpose().topLeftCorner( n_basis_hf, n_basis_hf ) ).transpose();

    std_mtx S_mix = C.transpose() * S * C;
    Eigen::SelfAdjointEigenSolver <std_mtx> FE;
    FE.compute( S_mix.bottomRightCorner( n_extra_fun, n_extra_fun ) );
    C.rightCols( n_extra_fun ) *= FE.eigenvectors();

    S_mix = C.transpose() * S * C;
    std_vec extra_norm = S_mix.diagonal().tail( n_extra_fun ).array().rsqrt();
    C.rightCols( n_extra_fun ) *= extra_norm.asDiagonal();
  };

  std_mtx T, V, PZ, T_mo, PZ_mo;
  ints::EngineSTV <real_> ( basis, basis, S, T, V, settings::nuclear_charge );
  ints::EngineUnified <real_> ( basis, basis, PZ, ints::IntegralType::Velocity );

  T = real_(2) * T / real_(3);
  T_mo  = C.transpose() * T * C;
  PZ_mo = C.transpose() * PZ * C;

  real_ target;
  target = T_mo.diagonal().head( settings::n_closed ).sum() + settings::rohf_coupling_f
         * T_mo.diagonal().segment( settings::n_closed, settings::n_open ).sum();

  target -= PZ_mo.leftCols( settings::n_closed ).cwiseAbs2().sum();
  target -= settings::rohf_coupling_f * PZ_mo.middleCols( settings::n_closed, settings::n_open ).cwiseAbs2().sum();
  return target;
};

real_ driver::Denominator( const basisSTO <real_> & basis, const std_mtx & C ) {
  std_mtx D1, D_closed, D_open;
  ints::EngineUnified <real_> ( basis, basis, D1, ints::IntegralType::Darwin1 );

  D_closed = C.leftCols( settings::n_closed ) * C.leftCols( settings::n_closed ).transpose();
  D_open   = C.block( 0, settings::n_closed, basis.n_fun, settings::n_open )
           * C.block( 0, settings::n_closed, basis.n_fun, settings::n_open ).transpose();

  return real_(2) * ( D_closed.cwiseProduct( D1 ).sum() + settings::rohf_coupling_f * D_open.cwiseProduct( D1 ).sum() );
};

real_ driver::NablaExact( const basisSTO <real_> & basis, const std_mtx & C ) {
  std_mtx S, T, V, PZ, T_mo, PZ_mo;
  ints::EngineSTV <real_> ( basis, basis, S, T, V, settings::nuclear_charge );
  ints::EngineUnified <real_> ( basis, basis, PZ, ints::IntegralType::Velocity );

  T = -real_(2) * T;
  T_mo  = C.transpose() * T * C;
  PZ_mo = C.transpose() * PZ * C;

  real_ p2;
  p2 = T_mo.diagonal().head( settings::n_closed ).sum() + settings::rohf_coupling_f
     * T_mo.diagonal().segment( settings::n_closed, settings::n_open ).sum();

  p2 += PZ_mo.topLeftCorner( settings::n_closed, settings::n_closed ).cwiseAbs2().sum();
  p2 += real_(2) * settings::rohf_coupling_f * PZ_mo.block( 0, settings::n_closed, settings::n_closed, settings::n_open ).cwiseAbs2().sum();
  p2 += settings::rohf_coupling_f * settings::rohf_coupling_f * settings::rohf_coupling_b
      * PZ_mo.block( settings::n_closed, settings::n_closed, settings::n_open, settings::n_open ).cwiseAbs2().sum();

  p2 *= real_(2);
  return p2;
};

real_ driver::NablaRI( const basisSTO <real_> & basis, const std_mtx & C ) {
  std_mtx PZ, PZ_mo;
  ints::EngineUnified <real_> ( basis, basis, PZ, ints::IntegralType::Velocity );
  PZ_mo = C.transpose() * PZ * C;

  real_ p2, a, b, c, d;
  a = -real_(3) * ( PZ_mo.leftCols( settings::n_closed ).cwiseAbs2().sum()
    - settings::rohf_coupling_f * PZ_mo.middleCols( settings::n_closed, settings::n_open ).cwiseAbs2().sum() );
  b = PZ_mo.topLeftCorner( settings::n_closed, settings::n_closed ).cwiseAbs2().sum();
  c = real_(2) * settings::rohf_coupling_f * PZ_mo.block( 0, settings::n_closed, settings::n_closed, settings::n_open ).cwiseAbs2().sum();
  d = settings::rohf_coupling_f * settings::rohf_coupling_f * settings::rohf_coupling_b
      * PZ_mo.block( settings::n_closed, settings::n_closed, settings::n_open, settings::n_open ).cwiseAbs2().sum();

  p2 = real_(2) * ( a + b + c + d );
  return p2;
};

real_ driver::DenominatorRI( const basisSTO <real_> & basis, const std_mtx & C_input ) {
  std_mtx C = C_input;
  ints::Integrals <real_> ERI;
  ints::EngineListERI( basis, ERI, settings::eri_threshold );

  std_mtx S, T, V, H;
  ints::EngineSTV <real_> ( basis, basis, S, T, V, settings::nuclear_charge );
  H = T + V;

  std_mtx D_closed, D_open;
  std_mtx J_closed, K_closed, F_closed;
  std_mtx J_open, K_open, F_open;
  std_mtx D_total, F_total;

  if( settings::n_closed > 0 ) {
    D_closed = C.leftCols( settings::n_closed ) * C.leftCols( settings::n_closed ).transpose();
    ints::FockFromList <real_> ( basis, ERI, D_closed, J_closed, K_closed );
  } else {
    D_closed.setZero( basis.n_fun, basis.n_fun );
    J_closed.setZero( basis.n_fun, basis.n_fun );
    K_closed.setZero( basis.n_fun, basis.n_fun );
  };

  if( settings::n_open > 0 ) {
    D_open = C.block( 0, settings::n_closed, basis.n_fun, settings::n_open )
           * C.block( 0, settings::n_closed, basis.n_fun, settings::n_open ).transpose();
    ints::FockFromList <real_> ( basis, ERI, D_open, J_open, K_open );
  } else {
    D_open.setZero( basis.n_fun, basis.n_fun );
    J_open.setZero( basis.n_fun, basis.n_fun );
    K_open.setZero( basis.n_fun, basis.n_fun );
  };

  D_total = real_(2) * D_closed;
  if( settings::n_open > 0 ) D_total += D_open / settings::n_open;
  scf::TotalFockOperator( H, C, S, J_closed, K_closed, J_open, K_open, F_closed, F_open, F_total );

  Eigen::SelfAdjointEigenSolver <std_mtx> FE;
  FE.compute( ( C.transpose() * F_total * C ).bottomRightCorner( settings::n_virtual, settings::n_virtual ) );
  C.rightCols( settings::n_virtual ) *= FE.eigenvectors();
  std_vec E = ( C.transpose() * F_total * C ).diagonal();

  std_mtx FC_mo, FO_mo;
  FC_mo = C.transpose() * F_closed * C;
  FO_mo = C.transpose() * F_open * C;

  std_mtx FC_xa, FC_xy, FC_xi, FC_ab;
  FC_xa = FC_mo.block( settings::n_closed, settings::n_closed + settings::n_open, settings::n_open, settings::n_virtual );
  FC_xy = FC_mo.block( settings::n_closed, settings::n_closed, settings::n_open, settings::n_open );
  FC_xi = FC_mo.block( settings::n_closed, 0, settings::n_open, settings::n_closed );
  FC_ab = FC_mo.bottomRightCorner( settings::n_virtual, settings::n_virtual );

  std_mtx FO_ij, FO_ai, FO_ab, FO_ix;
  FO_ij = FO_mo.topLeftCorner( settings::n_closed, settings::n_closed );
  FO_ai = FO_mo.block( settings::n_closed + settings::n_open, 0, settings::n_virtual, settings::n_closed );
  FO_ab = FO_mo.bottomRightCorner( settings::n_virtual, settings::n_virtual );
  FO_ix = FC_mo.block( 0, settings::n_closed, settings::n_closed, settings::n_open );

  std_mtx V_ix, V_ai, V_ax, DZ_ao, DZ_mo;
  ints::EngineUnified <real_> ( basis, basis, DZ_ao, ints::IntegralType::Velocity );

  DZ_mo = DZ_ao.triangularView <Eigen::StrictlyUpper> ();
  DZ_ao = DZ_mo + DZ_mo.transpose();
  DZ_mo = C.transpose() * DZ_ao * C;

  V_ix  = DZ_mo.block( settings::n_closed, 0, settings::n_open, settings::n_closed ).transpose();
  V_ai  = DZ_mo.block( settings::n_closed + settings::n_open, 0, settings::n_virtual, settings::n_closed );
  V_ax  = DZ_mo.block( settings::n_closed + settings::n_open, settings::n_closed, settings::n_virtual, settings::n_open );

  std_mtx DF_xi, DF_ai, DF_ax;
  cphf::FockDerivatives( V_ix, V_ai, V_ax, basis, ERI, C, DF_xi, DF_ai, DF_ax );

  std_mtx R_ix, R_ai, R_ax;
  if( settings::n_closed > 0 && settings::n_open > 0 )
    R_ix = DF_xi.transpose() + FO_ij * V_ix + FO_ai.transpose() * V_ax - ( FC_xa * V_ai ).transpose() + V_ix * FC_xy;
  if( settings::n_closed > 0 ) R_ai = DF_ai - FO_ab * V_ai - V_ax * FC_xi - ( V_ix * FC_xa ).transpose();
  if( settings::n_open > 0 ) R_ax = DF_ax - V_ai * FO_ix + FO_ai * V_ix - FC_ab * V_ax;

  for(int i=0; i<settings::n_closed; i++) for(int x=0; x<settings::n_open; x++)
    R_ix(i,x) -= ( E(i) + E(settings::n_closed + x) ) * V_ix(i,x);

  for(int i=0; i<settings::n_closed; i++) for(int a=0; a<settings::n_virtual; a++)
    R_ai(a,i) += ( E(settings::n_closed + settings::n_open + a) - E(i) ) * V_ai(a,i);

  for(int x=0; x<settings::n_open; x++) for(int a=0; a<settings::n_virtual; a++)
    R_ax(a,x) += ( E(settings::n_closed + settings::n_open + a) - E(settings::n_closed + x) ) * V_ax(a,x);

  real_ dtotal(0);
  if( settings::n_closed > 0 && settings::n_open > 0 ) dtotal += R_ix.cwiseProduct( V_ix ).sum();
  if( settings::n_closed > 0 ) dtotal += R_ai.cwiseProduct( V_ai ).sum();
  if( settings::n_open   > 0 ) dtotal += R_ax.cwiseProduct( V_ax ).sum();
  return -real_(6) * dtotal;
};

void driver::BetheDriver() {
  /*
   *   ---=== Basis set for Hartree-Fock calculation ===---
   */

  vector <int> n_prim_orbital;
  vector <int> n_principal_orbital;

  for(int l=0; l<settings::global_lmax; l++) {
    if( settings::orbital_momentum[l] != 0 ) {
      n_prim_orbital.push_back( settings::orbital_momentum[l] );
      n_principal_orbital.push_back( l + 1 );
    };
  };

  vector <vector <real_>> A( n_prim_orbital.size(), vector <real_>( settings::n_orbital_legendre ) );
  for(int l=0; l<settings::global_lmax; l++)
    if( settings::orbital_momentum[l] != 0 )
      for(int k=0; k<settings::n_orbital_legendre; k++)
        A[l][k] = settings::orbital_basis_params[l][k];

  basisSTO <real_> o_basis( A, n_prim_orbital, n_principal_orbital );
  settings::n_virtual = o_basis.n_fun - settings::n_closed - settings::n_open;
  cout << endl; o_basis.printDetails( "Orbital" );

  /*
   *   ---=== Identify the closed/open orbitals ===---
   */

  std_mtx L2_ao, L2_mo;
  std_vec L2_vec;
  ints::EngineUnified <real_> ( o_basis, o_basis, L2_ao, ints::IntegralType::Angular );
  L2_mo = scf::final_C.transpose() * L2_ao * scf::final_C;
  L2_vec = L2_mo.topLeftCorner( settings::n_closed + settings::n_open, settings::n_closed + settings::n_open ).diagonal();

  int l_orbital;
  bool check_aufbau = true;
  vector <int> orbital_types;

  for(uint i=0; i<L2_vec.size(); i++) {
    l_orbital = round( L2_vec(i) ).convert_to <int>();
    l_orbital = ( sqrt( 1 + 4 * l_orbital ) - 1 ) / 2;
    if( l_orbital != aufbau[i] ) check_aufbau = false;
    if( l_orbital == 0 ) orbital_types.push_back( l_orbital );
    if( l_orbital > 0 && orbital_types.back() != l_orbital )
      orbital_types.push_back( l_orbital );
  };

  if( check_aufbau ) cout << endl <<" Orbital occupations agree with the Aufbau principle." << endl;
    else { cout << endl << " Orbital occupations do not agree with the Aufbau principle." << endl; exit( EXIT_FAILURE ); };

  int occupied_momentum[settings::global_lmax_1] = { 0 };
  for(uint i=0; i<orbital_types.size(); i++) occupied_momentum[orbital_types[i]]++;

  cout << " Reference configuration composition:" << endl;
  for(int l=0; l<settings::global_lmax; l++) if( occupied_momentum[l] != 0 )
    cout << "  Number of " << settings::shell_names[l] << " occupied shells = " << occupied_momentum[l] << endl;

  /*
   *   ---=== Description of the occupied/open orbitals gradient ===---
   */

  real_ ri_without_extra, ri_with_extra;
  ri_without_extra = OrbitalsGradientTarget();

  uint n_extra_fun, n_extra_shl, n_basis_hf;
  n_basis_hf  = scf::final_C.rows();
  n_extra_shl = 0;
  n_extra_fun = 0;

  for(int l=0; l<settings::global_lmax; l++) {
    n_extra_shl += occupied_momentum[l];
    n_extra_fun += occupied_momentum[l] * ( l + l + 3 );
    for(int i=0; i<occupied_momentum[l]; i++) {
      settings::extra_cphf_n.push_back( l + 1 );
      settings::extra_cphf_l.push_back( l + 1 );
      settings::extra_cphf_a.push_back( settings::nuclear_charge / real_( i + 1 ) );
    };
  };

  cout << endl;
  cout << " Number of extra shells    = " << n_extra_shl << endl;
  cout << " Number of extra functions = " << n_extra_fun << endl;
  cout << endl;

  ri_with_extra = OrbitalsGradientTarget();

  cout << " RI error without extra functions   = " << ri_without_extra << endl;
  cout << " RI error including extra functions = " << ri_with_extra << endl;
  cout << endl;

  /*
   *   ---=== Form the extended basis (without optimization) ===---
   */

  for(uint i=0; i<settings::extra_cphf_n.size(); i++) {
    o_basis.max_l = max( o_basis.max_l, settings::extra_cphf_l[i] );
    o_basis.basis.push_back( orbitalSTO <real_> (
      settings::extra_cphf_n[i],
      settings::extra_cphf_l[i],
      settings::extra_cphf_a[i] ) );
  };

  o_basis.n_shl += n_extra_shl;
  o_basis.n_fun += n_extra_fun;

  std_mtx C, S;
  C = scf::final_C;

  C.conservativeResize( o_basis.n_fun, o_basis.n_fun );
  ints::EngineUnified <real_> ( o_basis, o_basis, S, ints::IntegralType::Overlap );

  if( n_extra_fun > 0 ) {
    C.bottomRightCorner( n_extra_fun, n_extra_fun ) = std_mtx::Identity( n_extra_fun, n_extra_fun );
    C.topRightCorner( n_basis_hf, n_extra_fun ) = -( ( S * C ).bottomLeftCorner( n_extra_fun, n_basis_hf )
      * C.transpose().topLeftCorner( n_basis_hf, n_basis_hf ) ).transpose();

    std_mtx S_mix = C.transpose() * S * C;
    Eigen::SelfAdjointEigenSolver <std_mtx> FE;
    FE.compute( S_mix.bottomRightCorner( n_extra_fun, n_extra_fun ) );
    C.rightCols( n_extra_fun ) *= FE.eigenvectors();

    S_mix = C.transpose() * S * C;
    std_vec extra_norm = S_mix.diagonal().tail( n_extra_fun ).array().rsqrt();
    C.rightCols( n_extra_fun ) *= extra_norm.asDiagonal();
  };

  settings::n_virtual = o_basis.n_fun - settings::n_closed - settings::n_open;
  o_basis.printDetails( "Extended orbital" );

  ints::EngineUnified <real_> ( o_basis, o_basis, S, ints::IntegralType::Overlap );

  /*
   *   ---=== Denominator and \nabla^2 expectation value ===---
   */

  real_ den_ex, den_ri, den_total, p2, p2_ri;
  p2    = NablaExact( o_basis, C );
  p2_ri = NablaRI( o_basis, C );

  den_ex  = Denominator( o_basis, C );
  den_ex *= -real_(2) * boost::math::constants::pi<real_>() * settings::nuclear_charge;
  den_ri  = DenominatorRI( o_basis, C );
  den_total = den_ri;

  cout << endl;
  cout << " Calculation of auxiliary quantities:" << endl;
  cout << "  Denominator - exact wavefunction = " << setw( settings::print_length + 10 ) << den_ex << endl;
  cout << "  Denominator - evaluation with RI = " << setw( settings::print_length + 10 ) << den_ri << endl;
  cout << "  Relative difference between them = " << setw( settings::print_length + 10 );
  cout << abs( real_(2) * ( den_ex - den_ri ) / ( den_ex + den_ri ) ) << endl << endl;
  cout << "  P**2 expectation value      = " << setw( settings::print_length + 10 ) << p2 << endl;
  cout << "  P**2 expectation value (RI) = " << setw( settings::print_length + 10 ) << p2_ri << endl;
  cout << "  Relative RI error in P**2   = " << setw( settings::print_length + 10 ) << ( p2 - p2_ri ) / p2 << endl;
  cout << endl;

  /*
   *   ---=== Generate grid for integration over momentum ===---
   */

  vector <real_> grid_t;
  vector <real_> grid_k;
  vector <real_> grid_w;

  real_ one, two;
  one = real_(1);
  two = real_(2);

  cout << " Calculation of Bethe logarithm:" << endl;
  cout << "  Size of the numerical momentum grid     = " << settings::n_points_bethe << endl;

  GaussLegendre <real_> gl( settings::n_points_bethe, real_(0), real_(1) );

  grid_t = gl.xk;
  grid_w = gl.wk;

  for( auto & t : grid_t ) grid_k.push_back( ( one / t / t - one ) / two );

  /*
   *   ---=== Set up extra functions for the response calculations ===---
   */

  vector <real_> scaling_factor;
  settings::n_extra_shl_opt = 0;

  for(int l=0; l<settings::global_lmax; l++) {
    settings::n_extra_shl_opt += occupied_momentum[l];
    for(int i=0; i<occupied_momentum[l]; i++) {
      settings::extra_cphf_n.insert( settings::extra_cphf_n.begin(), l + 1 );
      settings::extra_cphf_l.insert( settings::extra_cphf_l.begin(), l + 1 );
      settings::extra_cphf_a.insert( settings::extra_cphf_a.begin(), one );
      scaling_factor.insert( scaling_factor.begin(), one / pow( real_(2), i ) );
    };
  };

  cout << "  Number of extra shells for optimization = " << n_extra_shl << endl;

  /*
   *   ---=== Response calculations for small t (fitting) ===---
   */

  vector <real_> t_small;
  vector <real_> k_small;
  vector <real_> f_small;

  cout << endl;
  cout << " Calculation of auxiliary quantities:" << endl;
  cout << "  Small-t grid starting point = " << setw( settings::print_length + 10 ) << settings::grid_start_small << endl;
  cout << "  Small-t grid point spacing  = " << setw( settings::print_length + 10 ) << settings::grid_step_small << endl;
  cout << "  Number of small-t points    = " << setw( settings::print_length + 10 ) << settings::n_points_bethe_small << endl;

  for(int n=0; n<settings::n_points_bethe_small; n++) {
    t_small.push_back( settings::grid_start_small + n * settings::grid_step_small );
    k_small.push_back( ( one / t_small.back() / t_small.back() - one ) / two );
  };

  for(uint n=0; n<t_small.size(); n++) {
    settings::response_k = k_small[n];

    cout << endl << " Response calculations for momentum = " << settings::response_k << endl;
    for(int i=0; i<settings::n_extra_shl_opt; i++)
      settings::extra_cphf_a[i] = sqrt( settings::response_k + settings::response_k ) * scaling_factor[i];

    powell::OptimizeCPHF();
    real_ current_f = -cphf::EngineCPHF( true );
    current_f *= real_(3) * settings::response_k;
    current_f += p2 - real_(2) * den_total * t_small[n] * t_small[n];
    current_f  = current_f / t_small[n] / t_small[n] / t_small[n];
    f_small.push_back( current_f );

    cout << "  Total  : " << setw( settings::print_length + 10 ) << current_f << endl;
  };

  cout << endl << " Summary of the calculations:" << endl;
  for(uint n=0; n<grid_k.size(); n++) {
    cout << setw(5) << n;
    cout << setw( settings::print_length + 10 ) << t_small[n];
    cout << setw( settings::print_length + 10 ) << k_small[n];
    cout << setw( settings::print_length + 10 ) << f_small[n] << endl;
  };

/*
 *   ---=== Response calculations for other points ===---
 */

  vector <real_> f_response;
  for(uint n=0; n<grid_k.size(); n++) {
    settings::response_k = grid_k[n];

    cout << endl << " Response calculations for momentum = " << settings::response_k << endl;
    for(int i=0; i<settings::n_extra_shl_opt; i++)
      settings::extra_cphf_a[i] = sqrt( settings::response_k + settings::response_k ) * scaling_factor[i];

    powell::OptimizeCPHF();
    real_ current_f = -cphf::EngineCPHF( true );
    current_f *= real_(3) * settings::response_k;
    current_f += p2 - real_(2) * den_total * grid_t[n] * grid_t[n];
    current_f  = current_f / grid_t[n] / grid_t[n] / grid_t[n];
    f_response.push_back( current_f );

    cout << "  Total  : " << setw( settings::print_length + 10 ) << current_f << endl;
  };

  cout << endl << " Summary of the calculations:" << endl;
  for(uint n=0; n<grid_k.size(); n++) {
    cout << setw(5) << n;
    cout << setw( settings::print_length + 10 ) << grid_k[n];
    cout << setw( settings::print_length + 10 ) << grid_w[n];
    cout << setw( settings::print_length + 10 ) << f_response[n] << endl;
  };

  real_ numb(0);
  for(uint n=0; n<grid_k.size(); n++) numb += grid_w[n] * f_response[n];

  cout << endl;
  cout << " Total nominator   = " << setw( settings::print_length + 10 ) << -numb << endl;
  cout << " Total denominator = " << setw( settings::print_length + 10 ) << den_total << endl;
  cout << " Bethe logarithm   = " << setw( settings::print_length + 10 ) << -numb / den_total << endl;
};
