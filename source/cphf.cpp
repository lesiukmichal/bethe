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

#include "../headers/main.hpp"
#include "../headers/ints.hpp"
#include "../headers/scf.hpp"
#include "../headers/cphf.hpp"
#include "omp.h"

using namespace std;

std_mtx cphf::UO_ix_final;
std_mtx cphf::UC_ai_final;
std_mtx cphf::UO_ax_final;

void cphf::TrialMultiplication( const std_mtx & UO_ix, const std_mtx & UC_ai, const std_mtx & UO_ax,
                                const std_mtx & DF_xi, const std_mtx & DF_ai, const std_mtx & DF_ax,
                                const std_mtx & FC_xa, const std_mtx & FC_xy,
                                const std_mtx & FC_xi, const std_mtx & FC_ab,
                                const std_mtx & FO_ij, const std_mtx & FO_ai,
                                const std_mtx & FO_ab, const std_mtx & FO_ix, const std_vec & E,
                                std_mtx & R_ix, std_mtx & R_ai, std_mtx & R_ax ) {
  if( settings::n_closed > 0 && settings::n_open > 0 )
    R_ix = DF_xi.transpose() + FO_ij * UO_ix + FO_ai.transpose() * UO_ax - ( FC_xa * UC_ai ).transpose() + UO_ix * FC_xy;
  if( settings::n_closed > 0 ) R_ai = DF_ai - FO_ab * UC_ai - UO_ax * FC_xi - ( UO_ix * FC_xa ).transpose();
  if( settings::n_open > 0 ) R_ax = DF_ax - UC_ai * FO_ix + FO_ai * UO_ix - FC_ab * UO_ax;

  for(int i=0; i<settings::n_closed; i++) for(int x=0; x<settings::n_open; x++)
    R_ix(i,x) -= ( E(i) + E(settings::n_closed + x) - settings::response_k
              * ( settings::rohf_coupling_f - real_(1) ) ) * UO_ix(i,x);

  for(int i=0; i<settings::n_closed; i++) for(int a=0; a<settings::n_virtual; a++)
    R_ai(a,i) += ( E(settings::n_closed + settings::n_open + a) - E(i) + settings::response_k ) * UC_ai(a,i);

  for(int x=0; x<settings::n_open; x++) for(int a=0; a<settings::n_virtual; a++)
    R_ax(a,x) += ( E(settings::n_closed + settings::n_open + a) - E(settings::n_closed + x)
              +      settings::response_k * settings::rohf_coupling_f ) * UO_ax(a,x);
};

void cphf::FockDerivatives( const std_mtx & UO_ix, const std_mtx & UC_ai, const std_mtx & UO_ax,
                            const basisSTO <real_> & basis, const ints::Integrals <real_> & ERI,
                            const std_mtx & C, std_mtx & DF_xi, std_mtx & DF_ai, std_mtx & DF_ax ) {
  using namespace settings;
  std_mtx D_closed, J_closed, K_closed;
  std_mtx D_open, J_open, K_open;

  if( n_closed > 0 ) {
    D_closed  = C.rightCols( n_virtual ) * UC_ai * C.leftCols( n_closed ).transpose();
    if( n_open > 0 )
      D_closed -= C.leftCols( n_closed ) * UO_ix * C.block( 0, n_closed, basis.n_fun, n_open ).transpose();
    D_closed = D_closed + D_closed.transpose().eval();

    ints::FockFromList <real_> ( basis, ERI, D_closed, J_closed, K_closed );
  } else {
    J_closed.setZero( basis.n_fun, basis.n_fun );
    K_closed.setZero( basis.n_fun, basis.n_fun );
  };

  if( n_open > 0 ) {
    D_open  = C.leftCols( n_closed ) * UO_ix * C.block( 0, n_closed, basis.n_fun, n_open ).transpose();
    D_open += C.rightCols( n_virtual ) * UO_ax * C.block( 0, n_closed, basis.n_fun, n_open ).transpose();
    D_open  = D_open + D_open.transpose().eval();

    ints::FockFromList <real_> ( basis, ERI, D_open, J_open, K_open );
  } else {
    J_open.setZero( basis.n_fun, basis.n_fun );
    K_open.setZero( basis.n_fun, basis.n_fun );
  };

  std_mtx DF_closed_ao, DF_open_ao;
  real_ two = real_(2);

  DF_closed_ao = two * J_closed - K_closed + rohf_coupling_f * ( two * J_open - K_open );
  DF_open_ao   = rohf_coupling_f * ( two * J_closed - K_closed + rohf_coupling_f
               * ( two * rohf_coupling_a * J_open - rohf_coupling_b * K_open ) );

  std_mtx DF_closed_mo, DF_open_mo;
  DF_closed_mo = C.transpose() * DF_closed_ao * C;
  DF_open_mo   = C.transpose() * DF_open_ao * C;

  DF_xi = DF_open_mo.block( n_closed, 0, n_open, n_closed ) - DF_closed_mo.block( n_closed, 0, n_open, n_closed );
  DF_ai = DF_closed_mo.bottomLeftCorner( n_virtual, n_closed );
  DF_ax = DF_open_mo.block( n_closed + n_open, n_closed, n_virtual, n_open );
};

real_ cphf::EngineCPHF( const bool & print ) {
  std_mtx C;
  std_vec E;

  /*
   *   ---=== Set up the basis set for CPHF calculation ===---
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

  basisSTO <real_> cphf_basis( A, n_prim_orbital, n_principal_orbital );

  int n_extra_fun, n_extra_shl, n_basis_hf;
  n_basis_hf  = cphf_basis.n_fun;
  n_extra_shl = settings::extra_cphf_n.size();
  n_extra_fun = 0;
  for(int l=0; l<n_extra_shl; l++)
    n_extra_fun += 2 * settings::extra_cphf_l[l] + 1;

  if( print ) {
    cout << endl;
    cout << " Number of extra shells    = " << n_extra_shl << endl;
    cout << " Number of extra functions = " << n_extra_fun << endl;
    cout << endl;
  };

  for(int i=0; i<n_extra_shl; i++) {
    cphf_basis.max_l = max( cphf_basis.max_l, settings::extra_cphf_l[i] );
    cphf_basis.basis.push_back( orbitalSTO <real_> (
      settings::extra_cphf_n[i],
      settings::extra_cphf_l[i],
      settings::extra_cphf_a[i] ) );
  };

  cphf_basis.n_shl += n_extra_shl;
  cphf_basis.n_fun += n_extra_fun;
  settings::n_virtual = cphf_basis.n_fun - settings::n_closed - settings::n_open;
  if( print ) cphf_basis.printDetails( "CPHF" );

  /*
   *   ---=== Generate two-electron integrals ===---
   */

  std_mtx D_closed, D_open;
  std_mtx J_closed, K_closed, F_closed;
  std_mtx J_open, K_open, F_open;
  std_mtx D_total, F_total;

  ints::Integrals <real_> ERI;
  ints::EngineListERI( cphf_basis, ERI, settings::eri_threshold );

  /*
   *   ---=== Prepare the extra orbitals for CPHF calculation ===---
   */

  C = scf::final_C;
  E = scf::final_E;

  C.conservativeResize( cphf_basis.n_fun, cphf_basis.n_fun );
  E.conservativeResize( cphf_basis.n_fun );

  std_mtx S, T, V, H;
  ints::EngineSTV <real_> ( cphf_basis, cphf_basis, S, T, V, settings::nuclear_charge );
  H = T + V;

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

  if( settings::n_closed > 0 ) {
    D_closed = C.leftCols( settings::n_closed ) * C.leftCols( settings::n_closed ).transpose();
    ints::FockFromList <real_> ( cphf_basis, ERI, D_closed, J_closed, K_closed );
  } else {
    D_closed.setZero( cphf_basis.n_fun, cphf_basis.n_fun );
    J_closed.setZero( cphf_basis.n_fun, cphf_basis.n_fun );
    K_closed.setZero( cphf_basis.n_fun, cphf_basis.n_fun );
  };

  if( settings::n_open > 0 ) {
    D_open = C.block( 0, settings::n_closed, cphf_basis.n_fun, settings::n_open )
           * C.block( 0, settings::n_closed, cphf_basis.n_fun, settings::n_open ).transpose();
    ints::FockFromList <real_> ( cphf_basis, ERI, D_open, J_open, K_open );
  } else {
    D_open.setZero( cphf_basis.n_fun, cphf_basis.n_fun );
    J_open.setZero( cphf_basis.n_fun, cphf_basis.n_fun );
    K_open.setZero( cphf_basis.n_fun, cphf_basis.n_fun );
  };

  D_total = real_(2) * D_closed;
  if( settings::n_open > 0 ) D_total += D_open / real_( settings::n_open );
  scf::TotalFockOperator( H, C, S, J_closed, K_closed, J_open, K_open, F_closed, F_open, F_total );

  FE.compute( ( C.transpose() * F_total * C ).bottomRightCorner( settings::n_virtual, settings::n_virtual ) );
  C.rightCols( settings::n_virtual ) *= FE.eigenvectors();
  };

  /*
   *   ---=== Assemble the final Fock operators ===---
   */

  if( settings::n_closed > 0 ) {
    D_closed = C.leftCols( settings::n_closed ) * C.leftCols( settings::n_closed ).transpose();
    ints::FockFromList <real_> ( cphf_basis, ERI, D_closed, J_closed, K_closed );
  } else {
    D_closed.setZero( cphf_basis.n_fun, cphf_basis.n_fun );
    J_closed.setZero( cphf_basis.n_fun, cphf_basis.n_fun );
    K_closed.setZero( cphf_basis.n_fun, cphf_basis.n_fun );
  };

  if( settings::n_open > 0 ) {
    D_open = C.block( 0, settings::n_closed, cphf_basis.n_fun, settings::n_open )
           * C.block( 0, settings::n_closed, cphf_basis.n_fun, settings::n_open ).transpose();
    ints::FockFromList <real_> ( cphf_basis, ERI, D_open, J_open, K_open );
  } else {
    D_open.setZero( cphf_basis.n_fun, cphf_basis.n_fun );
    J_open.setZero( cphf_basis.n_fun, cphf_basis.n_fun );
    K_open.setZero( cphf_basis.n_fun, cphf_basis.n_fun );
  };

  D_total = real_(2) * D_closed;
  if( settings::n_open > 0 ) D_total += D_open / real_( settings::n_open );
  scf::TotalFockOperator( H, C, S, J_closed, K_closed, J_open, K_open, F_closed, F_open, F_total );

  if( print ) {
    cout << endl;
    cout << " Norm of the off-diagonal blocks of the Fock operator = ";
    cout << ( C.transpose() * F_total * C ).topRightCorner( settings::n_closed + settings::n_open, settings::n_virtual ).norm();
    cout << endl;
  };

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

  E = ( C.transpose() * F_total * C ).diagonal();

  /*
   *   ---=== Perturbation integrals, response amplitudes guess ===---
   */

  std_mtx V_xi, V_ai, V_ax, DZ_ao, DZ_mo;
  ints::EngineUnified <real_> ( cphf_basis, cphf_basis, DZ_ao, ints::IntegralType::Velocity );

  DZ_mo = DZ_ao.triangularView <Eigen::StrictlyUpper> ();
  DZ_ao = DZ_mo + DZ_mo.transpose();
  DZ_mo = C.transpose() * DZ_ao * C;

  V_xi  = DZ_mo.block( settings::n_closed, 0, settings::n_open, settings::n_closed );
  V_ai  = DZ_mo.block( settings::n_closed + settings::n_open, 0, settings::n_virtual, settings::n_closed );
  V_ax  = DZ_mo.block( settings::n_closed + settings::n_open, settings::n_closed, settings::n_virtual, settings::n_open );

  std_mtx UO_ix, UC_ai, UO_ax;
  std_mtx B_ix, B_ai, B_ax;
  std_mtx diagonal_ix, diagonal_ai, diagonal_ax;

  B_ix = ( real_(1) - settings::rohf_coupling_f ) * V_xi.transpose();
  B_ai = -V_ai;
  B_ax = -settings::rohf_coupling_f * V_ax;

  diagonal_ix.setZero( settings::n_closed, settings::n_open );
  diagonal_ai.setZero( settings::n_virtual, settings::n_closed );
  diagonal_ax.setZero( settings::n_virtual, settings::n_open );

  for(int x=0; x<settings::n_open; x++) for(int a=0; a<settings::n_virtual; a++)
    diagonal_ax(a,x) = E( settings::n_closed + settings::n_open + a ) - E( settings::n_closed + x ) + settings::response_k;

  for(int i=0; i<settings::n_closed; i++) {
    for(int x=0; x<settings::n_open; x++) diagonal_ix(i,x) = -E(i) - E( settings::n_closed + x ) + settings::response_k;
    for(int a=0; a<settings::n_virtual; a++) diagonal_ai(a,i) = E( settings::n_closed + settings::n_open + a ) - E(i) + settings::response_k;
  };

  if( UO_ix_final.size() > 0 ) {
    UO_ix = UO_ix_final;
    UC_ai = UC_ai_final;
    UO_ax = UO_ax_final;
  } else {
    UO_ix = B_ix.cwiseQuotient( diagonal_ix );
    UC_ai = B_ai.cwiseQuotient( diagonal_ai );
    UO_ax = B_ax.cwiseQuotient( diagonal_ax );
  };

  std_vec meta_diagonal( diagonal_ix.size() + diagonal_ai.size() + diagonal_ax.size() );
  meta_diagonal << Eigen::Map <std_vec> ( diagonal_ix.data(), diagonal_ix.size() ),
                   Eigen::Map <std_vec> ( diagonal_ai.data(), diagonal_ai.size() ),
                   Eigen::Map <std_vec> ( diagonal_ax.data(), diagonal_ax.size() );

  std_vec B( B_ix.size() + B_ai.size() + B_ax.size() );
  B << Eigen::Map <std_vec> ( B_ix.data(), B_ix.size() ),
       Eigen::Map <std_vec> ( B_ai.data(), B_ai.size() ),
       Eigen::Map <std_vec> ( B_ax.data(), B_ax.size() );

  /*
   *   ---=== Start CPHF calculations ===---
   */

  real_ a_coef, b_coef;
  real_ prop;
  bool did_cnv;
  int iter;

  std_vec X0, X1, AP, R0, R1, P0, P1;
  std_mtx DF_xi, DF_ai, DF_ax;
  std_mtx R_ix, R_ai, R_ax;

  X0.setZero( UO_ix.size() + UC_ai.size() + UO_ax.size() );
  X0 << Eigen::Map <std_vec> ( UO_ix.data(), UO_ix.size() ),
        Eigen::Map <std_vec> ( UC_ai.data(), UC_ai.size() ),
        Eigen::Map <std_vec> ( UO_ax.data(), UO_ax.size() );

  FockDerivatives( UO_ix, UC_ai, UO_ax, cphf_basis, ERI, C, DF_xi, DF_ai, DF_ax );
  TrialMultiplication( UO_ix, UC_ai, UO_ax,
                       DF_xi, DF_ai, DF_ax,
                       FC_xa, FC_xy, FC_xi, FC_ab,
                       FO_ij, FO_ai, FO_ab, FO_ix,
                       E, R_ix, R_ai, R_ax );

  R0.setZero( R_ix.size() + R_ai.size() + R_ax.size() );
  R0 << Eigen::Map <std_vec> ( R_ix.data(), R_ix.size() ),
        Eigen::Map <std_vec> ( R_ai.data(), R_ai.size() ),
        Eigen::Map <std_vec> ( R_ax.data(), R_ax.size() );
  R0 = B - R0;
  P0 = R0.cwiseQuotient( meta_diagonal );

  if( print ) {
    cout << endl << " CPHF calculations settings:" << endl;
    cout << "  Maximum number of iterations  : " << setw(4) << settings::cphf_maxit << endl;
    cout << "  Residual convergence threshold: " << setw(20) << settings::cphf_conv << endl;
    cout << "  Perturbation frequency        : " << setw(20) << settings::response_k << endl;

    cout << endl;
    cout << " ---------------------------------------------------------------------------------" << endl;
    cout << "  ###        RMS(error)                 MAX(error)              Response          " << endl;
    cout << " ---------------------------------------------------------------------------------" << endl;
  };

  did_cnv = false;
  for(iter=0; iter<settings::cphf_maxit; iter++) {
    UO_ix = Eigen::Map <std_mtx> ( P0.head( UO_ix.size() ).data(), UO_ix.rows(), UO_ix.cols() );
    UC_ai = Eigen::Map <std_mtx> ( P0.segment( UO_ix.size(), UC_ai.size() ).data(), UC_ai.rows(), UC_ai.cols() );
    UO_ax = Eigen::Map <std_mtx> ( P0.tail( UO_ax.size() ).data(), UO_ax.rows(), UO_ax.cols() );

    FockDerivatives( UO_ix, UC_ai, UO_ax, cphf_basis, ERI, C, DF_xi, DF_ai, DF_ax );
    TrialMultiplication( UO_ix, UC_ai, UO_ax,
                         DF_xi, DF_ai, DF_ax,
                         FC_xa, FC_xy, FC_xi, FC_ab,
                         FO_ij, FO_ai, FO_ab, FO_ix,
                         E, R_ix, R_ai, R_ax );

    AP.setZero( R_ix.size() + R_ai.size() + R_ax.size() );
    AP << Eigen::Map <std_vec> ( R_ix.data(), R_ix.size() ),
          Eigen::Map <std_vec> ( R_ai.data(), R_ai.size() ),
          Eigen::Map <std_vec> ( R_ax.data(), R_ax.size() );

    a_coef  = R0.cwiseProduct( R0.cwiseQuotient( meta_diagonal ) ).sum();
    a_coef /= P0.cwiseProduct( AP ).sum();

    X1 = X0 + a_coef * P0;
    R1 = R0 - a_coef * AP;
    prop = -real_(4) * X1.dot(B);

    if( print ) {
      cout << setw(5) << iter + 1 << setw( settings::print_length + 10 ) << R1.norm();
      cout << setw( settings::print_length + 10 ) << R1.cwiseAbs().maxCoeff();
      cout << setw( settings::print_length + 10 ) << prop << endl;
    };

    if( R1.cwiseAbs().maxCoeff() < settings::cphf_conv ) { did_cnv = true; break; }

    b_coef  = R1.cwiseProduct( R1.cwiseQuotient( meta_diagonal ) ).sum();
    b_coef /= R0.cwiseProduct( R0.cwiseQuotient( meta_diagonal ) ).sum();

    P1 = R1.cwiseQuotient( meta_diagonal ) + b_coef * P0;

    P0 = P1;
    X0 = X1;
    R0 = R1;
  };
  if( print ) cout << " ---------------------------------------------------------------------------------" << endl;

  if( !did_cnv ) return real_(0); /* exit( EXIT_FAILURE ); */

  if( print ) {
    cout << "  CPHF solution found within " << setw(4) << iter << " iterations." << endl;
    cout << endl;
  };

  UO_ix = Eigen::Map <std_mtx> ( X1.head( UO_ix.size() ).data(), UO_ix.rows(), UO_ix.cols() );
  UC_ai = Eigen::Map <std_mtx> ( X1.segment( UO_ix.size(), UC_ai.size() ).data(), UC_ai.rows(), UC_ai.cols() );
  UO_ax = Eigen::Map <std_mtx> ( X1.tail( UO_ax.size() ).data(), UO_ax.rows(), UO_ax.cols() );

  /*
   *   ---=== Evaluate the propagator ===---
   */

  prop = real_(2) * UO_ix.cwiseProduct( B_ix ).sum()
       - real_(2) * UC_ai.cwiseProduct( B_ai ).sum()
       - real_(2) * UO_ax.cwiseProduct( B_ax ).sum();

  /* prop = -real_(2) * X1.dot(B); */
  if( print ) cout << "  <<A;B>>: " << setw( settings::print_length + 10 ) << prop << endl;

  UO_ix_final = UO_ix;
  UC_ai_final = UC_ai;
  UO_ax_final = UO_ax;
  return prop;
};
