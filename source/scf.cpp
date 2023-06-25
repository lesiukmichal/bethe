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
#include "omp.h"

using namespace std;

std_mtx scf::final_C;
std_vec scf::final_E;

void scf::TotalFockOperator( const std_mtx & H, const std_mtx & C, const std_mtx & S,
                             const std_mtx & J_closed, const std_mtx & K_closed,
                             const std_mtx & J_open, const std_mtx & K_open,
                             std_mtx & F_closed, std_mtx & F_open, std_mtx & F_total ) {
  using namespace settings;
  real_ two = real_(2);

  F_closed = H + two * J_closed - K_closed + rohf_coupling_f * ( two * J_open - K_open );
  F_open   = rohf_coupling_f * ( H + two * J_closed - K_closed + rohf_coupling_f
           * ( two * rohf_coupling_a * J_open - rohf_coupling_b * K_open ) );

  std_mtx F_closed_mo, F_open_mo, F_total_mo;
  F_closed_mo = C.transpose() * F_closed * C;
  F_open_mo = C.transpose() * F_open * C;

  F_total_mo = F_open_mo;
  F_total_mo.topLeftCorner( n_closed, n_closed ) = F_closed_mo.topLeftCorner( n_closed, n_closed );
  F_total_mo.topRightCorner( n_closed, n_virtual ) = F_closed_mo.topRightCorner( n_closed, n_virtual );
  F_total_mo.bottomLeftCorner( n_virtual, n_closed ) = F_closed_mo.topRightCorner( n_closed, n_virtual ).transpose();
  F_total_mo.bottomRightCorner( n_virtual, n_virtual ) += F_closed_mo.bottomRightCorner( n_virtual, n_virtual );

  if( n_open > 0 ) {
    F_total_mo.block( n_closed, 0, n_open, n_closed ) -= F_closed_mo.block( n_closed, 0, n_open, n_closed );
    F_total_mo.block( 0, n_closed, n_closed, n_open ) = F_total_mo.block( n_closed, 0, n_open, n_closed ).eval().transpose();
  };

  F_total = S * C * F_total_mo * C.transpose() * S;
};

real_ scf::EngineSCF( const bool & print ) {
  std_mtx S, H;
  std_mtx C;
  std_vec E;

  /*
   *   ---=== Set up the basis set for Hartree-Fock calculation ===---
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
  if( print ) o_basis.printDetails( "Orbital" );

  /*
   *   ---=== Generate one-electron integrals ===---
   */

  std_mtx T, V, Q, DZ;
  std_mtx D_closed, D_open;
  std_mtx J_closed, K_closed, F_closed;
  std_mtx J_open, K_open, F_open;
  std_mtx D_total, F_total;

  ints::EngineSTV <real_> ( o_basis, o_basis, S, T, V, settings::nuclear_charge );
  H = T + V; Q = ints::generateOrtho( S, settings::linear_dependent, print );

  /*
   *   ---=== Generate two-electron integrals ===---
   */

  ints::Integrals <real_> ERI;
  ints::EngineListERI( o_basis, ERI, settings::eri_threshold );
  /* ERI.print(); */

  /*
   *   ---=== Initial guess for the orbitals ===---
   */
  auto start = std::chrono::steady_clock::now();
  Eigen::SelfAdjointEigenSolver <std_mtx> FE;

  D_closed.setZero( o_basis.n_fun, o_basis.n_fun );
  J_closed.setZero( o_basis.n_fun, o_basis.n_fun );
  K_closed.setZero( o_basis.n_fun, o_basis.n_fun );

  D_open.setZero( o_basis.n_fun, o_basis.n_fun );
  J_open.setZero( o_basis.n_fun, o_basis.n_fun );
  K_open.setZero( o_basis.n_fun, o_basis.n_fun );

  if( final_C.size() > 0 ) {
    C = final_C;
    E = final_E;
  } else {
    FE.compute( Q.transpose() * H * Q );
    C = Q * FE.eigenvectors();
    E = FE.eigenvalues();

    TotalFockOperator( H, C, S, J_closed, K_closed, J_open, K_open, F_closed, F_open, F_total );

    FE.compute( Q.transpose() * F_total * Q );
    C = Q * FE.eigenvectors();
    E = FE.eigenvalues();
  };

  /*
   *   ---=== Start SCF iterations ===---
   */

  real_ e_hf, e_hf_0 = real_(0);
  real_ d_err;

  vector <std_mtx> store_f_diis;
  vector <std_mtx> store_d_diis;
  vector <std_mtx> store_e_diis;

  if( print ) {
    cout << endl << " SCF calculations settings:" << endl;
    cout << "  Maximum number of iterations : " << setw(4) << settings::scf_maxit << endl;
    cout << "  Size of the DIIS subspace    : " << setw(4) << settings::n_diis << endl;
    cout << "  Energy  convergence threshold: " << setw(20) << settings::scf_e_conv << endl;
    cout << "  Density convergence threshold: " << setw(20) << settings::scf_d_conv << endl;

    cout << endl;
    cout << " ---------------------------------------------------------------------------------" << endl;
    cout << "  ###           Energy                 Energy diff.            Density diff.      " << endl;
    cout << " ---------------------------------------------------------------------------------" << endl;
  };

  bool did_cnv = false;
  for(int iter=0; iter<settings::scf_maxit; iter++) {
    if( settings::n_closed > 0 ) {
      D_closed = C.leftCols( settings::n_closed ) * C.leftCols( settings::n_closed ).transpose();
      ints::FockFromList <real_> ( o_basis, ERI, D_closed, J_closed, K_closed );
    } else {
      D_closed.setZero( o_basis.n_fun, o_basis.n_fun );
      J_closed.setZero( o_basis.n_fun, o_basis.n_fun );
      K_closed.setZero( o_basis.n_fun, o_basis.n_fun );
    };

    if( settings::n_open > 0 ) {
      D_open = C.block( 0, settings::n_closed, o_basis.n_fun, settings::n_open )
             * C.block( 0, settings::n_closed, o_basis.n_fun, settings::n_open ).transpose();
      ints::FockFromList <real_> ( o_basis, ERI, D_open, J_open, K_open );
    } else {
      D_open.setZero( o_basis.n_fun, o_basis.n_fun );
      J_open.setZero( o_basis.n_fun, o_basis.n_fun );
      K_open.setZero( o_basis.n_fun, o_basis.n_fun );
    };

    D_total = real_(2) * D_closed;
    if( settings::n_open > 0 ) D_total += D_open / real_( settings::n_open );
    TotalFockOperator( H, C, S, J_closed, K_closed, J_open, K_open, F_closed, F_open, F_total );

    e_hf  = D_closed.cwiseProduct( H + F_closed ).sum()
          + D_open.cwiseProduct( settings::rohf_coupling_f * H + F_open ).sum();
    if( iter > 0 ) d_err = ( D_total - store_d_diis.back() ).cwiseAbs().maxCoeff();
      else d_err = real_(0);

    if( print ) {
      cout << std::setw(5) << iter;
      cout << std::setw( settings::print_length + 10 ) << e_hf;
      cout << std::setw( settings::print_length + 10 ) << abs( e_hf - e_hf_0 );
      if( iter > 0 ) cout << std::setw( settings::print_length + 10 ) << d_err;
      cout << std::endl;
    };

    if( abs( e_hf - e_hf_0 ) < settings::scf_e_conv
      &&  d_err < settings::scf_d_conv ) { did_cnv = true; break; };

    if( iter < settings::n_diis_turn_on ) {
      store_f_diis.clear();
      store_d_diis.clear();
      store_e_diis.clear();
    };

    store_f_diis.push_back( F_total );
    store_d_diis.push_back( D_total );
    store_e_diis.push_back( ErrorDIIS( F_total, D_total, S ) );

    if( store_d_diis.size() > (uint)settings::n_diis ) store_d_diis.erase( store_d_diis.begin() );
    if( store_f_diis.size() > (uint)settings::n_diis ) store_f_diis.erase( store_f_diis.begin() );
    if( store_e_diis.size() > (uint)settings::n_diis ) store_e_diis.erase( store_e_diis.begin() );

    SolveDIIS( store_f_diis, store_d_diis, store_e_diis, F_total, D_total );

    FE.compute( Q.transpose() * F_total * Q );
    C = Q * FE.eigenvectors();
    E = FE.eigenvalues();

    e_hf_0 = e_hf;
  };

  if( !did_cnv ) return real_(0); /* exit( EXIT_FAILURE ); */
  if( print ) cout << " ---------------------------------------------------------------------------------" << endl;

  final_C = C;
  final_E = E;

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration <double> total_time = end - start;
  if( print ) cout << endl << " SCF computing time: " << (int) total_time.count() << " s." << endl;
  return e_hf;
};

void scf::SolveDIIS( const vector<std_mtx> & store_f_diis,
                     const vector<std_mtx> & store_d_diis,
                     const vector<std_mtx> & store_e_diis,
                     std_mtx & F, std_mtx & D ) {
  int n_sub = store_f_diis.size();
  std_mtx A_DIIS;
  std_vec B_DIIS;

  A_DIIS.setZero( n_sub + 1, n_sub + 1 );
  B_DIIS.setZero( n_sub + 1 );
  B_DIIS( n_sub ) = -1;

  for(int i=0; i<n_sub; i++) {
    A_DIIS( i, n_sub ) = -1.0;
    A_DIIS( n_sub, i ) = -1.0;
    for(int j=0; j<=i; j++) {
      A_DIIS( i, j ) = store_e_diis[i].cwiseProduct( store_e_diis[j] ).sum();
      A_DIIS( j, i ) = A_DIIS( i, j );
    };
  };

  Eigen::ColPivHouseholderQR <std_mtx> lin_sol( A_DIIS );
  std_vec X_DIIS = lin_sol.solve( B_DIIS );

  F.setZero( F.rows(), F.cols() );
  D.setZero( D.rows(), D.cols() );

  for(int i=0; i<n_sub; i++) {
    F += store_f_diis[i] * X_DIIS(i);
    D += store_d_diis[i] * X_DIIS(i);
  };
};
