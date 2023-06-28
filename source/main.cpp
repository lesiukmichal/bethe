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

#include "../mINI-master/src/mini/ini.h"

#include "../headers/main.hpp"
#include "../headers/ints.hpp"
#include "../headers/scf.hpp"
#include "../headers/cphf.hpp"
#include "../headers/powell.hpp"
#include "../headers/bethe.hpp"
#include "omp.h"

using namespace std;

#include "globals.cpp"

int main(int argc, char* argv[]) {
  auto start = std::chrono::steady_clock::now();

  cout << scientific;
  cout << setprecision( settings::print_length );

  if( argc > 1 ) {
    cout << " Name of the input file: " << argv[1] << endl;
  } else {
    cout << " Input file name not provided! " << endl;
    std::exit( EXIT_FAILURE );
  };

  Initialize();

  /*
   *   ---=== Read the input file ===---
   */

  mINI::INIFile ini_file( argv[1] );
  mINI::INIStructure ini;
  ini_file.read(ini);

  if( ini["job"].has( "mode" ) )
    settings::job_mode = stoi( ini["job"]["mode"] );
      else settings::job_mode = -1;

  if( ini["parallel"].has( "eigen" ) )
    settings::num_threads_eigen = stoi( ini["parallel"]["eigen"] );
      else settings::num_threads_eigen = 1;

  if( ini["parallel"].has( "openmp" ) )
    settings::num_threads_openmp = stoi( ini["parallel"]["openmp"] );
      else settings::num_threads_openmp = 1;

  omp_set_num_threads( settings::num_threads_openmp );
  Eigen::setNbThreads( settings::num_threads_eigen  );

  cout << endl;
  cout << " Number of Eigen  threads:" <<  setw(6) << settings::num_threads_eigen << endl;
  cout << " Number of OpenMP threads:" <<  setw(6) << settings::num_threads_openmp << endl;
  cout << endl;

  string tmp;
  if( ini["orbital_basis"].has( "n_orbital_legendre" ) )
    settings::n_orbital_legendre = stoi( ini["orbital_basis"]["n_orbital_legendre"] );
      else settings::n_orbital_legendre = 2;

  assert( ini["orbital_basis"].has( "orbital_momentum" ) );
  tmp = ini["orbital_basis"]["orbital_momentum"];
  istringstream is( tmp );
  uint index = 0; while( is >> settings::orbital_momentum[index] ) { index++; };

  assert( ini["system"].has( "nuclear_charge" ) );
  settings::nuclear_charge = real_( ini["system"]["nuclear_charge"] );

  assert( ini["system"].has( "nuclear_charge" ) );
  settings::nuclear_charge = real_( ini["system"]["nuclear_charge"] );

  assert( ini["system"].has( "n_electrons" ) );
  settings::n_electrons = stoi( ini["system"]["n_electrons"] );

  assert( ini["system"].has( "n_closed" ) );
  settings::n_closed = stoi( ini["system"]["n_closed"] );

  assert( ini["system"].has( "n_open" ) );
  settings::n_open = stoi( ini["system"]["n_open"] );

  if( ini["thresh"].has( "linear_dependent" ) )
    settings::linear_dependent = real_( ini["thresh"]["linear_dependent"] );
      else settings::linear_dependent = real_("1.0e-16");

  if( ini["thresh"].has( "eri_threshold" ) )
    settings::eri_threshold = real_( ini["thresh"]["eri_threshold"] );
      else settings::eri_threshold = real_("1.0e-30");

  if( ini["scf"].has( "maxit" ) )
    settings::scf_maxit = stoi( ini["scf"]["maxit"] );
      else settings::scf_maxit = 100;

  if( ini["scf"].has( "n_diis" ) )
    settings::n_diis = stoi( ini["scf"]["n_diis"] );
      else settings::n_diis = 100;

  if( ini["scf"].has( "n_diis_turn_on" ) )
    settings::n_diis_turn_on = stoi( ini["scf"]["n_diis_turn_on"] );
      else settings::n_diis_turn_on = 0;

  if( ini["scf"].has( "e_conv" ) )
    settings::scf_e_conv = real_( ini["scf"]["e_conv"] );
      else settings::scf_e_conv = real_(1.0e-16);

  if( ini["scf"].has( "d_conv" ) )
    settings::scf_d_conv = real_( ini["scf"]["d_conv"] );
      else settings::scf_d_conv = real_(1.0e-08);

  assert( ini["scf"].has( "rohf_coupling_f" ) );
  settings::rohf_coupling_f = real_( ini["scf"]["rohf_coupling_f"] );

  assert( ini["scf"].has( "rohf_coupling_a" ) );
  settings::rohf_coupling_a = real_( ini["scf"]["rohf_coupling_a"] );

  assert( ini["scf"].has( "rohf_coupling_b" ) );
  settings::rohf_coupling_b = real_( ini["scf"]["rohf_coupling_b"] );

  for(uint n=0; n<settings::global_lmax_1; n++) {
    if( settings::orbital_momentum[n] == 0 ) continue ;
    tmp = "basis_params_" + string( 1, settings::shell_names[n] );
    assert( ini["orbital_basis"].has( tmp ) );
    is.clear(); is.str( ini["orbital_basis"][tmp] );
    index = 0; while( is >> settings::orbital_basis_params[n][index] ) { index++; };
  };

  cout << " Input parameters:" << endl;
  cout << "  Nuclear charge            = " << settings::nuclear_charge << endl;
  cout << "  Number of electrons       = " << settings::n_electrons << endl;
  cout << "  Number of closed shells   = " << settings::n_closed << endl;
  cout << "  Number of open shells     = " << settings::n_open << endl;
  cout << endl;
  cout << "  Linear dependency thresh. = " << settings::linear_dependent << endl;
  cout << "  Integrals threshold       = " << settings::eri_threshold << endl;
  cout << endl;
  cout << "  Maximum SCF iterations    = " << settings::scf_maxit << endl;
  cout << "  Size of the DIIS subspace = " << settings::n_diis << endl;
  cout << "  SCF energy convergence    = " << settings::scf_e_conv << endl;
  cout << "  SCF density convergence   = " << settings::scf_d_conv << endl;
  cout << endl;
  cout << "  ROHF coupling parameter f = " << settings::rohf_coupling_f << endl;
  cout << "  ROHF coupling parameter a = " << settings::rohf_coupling_a << endl;
  cout << "  ROHF coupling parameter b = " << settings::rohf_coupling_b << endl;

  cout << endl;
  cout << "  Legendre expansion length = " << settings::n_orbital_legendre << endl;
  cout << "  Basis set composition     = ";
  for(uint n=0; n<settings::global_lmax_1; n++) {
    if( settings::orbital_momentum[n] == 0 ) continue ;
    cout << settings::shell_names[n] << settings::orbital_momentum[n];
  };

  cout << endl << "  Basis set tempering parameters:" << endl;
  for(int n=0; n<settings::global_lmax_1; n++) {
    if( settings::orbital_momentum[n] == 0 ) continue ;
    for(int k=0; k<settings::n_orbital_legendre; k++) cout << "  " << settings::orbital_basis_params[n][k] << " ";
    cout << endl;
  };
  cout << endl;

  /*
   *   ---=== SCF calculation/optimization ===---
   */

  if( ini["powell"].has( "mx_macro" ) )
    powell::mx_macro = stoi( ini["powell"]["mx_macro"] );
      else powell::mx_macro = 100;

  if( ini["powell"].has( "mx_powell" ) )
    powell::mx_powell = stoi( ini["powell"]["mx_powell"] );
      else powell::mx_powell = 100;

  if( ini["powell"].has( "mx_search" ) )
    powell::mx_search = stoi( ini["powell"]["mx_search"] );
      else powell::mx_search = 100;

  if( ini["powell"].has( "mx_gold" ) )
    powell::mx_gold = stoi( ini["powell"]["mx_gold"] );
      else powell::mx_gold = 100;

  if( ini["powell"].has( "stop_powell" ) )
    powell::stop_powell = real_( ini["powell"]["stop_powell"] );
      else powell::stop_powell = real_("1.0e-12");

  if( ini["powell"].has( "stop_gold" ) )
    powell::stop_gold = real_( ini["powell"]["stop_gold"] );
      else powell::stop_gold = real_("1.0e-12");

  cout << " Powell algorithm settings:" << endl;
  cout << "  Max. number of macroiterations    = " << powell::mx_macro << endl;
  cout << "  Max. number of Powell iterations  = " << powell::mx_powell << endl;
  cout << "  Max. number of bracketing steps   = " << powell::mx_search << endl;
  cout << "  Max. number of golden bisections  = " << powell::mx_gold << endl;
  cout << "  Powell method stopping condition  = " << powell::stop_powell << endl;
  cout << "  Golden section stopping condition = " << powell::stop_gold << endl;
  cout << endl;

  if( settings::job_mode == 0 ) powell::OptimizeBasis();
  scf::EngineSCF( true );

  /*
   *   ---=== Calculate the Bethe logarithm ===---
   */

  if( settings::job_mode == 1 || settings::job_mode == 2 ) {
    if( ini["cphf"].has( "maxit" ) )
      settings::cphf_maxit = stoi( ini["cphf"]["maxit"] );
        else settings::cphf_maxit = 100;

    if( ini["cphf"].has( "conv" ) )
      settings::cphf_conv = real_( ini["cphf"]["conv"] );
        else settings::cphf_conv = real_(1.0e-08);

    if( ini["bethe"].has( "n_grid" ) )
      settings::n_points_bethe = stoi( ini["bethe"]["n_grid"] );
        else settings::n_points_bethe = 50;

    if( ini["bethe"].has( "n_points_bethe_small" ) )
      settings::n_points_bethe_small = stoi( ini["bethe"]["n_points_bethe_small"] );
        else settings::n_points_bethe_small = 25;

    if( ini["bethe"].has( "grid_start_small" ) )
      settings::grid_start_small = real_( ini["bethe"]["grid_start_small"] );
        else settings::grid_start_small = real_(0.001);

    if( ini["bethe"].has( "grid_step_small" ) )
      settings::grid_step_small = real_( ini["bethe"]["grid_step_small"] );
        else settings::grid_step_small = real_(0.001);

    if( ini["bethe"].has( "k_single_shot" ) )
      settings::k_single_shot = real_( ini["bethe"]["k_single_shot"] );
        else settings::k_single_shot = real_(0.0);

    driver::BetheDriver();
  };

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration <double> total_time = end - start;
  cout << endl << " Total computing time: " << (int) total_time.count() << " s." << endl;
};

void Initialize() {
  for(uint n=0; n<coupling::n_binomial_max; n++)
  for(uint k=0; k<coupling::n_binomial_max; k++)
    coupling::binomial_table[n][k] = coupling::binomial_recursive <real_> ( n, k );

  for(int l1=0; l1<=settings::global_lmax; l1++)
  for(int l2=0; l2<=settings::global_lmax; l2++)
  for(int l3=0; l3<=settings::global_lmax; l3++)
  for(int m1=-l1; m1<=l1; m1++) for(int m2=-l2; m2<=l2; m2++)
    coupling::wigner3j[l1][l2][l3][m1+l1][m2+l2] =
    coupling::Wigner_3J_explicit <real_> ( l1, l2, l3, m1, m2, -m1 -m2 );
};
