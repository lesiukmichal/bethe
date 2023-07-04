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
#include "../headers/powell.hpp"
#include "omp.h"

using namespace std;

void powell::OptimizeCPHF() {
  powell::nvar = settings::n_extra_shl_opt;
  powell::params = settings::extra_cphf_a.data();
  powell::PowellSingleShell();

  cout << endl << "  Current basis set parameters:" << endl;
  for(int n=0; n<settings::n_extra_shl_opt; n++) cout << "  " << settings::extra_cphf_a[n] << " ";
  cout << endl;
};

void powell::OptimizeBasis() {
  real_ e_now = real_(0);
  vector <real_> energy_history;
  energy_history.push_back( e_now );

  powell::nvar = settings::n_orbital_legendre;
  for(int n_macro=0; n_macro<powell::mx_macro; n_macro++) {

    for(int l=0; l<settings::global_lmax_1; l++) {
      if( settings::orbital_momentum[l] == 0 ) continue ;
      cout << " ---------------------------------------------------------------------------------" << endl;
      cout << " --== MACROITERATION " << setw(4) << n_macro << " | Angular momentum: " << setw(2) << l << endl;
      cout << " ---------------------------------------------------------------------------------" << endl;
      cout << endl;

      powell::params = settings::orbital_basis_params[l];
      e_now = powell::PowellSingleShell();

      cout << endl << "  Current basis set tempering parameters:" << endl;
      for(int n=0; n<settings::global_lmax_1; n++) {
        if( settings::orbital_momentum[n] == 0 ) continue ;
        for(int k=0; k<settings::n_orbital_legendre; k++) cout << "  " << settings::orbital_basis_params[n][k] << " ";
        cout << endl;
      };
      cout << endl;
    };

    if( real_(2) * abs( e_now - energy_history.back() ) / abs( e_now + energy_history.back() ) < stop_powell ) {
      energy_history.push_back( e_now ); break ;
    }; energy_history.push_back( e_now );
  };

  cout << "  Macroiteration history:" << endl;
  for(uint n=0; n<energy_history.size(); n++) {
    cout << setw(4) << n;
    cout << setw( settings::print_length + 10 ) << energy_history[n];
    if( n > 0 ) cout << setw( settings::print_length + 10 ) << energy_history[n] - energy_history[n-1];
    cout << endl;
  };
  cout << " ---------------------------------------------------------------------------------" << endl;
  cout << endl;
};

real_ powell::PowellSingleShell() {
  real_ e0, e1, emax;
  real_ f0, fn, fe;
  real_ e_new, e_old;
  int rmax;

  std_vec xi, t0, df;
  xi.setZero( nvar );
  t0.setZero( nvar );
  df.setZero( nvar + 1 );

  std_mtx uk, pk;
  uk.setIdentity( nvar, nvar );
  pk.setZero( nvar + 1, nvar + 1 );
  for(int k=0; k<nvar; k++) xi(k) = params[k];

  e_old = 0.0;
  for(int n_powell=0; n_powell<mx_powell; n_powell++) {
    pk.col(0).head(nvar) = xi.transpose();
    e_new = powell::Energy( xi );
    if( real_(2) * abs( e_old - e_new ) / abs( e_old + e_new ) < stop_powell ) break ;
    e_old = e_new;

    for(int n_dir=0; n_dir<nvar; n_dir++) {
      t0 = pk.col(n_dir).head(nvar).transpose();
      e0 = powell::Energy( t0 );
      powell::LineSearch( t0, uk.row(n_dir).transpose(), nvar );
      e1 = powell::Energy( t0 );

      cout << "  Line search in the direction " << setw(3) << n_dir << ":" << endl;
      cout << "  Starting energy = " << e0 << endl;
      cout << "  Final energy  = " << e1 << endl;
      df[n_dir] = e1 - e0;
      pk.col(n_dir+1).head(nvar) = t0.transpose();

      cout << "  Parameters: " << endl;
      for(int k=0; k<nvar; k++) cout << "   " << k << ": " << t0(k) << endl;
    };

    emax = abs( df[0] );
    rmax = 0;

    for(int n_dir=0; n_dir<nvar; n_dir++) {
      if( abs( df[n_dir] ) > emax ) {
        emax = abs( df[n_dir] );
        rmax = n_dir;
      };
    };
    cout << "  Largest energy change (direction " << setw(3) << rmax << ") = " << emax << endl;
    cout << endl;

    t0 = pk.col(0).head(nvar).transpose();
    f0 = powell::Energy( t0 );
    t0 = pk.col(nvar).head(nvar).transpose();
    fn = powell::Energy( t0 );
    t0 = real_(2) * pk.col(nvar).head(nvar).transpose() - pk.col(0).head(nvar).transpose();
    fe = powell::Energy( t0 );

    if( fe >= f0 || real_(2) * ( f0 - real_(2) * fn + fe ) * pow( f0 - fn - emax, 2 ) >= emax * pow( f0 - fe, 2 ) ) {
      xi = pk.col(nvar).head(nvar).transpose();
      continue ;
    };

    uk.row(rmax) = pk.col(nvar).head(nvar).transpose() - pk.col(0).head(nvar).transpose();
    t0 = pk.col(0).head(nvar).transpose();
    powell::LineSearch( t0, uk.row(rmax).transpose(), nvar );
    xi = t0;
  };
  return e_new;
};

void powell::LineSearch( std_vec & x0, const std_vec dir, const int & nvar ) {
  real_ old_tar, tar, gr, stp;
  real_ a, b, c, d, fc, fd;

  stp = real_("1.0e-04");
  old_tar = powell::Energy( x0 );
  cout << "   Stage 1: Bracketing." << endl;

  for( int i=1; i<=mx_search; i++ ) {
    x0 += stp * dir; tar = powell::Energy( x0 ); x0 -= stp * dir;
    cout << "  " << setw(4) << i << setw( settings::print_length + 10 ) << stp << setw( settings::print_length + 10 ) << tar << endl;
    if( old_tar <= tar ) break;
    old_tar = tar; stp = real_(2) * stp;
  };

  b   = stp;
  stp = real_("-1.0e-04");
  old_tar = powell::Energy( x0 );

  for( int i=1; i<=mx_search; i++ ) {
    x0 += stp * dir; tar = powell::Energy( x0 ); x0 -= stp * dir;
    cout << "  " << setw(4) << i << setw( settings::print_length + 10 ) << stp << setw( settings::print_length + 10 ) << tar << endl;
    if( old_tar <= tar ) break;
    old_tar = tar; stp = real_(2) * stp;
  };

  a  = stp;
  gr = ( real_(1) + sqrt( real_(5) ) ) / real_(2);
  c  = b - ( b - a ) / gr;
  d  = a + ( b - a ) / gr;

  cout << "   Stage 2: Golden section search." << endl;
  for( int i=1; i<=mx_gold; i++ ) {
    x0 += c * dir; fc = powell::Energy( x0 ); x0 -= c * dir;
    x0 += d * dir; fd = powell::Energy( x0 ); x0 -= d * dir;

    cout << "  " << setw( settings::print_length + 10 ) << c << setw( settings::print_length + 10 ) << d << setw( settings::print_length + 10 ) << fc << setw( settings::print_length + 10 ) << fd << endl;

    if( fc < fd ) { b = d; } else { a = c; };
    c  = b - ( b - a ) / gr;
    d  = a + ( b - a ) / gr;
    if( real_(2) * abs( fc - fd ) / abs( fc + fd ) < stop_gold ) break ;
  };

  stp = ( a + b ) / real_(2);
  cout << "   Minimum found = " << setw( settings::print_length + 10 ) << stp << endl;
  x0 += stp * dir;
};
