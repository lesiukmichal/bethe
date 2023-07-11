#ifndef POWELL_H
#define POWELL_H

#include <string>
#include <vector>

namespace powell {
  extern int mx_macro;
  extern int mx_powell;
  extern int mx_search;
  extern int mx_gold;

  extern real_ stop_powell;
  extern real_ stop_gold;

  extern real_* params;
  extern int nvar;

  inline real_ Energy( const std_vec & x0 ) {
    for(uint k=0; k<x0.size(); k++) {
      if( x0[k] < real_(0) ) return real_(0);
      params[k] = x0[k];
    };

    if( settings::job_mode == 0 ) return scf::EngineSCF();
    if( settings::job_mode == 1 || settings::job_mode == 2 ) return cphf::EngineCPHF();
    return real_(0);
  };

  real_ PowellSingleShell();
  void LineSearch( std_vec & x0, const std_vec dir, const int & nvar );
  void LineSearchGolden( std_vec & x0, const std_vec dir, const int & nvar );
  void OptimizeBasis();
  void OptimizeCPHF();
};

#endif
