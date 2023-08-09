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

/** @brief Calculate the target function which is being optimized.
 *  @param[in] x0: current values of the parameters
 */
  inline real_ Energy( const std_vec & x0 ) {
    for(uint k=0; k<x0.size(); k++) {
      if( x0[k] < real_(0) ) return real_(0);
      params[k] = x0[k];
    };

    if( settings::job_mode == 0 ) return scf::EngineSCF();
    if( settings::job_mode == 1 || settings::job_mode == 2 ) return cphf::EngineCPHF();
    return real_(0);
  };

/** @brief Optimization for a single basis set shell.
 */
  real_ PowellSingleShell();

/** @brief Line search for a minimum in a given direction:
 *         parabolic interpolation method.
 *  @param[in,out] x0: current values of the parameters
 *  @param[in] dir: direction in which we search
 *  @param[in] nvar: number of variables
 */
  void LineSearch( std_vec & x0, const std_vec dir, const int & nvar );

/** @brief Line search for a minimum in a given direction:
*          brute-force golden section search.
 *  @param[in,out] x0: current values of the parameters
 *  @param[in] dir: direction in which we search
 *  @param[in] nvar: number of variables
 */
  void LineSearchGolden( std_vec & x0, const std_vec dir, const int & nvar );

/** @brief Optimization of the whole basis.
 */
  void OptimizeBasis();

/** @brief Optimization of the basis for the CPHF response function.
 */
  void OptimizeCPHF();
};

#endif
