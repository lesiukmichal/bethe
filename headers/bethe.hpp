/** @file bethe.hpp
 *  @brief Bethe logarithm implementation.
 *  @author Michal Lesiuk
 */
#ifndef DENOM_H
#define DENOM_H

#include <string>
#include <vector>

namespace driver {
  /* standard Madelung rule: 1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, 4d, 5p */
  static constexpr int aufbau[] = { 0,0,1,1,1,0,1,1,1,0,2,2,2,2,2,1,1,1,0,2,2,2,2,2,1,1,1 };

/** @brief Calculate the denominator using the "exact" formalism,
 *         i.e. Eq. (17) in the paper.
 *  @param[in] basis: basis set class
 *  @param[in] C: Hartree-Fock orbitals coefficients matrix
 */
  real_ Denominator( const basisSTO <real_> & basis, const std_mtx & C );

/** @brief Calculate the denominator using the RI formalism,
 *         i.e. Eqs. (26) and (27) in the paper.
 *  @param[in] basis: basis set class
 *  @param[in] C: Hartree-Fock orbitals coefficients matrix
 */
  real_ DenominatorRI( const basisSTO <real_> & basis, const std_mtx & C );

/** @brief Calculate the \f$ \nabla^2 \f$ expectation value on the mean-field
 *         wavefunction using the exact formalism (Slater-Condon rules).
 *  @param[in] basis: basis set class
 *  @param[in] C: Hartree-Fock orbitals coefficients matrix
 */
  real_ NablaExact( const basisSTO <real_> & basis, const std_mtx & C );

/** @brief Calculate the \f$ \nabla^2 \f$ expectation value on the mean-field
 *         wavefunction using the RI formalism to approximate two-electron terms.
 *  @param[in] basis: basis set class
 *  @param[in] C: Hartree-Fock orbitals coefficients matrix
 */
  real_ NablaRI( const basisSTO <real_> & basis, const std_mtx & C );

/** @brief Calculate the average error of the RI approximation in the
 *         action of the \f$ \nabla \f$ operator on the occupied orbitals.
 */
  real_ OrbitalsGradientTarget();

/** @brief Main driver for calculation of the Bethe logarithm.
 */
  void BetheDriver();

/** @brief This class implements the Gauss-Legendre numerical quadrature -
 *         the nodes and weights for quadrature order \f$ n=5,10,15,\ldots 100 \f$
 *         are stored on the disk.
 *         constructor.
 */
  template <typename T>
  class GaussLegendre {
  public:
    int n;
    T a, b;

    std::vector <T> xk;
    std::vector <T> wk;

/** @brief The custom constructor.
 *  @param[in] n: the order of the quadrature
 *  @param[in] a: start of the integration interval
 *  @param[in] b: end of the integration interval
 */
    GaussLegendre( const int & n_, const T a_, const T b_ ) : n(n_), a(a_), b(b_) {
      xk.reserve( n );
      wk.reserve( n );

      std::string file_name = "private/gauss-legendre-" + std::to_string(n) + ".dat";
      std::ifstream file( file_name ); assert( file.good() );

      real_ x, w;
      while( file >> x >> w ) {
        xk.push_back( x );
        wk.push_back( w );

        T A, B;
        A = ( b - a ) / real_(2);
        B = ( b - a ) / real_(2);

        xk.back()  = A * xk.back() + B;
        wk.back() *= A;
      };
    };
  };
};

#endif
