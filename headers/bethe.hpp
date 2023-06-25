#ifndef DENOM_H
#define DENOM_H

#include <string>
#include <vector>

namespace driver {
  /* standard Madelung rule: 1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, 4d, 5p */
  static constexpr int aufbau[] = { 0,0,1,1,1,0,1,1,1,0,2,2,2,2,2,1,1,1,0,2,2,2,2,2,1,1,1 };

  real_ Denominator( const basisSTO <real_> & basis, const std_mtx & C );
  real_ DenominatorRI( const basisSTO <real_> & basis, const std_mtx & C );


  real_ NablaExact( const basisSTO <real_> & basis, const std_mtx & C );
  real_ NablaRI( const basisSTO <real_> & basis, const std_mtx & C );

  real_ OrbitalsGradientTarget();

  void BetheDriver();

  template <typename T>
  class GaussLegendre {
  public:
    int n;
    T a, b;

    std::vector <T> xk;
    std::vector <T> wk;

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
