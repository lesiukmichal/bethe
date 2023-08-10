/** @file ints.hpp
 *  @brief Various one- and two-electron integrals.
 *  @author Michal Lesiuk
 */
#ifndef INTS_H
#define INTS_H

#include <string>
#include <vector>

#include "boost/math/constants/constants.hpp"
#include "boost/multi_array.hpp"
#include "omp.h"

template <typename T> class orbitalSTO;
template <typename T> class basisSTO;

namespace coupling {
  static constexpr int n_binomial_max = 100;
  extern real_ binomial_table[n_binomial_max][n_binomial_max];

  extern boost::multi_array<real_, 5> wigner3j;

/** @brief Binomial coefficients from a look-up table.
 *  @param[in] n: upper index
 *  @param[in] k: lower index
 */
  template <typename T>
  inline T binomial( const int & n, const int & k ) {
    assert( n < n_binomial_max && k < n_binomial_max );
    return T( binomial_table[n][k] );
  };

/** @brief Binomial coefficients - recursive calculation.
 *  @param[in] n: upper index
 *  @param[in] k: lower index
 */
  template <typename T>
  T binomial_recursive( const int & n, const int & k ) {
    if( n < k ) return binomial_recursive <T> (k,n);
    if( k == 0 ) return T(1);
    return binomial_recursive <T> ( n - 1, k - 1 ) * T(n) / T(k);
  };

/** @brief Checking if the integer is odd.
 *  @param[in] val: the input integer
 */
  inline bool is_odd ( int val ) { return ((val) & 1); };

/** @brief Checking if the integer is even.
 *  @param[in] val: the input integer
 */
  inline bool is_even( int val ) { return (!(is_odd(val))); };

/** @brief Checking if a triplet of integers satisfies the
 *         triangle inequality condition.
 *  @param[in] two_ja: first integer
 *  @param[in] two_jb: second integer
 *  @param[in] two_jc: third integer
 */
  inline bool triangle_selection_fails( int two_ja, int two_jb, int two_jc ) {
    return ( ( two_jb < abs( two_ja - two_jc ) ) || ( two_jb > two_ja + two_jc ) );
  };

/** @brief Checking if the Clebsch–Gordan coefficient is non-zero.
 *  @param[in] two_ja: first upper index
 *  @param[in] two_jb: second upper index
 *  @param[in] two_jc: third upper index
 *  @param[in] two_ma: first lower index
 *  @param[in] two_mb: second lower index
 *  @param[in] two_mc: third lower index
 */
  inline bool m_selection_fails( int two_ja, int two_jb, int two_jc,
                          int two_ma, int two_mb, int two_mc ) {
    return ( abs(two_ma) > two_ja    || abs(two_mb) > two_jb
          || abs(two_mc) > two_jc    || is_odd(two_ja + two_ma)
          || is_odd(two_jb + two_mb) || is_odd(two_jc + two_mc)
          || (two_ma + two_mb + two_mc) != 0 );
  };

/** @brief Calculation of the Clebsch–Gordan coefficients using explicit formulas.
 *  @param[in] two_ja: first upper index
 *  @param[in] two_jb: second upper index
 *  @param[in] two_jc: third upper index
 *  @param[in] two_ma: first lower index
 *  @param[in] two_mb: second lower index
 *  @param[in] two_mc: third lower index
 */
  template <typename T>
  T Coupling_3J( int two_ja, int two_jb, int two_jc,
                 int two_ma, int two_mb, int two_mc ) {
  if( two_ja < 0 || two_jb < 0 || two_jc < 0 ) return T(0);
  if( triangle_selection_fails( two_ja, two_jb, two_jc )
   || m_selection_fails( two_ja, two_jb, two_jc, two_ma, two_mb, two_mc ) ) return T(0);

    int jca  = (-two_ja + two_jb + two_jc) / 2;
    int jcb  = ( two_ja - two_jb + two_jc) / 2;
    int jcc  = ( two_ja + two_jb - two_jc) / 2;
    int jmma = ( two_ja - two_ma) / 2;
    int jmmb = ( two_jb - two_mb) / 2;
    int jmmc = ( two_jc - two_mc) / 2;
    int jpma = ( two_ja + two_ma) / 2;
    int jpmb = ( two_jb + two_mb) / 2;
    int jpmc = ( two_jc + two_mc) / 2;
    int jsum = ( two_ja + two_jb + two_jc) / 2;
    int kmin = std::max( std::max( 0, jpmb-jmmc ), jmma - jpmc );
    int kmax = std::min( std::min( jcc, jmma ), jpmb );
    int k;
    int sign = is_odd( kmin - jpma + jmmb ) ? -1 : 1;
    T sum_pos = T(0);
    T sum_neg = T(0);
    T norm, term;
    T bc1, bc2, bc3, bcn1, bcn2, bcd1, bcd2, bcd3, bcd4;

    bcn1 = binomial <T> ( two_ja, jcc );
    bcn2 = binomial <T> ( two_jb, jcc );
    bcd1 = binomial <T> ( jsum+1, jcc );
    bcd2 = binomial <T> ( two_ja, jmma );
    bcd3 = binomial <T> ( two_jb, jmmb );
    bcd4 = binomial <T> ( two_jc, jpmc );

    norm = sqrt( bcn1 * bcn2 ) / sqrt( bcd1 * bcd2 * bcd3 * bcd4 * T(two_jc+1) );

    for (k = kmin; k <= kmax; k++) {
      bc1 = binomial <T> ( jcc, k );
      bc2 = binomial <T> ( jcb, jmma-k );
      bc3 = binomial <T> ( jca, jpmb-k );

      term = bc1 * bc2 * bc3;
      if( sign < 0 ) sum_neg += norm * term;
        else sum_pos += norm * term;
      sign = -sign;
    };
    return sum_pos - sum_neg;
  };

/** @brief Calculation of the Wigner 3J symbols using explicit formulas.
 *  @param[in] two_ja: first upper index
 *  @param[in] two_jb: second upper index
 *  @param[in] two_jc: third upper index
 *  @param[in] two_ma: first lower index
 *  @param[in] two_mb: second lower index
 *  @param[in] two_mc: third lower index
 */
  template <typename T>
  T Wigner_3J_explicit( int ja, int jb, int jc,
                        int ma, int mb, int mc ) {
    int two_ja, two_jb, two_jc;
    int two_ma, two_mb, two_mc;
    two_ja = 2 * ja;
    two_jb = 2 * jb;
    two_jc = 2 * jc;
    two_ma = 2 * ma;
    two_mb = 2 * mb;
    two_mc = 2 * mc;

    return Coupling_3J <T> ( two_ja, two_jb, two_jc,
                             two_ma, two_mb, two_mc );
  };

/** @brief Calculation of the Wigner 3J symbols using look-up tables.
 *  @param[in] two_ja: first upper index
 *  @param[in] two_jb: second upper index
 *  @param[in] two_jc: third upper index
 *  @param[in] two_ma: first lower index
 *  @param[in] two_mb: second lower index
 *  @param[in] two_mc: third lower index
 */
  template <typename T>
  T Wigner_3J( int ja, int jb, int jc,
               int ma, int mb, int mc ) {
    //    Wigner 3J symbol
    //         | ja jb jc |
    //         | ma mb mc |
    if( ma + mb != -mc ) return T(0);
    return T( coupling::wigner3j[ja][jb][jc][ja+ma][jb+mb] );
  };
};

namespace ints {
/** @brief A class storing the basic two-electron integrals.
 */
  template <typename T>
  class Integrals {
  public:
    std::vector <usint> i;
    std::vector <usint> j;
    std::vector <usint> k;
    std::vector <usint> l;

    std::vector <T> eri;

    inline void clear() {
      i.clear(); j.clear();
      k.clear(); l.clear();
      eri.clear();
    };

    inline void print() const {
      for(ulint n=0; n<eri.size(); n++) {
        std::cout << std::setw(5) << i[n];
        std::cout << std::setw(5) << j[n];
        std::cout << std::setw(5) << k[n];
        std::cout << std::setw(5) << l[n];
        std::cout << std::setw( settings::print_length + 10 ) << eri[n];
        std::cout << std::endl;
      };
    };
  };

/** @brief Given an overlap matrix, generates the transformation
 *         matrix to the orthonormal basis using the Lowdin
 *         symmetric orthogonalization.
 *  @param[in] S: the overlap matrix
 *  @param[in] overlap_threshold: threshold for dropping linearly dependent functions
 *  @param[in] print: if true, print additional information
 */
  template <typename T>
  inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    generateOrtho( const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & S,
                   const T & overlap_threshold, const bool & print = true ) {
    Eigen::SelfAdjointEigenSolver <Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> SE;
    SE.compute( S );

    int i_sing = 0;
    for(int i=0; i<S.rows(); i++) if( SE.eigenvalues()(i) < overlap_threshold ) i_sing++;

    if( print ) {
      std::cout << std::endl;
      std::cout << "  There are " << std::setw(3) << i_sing << " eigenvalues below ";
      std::cout << std::setw( 50 ) << overlap_threshold << "." << std::endl;
      std::cout << "  The lowest eigenvalue of the S matrix = " << std::setw(20) << SE.eigenvalues()(0) << std::endl;
    };

    int n_bas_eff = S.rows() - i_sing;
    Eigen::Matrix<T, Eigen::Dynamic, 1> tmp_eig_s = SE.eigenvalues().tail( n_bas_eff );
    tmp_eig_s = tmp_eig_s.array().rsqrt();
    return SE.eigenvectors().block( 0, i_sing, S.rows(), n_bas_eff ) * tmp_eig_s.asDiagonal();
  };

/** @brief Calculates the Legendre polynomials with the order
 *         \f$ n\in[0,nmax] \f$ using recursion relations
 *  @param[in] nmax: the maximum polynomial order
 *  @param[in] x: the real argument of the polynomial
 *  @param[out] pn: output as std::vector <T>
 */
  template <typename T>
  inline void LegendreP( const int & nmax, const T & x,
                         std::vector <T> & pn ) {
    pn.resize( nmax + 1 );
    pn[0] = T(1);
    pn[1] = x;

    int n1;
    for(int n=1; n<nmax; n++) {
      n1 = n + 1;
      pn[n1] = ( T( n + n1 ) * x * pn[n] - T(n) * pn[n-1] ) / T(n1);
    };
  };

/** @brief Calculates the basic one-electron integrals between a pair
 *         of STOs shells.
 *  @param[in] s1: bra shell class
 *  @param[in] s2: ket shell class
 *  @param[out] S_12: overlap integrals
 *  @param[out] T_12: overlap integrals
 *  @param[out] V_12: overlap integrals
 *  @param[in] Z: the nuclear charge
 */
  template <typename T>
  inline void IntegralsSTV( const orbitalSTO <T> & s1,
                            const orbitalSTO <T> & s2,
                            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & S_12,
                            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & T_12,
                            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & V_12,
                            const T & Z ) {
    S_12.setConstant( s1.shell_size, s2.shell_size, T(0) );
    T_12.setConstant( s1.shell_size, s2.shell_size, T(0) );
    V_12.setConstant( s1.shell_size, s2.shell_size, T(0) );

    if( s1.l != s2.l ) return ;
    T radial_s, radial_v, radial_t;
    T radial_t_1, radial_t_2, radial_t_3;

    radial_t_3 = tgamma( T( s1.n + s2.n + 1 ) ) / pow( s1.a + s2.a, s1.n + s2.n + 1 );
    radial_s   = radial_t_3 * s1.getNorm() * s2.getNorm();
    S_12.diagonal().setConstant( radial_s );

    radial_v = -Z * ( s1.a + s2.a ) / T( s1.n + s2.n ) * radial_s;
    V_12.diagonal().setConstant( radial_v );

    T fac1, fac2, fac3;
    fac1 = -T( s2.n * ( s2.n - 1 ) - s2.l * ( s2.l + 1 ) ) / T(2);
    fac2 = s2.a * T( s2.n );
    fac3 = -s2.a * s2.a / T(2);

    if( s2.n != s2.l + 1 )
      radial_t_1 = tgamma( T( s1.n + s2.n - 1 ) ) / pow( s1.a + s2.a, s1.n + s2.n - 1 );
        else radial_t_1 = T(0);
    radial_t_2 = tgamma( T( s1.n + s2.n ) ) / pow( s1.a + s2.a, s1.n + s2.n );

    radial_t  = fac1 * radial_t_1 + fac2 * radial_t_2 + fac3 * radial_t_3;
    radial_t *= s1.getNorm() * s2.getNorm();
    T_12.diagonal().setConstant( radial_t );
  };

/** @brief Calculates the basic one-electron integrals for
 *         the whole atomic basis supplied at the input.
 *  @param[in] basis_1: bra basis class
 *  @param[in] basis_2: ket basis class
 *  @param[out] SM: overlap integrals matrix
 *  @param[out] TM: overlap integrals matrix
 *  @param[out] VM: overlap integrals matrix
 *  @param[in] Z: the nuclear charge
 */
  template <typename T>
  inline void EngineSTV( const basisSTO <T> & basis_1, const basisSTO <T> & basis_2,
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & SM,
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & TM,
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & VM,
                         const T & Z ) {

    SM.setZero( basis_1.n_fun, basis_2.n_fun );
    TM.setZero( basis_1.n_fun, basis_2.n_fun );
    VM.setZero( basis_1.n_fun, basis_2.n_fun );

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> S_shl;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> T_shl;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> V_shl;

    int n_pos_1, n_pos_2;

    n_pos_1 = 0;
    for(int i_shl=0; i_shl<basis_1.n_shl; i_shl++) {
      n_pos_2 = 0;

      for(int j_shl=0; j_shl<basis_2.n_shl; j_shl++) {
        IntegralsSTV <T> ( basis_1.basis[i_shl], basis_2.basis[j_shl], S_shl, T_shl, V_shl, Z );

        SM.block( n_pos_1, n_pos_2, basis_1.basis[i_shl].shell_size, basis_2.basis[j_shl].shell_size ) = S_shl;
        TM.block( n_pos_1, n_pos_2, basis_1.basis[i_shl].shell_size, basis_2.basis[j_shl].shell_size ) = T_shl;
        VM.block( n_pos_1, n_pos_2, basis_1.basis[i_shl].shell_size, basis_2.basis[j_shl].shell_size ) = V_shl;

        n_pos_2 += basis_2.basis[j_shl].shell_size;
      };
      n_pos_1 += basis_1.basis[i_shl].shell_size;
    };
  };

/** @brief Auxiliary integrals \f$ \int_0^\infty dr\;r^n\,e^{-ar} \f$
 *         for \f$ n\in[0,nmax] \f$ calculated using recursion relation.
 *  @param[in] a: the value of the exponent
 *  @param[in] n_max: the maximum power in the integrand
 *  @param[out] rn: output as std::vector <T>
 */
  template <typename T>
  inline void Rn( const T & a,
                  const int & n_max,
                  std::vector <T> & rn ) {
    rn[0] = T(1) / a;
    for(int n=1; n<=n_max; n++) rn[n] = rn[n-1] * T(n) / a;
  };

/** @brief Auxiliary radial integrals needed for calculation of the
 *         atomic two-electron Coulomb integrals.
 *  @param[in] m_max: the maximum power of \f$ r_1 \f$ in the integrand
 *  @param[in] n_max: the maximum power of \f$ r_2 \f$ in the integrand
 *  @param[in] a: the value of the exponent of the first electron
 *  @param[in] b: the value of the exponent of the second electron
 *  @param[out] kmn: output as std::vector <std::vector<T>>
 */
  template <typename T>
  void Kmn( const int & m_max, const int & n_max,
            const T & a, const T & b,
            std::vector <std::vector<T>> & kmn ) {

    std::vector <T> ln( m_max + n_max + 1 );
    Rn( T(a + b), m_max + n_max, ln );

    for(int m=0; m<=m_max; m++) {
      kmn[m][0] = ln[m] / b;
      for(int n=1; n<=n_max; n++) kmn[m][n] = ( kmn[m][n-1] * T(n) + ln[m+n] ) / b;
    };
  };

/** @brief Auxiliary radial integrals needed for calculation of the
 *         atomic two-electron Coulomb integrals - a special case
 *         required by the 1p orbitals.
 *  @param[in] m_max: the maximum power of \f$ r_1 \f$ in the integrand
 *  @param[in] a: the value of the exponent of the first electron
 *  @param[in] b: the value of the exponent of the second electron
 *  @param[out] km: output as std::vector<T>
 */
  template <typename T>
  void KmSpecial( const int & m_max, const T & a, const T & b,
                  std::vector <T> & km ) {

    std::vector <T> ln( m_max + 1 );
    Rn( T(a + b), m_max, ln );

    km[0] = log( ( a + b ) / b ) / a;
    for(int m=1; m<=m_max; m++) km[m] = ( T(m) * km[m-1] - ln[m-1] ) / a;
  };

/** @brief Calculates the basic two-electron Coulomb integrals between
 *         a quartet of STOs shells, (11|22) notation.
 *  @param[in] s1: bra first shell class
 *  @param[in] s2: bra second shell class
 *  @param[in] s3: ket first shell class
 *  @param[in] s4: ket second shell class
 *  @param[out] shell: packed integrals as std::vector <T>
 */
  template <typename T>
  void IntegralsERI( const orbitalSTO <T> & s1,
                     const orbitalSTO <T> & s2,
                     const orbitalSTO <T> & s3,
                     const orbitalSTO <T> & s4,
                     std::vector <T> & shell ) {
    T a12, a34, fac1, fac2;
    int n12, n34;
    int l1max, l2max, l1min, l2min;
    int lmax, lmin, lmax1, kappa;
    int mm1, mm2;

    a12 = s1.a + s2.a;
    a34 = s3.a + s4.a;
    n12 = s1.n + s2.n;
    n34 = s3.n + s4.n;

    l1min = std::abs( s1.l - s2.l );
    l2min = std::abs( s3.l - s4.l );
    l1max = s1.l + s2.l;
    l2max = s3.l + s4.l;

    lmax  = std::min( l1max, l2max );
    lmin  = std::max( l1min, l2min );
    lmax1 = lmax + 1;

    kappa = s1.shell_size * s2.shell_size * s3.shell_size * s4.shell_size;
    shell.resize( kappa );

    fac1  = T(kappa);
    fac1  = sqrt( fac1 );
    fac2  = s1.getNorm() * s2.getNorm() * s3.getNorm() * s4.getNorm();
    fac2  = fac2 * fac1;

    int m_max_1, n_max_1;
    m_max_1 = n12 + lmax;
    n_max_1 = n34 - 1;

    std::vector <std::vector <T>> i1( m_max_1 + 1, std::vector <T> ( n_max_1 + 1 ) );
    Kmn( m_max_1, n_max_1, a12, a34, i1 );

    int m_max_2, n_max_2;
    m_max_2 = n34 + lmax;
    n_max_2 = n12 - 1;

    std::vector <std::vector <T>> i2( m_max_2 + 1, std::vector <T> ( n_max_2 + 1 ) );
    Kmn( m_max_2, n_max_2, a34, a12, i2 );

    std::vector <T> cg1( lmax1 );
    std::vector <T> cg2( lmax1 );
    for(int l=lmin; l<=lmax; l++)
      cg2[l] = coupling::Wigner_3J <T> ( s1.l, l, s2.l, 0, 0, 0 )
             * coupling::Wigner_3J <T> ( s3.l, l, s4.l, 0, 0, 0 );

    std::vector <T> k_special_1;
    std::vector <T> k_special_2;

    bool special = false;
    if( n34 == lmax ) {
      special = true;
      k_special_1.resize( n12 + lmax + 1 );
      KmSpecial <T> ( n12 + lmax, a12, a34, k_special_1 );
    };

    if( n12 == lmax ) {
      special = true;
      k_special_2.resize( n34 + lmax + 1 );
      KmSpecial <T> ( n34 + lmax, a34, a12, k_special_2 );
    };

    int ii = 0;
    for(int m1=-s1.l; m1<=s1.l; m1++) for(int m2=-s2.l; m2<=s2.l; m2++)
    for(int m3=-s3.l; m3<=s3.l; m3++) for(int m4=-s4.l; m4<=s4.l; m4++) {
      shell[ii] = T(0);

      mm1 = m2 - m1;
      mm2 = m4 - m3;
      if( mm1 != -mm2 ) { ii++; continue; };

      std::fill( cg1.begin(), cg1.end(), T(0) );
      for(int l=lmin; l<=lmax; l++) for(int m=-l; m<=l; m++)
        cg1[l] += T( std::pow( -1, m ) ) * coupling::Wigner_3J <T> ( s1.l, l, s2.l, -m1, -m, m2 )
                                         * coupling::Wigner_3J <T> ( s3.l, l, s4.l, -m3, +m, m4 );

      for(int l=lmin; l<lmax; l++) shell[ii] += cg1[l] * cg2[l] * ( i1[n12+l][n34-l-1] + i2[n34+l][n12-l-1] );

      if( special ) {
        if( n34 == lmax && n12 == lmax ) shell[ii] += cg1[lmax] * cg2[lmax] * ( k_special_1[n12+lmax] + k_special_2[n34+lmax] );
          else if( n34 == lmax ) shell[ii] += cg1[lmax] * cg2[lmax] * ( k_special_1[n12+lmax] + i2[n34+lmax][n12-lmax-1] );
          else if( n12 == lmax ) shell[ii] += cg1[lmax] * cg2[lmax] * ( i1[n12+lmax][n34-lmax-1] + k_special_2[n34+lmax] );
      } else shell[ii] += cg1[lmax] * cg2[lmax] * ( i1[n12+lmax][n34-lmax-1] + i2[n34+lmax][n12-lmax-1] );

      shell[ii] = shell[ii] * fac2 * T( std::pow( -1, std::abs( m1 + m3 ) ) );
      ii++;
    };
  };

/** @brief Calculates the two-electron part of the "skeleton" Coulomb
 *         and exchange matrices using the supplied density matrix
 *         and two-electron integrals.
 *  @param[in] basis: basis set class
 *  @param[in] ERI: integrals class
 *  @param[in] D: density matrix
 *  @param[out] J: Coulomb matrix
 *  @param[out] K: exchange matrix
 */
  template <typename T>
  void FockFromList( const basisSTO <T> & basis, const Integrals <T> & ERI,
                     const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & D,
                     Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & J,
                     Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & K ) {
    J.setZero( basis.n_fun, basis.n_fun );
    K.setZero( basis.n_fun, basis.n_fun );

    std::vector <std_mtx> J_parallel( settings::num_threads_openmp );
    std::vector <std_mtx> K_parallel( settings::num_threads_openmp );
    for(int i=0; i<settings::num_threads_openmp; i++) {
      J_parallel[i].setZero( basis.n_fun, basis.n_fun );
      K_parallel[i].setZero( basis.n_fun, basis.n_fun );
    };

    #pragma omp parallel default(none) \
     shared( ERI, D, J_parallel, K_parallel )
    {
    const uint id = omp_get_thread_num();
    int ind_1, ind_2;
    int ind_3, ind_4;
    T val;

    #pragma omp for schedule(static)
    for(ulint i=0; i<ERI.eri.size(); i++) {
      ind_1 = ERI.i[i];
      ind_2 = ERI.j[i];
      ind_3 = ERI.k[i];
      ind_4 = ERI.l[i];
      val = ERI.eri[i];

      J_parallel[id]( ind_1, ind_2 ) += val * D( ind_3, ind_4 );
      K_parallel[id]( ind_1, ind_4 ) += val * D( ind_2, ind_3 );
    };
    };

    for(int i=0; i<settings::num_threads_openmp; i++) {
      J += J_parallel[i]; K += K_parallel[i];
    };
  };

/** @brief Calculates the list of two-electron integrals from
 *         a given basis set.
 *  @param[in] basis: basis set class
 *  @param[out] ERI: integrals class
 *  @param[in] thresh_integrals: cutoff for dropping small integrals
 */
 template <typename T>
  void EngineListERI( const basisSTO <T> & basis,
                      Integrals <T> & ERI,
                      const T & thresh_integrals ) {
    std::vector <int> shell_map(1,0);
    for(int i_shl=1; i_shl<basis.n_shl; i_shl++)
      shell_map.push_back( shell_map.back() + basis.basis[i_shl-1].shell_size );

    #pragma omp parallel default(none) \
     shared( ERI, shell_map, basis, thresh_integrals )
    {
    int n_pos_1, n_pos_2;
    int n_pos_3, n_pos_4;
    int l1, l2, l3, l4;

    int ind_1, ind_2;
    int ind_3, ind_4;
    int ii; T val;

    std::vector <T> shell;

    #pragma omp for schedule(dynamic,1) collapse(2)
    for(int i_shl=0; i_shl<basis.n_shl; i_shl++)
    for(int j_shl=0; j_shl<basis.n_shl; j_shl++) {
      l1 = basis.basis[i_shl].l;
      l2 = basis.basis[j_shl].l;

      n_pos_1 = shell_map[i_shl];
      n_pos_2 = shell_map[j_shl];

      for(int k_shl=0; k_shl<basis.n_shl; k_shl++) {
        l3 = basis.basis[k_shl].l;
        n_pos_3 = shell_map[k_shl];

        for(int l_shl=0; l_shl<basis.n_shl; l_shl++) {
          l4 = basis.basis[l_shl].l;
          n_pos_4 = shell_map[l_shl];

          IntegralsERI <T> ( basis.basis[i_shl], basis.basis[j_shl],
                             basis.basis[k_shl], basis.basis[l_shl], shell );

          ii = 0;
          for(int m1=0; m1<=l1+l1; m1++) for(int m2=0; m2<=l2+l2; m2++)
          for(int m3=0; m3<=l3+l3; m3++) for(int m4=0; m4<=l4+l4; m4++) {
            val = shell[ii];

            if( abs( val ) < thresh_integrals ) { ii++; continue; };

            ind_1 = n_pos_1 + m1;
            ind_2 = n_pos_2 + m2;
            ind_3 = n_pos_3 + m3;
            ind_4 = n_pos_4 + m4;

            #pragma omp critical
            {
              ERI.i.push_back( ind_1 );
              ERI.j.push_back( ind_2 );
              ERI.k.push_back( ind_3 );
              ERI.l.push_back( ind_4 );
              ERI.eri.push_back( val );
            };
            ii++;
          };

        };
      };
    };
    };
  };

/** @brief Calculates the dipole moment integrals between a pair
 *         of STOs shells.
 *  @param[in] s1: bra shell class
 *  @param[in] s2: ket shell class
 *  @param[out] DZ_shl: overlap integrals
 */
  template <typename T>
  inline void IntegralsDipoleZ( const orbitalSTO <T> & s1,
                                const orbitalSTO <T> & s2,
                                Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & DZ_shl ) {
    int m_limit = std::min( s1.l, s2.l );
    DZ_shl.setConstant( s1.shell_size, s2.shell_size, T(0) );

    T radial = tgamma( T( s1.n + s2.n + 2 ) ) / pow( s1.a + s2.a, s1.n + s2.n + 2 );
    radial  *= sqrt( T( s1.shell_size * s2.shell_size ) ) * s1.getNorm() * s2.getNorm();

    for(int m1=-m_limit; m1<=m_limit; m1++)
      DZ_shl( m1 + s1.l, m1 + s2.l ) = pow( -1, m1 ) * radial
                                     * coupling::Wigner_3J <T> ( s1.l, 1, s2.l,  0 , 0,  0  )
                                     * coupling::Wigner_3J <T> ( s1.l, 1, s2.l, -m1, 0, +m1 );
  };

/** @brief Calculates the dipole velocity integrals between a pair
 *         of STOs shells.
 *  @param[in] s1: bra shell class
 *  @param[in] s2: ket shell class
 *  @param[out] PZ_shl: overlap integrals
 */
  template <typename T>
  inline void IntegralsVelocityZ( const orbitalSTO <T> & s1,
                                  const orbitalSTO <T> & s2,
                                  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & PZ_shl ) {
    int m_limit = std::min( s1.l, s2.l );
    PZ_shl.setConstant( s1.shell_size, s2.shell_size, T(0) );

    if( std::abs( s1.l - s2.l ) != 1 ) return ;

    T radial, factor;
    radial  = tgamma( T( s1.n + s2.n ) ) / pow( s1.a + s2.a, s1.n + s2.n + 1 );
    radial *= s1.getNorm() * s2.getNorm() / sqrt( T( s2.shell_size ) );

    if( s1.l == s2.l + 1 ) {
      for(int m2=-m_limit; m2<=m_limit; m2++) {
        factor = sqrt( T( s2.l + m2 + 1 ) * T( s2.l - m2 + 1 ) / T( s2.l + s2.l + 3 ) )
               * ( T( s2.n - s2.l - 1 ) * ( s1.a + s2.a ) - s2.a * T( s1.n + s2.n ) );
        PZ_shl( m2 + s1.l, m2 + s2.l ) = factor * radial;
      };
    } else if( s1.l == s2.l - 1 ) {
      for(int m2=-m_limit; m2<=m_limit; m2++) {
        factor = sqrt( T( s2.l + m2 ) * T( s2.l - m2 ) / T( s2.l + s2.l - 1 ) )
               * ( T( s2.n + s2.l ) * ( s1.a + s2.a ) - s2.a * T( s1.n + s2.n ) );
        PZ_shl( m2 + s1.l, m2 + s2.l ) = factor * radial;
      };
    };
  };

/** @brief Calculates the one-electron Darwin integrals between a pair
 *         of STOs shells.
 *  @param[in] s1: bra shell class
 *  @param[in] s2: ket shell class
 *  @param[out] D1_shl: overlap integrals
 */
  template <typename T>
  inline void IntegralsDarwin1( const orbitalSTO <T> & s1,
                                const orbitalSTO <T> & s2,
                                Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & D1_shl ) {
    D1_shl.setConstant( s1.shell_size, s2.shell_size, T(0) );

    if( s1.n + s2.n != 2 ) return ;
    D1_shl( s1.l, s2.l )  = s1.getNorm() * s2.getNorm();
    D1_shl( s1.l, s2.l ) *= sqrt( T( s1.shell_size * s2.shell_size ) );
    D1_shl( s1.l, s2.l ) /= ( T(4) * boost::math::constants::pi<T>() );
  };

/** @brief Calculates the overlap integrals between a pair of STOs shells.
 *  @param[in] s1: bra shell class
 *  @param[in] s2: ket shell class
 *  @param[out] S_shl: overlap integrals
 */
  template <typename T>
  inline void IntegralsOverlap( const orbitalSTO <T> & s1,
                                const orbitalSTO <T> & s2,
                                Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & S_shl ) {
    S_shl.setConstant( s1.shell_size, s2.shell_size, T(0) );
    if( s1.l != s2.l ) return ;
    T radial_s = s1.getNorm() * s2.getNorm();
    radial_s *= tgamma( T( s1.n + s2.n + 1 ) ) / pow( s1.a + s2.a, s1.n + s2.n + 1 );
    S_shl.diagonal().setConstant( radial_s );
  };

/** @brief Calculates the angular momentum operator integrals between
 *         a pair of STOs shells.
 *  @param[in] s1: bra shell class
 *  @param[in] s2: ket shell class
 *  @param[out] L2_shl: overlap integrals
 */
  template <typename T>
  inline void IntegralsAngular( const orbitalSTO <T> & s1,
                                const orbitalSTO <T> & s2,
                                Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & L2_shl ) {
    L2_shl.setConstant( s1.shell_size, s2.shell_size, T(0) );
    if( s1.l != s2.l ) return ;
    T radial_s = s1.getNorm() * s2.getNorm() * s2.l * ( s2.l + 1 );
    radial_s *= tgamma( T( s1.n + s2.n + 1 ) ) / pow( s1.a + s2.a, s1.n + s2.n + 1 );
    L2_shl.diagonal().setConstant( radial_s );
  };

/** @brief Calculates the quadruple moment integrals (trace) between
 *         a pair of STOs shells.
 *  @param[in] s1: bra shell class
 *  @param[in] s2: ket shell class
 *  @param[out] R_shl: overlap integrals
 */
  template <typename T>
  inline void IntegralsMeanR( const orbitalSTO <T> & s1,
                              const orbitalSTO <T> & s2,
                              Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & R_shl ) {
    R_shl.setConstant( s1.shell_size, s2.shell_size, T(0) );
    if( s1.l != s2.l ) return ;
    T radial_r = s1.getNorm() * s2.getNorm();
    radial_r *= tgamma( T( s1.n + s2.n + 2 ) ) / pow( s1.a + s2.a, s1.n + s2.n + 2 );
    R_shl.diagonal().setConstant( radial_r );
  };

  enum class IntegralType { Overlap = 0, Dipole = 1, Velocity = 2, Darwin1 = 3, Angular = 4, MeanR = 5 };

/** @brief A wrapper for calculation of various one-electron integrals
 *         in the STOs basis.
 *  @param[in] basis_1: bra basis class
 *  @param[in] basis_2: ket basis class
 *  @param[out] I: integrals matrix
 *  @param[in] label: the integrals label (see the enum above)
 */
  template <typename T>
  inline void EngineUnified( const basisSTO <T> & basis_1, const basisSTO <T> & basis_2,
                             Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> & I,
                             const IntegralType & label ) {
    I.setZero( basis_1.n_fun, basis_2.n_fun );
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> I_shl;
    int n_pos_1, n_pos_2;

    n_pos_1 = 0;
    for(int i_shl=0; i_shl<basis_1.n_shl; i_shl++) {
      n_pos_2 = 0;

      for(int j_shl=0; j_shl<basis_2.n_shl; j_shl++) {
        switch ( label ) {
          case IntegralType::Overlap:
            IntegralsOverlap <T> ( basis_1.basis[i_shl], basis_2.basis[j_shl], I_shl ); break ;
          case IntegralType::Dipole:
            IntegralsDipoleZ <T> ( basis_1.basis[i_shl], basis_2.basis[j_shl], I_shl ); break ;
          case IntegralType::Velocity:
            IntegralsVelocityZ <T> ( basis_1.basis[i_shl], basis_2.basis[j_shl], I_shl ); break ;
          case IntegralType::Darwin1:
            IntegralsDarwin1 <T> ( basis_1.basis[i_shl], basis_2.basis[j_shl], I_shl ); break ;
          case IntegralType::Angular:
            IntegralsAngular <T> ( basis_1.basis[i_shl], basis_2.basis[j_shl], I_shl ); break ;
          case IntegralType::MeanR:
            IntegralsMeanR <T> ( basis_1.basis[i_shl], basis_2.basis[j_shl], I_shl ); break ;
          default:
            std::cout << " Invalid integral label! Emergency halt." << std::endl << std::endl;
            std::exit( EXIT_FAILURE );
        };

        I.block( n_pos_1, n_pos_2, basis_1.basis[i_shl].shell_size, basis_2.basis[j_shl].shell_size ) = I_shl;

        n_pos_2 += basis_2.basis[j_shl].shell_size;
      };
      n_pos_1 += basis_1.basis[i_shl].shell_size;
    };
  };
};

#endif
