#ifndef MAIN_H
#define MAIN_H

#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

using usint = unsigned short int;
using  sint = signed short int;
using  lint = long long int;
using  uint = unsigned int;
using ulint = unsigned long long int;

//#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/eigen.hpp>

using real_ = boost::multiprecision::float128;
//using real_ = boost::multiprecision::cpp_dec_float_50;
//using real_ = double;
using std_mtx  = Eigen::Matrix<real_, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using std_vec  = Eigen::Matrix<real_, Eigen::Dynamic, +1>;

namespace settings {
  extern int num_threads_eigen;
  extern int num_threads_openmp;

  static constexpr int global_lmax = 6;
  static constexpr int global_kmax = 4;
  static constexpr int global_lmax_1 = global_lmax + 1;
  static constexpr int global_mmax = global_lmax + global_lmax + 1;
  static constexpr int global_mmax_1 = global_mmax + 1;
  static constexpr int print_length = 15;

  extern int orbital_momentum[global_lmax_1];
  extern int n_orbital_legendre;
  extern real_ orbital_basis_params[global_lmax_1][global_kmax];

  extern real_ nuclear_charge;

  extern int n_closed;
  extern int n_open;
  extern int n_virtual;
  extern int n_electrons;

  extern real_ linear_dependent;
  extern real_ eri_threshold;

  extern int scf_maxit;
  extern int n_diis;
  extern int n_diis_turn_on;
  extern real_ scf_e_conv;
  extern real_ scf_d_conv;

  extern real_ rohf_coupling_f;
  extern real_ rohf_coupling_a;
  extern real_ rohf_coupling_b;

  extern int job_mode;

  static constexpr char shell_names[] = "SPDFGHIKLMN";

  extern std::vector <int> extra_cphf_n;
  extern std::vector <int> extra_cphf_l;
  extern std::vector <real_> extra_cphf_a;
  extern int n_extra_shl_opt;

  extern real_ response_k;
  extern int n_points_bethe;

  extern int cphf_maxit;
  extern real_ cphf_conv;

  extern int n_points_bethe_small;
  extern real_ grid_start_small;
  extern real_ grid_step_small;

  extern real_ k_single_shot;
};

#include "../headers/ints.hpp"

template <typename T>
class orbitalSTO {
  public:
    int n;
    int l;
    int shell_size;

    T a;

    orbitalSTO() : n(1), l(0), a(T(1)) {};
    orbitalSTO( const int & n_,
                const int & l_,
                const T & a_ ) : n(n_), l(l_), a(a_) {
      assert( l >=0 && n >= l && a > T(0) );
      setNorm();
      shell_size = l + l + 1;
    };

    inline void setNorm() {
      int n21 = n + n + 1;
      norm = pow( T(2) * a, T(n21) / T(2) );
      norm = norm / sqrt( tgamma( T(n21) ) );
    };

    inline T getNorm() const {
      return norm;
    };

  private:
    T norm;
};

template <typename T>
class basisSTO {
  public:
    int n_shl;
    int n_fun;
    int max_l;

    std::vector <orbitalSTO <T>> basis;

    basisSTO();

    basisSTO( const std::vector <std::vector <T>> & A,
              const std::vector <int> & n_prim,
              const std::vector <int> & n_principal ) {
      initializeTempering( A, n_prim, n_principal );
    };

    void initializeTempering( const std::vector <std::vector <T>> & A,
                              const std::vector <int> & n_prim,
                              const std::vector <int> & n_principal ) {
      n_shl = 0; n_fun = 0;
      max_l = n_prim.size() - 1;
      T lna; std::vector <T> pn;

      basis.clear();
      for(int l=0; l<=max_l; l++) for(int j=0; j<n_prim[l]; j++) {
        n_shl++; n_fun += l + l + 1;

        ints::LegendreP( A[0].size(), T(j+j) / T( std::max( n_prim[l] - 1, 1 ) ), pn );
        lna = T(0); for(uint k=0; k<std::min( (uint)A[0].size(), (uint)n_prim[l] ); k++) lna += A[l][k] * pn[k];
        basis.push_back( orbitalSTO <T> ( n_principal[l], l, exp(lna) ) );
      };
    };

    void printDetails( const std::string & name ) {
      std::cout << " " << name << " basis set data:" << std::endl;
      for(int n=0; n<n_shl; n++) {
        std::cout << "  " << basis[n].n << settings::shell_names[basis[n].l];
        std::cout << std::setw( settings::print_length + 10 ) << basis[n].a << std::endl;
      };
    };
};

void Initialize();

#endif
