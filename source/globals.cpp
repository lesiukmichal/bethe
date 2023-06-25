int settings::num_threads_eigen;
int settings::num_threads_openmp;
int settings::job_mode;

int settings::orbital_momentum[settings::global_lmax_1] = {0};
int settings::n_orbital_legendre;
real_ settings::orbital_basis_params
  [settings::global_lmax_1][settings::global_kmax] = {0};

real_ settings::nuclear_charge;

int settings::n_closed;
int settings::n_open;
int settings::n_virtual;
int settings::n_electrons;

real_ settings::linear_dependent;
real_ settings::eri_threshold;

int settings::scf_maxit;
int settings::n_diis;
int settings::n_diis_turn_on;
real_ settings::scf_e_conv;
real_ settings::scf_d_conv;

real_ settings::rohf_coupling_f;
real_ settings::rohf_coupling_a;
real_ settings::rohf_coupling_b;

int powell::mx_macro;
int powell::mx_powell;
int powell::mx_search;
int powell::mx_gold;

real_ powell::stop_powell;
real_ powell::stop_gold;

int powell::nvar;
real_* powell::params;

std::vector <int> settings::extra_cphf_n;
std::vector <int> settings::extra_cphf_l;
std::vector <real_> settings::extra_cphf_a;
int settings::n_extra_shl_opt;

real_ settings::response_k;
int settings::n_points_bethe;

int settings::cphf_maxit;
real_ settings::cphf_conv;

int settings::n_points_bethe_small;
real_ settings::grid_start_small;
real_ settings::grid_step_small;

real_ coupling::binomial_table[coupling::n_binomial_max][coupling::n_binomial_max];

boost::multi_array<real_, 5> coupling::wigner3j(
  boost::extents[settings::global_lmax_1][settings::global_lmax_1][settings::global_lmax_1][settings::global_mmax_1][settings::global_mmax_1] );
