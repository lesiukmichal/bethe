#ifndef SCF_H
#define SCF_H

#include <string>
#include <vector>

namespace scf {
  extern std_mtx final_C;
  extern std_vec final_E;

/** @brief The main driver for solving the restricted open-shell
 *         Hartree-Fock equations, Eqs. (9)-(13) in the paper.
 *  @param[in] print: turn on additional output (default = false)
 *
 *  Details are in a classic paper:\n
 *  C. C. J. Roothaan, Rev. Mod. Phys. 32, 179 (1960).
 */
  real_ EngineSCF( const bool & print = false );

/** @brief Calculate the DIIS error matrix (in AO basis).
 *  @param[in] F: Fock matrix
 *  @param[in] D: density matrix
 *  @param[in] S: overlap matrix
 */
  inline std_mtx ErrorDIIS( const std_mtx & F,
                            const std_mtx & D,
                            const std_mtx & S ) {
    return F * D * S - S * D * F;
  };

/** @brief Construct the total closed and open Fock operators,
 *         Eqs. (10) and (11) in the paper.
 *  @param[in] H: one-electron integrals matrix (bare nuclei Hamiltonian)
 *  @param[in] C: Hartree-Fock orbitals coefficients matrix
 *  @param[in] S: overlap matrix (AO basis)
 *  @param[in] J_closed: closed Coulomb operator, Eq. (12)
 *  @param[in] K_closed: closed exchange operator, Eq. (13)
 *  @param[in] J_open: open Coulomb operator
 *  @param[in] K_open: open exchange operator
 *  @param[out] F_closed: closed Fock operator, Eq. (10)
 *  @param[out] F_open: open Fock operator, Eq. (11)
 *  @param[out] F_total: total Fock operator, Eq. (9)
 */
  void TotalFockOperator( const std_mtx & H, const std_mtx & C, const std_mtx & S,
                          const std_mtx & J_closed, const std_mtx & K_closed,
                          const std_mtx & J_open, const std_mtx & K_open,
                          std_mtx & F_closed, std_mtx & F_open, std_mtx & F_total );

/** @brief Solve the DIIS equations to obtain the next guess for
 *         the density matrix (convergence acceleration)
 *  @param[in] store_f_diis: previous Fock matrices
 *  @param[in] store_d_diis: previous density matrices
 *  @param[in] store_e_diis: previous DIIS error matrices
 *  @param[out] F: next Fock matrix
 *  @param[out] D: next density matrix
 */
  void SolveDIIS( const std::vector <std_mtx> & store_f_diis,
                  const std::vector <std_mtx> & store_d_diis,
                  const std::vector <std_mtx> & store_e_diis,
                  std_mtx & F, std_mtx & D );
};

#endif
