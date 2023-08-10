/** @file cphf.hpp
 *  @brief Response equations solver.
 *  @author Michal Lesiuk
 */
#ifndef CPHF_H
#define CPHF_H

#include <string>
#include <vector>

namespace cphf {
  extern std_mtx UO_ix_final;
  extern std_mtx UC_ai_final;
  extern std_mtx UO_ax_final;

/** @brief The main driver for solving the coupled-perturbed
 *         Hartree-Fock equations, Eqs. (20)-(22) in the paper.
 *  @param[in] print: turn on additional output (default = false)
 */
  real_ EngineCPHF( const bool & print = false );

/** @brief Calculates the derivatives of the Fock matrices, defined
 *         by Eqs. (21)-(22) in the paper.
 *  @param[in] UO_ix: CPHF coefficients, closed-open block
 *  @param[in] UO_ai: CPHF coefficients, virtual-closed block
 *  @param[in] UO_ax: CPHF coefficients, virtual-open block
 *  @param[in] basis: basis set class
 *  @param[in] ERI: two-electron integrals class
 *  @param[in] C: Hartree-Fock orbitals coefficients matrix
 *  @param[out] DF_xi: Fock matrix derivative, open-closed block
 *  @param[out] DF_ai: Fock matrix derivative, virtual-closed block
 *  @param[out] DF_ax: Fock matrix derivative, virtual-open block
 */
  void FockDerivatives( const std_mtx & UO_ix,
                        const std_mtx & UC_ai,
                        const std_mtx & UO_ax,
                        const basisSTO <real_> & basis,
                        const ints::Integrals <real_> & ERI,
                        const std_mtx & C, std_mtx & DF_xi,
                        std_mtx & DF_ai, std_mtx & DF_ax );

/** @brief Essentially, calculates the left-hand side of Eq. (20) for given
 *         values of the CPHF expansion coefficients.
 *  @param[in] UO_ix: CPHF coefficients, closed-open block
 *  @param[in] UO_ai: CPHF coefficients, virtual-closed block
 *  @param[in] UO_ax: CPHF coefficients, virtual-open block
 *  @param[in] DF_xi: Fock matrix derivative, open-closed block
 *  @param[in] DF_ai: Fock matrix derivative, virtual-closed block
 *  @param[in] DF_ax: Fock matrix derivative, virtual-open block
 *  @param[in] FC_xa: closed Fock matrix, open-virtual block
 *  @param[in] FC_xy: closed Fock matrix, open-open block
 *  @param[in] FC_xi: closed Fock matrix, open-closed block
 *  @param[in] FC_ab: closed Fock matrix, virtual-virtual block
 *  @param[in] FO_ij: open Fock matrix, closed-closed block
 *  @param[in] FO_ai: open Fock matrix, virtual-closed block
 *  @param[in] FO_ab: open Fock matrix, virtual-virtual block
 *  @param[in] FO_ix: open Fock matrix, closed-openblock
 *  @param[in] E: orbital energies vector
 *  @param[out] R_ix: result, closed-open block
 *  @param[out] R_ai: result, virtual-closed block
 *  @param[out] R_ax: result, virtual-open block
 */
  void TrialMultiplication( const std_mtx & UO_ix, const std_mtx & UC_ai, const std_mtx & UO_ax,
                            const std_mtx & DF_xi, const std_mtx & DF_ai, const std_mtx & DF_ax,
                            const std_mtx & FC_xa, const std_mtx & FC_xy,
                            const std_mtx & FC_xi, const std_mtx & FC_ab,
                            const std_mtx & FO_ij, const std_mtx & FO_ai,
                            const std_mtx & FO_ab, const std_mtx & FO_ix, const std_vec & E,
                            std_mtx & R_ix, std_mtx & R_ai, std_mtx & R_ax );
};

#endif
