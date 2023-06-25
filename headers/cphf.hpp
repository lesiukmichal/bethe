#ifndef CPHF_H
#define CPHF_H

#include <string>
#include <vector>

namespace cphf {
  extern std_mtx UO_ix_final;
  extern std_mtx UC_ai_final;
  extern std_mtx UO_ax_final;

  real_ EngineCPHF( const bool & print = false );

  void FockDerivatives( const std_mtx & UO_ix,
                        const std_mtx & UC_ai,
                        const std_mtx & UO_ax,
                        const basisSTO <real_> & basis,
                        const ints::Integrals <real_> & ERI,
                        const std_mtx & C, std_mtx & DF_xi,
                        std_mtx & DF_ai, std_mtx & DF_ax );

  void TrialMultiplication( const std_mtx & UO_ix, const std_mtx & UC_ai, const std_mtx & UO_ax,
                            const std_mtx & DF_xi, const std_mtx & DF_ai, const std_mtx & DF_ax,
                            const std_mtx & FC_xa, const std_mtx & FC_xy,
                            const std_mtx & FC_xi, const std_mtx & FC_ab,
                            const std_mtx & FO_ij, const std_mtx & FO_ai,
                            const std_mtx & FO_ab, const std_mtx & FO_ix, const std_vec & E,
                            std_mtx & R_ix, std_mtx & R_ai, std_mtx & R_ax );





};

#endif
