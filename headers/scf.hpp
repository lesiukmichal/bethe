#ifndef SCF_H
#define SCF_H

#include <string>
#include <vector>

namespace scf {
  extern std_mtx final_C;
  extern std_vec final_E;

  real_ EngineSCF( const bool & print = false );

  inline std_mtx ErrorDIIS( const std_mtx & F,
                            const std_mtx & D,
                            const std_mtx & S ) {
    return F * D * S - S * D * F;
  };

  void TotalFockOperator( const std_mtx & H, const std_mtx & C, const std_mtx & S,
                          const std_mtx & J_closed, const std_mtx & K_closed,
                          const std_mtx & J_open, const std_mtx & K_open,
                          std_mtx & F_closed, std_mtx & F_open, std_mtx & F_total );

  void SolveDIIS( const std::vector <std_mtx> & store_f_diis,
                  const std::vector <std_mtx> & store_d_diis,
                  const std::vector <std_mtx> & store_e_diis,
                  std_mtx & F, std_mtx & D );
};

#endif
