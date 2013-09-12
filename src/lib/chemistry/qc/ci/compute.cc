#include "chemistry/qc/ci/ci.h"
#include "util/misc/consumableresources.h"
#include <stdexcept>
#include <cassert>

//#define MPQC_PROFILE_ENABLE
#include "mpqc/utility/profile.hpp"

#include <chemistry/qc/nbody/nbwfn.h>

#include "mpqc/ci/integrals.hpp"
#include "mpqc/ci/full/ci.hpp"
//#include "mpqc/ci/restricted.hpp"
#include "mpqc/ci/direct.hpp"
#include "mpqc/utility/string.hpp"

#include <memory>

std::vector<double> sc::CI::compute(const Ref<RefWavefunction> &wfn,
                                    const mpqc::ci::Config& config) {

  using mpqc::Vector;
  using mpqc::Matrix;

  ExEnv::out0() << indent << "Beginning CI\n";

  Matrix C = Matrix(wfn->orbs()->coefs()).transpose();

  Vector h;
  Matrix V;

  // compute molecular integrals
  {
    range occ_orbs(0, config.core + config.orbitals);
    range act_orbs(config.core, config.core + config.orbitals);
    range core_orbs(0, config.core);
    range ao(0, C.cols());
    
    assert(act_orbs.size() <= C.rows());
    
    const auto &basis = wfn->basis();
    
    // need all h elements in order to compute core constribution to the energy
    mpqc::ci::integrals(basis, C(occ_orbs, ao), wfn->integral()->hcore(), h);
    mpqc::ci::integrals(basis, C(act_orbs, ao), wfn->integral()->electron_repulsion(), V);

    if (config.core != 0) {
      // only need integrals with 2 core indices here ... this is too much work TODO fix this
      Matrix V_occ;
      mpqc::ci::integrals(basis, C(occ_orbs, ao),
                          wfn->integral()->electron_repulsion(), V_occ);

      // we need core contribution to the energy
      // E = \sum_i 2 h(i,i) + V(i,i;i,i) + \sum_i<j 4 V(i,i;j,j) - 2 V(i,j;i,j)
      config.e_core = 0.0;
      foreach(auto i, core_orbs) {
        const size_t i_offset = i*(i+1)/2;
        const size_t ii = i_offset + i;
        double E_core_i = 2 * h(ii);
        E_core_i += V_occ(ii,ii);
        foreach(auto j, range(0, i)) {
          const size_t ij = i_offset + j;
          const size_t j_offset = j*(j+1)/2;
          const size_t jj = j_offset + j;
          E_core_i += 4 * V_occ(ii, jj) - 2 * V_occ(ij, ij);
        }

        config.e_core += E_core_i;
      }
      //ExEnv::out0() << indent << scprintf("e_core = %20.15lf",config.e_core) << std::endl;
      // after this h only needed in active orbitals
      {
        const size_t nact = act_orbs.size();
        Vector h_act(nact*(nact+1)/2);
        size_t ij_act = 0;
        foreach(auto i, act_orbs) {
          const size_t i_offset = i*(i+1)/2;
          foreach(auto j, range(act_orbs.front(), i+1)) {
            h_act(ij_act) = h(i_offset + j);
            ++ij_act;
          }
        }

        h = h_act;
      }

      // h needs to include core contributions. For closed-shell case:
      // h'(p,q) = h(p,q) + \sum_i 2 V(p,q;i,i) - V(p,i;q,i)
      foreach(auto p, act_orbs) {
        const size_t p_offset = p*(p+1)/2;
        foreach(auto q, range(act_orbs.front(), p+1)) {
          const size_t pq = p_offset + q;
          const size_t q_offset = q*(q+1)/2;

          double h_pq_core = 0.0;
          foreach(auto i, core_orbs) {
            const size_t pi = p_offset + i;
            const size_t qi = q_offset + i;
            const size_t ii = i*(i+1)/2 + i;

            h_pq_core += 2 * V_occ(pq, ii) - V_occ(pi, qi);
          }

          const size_t p_act = p-config.core;
          const size_t q_act = q-config.core;
          const size_t pq_act = p_act*(p_act+1)/2 + q_act;
          h(pq_act) += h_pq_core;
        }
      }

    } // nfzc > 0

  }

  mpqc::MPI::Comm comm = mpqc::MPI::Comm::World();

  std::string fname =
      ConsumableResources::get_default_instance()->disk_location() +
      SCFormIO::fileext_to_filename(".h5");

  std::auto_ptr<mpqc::File> file;
  file.reset(new mpqc::File(fname + "." + mpqc::string_cast(comm.rank())));

  std::vector<double> E;

  if (config.rank > 0) {
      //mpqc::ci::RestrictedCI ci(config, comm, file->group());
      //E = mpqc::ci::direct(ci, h, V);
  } else {
      mpqc::ci::FullCI ci(config, comm, file->group());
      E = mpqc::ci::direct(ci, h, V);
  }

  // delete file;

  return E;
}
