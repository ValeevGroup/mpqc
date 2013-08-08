#include "chemistry/qc/ci/ci.h"
#include "util/misc/consumableresources.h"
#include <stdexcept>

//#define MPQC_PROFILE_ENABLE
#include "mpqc/utility/profile.hpp"

#include <chemistry/qc/nbody/nbwfn.h>
#include "mpqc/ci/integrals.hpp"
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
    range mo(config.core, config.core + config.orbitals);
    range ao(0, C.cols());
    
    if (mo.size() > C.rows())
      throw std::runtime_error("Number of orbitals too great");
    
    const auto &basis = wfn->basis();
    
    mpqc::ci::integrals(basis, C(mo, ao), wfn->integral()->hcore(), h);
    mpqc::ci::integrals(basis, C(mo, ao), wfn->integral()->electron_repulsion(), V);
  }

  mpqc::MPI::Comm comm = mpqc::MPI::Comm::World();

  std::string fname =
      ConsumableResources::get_default_instance()->disk_location() +
      SCFormIO::fileext_to_filename(".h5");

  std::auto_ptr<mpqc::File> file;
  file.reset(new mpqc::File(fname + "." + mpqc::string_cast(comm.rank())));

  std::vector<double> E;

  if (config.level > 0) {
      mpqc::ci::CI<mpqc::ci::Truncated> ci(config, comm, file->group());
      E = mpqc::ci::direct(ci, h, V);
  } else {
      mpqc::ci::CI<mpqc::ci::Full> ci(config, comm, file->group());
      E = mpqc::ci::direct(ci, h, V);
  }

  // delete file;

  return E;
}
