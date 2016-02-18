//
// scfvector.cc --- implementation of SCF::compute_vector
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include <sstream>

#include <util/misc/regtime.h>
#include <util/misc/formio.h>

#include <util/state/state_bin.h>
#include <util/group/mstate.h>

#include <math/scmat/offset.h>
#include <math/scmat/blocked.h>
#include <math/scmat/blkiter.h>

#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>

#include <chemistry/qc/basis/symmint.h>

#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/scf/scfops.h>
#include <chemistry/qc/scf/scflocal.h>

#include <errno.h>

#ifdef MPQC_NEW_FEATURES
#include <chemistry/qc/scf/iter_logger.h>
#include <util/misc/xmlwriter.h>
#endif

using namespace std;
using namespace sc;

#undef GENERALIZED_EIGENSOLVER

namespace sc {

///////////////////////////////////////////////////////////////////////////

extern "C" {
  void
  dsygv_(int *ITYPE, const char *JOBZ, const char *UPLO,
         int *N, double *A, int *LDA, double *B, int *LDB,
         double *W, double *WORK, int *LWORK, int *INFO);
}

void
SCF::savestate_to_file(const std::string &filename)
{
  std::string filename_to_delete = previous_savestate_file_;
  std::string filename_to_use;
  if (scf_grp_->me() == 0) {
    filename_to_use = filename;
    previous_savestate_file_ = filename;
  }
  else {
    filename_to_use = "/dev/null";
  }
  StateOutBin so(filename_to_use.c_str());
  save_state(this,so);
  so.close();
  if (filename_to_delete.size() > 0) {
    if (unlink(filename_to_delete.c_str())) {
      int unlink_errno = errno;
      ExEnv::out0() << indent
                    << "WARNING: SCF::compute_vector(): "
                    << "unlink of temporary checkpoint file"
                    << endl
                    << indent
                    << "         \"" << filename_to_delete << "\" "
                    << "failed with error: "
                    << strerror(unlink_errno)
                    << endl;
    }
  }
}

void
SCF::savestate_iter(int iter)
{
  char *ckptfile=0, *oldckptfile=0;
  const char *devnull=0;
  const char *filename=0;

  bool savestate = if_to_checkpoint();
  int savestate_freq = checkpoint_freq();

  if (savestate && ( (iter+1)%savestate_freq==0) ) {
    ostringstream sstr;
    const char *filename = checkpoint_file();
    sstr << filename << "." << iter+1 << ".tmp";
    savestate_to_file(sstr.str());
  }
}

void
SCF::iter_print(int iter,
                double energy,
                double delta,
                double walltime,
                std::ostream& os) {
  os << sc::indent
     << scprintf("iter %5d energy = %20.12f delta = %8.3e  (%8.2f sec)",
                 iter+1, energy, delta, walltime)
     << std::endl;
}

double
SCF::compute_vector(double& eelec, double nucrep)
{
  Timer tim("vector");
  int i;

  // reinitialize the extrapolation object
  extrap_->reinitialize(initial_extrap_data());

  // create level shifter
  LevelShift *level_shift = new LevelShift(this);
  level_shift->reference();

  // calculate the core Hamiltonian
  hcore_ = core_hamiltonian();

  // add density independent contributions to Hcore
  accumdih_->accum(hcore_);

  // set up subclass for vector calculation
  init_vector();

  RefDiagSCMatrix evals(oso_dimension(), basis_matrixkit());

  double delta = 1.0;
  int iter, iter_since_reset = 0;
  double accuracy = 1.0;

  ExEnv::out0() << indent
                << "Beginning iterations.  Basis is "
                << basis()->label() << '.' << std::endl;
  for (iter=0; iter < maxiter_; iter++, iter_since_reset++) {

#ifdef MPQC_NEW_FEATURES
    if(iter_log_.nonnull()) iter_log_->new_iteration();
#endif // MPQC_NEW_FEATURES

    const double wall_time_start = RegionTimer::get_wall_time();

    // form the density from the current vector
    tim.enter("density");
    delta = new_density();
    tim.exit("density");

    // reset the density from time to time
    if (iter_since_reset && !(iter_since_reset%dens_reset_freq_)) {
      reset_density();
      iter_since_reset = 0;
    }

    // form the AO basis fock matrix & add density dependant H
    tim.enter("fock");
    double base_accuracy = delta;
    if (base_accuracy < desired_value_accuracy())
      base_accuracy = desired_value_accuracy();
    double new_accuracy = 0.01 * base_accuracy;
    if (new_accuracy > 0.0001) new_accuracy = 0.0001;
    if (iter == 0) accuracy = new_accuracy;
    else if (new_accuracy < accuracy) {
      accuracy = new_accuracy/10.0;
      if (iter_since_reset > 0) {
        reset_density();
        iter_since_reset = 0;
      }
    }
    ExEnv::out0() << indent << "accuracy = " << accuracy << " new_accuracy = " << new_accuracy << std::endl;
    ao_fock(accuracy);
    tim.exit("fock");

    // calculate the electronic energy
    eelec = scf_energy();
    double eother = 0.0;
    if (accumddh_) eother = accumddh_->e();

#ifdef MPQC_NEW_FEATURES
    if(iter_log_.nonnull()) {
      using boost::property_tree::ptree;
      iter_log_->log_iter_misc([eelec,eother,nucrep,delta](ptree& parent, const XMLWriter& writer) {
        parent.put("energy", eelec+eother+nucrep);
        parent.put("delta", delta);
      });
    }
#endif // MPQC_NEW_FEATURES

    if(fake_scf_convergence_after_fock_build_ ||
        (fake_scf_convergence_after_n_iter_ > 0 && iter+1 >= fake_scf_convergence_after_n_iter_)
    ) {
      delta = 0.0;
      accuracy = 0.0;
      ExEnv::out0()
          << indent << "#=================== WARNING =====================#" << endl
          << indent << "# delta is being artificially set to 0.0 in order #" << endl
          << indent << "# to emulate convergence for testing purposes.    #" << endl
          << indent << "# Energy and density are not actually converged   #" << endl
          << indent << "# and ARE NOT CORRECT!!!  If you are seeing this  #" << endl
          << indent << "# message and did not expect to see it, please    #" << endl
          << indent << "# contact a developer immediately, as something   #" << endl
          << indent << "# has gone seriously wrong.  This feature is only #" << endl
          << indent << "# for testing and benchmarking purposes.  If you  #" << endl
          << indent << "# are doing real chemistry, or anything other     #" << endl
          << indent << "# than benchmarking SCF code, you should NEVER    #" << endl
          << indent << "# see this message.  DO NOT USE THE ENERGY OUTPUT #" << endl
          << indent << "# BELOW OR ANYTHING ELSE IN THIS FILE FOR         #" << endl
          << indent << "# ANYTHING OTHER THAN BENCHMARKING!!!!!!!         #" << endl
          << indent << "#=================================================#" << endl;
    }

    // check convergence
    if (delta < desired_value_accuracy()
        && accuracy < desired_value_accuracy()
        && iter+1 >= miniter_ ) {
      const double wall_time_end = RegionTimer::get_wall_time();
      iter_print(iter,
                 eelec+eother+nucrep, delta, wall_time_end - wall_time_start,
                 ExEnv::out0());
      break;
    }

    // now extrapolate the fock matrix
    tim.enter("extrap");
    Ref<SCExtrapData> data = extrap_data();
    Ref<SCExtrapError> error = extrap_error();
    extrap_->extrapolate(data,error);
    data=0;
    error=0;
    tim.exit("extrap");

#ifdef GENERALIZED_EIGENSOLVER
    // Get the fock matrix and overlap in the SO basis.  The fock matrix
    // used here works for CLOSED SHELL ONLY.
    RefSymmSCMatrix bfmatref = fock(0);
    RefSymmSCMatrix bsmatref = overlap();
    BlockedSymmSCMatrix *bfmat
      = dynamic_cast<BlockedSymmSCMatrix*>(bfmatref.pointer());
    BlockedSymmSCMatrix *bsmat
      = dynamic_cast<BlockedSymmSCMatrix*>(bsmatref.pointer());
    BlockedDiagSCMatrix *bevals
      = dynamic_cast<BlockedDiagSCMatrix*>(evals.pointer());
    BlockedSCMatrix *bvec
      = dynamic_cast<BlockedSCMatrix*>(oso_scf_vector_.pointer());

    ExEnv::out0() << indent
                  << "solving generalized eigenvalue problem" << endl;

    for (int iblock=0; iblock<bfmat->nblocks(); iblock++) {
      RefSymmSCMatrix fmat = bfmat->block(iblock);
      RefSymmSCMatrix smat = bsmat->block(iblock);
      RefDiagSCMatrix evalblock = bevals->block(iblock);
      RefSCMatrix oso_scf_vector_block = bvec->block(iblock);
      int nbasis = fmat.dim().n();
      int nbasis2 = nbasis*nbasis;

      if (!nbasis) continue;

      // Convert to the lapack storage format.
      double *fso = new double[nbasis2];
      double *sso = new double[nbasis2];
      int ij=0;
      for (int i=0; i<nbasis; i++) {
        for (int j=0; j<nbasis; j++,ij++) {
          fso[ij] = fmat(i,j);
          sso[ij] = smat(i,j);
        }
      }

      // solve generalized eigenvalue problem with DSYGV
      int itype = 1;
      double *epsilon = new double[nbasis];
      int lwork = -1;
      int info;
      double optlwork;
      dsygv_(&itype,"V","U",&nbasis,fso,&nbasis,sso,&nbasis,
             epsilon,&optlwork,&lwork,&info);
      if (info) {
        ExEnv::outn() << "dsygv could not determine work size: info = "
                      << info << endl;
        abort();
      }
      lwork = (int)optlwork;
      double *work = new double[lwork];
      dsygv_(&itype,"V","U",&nbasis,fso,&nbasis,sso,&nbasis,
             epsilon,work,&lwork,&info);
      if (info) {
        ExEnv::outn() << "dsygv could not diagonalize matrix: info = "
                      << info << endl;
        abort();
      }
      double *z = fso; // the vector is placed in fso

      // make sure everyone agrees on the new arrays
      scf_grp_->bcast(z, nbasis2);
      scf_grp_->bcast(epsilon, nbasis);

      evalblock->assign(epsilon);
      oso_scf_vector_block->assign(z);

      // cleanup
      delete[] fso;
      delete[] work;
      delete[] sso;
      delete[] epsilon;
    }

    oso_scf_vector_ = (oso_scf_vector_ * so_to_orthog_so_inverse()).t();
#else
    // diagonalize effective MO fock to get MO vector
    tim.enter("evals");
    RefSCMatrix nvector(oso_dimension(),oso_dimension(),basis_matrixkit());

    RefSymmSCMatrix eff = effective_fock();

    // level shift effective fock
    level_shift->set_shift(level_shift_);
    eff.element_op(level_shift);

    if (debug_>1) {
      eff.print("effective 1 body hamiltonian in current mo basis");
    }

    // transform eff to the oso basis to diagonalize it
    RefSymmSCMatrix oso_eff(oso_dimension(), basis_matrixkit());
    oso_eff.assign(0.0);
    oso_eff.accumulate_transform(oso_scf_vector_,eff);
    eff = 0;
    oso_eff.diagonalize(evals, oso_scf_vector_);

    tim.exit("evals");

    if (debug_>0 && level_shift_ != 0.0) {
      evals.print("level shifted scf eigenvalues");
    }

    // now un-level shift eigenvalues
    level_shift->set_shift(-level_shift_);
    evals.element_op(level_shift);
#endif

    current_evals_ = evals;

    if (debug_>0) {
      evals.print("scf eigenvalues");
    }

#ifdef MPQC_NEW_FEATURES
    if (iter_log_.nonnull()){
      iter_log_->log_evals(evals.copy());

      if (iter_log_->log_coeffs_enabled()){
        // Get the mospace_ object from the OneBodyWavefunction
        Ref<PetiteList> pl = integral()->petite_list();
        RefSCMatrix aocoefs_full = pl->evecs_to_AO_basis((oso_scf_vector_.t() * so_to_orthog_so()).t());
        Ref<OrbitalSpace> mospace = new OrbitalSpace("p", "energy-ordered MOs to evaluate",
                                                     aocoefs_full, basis(), integral(),
                                                     evals,
                                                     0, 0,
                                                     OrbitalSpace::energy);
        iter_log_->log_coeffs(
            mospace->coefs_nb()
            //(oso_scf_vector_.t() * so_to_orthog_so() * (integral()->petite_list()->aotoso()).t()).t()
        );
      }
    }
#endif // MPQC_NEW_FEATURES

    if (reset_occ_)
      set_occupations(evals);

    if (debug_>1) {
      oso_scf_vector_.print("OSO basis scf vector");
      (oso_scf_vector_.t()*oso_scf_vector_).print(
        "vOSO.t()*vOSO",ExEnv::out0(),14);
    }

    savestate_iter(iter);

    const double wall_time_end = RegionTimer::get_wall_time();
    iter_print(iter,
               eelec+eother+nucrep, delta, wall_time_end - wall_time_start,
               ExEnv::out0());
  }

  eigenvalues_ = evals;
  eigenvalues_.computed() = 1;
  eigenvalues_.set_actual_accuracy(accuracy<delta?delta:accuracy);

  // search for HOMO and LUMO
  // first convert evals to something we can deal with easily
  BlockedDiagSCMatrix *evalsb = require_dynamic_cast<BlockedDiagSCMatrix*>(evals,
                                                 "SCF::compute_vector");

  CharacterTable ct = molecule()->point_group()->char_table();

  struct Orbital {

  };

  int homo_ir=0, lumo_ir=0;
  int homo_mo = -1, lumo_mo = -1;
  double homo=-1e99, lumo=1e99;
  for (i=0; i < oso_dimension()->blocks()->nblock(); i++) {
    int nf=oso_dimension()->blocks()->size(i);
    if (nf) {
      double *vals = new double[nf];
      evalsb->block(i)->convert(vals);

      for (int mo=0; mo < nf; mo++) {
        if (occupation(i, mo) > 0.0) {
          if (vals[mo] > homo) {
            homo = vals[mo];
            homo_ir = i;
            homo_mo = mo;
          }
        } else {
          if (vals[mo] < lumo) {
            lumo = vals[mo];
            lumo_ir = i;
            lumo_mo = mo;
          }
        }
      }

      if (print_all_evals_ || print_occ_evals_) {
        ExEnv::out0() << endl
             << indent << ct.gamma(i).symbol() << endl << incindent;
        for (int m=0; m < nf; m++) {
          if (occupation(i,m) < 1e-8 && !print_all_evals_)
            break;
          ExEnv::out0() << indent
               << scprintf("%5d %10.5f %10.5f", m+1, vals[m], occupation(i,m))
               << endl;
        }
        ExEnv::out0() << decindent;
      }

      delete[] vals;
    }
  }

  if (homo_mo >= 0) {
    ExEnv::out0() << endl << indent
         << scprintf("HOMO is %5d %3s = %10.6f",
                     homo_mo+1,
                     ct.gamma(homo_ir).symbol(),
                     homo)
         << endl;
  }
  if (lumo_mo >= 0) {
    ExEnv::out0() << indent
         << scprintf("LUMO is %5d %3s = %10.6f",
                     lumo_mo+1,
                     ct.gamma(lumo_ir).symbol(),
                     lumo)
         << endl;
  }

  // free up evals
  evals = 0;

  oso_eigenvectors_ = oso_scf_vector_;
  oso_eigenvectors_.computed() = 1;
  oso_eigenvectors_.set_actual_accuracy(delta);
  // Checkpoint wavefunction, if needed, so that if converged
  // on the last iteration then the wavefunction is marked as computed
  if (if_to_checkpoint()) {
    const char *checkpoint_filename = checkpoint_file();
    std::string state_filename = checkpoint_filename;
    state_filename += ".tmp";
    savestate_to_file(state_filename);
  }

  // now clean up
  done_vector();
  hcore_ = 0;

  level_shift->dereference();
  delete level_shift;

  tim.exit("vector");
  //tim.print();

  return delta;
}

////////////////////////////////////////////////////////////////////////////

class ExtrapErrorOp : public BlockedSCElementOp {
  private:
    SCF *scf_;

  public:
    ExtrapErrorOp(SCF *s) : scf_(s) {}
    ~ExtrapErrorOp() {}

    int has_side_effects() { return 1; }

    void process(SCMatrixBlockIter& bi) {
      int ir=current_block();

      for (bi.reset(); bi; bi++) {
        int i=bi.i();
        int j=bi.j();
        if (scf_->occupation(ir,i) == scf_->occupation(ir,j))
          bi.set(0.0);
      }
    }
};

Ref<SCExtrapError>
SCF::extrap_error()
{
  RefSymmSCMatrix mofock = effective_fock();

  Ref<SCElementOp> op = new ExtrapErrorOp(this);
  mofock.element_op(op);

  RefSymmSCMatrix aoerror(so_dimension(), basis_matrixkit());
  aoerror.assign(0.0);
  aoerror.accumulate_transform(so_to_orthog_so().t()*oso_scf_vector_, mofock);
  mofock=0;

  Ref<SCExtrapError> error = new SymmSCMatrixSCExtrapError(aoerror);
  aoerror=0;

  return error;
}

/////////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
