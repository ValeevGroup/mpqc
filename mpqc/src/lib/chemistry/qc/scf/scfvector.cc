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

#include <util/misc/timer.h>
#include <util/misc/formio.h>

#include <math/scmat/offset.h>
#include <math/scmat/blocked.h>
#include <math/scmat/blkiter.h>

#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>

#include <chemistry/qc/basis/symmint.h>

#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/scf/scfops.h>
#include <chemistry/qc/scf/scflocal.h>

///////////////////////////////////////////////////////////////////////////

double
SCF::compute_vector(double& eelec)
{
  tim_enter("vector");
  int i;

  // reinitialize the extrapolation object
  extrap_->reinitialize();
  
  // create level shifter
  LevelShift *level_shift = new LevelShift(this);
  level_shift->reference();
  
  // calculate the core Hamiltonian
  hcore_ = core_hamiltonian();

  // add density independant contributions to Hcore
  accumdih_->accum(hcore_);

  // set up subclass for vector calculation
  init_vector();
  
  // calculate the nuclear repulsion energy
  double nucrep = molecule()->nuclear_repulsion_energy();
  ExEnv::out() << node0 << indent
       << scprintf("nuclear repulsion energy = %15.10f", nucrep)
       << endl << endl;

  RefDiagSCMatrix evals(oso_dimension(), basis_matrixkit());

  double delta = 1.0;
  int iter, iter_since_reset = 0;
  double accuracy = 1.0;
  for (iter=0; iter < maxiter_; iter++, iter_since_reset++) {
    // form the density from the current vector 
    tim_enter("density");
    delta = new_density();
    tim_exit("density");
    
    // check convergence
    if (delta < desired_value_accuracy()
        && accuracy < desired_value_accuracy()) break;

    // reset the density from time to time
    if (iter_since_reset && !(iter_since_reset%dens_reset_freq_)) {
      reset_density();
      iter_since_reset = 0;
    }
      
    // form the AO basis fock matrix & add density dependant H
    tim_enter("fock");
    double base_accuracy = delta;
    if (base_accuracy < desired_value_accuracy())
      base_accuracy = desired_value_accuracy();
    double new_accuracy = 0.01 * base_accuracy;
    if (new_accuracy > 0.001) new_accuracy = 0.001;
    if (iter == 0) accuracy = new_accuracy;
    else if (new_accuracy < accuracy) {
      accuracy = new_accuracy/10.0;
      if (iter_since_reset > 0) {
        reset_density();
        iter_since_reset = 0;
      }
    }
    ao_fock(accuracy);
    tim_exit("fock");

    // calculate the electronic energy
    eelec = scf_energy();
    double eother = 0.0;
    if (accumddh_.nonnull()) eother = accumddh_->e();
    ExEnv::out() << node0 << indent
         << scprintf("iter %5d energy = %15.10f delta = %10.5e",
                     iter+1, eelec+eother+nucrep, delta)
         << endl;

    // now extrapolate the fock matrix
    tim_enter("extrap");
    RefSCExtrapData data = extrap_data();
    RefSCExtrapError error = extrap_error();
    extrap_->extrapolate(data,error);
    data=0;
    error=0;
    tim_exit("extrap");

    // diagonalize effective MO fock to get MO vector
    tim_enter("evals");
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

    tim_exit("evals");

    if (debug_>0 && level_shift_ != 0.0) {
      evals.print("level shifted scf eigenvalues");
    }

    // now un-level shift eigenvalues
    level_shift->set_shift(-level_shift_);
    evals.element_op(level_shift);

    if (debug_>0) {
      evals.print("scf eigenvalues");
    }
    
    if (reset_occ_)
      set_occupations(evals);

    eff=0;
    
    if (debug_>1) {
      oso_scf_vector_.print("orthogonalized SO basis scf vector");
    }
  }
      
  eigenvalues_ = evals;
  eigenvalues_.computed() = 1;
  eigenvalues_.set_actual_accuracy(accuracy<delta?delta:accuracy);

  // search for HOMO and LUMO
  // first convert evals to something we can deal with easily
  BlockedDiagSCMatrix *evalsb = BlockedDiagSCMatrix::require_castdown(evals,
                                                 "SCF::compute_vector");
  
  CharacterTable ct = molecule()->point_group()->char_table();
  
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
        ExEnv::out() << node0 << endl
             << indent << ct.gamma(i).symbol() << endl << incindent;
        for (int m=0; m < nf; m++) {
          if (occupation(i,m) < 1e-8 && !print_all_evals_)
            break;
          ExEnv::out() << node0 << indent
               << scprintf("%5d %10.5f %10.5f", m+1, vals[m], occupation(i,m))
               << endl;
        }
        ExEnv::out() << node0 << decindent;
      }

      delete[] vals;
    }
  }

  if (homo_mo >= 0) {
    ExEnv::out() << node0 << endl << indent
         << scprintf("HOMO is %5d %3s = %10.6f",
                     homo_mo+1, 
                     ct.gamma(homo_ir).symbol(),
                     homo)
         << endl;
  }
  if (lumo_mo >= 0) {
    ExEnv::out() << node0 << indent
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

  // now clean up
  done_vector();
  hcore_ = 0;

  level_shift->dereference();
  delete level_shift;

  tim_exit("vector");
  //tim_print(0);

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

RefSCExtrapError
SCF::extrap_error()
{
  RefSymmSCMatrix mofock = effective_fock();
  
  RefSCElementOp op = new ExtrapErrorOp(this);
  mofock.element_op(op);
  
  RefSymmSCMatrix aoerror(so_dimension(), basis_matrixkit());
  aoerror.assign(0.0);
  aoerror.accumulate_transform(so_to_orthog_so().t()*oso_scf_vector_, mofock);
  mofock=0;

  RefSCExtrapError error = new SymmSCMatrixSCExtrapError(aoerror);
  aoerror=0;

  return error;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
