//
// solvent.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/timer.h>
#include <util/misc/formio.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/wfn/solvent.h>

#define CLASSNAME BEMSolventH
#define VERSION 1
#define PARENTS public AccumH
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
BEMSolventH::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = AccumH::_castdown(cd);
  return do_castdowns(casts,cd);
}

BEMSolventH::BEMSolventH(const RefKeyVal&keyval):
  AccumH(keyval)
{
  charge_positions_ = 0;
  normals_ = 0;
  efield_dot_normals_ = 0;
  charges_ = 0;
  charges_n_ = 0;
  solvent_ = keyval->describedclassvalue("solvent");
  gamma_ = keyval->doublevalue("gamma");
  if (keyval->error() != KeyVal::OK) {
      RefUnits npm = new Units("dyne/cm");
      gamma_ = 72.75 * npm->to_atomic_units();
    }
}

BEMSolventH::BEMSolventH(StateIn&s):
  AccumH(s)
  maybe_SavableState(s)
{
  charge_positions_ = 0;
  normals_ = 0;
  efield_dot_normals_ = 0;
  charges_ = 0;
  charges_n_ = 0;
  escalar_ = 0;

  wfn_.restore_state(s);
  //solvent_.restore_state(s);
  abort();
}

BEMSolventH::~BEMSolventH()
{
  // just in case
  done();
}

void
BEMSolventH::save_data_state(StateOut&s)
{
  AccumH::save_data_state(s);

  wfn_.save_state(s);
  //solvent_.save_state(s);
  abort();
}

void
BEMSolventH::init(const RefWavefunction& wfn)
{
  tim_enter("solvent");
  tim_enter("init");
  wfn_ = wfn;
  // just in case
  done();
  solvent_->init();
  charge_positions_ = solvent_->alloc_charge_positions();
  normals_ = solvent_->alloc_normals();
  efield_dot_normals_ = solvent_->alloc_efield_dot_normals();
  charges_ = solvent_->alloc_charges();
  charges_n_ = solvent_->alloc_charges();

  // get the positions of the charges
  solvent_->charge_positions(charge_positions_);

  // get the surface normals
  solvent_->normals(normals_);
  tim_exit("init");
  tim_exit("solvent");
}

void
BEMSolventH::accum(const RefSymmSCMatrix& h)
{
  tim_enter("solvent");
  tim_enter("accum");
  int i,j;

  //// compute the polarization charges

  // compute the e-field at each point and dot with normals
  tim_enter("efield");
  int ncharge = solvent_->ncharge();
  RefEfieldDotVectorData efdn_dat = new EfieldDotVectorData;
  RefOneBodyInt efdn = wfn_->integral()->efield_dot_vector(efdn_dat);
  RefSCElementOp efdn_op = new OneBodyIntOp(efdn);
  RefSymmSCMatrix ao_density = wfn_->ao_density()->copy();
  RefSymmSCMatrix efdn_mat(ao_density->dim(), ao_density->kit());
  // for the scalar products, scale the density's off-diagonals by two
  ao_density->scale(2.0);
  ao_density->scale_diagonal(0.5);
  RefSCElementScalarProduct sp = new SCElementScalarProduct;
  RefSCElementOp2 generic_sp(sp.pointer());
  for (i=0; i<ncharge; i++) {
      efdn_dat->set_position(charge_positions_[i]);
      efdn_dat->set_vector(normals_[i]);
      efdn->reinitialize();
      efdn_mat->assign(0.0);
      efdn_mat->element_op(efdn_op);
      sp->init();
      efdn_mat->element_op(generic_sp, ao_density);
      efield_dot_normals_[i] = sp->result();
    }
  RefSCDimension aodim = ao_density.dim();
  RefSCMatrixKit aokit = ao_density.kit();
  ao_density = 0;
  efdn_mat = 0;
  tim_exit("efield");

  // compute a new set of charges
  tim_enter("charges");
  // electron contrib
  solvent_->compute_charges(efield_dot_normals_, charges_);
  double qeenc = solvent_->computed_enclosed_charge();
  // nuclear contrib
  for (i=0; i<ncharge; i++) {
      double nuc_efield[3];
      wfn_->molecule()->nuclear_efield(charge_positions_[i], nuc_efield);
      double tmp = 0.0;
      for (j=0; j<3; j++) {
          tmp += nuc_efield[j] * normals_[i][j];
        }
      efield_dot_normals_[i] = tmp;
    }
  solvent_->compute_charges(efield_dot_normals_, charges_n_);
  double qnenc = solvent_->computed_enclosed_charge();
  tim_exit("charges");

  // normalize the charges
  // e and n are independently normalized since the nature of the
  // errors in e and n are different: n error is just numerical and
  // e error is numerical plus diffuseness of electron distribution
  tim_enter("norm");
  // electron contrib
  solvent_->normalize_charge(-wfn_->nelectron(), charges_);
  // nuclear contrib
  solvent_->normalize_charge(wfn_->molecule()->nuclear_charge(), charges_n_);
  // sum the nuclear and electron contrib
  for (i=0; i<ncharge; i++) charges_[i] += charges_n_[i];
  tim_exit("norm");

  //// compute scalar contributions
  double A = solvent_->area();

  // the cavitation energy
  ecavitation_ = A * gamma_;

  // the dispersion energy
  double edisprep = 0.0;

  // compute the nuclear-surface interaction energy
  tim_enter("n-s");
  double enucsurf
      = solvent_->nuclear_interaction_energy(charge_positions_, charges_);
  tim_exit("n-s");

  //// compute one body contributions

  // compute the electron-surface interaction matrix elements
  tim_enter("e-s");
  RefPointChargeData pc_dat = new PointChargeData(ncharge,
                                                  charge_positions_, charges_);
  RefOneBodyInt pc = wfn_->integral()->point_charge(pc_dat);
  RefSCElementOp pc_op = new OneBodyIntOp(pc);

  // compute matrix elements in the ao basis
  RefSymmSCMatrix h_ao(aodim, aokit);
  h_ao.assign(0.0);
  h_ao.element_op(pc_op);
  // transform to the so basis and add to h
  RefSymmSCMatrix h_so = wfn_->integral()->petite_list()->to_SO_basis(h_ao);
  // why?
  h_so.scale(2.0);
  h->accumulate(h_so);
  // compute the contribution to the energy
  sp->init();
  RefSymmSCMatrix so_density = wfn_->density()->copy();
  // for the scalar products, scale the density's off-diagonals by two
  //so_density->scale(2.0);
  so_density->scale_diagonal(0.5);
  h_so->element_op(generic_sp, so_density);
  double eelecsurf = sp->result();
  tim_exit("e-s");

  // compute the surface-surface interaction energy
  double esurfsurf = -0.5*(eelecsurf+enucsurf);
  // (this can also be computed as below, but is much more expensive)
  //tim_enter("s-s");
  //double esurfsurf;
  //esurfsurf = solvent_->self_interaction_energy(charge_positions_, charges_);
  //tim_exit("s-s");

  escalar_ = enucsurf + esurfsurf + ecavitation_ + edisprep;

  cout << incindent;
  cout << node0 << indent
       << "Solvent: "
       << scprintf("q(e-enc)=%12.10f q(n-enc)=%12.10f", qeenc, qnenc)
       << endl;
  cout <<  incindent;
  cout << node0 << indent
       << scprintf("E(c)=%10.8f ", ecavitation_)
       << scprintf("E(disp-rep)=%10.8f", edisprep)
       << endl;
  cout << node0 << indent
       << scprintf("E(n-s)=%10.8f ", enucsurf)
       << scprintf("E(e-s)=%10.8f", eelecsurf)
       << scprintf("E(s-s)=%10.8f ", esurfsurf)
       << endl;
  cout << decindent;
  cout << decindent;

  tim_exit("accum");
  tim_exit("solvent");
}

void
BEMSolventH::done()
{
  solvent_->free_normals(normals_);
  normals_ = 0;
  solvent_->free_efield_dot_normals(efield_dot_normals_);
  efield_dot_normals_ = 0;
  solvent_->free_charges(charges_);
  solvent_->free_charges(charges_n_);
  charges_ = 0;
  charges_n_ = 0;
  solvent_->free_charge_positions(charge_positions_);
  charge_positions_ = 0;
  solvent_->done();
}

double
BEMSolventH::e()
{
  return escalar_;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
