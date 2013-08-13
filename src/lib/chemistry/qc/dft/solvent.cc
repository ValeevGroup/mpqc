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

#include <util/misc/regtime.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/dft/solvent.h>

#include <math/isosurf/volume.h>
#include <chemistry/qc/dft/integrator.h>
#include <chemistry/qc/dft/functional.h>

#include <iomanip>

using namespace std;
using namespace sc;

namespace sc {

//. The \clsnm{NElFunctional} computes the number of electrons.
//. It is primarily for testing the integrator.
class NElInShapeFunctional: public DenFunctional {
  private:
    Ref<Volume> vol_;
    double isoval_;
  public:
    NElInShapeFunctional(const Ref<Volume> &, double);
    ~NElInShapeFunctional();

    void point(const PointInputData&, PointOutputData&);
};

/////////////////////////////////////////////////////////////////////////////
// NElFunctional

static ClassDesc NElInShapeFunctional_cd(
  typeid(NElInShapeFunctional),"NElInShapeFunctional",1,"public DenFunctional",
  0, 0, 0);

NElInShapeFunctional::NElInShapeFunctional(const Ref<Volume>& vol,
                                           double isoval)
{
  vol_ = vol;
  isoval_ = isoval;
}

NElInShapeFunctional::~NElInShapeFunctional()
{
}

void
NElInShapeFunctional::point(const PointInputData &id,
                            PointOutputData &od)
{
  vol_->set_x(id.r);
  if (vol_->value() <= isoval_) {
      od.energy = id.a.rho + id.b.rho;
    }
  else {
      od.energy = 0.0;
    }
}

/////////////////////////////////////////////////////////////////////////////

static ClassDesc BEMSolventH_cd(
  typeid(BEMSolventH),"BEMSolventH",1,"public AccumH",
  0, create<BEMSolventH>, create<BEMSolventH>);

BEMSolventH::BEMSolventH(const Ref<KeyVal>&keyval):
  AccumH(keyval)
{
  charge_positions_ = 0;
  normals_ = 0;
  efield_dot_normals_ = 0;
  charges_ = 0;
  charges_n_ = 0;
  solvent_ << keyval->describedclassvalue("solvent");
  gamma_ = keyval->doublevalue("gamma");
  if (keyval->error() != KeyVal::OK) {
      Ref<Units> npm = new Units("dyne/cm");
      gamma_ = 72.75 * npm->to_atomic_units();
    }
  // If onebody add a term to the one body hamiltonian, h.
  // Otherwise the energy contribution is scalar.
  onebody_ = keyval->booleanvalue("onebody");
  if (keyval->error() != KeyVal::OK) onebody_ = 1;
  // Normalize the charges if normalize_q is set.
  normalize_q_ = keyval->booleanvalue("normalize_q");
  if (keyval->error() != KeyVal::OK) normalize_q_ = 1;
  // Compute separately contributes to the energy from surfaces
  // charges induced by the nuclear and electronic charge densities.
  separate_surf_charges_ = keyval->booleanvalue("separate_surf_charges");
  if (keyval->error() != KeyVal::OK) separate_surf_charges_ = 0;
  // The Cammi-Tomasi Y term is set equal to the J term (as it formally is).
  y_equals_j_ = keyval->booleanvalue("y_equals_j");
  if (keyval->error() != KeyVal::OK) y_equals_j_ = 0;
  // As a test, integrate the number of electrons inside the surface.
  integrate_nelectron_ = keyval->booleanvalue("integrate_nelectron");
  if (keyval->error() != KeyVal::OK) integrate_nelectron_ = 0;
}

BEMSolventH::BEMSolventH(StateIn&s):
  SavableState(s),
  AccumH(s)
{
  charge_positions_ = 0;
  normals_ = 0;
  efield_dot_normals_ = 0;
  charges_ = 0;
  charges_n_ = 0;
  escalar_ = 0;

  wfn_ << SavableState::restore_state(s);
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

  SavableState::save_state(wfn_.pointer(),s);
  //solvent_.save_state(s);
  abort();
}

void
BEMSolventH::init(const Ref<Wavefunction>& wfn)
{
  Timer tim("solvent");
  tim.enter("init");
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

  if (integrate_nelectron_) {
      Ref<DenIntegrator> integrator = new RadialAngularIntegrator();
      Ref<DenFunctional> functional
          = new NElInShapeFunctional(solvent_->surface()->volume_object(),
                                     solvent_->surface()->isovalue());
      integrator->init(wfn_);
      integrator->integrate(functional);
      integrator->done();
      ExEnv::out0() << indent
           << scprintf("N(e) in isosurf = %12.8f", integrator->value())
           << endl;
    }

  edisprep_ = solvent_->disprep();

  tim.exit("init");
  tim.exit("solvent");
}

// This adds J + X to h, where J and X are the matrices defined
// by Canni and Tomasi, J Comp Chem, 16(12), 1457, 1995.
// The resulting SCF free energy expression is
//    G = 1/2TrP[h' + F'] + Une + Unn + Vnn
//        -1/2(Uee+Uen+Une+Unn)
// which in the Canni-Tomasi notation is
//      = 1/2TrP[h+1/2(X+J+Y+G)] + Vnn + 1/2Unn
// which is identical to the Canni-Tomasi energy expression.
// My Fock matrix is
//       F' = h + J + X + G
// while the Canni-Tomasi Fock matrix is F' = h + 1/2(J+Y) + X + G.
// However, since J = Y formally, (assuming no numerical errors
// and all charge is enclosed, Canni-Tomasi use F' = h + J + X + G
// to get better numerical results.
//
//   If the y_equals_j option is true, the energy expression used
// here is G = 1/2TrP[h+1/2(X+2J+G)] + Vnn + 1/2Unn, however, THIS
// IS NOT RECOMMENDED.
void
BEMSolventH::accum(const RefSymmSCMatrix& h)
{
  Timer tim("solvent");
  tim.enter("accum");
  int i,j;

  //// compute the polarization charges

  // compute the e-field at each point and dot with normals
  tim.enter("efield");
  int ncharge = solvent_->ncharge();
  Ref<EfieldDotVectorData> efdn_dat = new EfieldDotVectorData;
  Ref<OneBodyInt> efdn = wfn_->integral()->efield_dot_vector(efdn_dat);
  Ref<SCElementOp> efdn_op = new OneBodyIntOp(efdn);
  RefSymmSCMatrix ao_density = wfn_->ao_density()->copy();
  RefSymmSCMatrix efdn_mat(ao_density->dim(), ao_density->kit());
  // for the scalar products, scale the density's off-diagonals by two
  ao_density->scale(2.0);
  ao_density->scale_diagonal(0.5);
  Ref<SCElementScalarProduct> sp = new SCElementScalarProduct;
  Ref<SCElementOp2> generic_sp(sp.pointer());
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
  Ref<SCMatrixKit> aokit = ao_density.kit();
  ao_density = 0;
  efdn_mat = 0;
  tim.exit("efield");

  // compute a new set of charges
  tim.enter("charges");
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
  tim.exit("charges");

  // normalize the charges
  // e and n are independently normalized since the nature of the
  // errors in e and n are different: n error is just numerical and
  // e error is numerical plus diffuseness of electron distribution
  if (normalize_q_) {
      tim.enter("norm");
      // electron contrib
      solvent_->normalize_charge(-wfn_->nelectron(), charges_);
      // nuclear contrib
      solvent_->normalize_charge(wfn_->molecule()->total_charge(),
                                 charges_n_);
      tim.exit("norm");
    }
  // sum the nuclear and electron contrib
  for (i=0; i<ncharge; i++) charges_[i] += charges_n_[i];

  //// compute scalar contributions
  double A = solvent_->area();

  // the cavitation energy
  ecavitation_ = A * gamma_;

  // compute the nuclear-surface interaction energy
  tim.enter("n-s");
  enucsurf_
      = solvent_->nuclear_interaction_energy(charge_positions_, charges_);
  tim.exit("n-s");

  double enqn = 0.0, enqe = 0.0;
  if (y_equals_j_ || separate_surf_charges_) {
      tim.enter("n-qn");
      enqn = solvent_->nuclear_interaction_energy(charge_positions_,
                                                  charges_n_);
      enqe = enucsurf_ - enqn;
      tim.exit("n-qn");
    }

  //// compute one body contributions

  // compute the electron-surface interaction matrix elements
  tim.enter("e-s");
  Ref<PointChargeData> pc_dat = new PointChargeData(ncharge,
                                                  charge_positions_, charges_);
  Ref<OneBodyInt> pc = wfn_->integral()->point_charge(pc_dat);
  Ref<SCElementOp> pc_op = new OneBodyIntOp(pc);

  // compute matrix elements in the ao basis
  RefSymmSCMatrix h_ao(aodim, aokit);
  h_ao.assign(0.0);
  h_ao.element_op(pc_op);
  // transform to the so basis and add to h
  RefSymmSCMatrix h_so = wfn_->integral()->petite_list()->to_SO_basis(h_ao);
  if (onebody_) h->accumulate(h_so);
  // compute the contribution to the energy
  sp->init();
  RefSymmSCMatrix so_density = wfn_->density()->copy();
  // for the scalar products, scale the density's off-diagonals by two
  so_density->scale(2.0);
  so_density->scale_diagonal(0.5);
  h_so->element_op(generic_sp, so_density);
  eelecsurf_ = sp->result();
  tim.exit("e-s");

  double eeqn = 0.0, eeqe = 0.0;
  if (y_equals_j_ || separate_surf_charges_) {
      tim.enter("e-qn");
      pc_dat = new PointChargeData(ncharge, charge_positions_, charges_n_);
      pc = wfn_->integral()->point_charge(pc_dat);
      pc_op = new OneBodyIntOp(pc);

      // compute matrix elements in the ao basis
      h_ao.assign(0.0);
      h_ao.element_op(pc_op);
      // transform to the so basis
      h_so = wfn_->integral()->petite_list()->to_SO_basis(h_ao);
      // compute the contribution to the energy
      sp->init();
      h_so->element_op(generic_sp, so_density);
      eeqn = sp->result();
      eeqe = eelecsurf_ - eeqn;
      tim.exit("e-qn");
    }

  if (y_equals_j_) {
      // Remove the y term (enqe) and add the j term (eeqn).  Formally,
      // they are equal, but they are not because some e-density is outside
      // the surface and because of the numerical approximations.
      enucsurf_ += eeqn - enqe;
    }

  // compute the surface-surface interaction energy
  esurfsurf_ = -0.5*(eelecsurf_+enucsurf_);
  // (this can also be computed as below, but is much more expensive)
  //tim.enter("s-s");
  //double esurfsurf_;
  //esurfsurf_ = solvent_->self_interaction_energy(charge_positions_, charges_);
  //tim.exit("s-s");

  escalar_ = enucsurf_ + esurfsurf_ + ecavitation_ + edisprep_;
  // NOTE: SCF currently only adds h_so to the Fock matrix
  // so a term is missing in the energy.  This term is added here
  // and when SCF is fixed, should no longer be included.
  if (onebody_) escalar_ += 0.5 * eelecsurf_;

  if (!onebody_) escalar_ += eelecsurf_;

  ExEnv::out0() << incindent;
  ExEnv::out0() << indent
       << "Solvent: "
       << scprintf("q(e-enc)=%12.10f q(n-enc)=%12.10f", qeenc, qnenc)
       << endl;
  ExEnv::out0() << incindent;
  if (separate_surf_charges_) {
      ExEnv::out0() << indent
           << scprintf("E(n-qn)=%10.8f ", enqn)
           << scprintf("E(n-qe)=%10.8f", enqe)
           << endl;
      ExEnv::out0() << indent
           << scprintf("E(e-qn)=%10.8f ", eeqn)
           << scprintf("E(e-qe)=%10.8f", eeqe)
           << endl;
      //ExEnv::out0() << indent
      //     << scprintf("DG = %12.8f ", 0.5*627.51*(enqn+enqe+eeqn+eeqe))
      //     << scprintf("DG(Y=J) = %12.8f", 0.5*627.51*(enqn+2*eeqn+eeqe))
      //     << endl;
    }
  ExEnv::out0() << indent
       << scprintf("E(c)=%10.8f ", ecavitation_)
       << scprintf("E(disp-rep)=%10.8f", edisprep_)
       << endl;
  ExEnv::out0() << indent
       << scprintf("E(n-s)=%10.8f ", enucsurf_)
       << scprintf("E(e-s)=%10.8f ", eelecsurf_)
       << scprintf("E(s-s)=%10.8f ", esurfsurf_)
       << endl;
  ExEnv::out0() << decindent;
  ExEnv::out0() << decindent;

  tim.exit("accum");
  tim.exit("solvent");
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

void
BEMSolventH::print_summary()
{
  Ref<Units> unit = new Units("kcal/mol");
  ExEnv::out0() << endl;
  ExEnv::out0() << "Summary of solvation calculation:" << endl;
  ExEnv::out0() << "_______________________________________________" << endl;
  ExEnv::out0() << endl;
  ExEnv::out0().setf(ios::scientific,ios::floatfield); // use scientific format
  ExEnv::out0().precision(5);
  ExEnv::out0() << indent << "E(nuc-surf):              " 
       << setw(12) << setfill(' ')
       << enucsurf_*unit->from_atomic_units() << " kcal/mol" << endl; 
  ExEnv::out0() << indent << "E(elec-surf):             " 
       << setw(12) << setfill(' ')
       << eelecsurf_*unit->from_atomic_units() << " kcal/mol" << endl; 
  ExEnv::out0() << indent << "E(surf-surf):             " 
       << setw(12) << setfill(' ')
       << esurfsurf_*unit->from_atomic_units() << " kcal/mol" << endl; 
  ExEnv::out0() << indent << "Electrostatic energy:     " 
       << setw(12) << setfill(' ')
       << (enucsurf_+eelecsurf_+esurfsurf_)*unit->from_atomic_units()
       << " kcal/mol" << endl; 
  ExEnv::out0() << "_______________________________________________" << endl;
  ExEnv::out0() << endl;
  ExEnv::out0() << indent << "E(cav):                   " 
       << setw(12) << setfill(' ')
       << ecavitation_*unit->from_atomic_units() << " kcal/mol" << endl; 
  ExEnv::out0() << indent << "E(disp):                  " 
       << setw(12) << setfill(' ')
       << solvent_->disp()*unit->from_atomic_units() << " kcal/mol" << endl; 
  ExEnv::out0() << indent << "E(rep):                   " 
       << setw(12) << setfill(' ')
       << solvent_->rep()*unit->from_atomic_units() << " kcal/mol" << endl; 
  ExEnv::out0() << indent << "Non-electrostatic energy: "
       << setw(12) << setfill(' ')
       << (ecavitation_+solvent_->disp()+solvent_->rep())
          *unit->from_atomic_units() << " kcal/mol" << endl; 
  ExEnv::out0() << "_______________________________________________" << endl;

}

double
BEMSolventH::e()
{
  return escalar_;
}

/////////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
