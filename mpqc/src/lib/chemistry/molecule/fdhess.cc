//
// fdhess.cc
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

#include <stdlib.h>
#include <sys/stat.h>

#include <util/class/scexception.h>
#include <util/misc/formio.h>
#include <util/misc/timer.h>
#include <util/state/stateio.h>
#include <util/state/state_bin.h>
#include <util/group/mstate.h>
#include <util/keyval/keyval.h>
#include <math/scmat/blocked.h>
#include <math/symmetry/corrtab.h>
#include <chemistry/molecule/fdhess.h>

using namespace std;
using namespace sc;

#define DEFAULT_CHECKPOINT  1
#define DEFAULT_RESTART     1

/////////////////////////////////////////////////////////////////
// FinDispMolecularHessian

static ClassDesc FinDispMolecularHessian_cd(
  typeid(FinDispMolecularHessian),"FinDispMolecularHessian",1,"public MolecularHessian",
  0, create<FinDispMolecularHessian>, create<FinDispMolecularHessian>);

FinDispMolecularHessian::FinDispMolecularHessian(const Ref<MolecularEnergy> &e):
  mole_(e)
{
  only_totally_symmetric_ = 0;
  eliminate_cubic_terms_ = 1;
  do_null_displacement_ = 1;
  disp_ = 1.0e-2;
  ndisp_ = 0;
  debug_ = 0;
  gradients_ = 0;
  accuracy_ = disp_/1000;
  restart_ = DEFAULT_RESTART;
  checkpoint_ = DEFAULT_CHECKPOINT;
  checkpoint_file_ = 0;
  restart_file_ = 0;
  checkpoint_file_ = SCFormIO::fileext_to_filename(".ckpt.hess");
  restart_file_ = SCFormIO::fileext_to_filename(".ckpt.hess");
}

FinDispMolecularHessian::FinDispMolecularHessian(const Ref<KeyVal>&keyval):
  MolecularHessian(keyval)
{
  mole_ << keyval->describedclassvalue("energy");

  debug_ = keyval->booleanvalue("debug");

  displacement_point_group_ << keyval->describedclassvalue("point_group");

  disp_ = keyval->doublevalue("displacement",KeyValValuedouble(1.0e-2));

  KeyValValueboolean def_restart(DEFAULT_RESTART);
  KeyValValueboolean def_checkpoint(DEFAULT_CHECKPOINT);
  KeyValValueboolean truevalue(1);
  KeyValValueboolean falsevalue(0);

  restart_ = keyval->booleanvalue("restart", def_restart);
  char *def_file = SCFormIO::fileext_to_filename(".ckpt.hess");
  KeyValValueString def_restart_file(def_file,KeyValValueString::Steal);
  restart_file_ = keyval->pcharvalue("restart_file", def_restart_file);
  checkpoint_ = keyval->booleanvalue("checkpoint", def_checkpoint);
  checkpoint_file_ = keyval->pcharvalue("checkpoint_file", def_restart_file);
  only_totally_symmetric_ = keyval->booleanvalue("only_totally_symmetric",
                                                 falsevalue);
  eliminate_cubic_terms_ = keyval->booleanvalue("eliminate_cubic_terms",
                                                truevalue);

  do_null_displacement_ = keyval->booleanvalue("do_null_displacement",
                                               truevalue);

  accuracy_ = keyval->doublevalue("gradient_accuracy",
                                  KeyValValuedouble(disp_/1000));

  gradients_ = 0;
  ndisp_ = 0;
}

FinDispMolecularHessian::FinDispMolecularHessian(StateIn&s):
  SavableState(s),
  MolecularHessian(s)
{
  mole_ << SavableState::restore_state(s);
  s.get(checkpoint_);
  s.get(debug_);
  s.get(accuracy_);
  s.getstring(checkpoint_file_);
  s.getstring(restart_file_);

  gradients_ = 0;
  restore_displacements(s);
}

FinDispMolecularHessian::~FinDispMolecularHessian()
{
  delete[] gradients_;
  delete[] checkpoint_file_;
  delete[] restart_file_;
}

void
FinDispMolecularHessian::save_data_state(StateOut&s)
{
  MolecularHessian::save_data_state(s);

  SavableState::save_state(mole_.pointer(),s);
  s.put(checkpoint_);
  s.put(debug_);
  s.put(accuracy_);
  s.putstring(checkpoint_file_);
  s.putstring(restart_file_);

  checkpoint_displacements(s);
}

void
FinDispMolecularHessian::init()
{
  if (mole_.null()) return;

  mol_ = mole_->molecule();

  if (displacement_point_group_.null()) {
    displacement_point_group_
      = new PointGroup(*mol_->point_group().pointer());
    }

  nirrep_ = displacement_point_group_->char_table().nirrep();

  original_point_group_ = mol_->point_group();
  original_geometry_ = matrixkit()->vector(d3natom());

  int i, coor;
  for (i=0, coor=0; i<mol_->natom(); i++) {
      for (int j=0; j<3; j++, coor++) {
          original_geometry_(coor) = mol_->r(i,j);
        }
    }

  ndisp_ = 0;
  symbasis_ = cartesian_to_symmetry(mol_,
                                      displacement_point_group_,
                                      matrixkit());
  delete[] gradients_;
  gradients_ = new RefSCVector[ndisplace()];
}

void
FinDispMolecularHessian::set_energy(const Ref<MolecularEnergy> &energy)
{
  mole_ = energy;
}

MolecularEnergy*
FinDispMolecularHessian::energy() const
{
  return mole_.pointer();
}

void
FinDispMolecularHessian::restart()
{
  int statresult, statsize;
  Ref<MessageGrp> grp = MessageGrp::get_default_messagegrp();
  if (grp->me() == 0) {
    struct stat sb;
    statresult = stat(restart_file_,&sb);
    statsize = (statresult==0) ? sb.st_size : 0;
    }
  grp->bcast(statsize);
  if (statsize) {
    BcastStateInBin si(grp,restart_file_);
    restore_displacements(si);
    mol_ = mole_->molecule();

    if (ndisplacements_done() >= ndisplace()) {
        restart_=0;
        return;
      }
    }

  if (ndisp_) {
    int irrep, index;
    double coef;
    get_disp(ndisplacements_done(), irrep, index, coef);
    if (irrep != 0 && index != 0) {
        displace(ndisplacements_done());
        mole_->symmetry_changed();
      }
    }
  else {
    init();
    }
  restart_ = 0;
}

void
FinDispMolecularHessian::restore_displacements(StateIn& s)
{
  int i;
  displacement_point_group_ << SavableState::restore_state(s);
  original_point_group_ << SavableState::restore_state(s);
  original_geometry_ = matrixkit()->vector(d3natom());
  original_geometry_.restore(s);

  s.get(disp_);
  s.get(ndisp_);
  s.get(nirrep_);
  s.get(only_totally_symmetric_);
  s.get(eliminate_cubic_terms_);
  s.get(do_null_displacement_);

  if (ndisp_) {
    RefSCDimension symrow, symcol;
    symrow << SavableState::restore_state(s);
    symcol << SavableState::restore_state(s);
    Ref<SCMatrixKit> symkit = new BlockedSCMatrixKit(matrixkit());
    symbasis_ = symkit->matrix(symrow,symcol);
    symbasis_.restore(s);

    delete[] gradients_;
    gradients_ = new RefSCVector[ndisplace()];
    for (i=0; i < ndisp_; i++) {
      int ndisp;
      s.get(ndisp);
      RefSCDimension ddisp = new SCDimension(ndisp);
      gradients_[i] = matrixkit()->vector(ddisp);
      gradients_[i].restore(s);
      }
    }
}

void
FinDispMolecularHessian::checkpoint_displacements(StateOut& s)
{
  int i;
  SavableState::save_state(displacement_point_group_.pointer(),s);
  SavableState::save_state(original_point_group_.pointer(),s);
  original_geometry_.save(s);

  s.put(disp_);
  s.put(ndisp_);
  s.put(nirrep_);
  s.put(only_totally_symmetric_);
  s.put(eliminate_cubic_terms_);
  s.put(do_null_displacement_);

  if (ndisp_) {
    SavableState::save_state(symbasis_.rowdim().pointer(),s);
    SavableState::save_state(symbasis_.coldim().pointer(),s);
    symbasis_.save(s);

    for (i=0; i < ndisp_; i++) {
      s.put(gradients_[i].n());
      gradients_[i].save(s);
      }
    }
}

RefSCMatrix
FinDispMolecularHessian::displacements(int irrep) const
{
  BlockedSCMatrix *bsymbasis = dynamic_cast<BlockedSCMatrix*>(symbasis_.pointer());
  RefSCMatrix block = bsymbasis->block(irrep);
  if (block.null() || (only_totally_symmetric_ && irrep > 0)) {
    RefSCDimension zero = new SCDimension(0);
    block = matrixkit()->matrix(zero,zero);
    return block;
    }
  return block.t();
}

void
FinDispMolecularHessian::get_disp(int disp, int &irrep,
                                  int &index, double &coef)
{
  int disp_offset = 0;

  if (do_null_displacement_ && disp == 0) {
    irrep = 0;
    coef = 0.0;
    index = -1;
    return;
    }
  disp_offset++;
  // check for +ve totally symmetric displacements
  if (disp < disp_offset + displacements(0).ncol()) {
    irrep = 0;
    coef = 1.0;
    index = disp - disp_offset;
    return;
    }
  disp_offset += displacements(0).ncol();
  // check for -ve totally symmetric displacements
  if (eliminate_cubic_terms_) {
    if (disp < disp_offset + displacements(0).ncol()) {
      irrep = 0;
      coef = -1.0;
      index = disp - disp_offset;
      return;
      }
    disp_offset += displacements(0).ncol();
    }
  for (int i=1; i<nirrep_; i++) {
    if (disp < disp_offset + displacements(i).ncol()) {
      irrep = i;
      coef = 1.0;
      index = disp - disp_offset;
      return;
      }
    disp_offset += displacements(i).ncol();
    }
  throw ProgrammingError("bad displacement number",
                         __FILE__, __LINE__, class_desc());
}

int
FinDispMolecularHessian::ndisplace() const
{
  int ndisp = displacements(0).ncol();
  if (eliminate_cubic_terms_) {
    ndisp *= 2;
    }
  for (int i=1; i<nirrep_; i++) {
    ndisp += displacements(i).ncol();
    }
  if (do_null_displacement_) ndisp++;
  return ndisp;
}

void
FinDispMolecularHessian::displace(int disp)
{
  int irrep, index;
  double coef;
  get_disp(disp, irrep, index, coef);

  if (mole_.nonnull()) mole_->obsolete();

  for (int i=0, coor=0; i<mol_->natom(); i++) {
    for (int j=0; j<3; j++, coor++) {
      if (index >= 0) {
        mol_->r(i,j) = original_geometry_(coor)
                       + coef * disp_
                       * displacements(irrep)->get_element(coor,index);

        }
      else {
        mol_->r(i,j) = original_geometry_(coor);
        }
      }
    }

  if (irrep == 0) {
    mol_->set_point_group(original_point_group_);
    }
  else {
    Ref<PointGroup> oldpg = mol_->point_group();
    Ref<PointGroup> newpg = mol_->highest_point_group();
    CorrelationTable corrtab;
    if (corrtab.initialize_table(original_point_group_, newpg)) {
      // something went wrong so use c1 symmetry
      newpg = new PointGroup("c1");
      }
    if (!oldpg->equiv(newpg)) {
      mol_->set_point_group(newpg);
      mole_->symmetry_changed();
      }
    }

#ifdef DEBUG
  ExEnv::out0() << indent
       << "Displacement point group: " << endl
       << incindent << displacement_point_group_ << decindent;
  ExEnv::out0() << indent
       << "Displaced molecule: " << endl
       << incindent << mol_ << decindent;
#endif

  ExEnv::out0() << indent
       << "Displacement is "
       << displacement_point_group_->char_table().gamma(irrep).symbol()
       << " in " << displacement_point_group_->symbol()
       << ".  Using point group "
       << mol_->point_group()->symbol()
       << " for displaced molecule."
       << endl;
}

void
FinDispMolecularHessian::original_geometry()
{
  if (mole_.nonnull()) mole_->obsolete();

  for (int i=0, coor=0; i<mol_->natom(); i++) {
    for (int j=0; j<3; j++, coor++) {
      mol_->r(i,j) = original_geometry_(coor);
      }
    }

  if (!mol_->point_group()->equiv(original_point_group_)) {
    mol_->set_point_group(original_point_group_);
    mole_->symmetry_changed();
    }
}

void
FinDispMolecularHessian::set_gradient(int disp, const RefSCVector &grad)
{
  int irrep, index;
  double coef;
  get_disp(disp, irrep, index, coef);

  // transform the gradient into symmetrized coordinates
  gradients_[disp] = displacements(irrep).t() * grad;
  if (debug_) {
    grad.print("cartesian gradient");
    gradients_[disp].print("internal gradient");
    }

  ndisp_++;
}

RefSymmSCMatrix
FinDispMolecularHessian::compute_hessian_from_gradients()
{
  int i;

  RefSymmSCMatrix dhessian;

  RefSymmSCMatrix xhessian = matrixkit()->symmmatrix(d3natom());
  xhessian.assign(0.0);

  // start with the totally symmetric displacments
  int offset = 0;
  if (do_null_displacement_) offset++;
  RefSCMatrix dtrans = displacements(0);
  RefSCDimension ddim = dtrans.coldim();
  dhessian = matrixkit()->symmmatrix(ddim);
  for (i=0; i<ddim.n(); i++) {
    for (int j=0; j<=i; j++) {
      double hij = gradients_[i+offset](j) + gradients_[j+offset](i);
      double ncontrib = 2.0;
      if (do_null_displacement_) {
        hij -= gradients_[0](j) + gradients_[0](i);
        }
      if (eliminate_cubic_terms_) {
        hij -=   gradients_[i+ddim.n()+offset](j)
                 + gradients_[j+ddim.n()+offset](i);
        ncontrib += 2.0;
        if (do_null_displacement_) {
          hij += gradients_[0](j) + gradients_[0](i);
          }
        }
      hij /= ncontrib*disp_;
      dhessian(i,j) = hij;
      }
    }
  do_hess_for_irrep(0, dhessian, xhessian);

  offset += ddim.n();
  if (eliminate_cubic_terms_) offset += ddim.n();
  for (int irrep=1; irrep<nirrep_; irrep++) {
    dtrans = displacements(irrep);
    ddim = dtrans.coldim();
    if (ddim.n() == 0) continue;
    dhessian = matrixkit()->symmmatrix(ddim);
    for (i=0; i<ddim.n(); i++) {
      for (int j=0; j<=i; j++) {
        dhessian(i,j) = (gradients_[i+offset](j)
                         + gradients_[j+offset](i))
                        /(2.0*disp_);
        }
      }
    do_hess_for_irrep(irrep, dhessian, xhessian);
    offset += ddim.n();
    }

  if (debug_) {
    xhessian.print("xhessian");
    }

  return xhessian;
}

void
FinDispMolecularHessian::do_hess_for_irrep(int irrep,
                                        const RefSymmSCMatrix &dhessian,
                                        const RefSymmSCMatrix &xhessian)
{
  RefSCMatrix dtrans = displacements(irrep);
  RefSCDimension ddim = dtrans.coldim();
  if (ddim.n() == 0) return;
  if (debug_) {
    dhessian.print("dhessian");
    dtrans.print("dtrans");
    }
  xhessian.accumulate_transform(dtrans, dhessian);
}

RefSymmSCMatrix
FinDispMolecularHessian::cartesian_hessian()
{
  tim_enter("hessian");

  if (restart_) restart();
  else init();

  ExEnv::out0() << indent
       << "Computing molecular hessian from "
       << ndisplace() << " displacements:" << endl
       << indent << "Starting at displacement: "
       << ndisplacements_done() << endl;
  ExEnv::out0() << indent << "Hessian options: " << endl;
  ExEnv::out0() << indent << "  displacement: " << disp_
               << " bohr" << endl;
  ExEnv::out0() << indent << "  gradient_accuracy: "
               << accuracy_ << " au" << endl;
  ExEnv::out0() << indent << "  eliminate_cubic_terms: "
               << (eliminate_cubic_terms_==0?"no":"yes") << endl;
  ExEnv::out0() << indent << "  only_totally_symmetric: "
               << (only_totally_symmetric_==0?"no":"yes") << endl;

  for (int i=ndisplacements_done(); i<ndisplace(); i++) {
    // This produces side-effects in mol and may even change
    // its symmetry.
    ExEnv::out0() << endl << indent
         << "Beginning displacement " << i << ":" << endl;
    displace(i);

    mole_->obsolete();
    double original_accuracy;
    original_accuracy = mole_->desired_gradient_accuracy();
    if (accuracy_ > 0.0)
      mole_->set_desired_gradient_accuracy(accuracy_);
    else
      mole_->set_desired_gradient_accuracy(disp_/1000.0);
    RefSCVector gradv = mole_->get_cartesian_gradient();
    mole_->set_desired_gradient_accuracy(original_accuracy);
    set_gradient(i, gradv);

    if (checkpoint_) {
      const char *hessckptfile;
      if (MessageGrp::get_default_messagegrp()->me() == 0) {
        hessckptfile = checkpoint_file_;
        }
      else {
        hessckptfile = "/dev/null";
        }
      StateOutBin so(hessckptfile);
      checkpoint_displacements(so);
      }
    }
  original_geometry();
  RefSymmSCMatrix xhessian = compute_hessian_from_gradients();
  tim_exit("hessian");

  symbasis_ = 0;
  delete[] gradients_;
  gradients_ = 0;
  ndisp_ = 0;

  return xhessian;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
