//
// findisp.cc
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

#include <stdlib.h>
#include <cassert>
#include <sys/stat.h>

#include <util/misc/scexception.h>
#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <util/state/stateio.h>
#include <util/state/state_bin.h>
#include <util/state/state_text.h>
#include <util/group/mstate.h>
#include <util/keyval/keyval.h>
#include <math/scmat/blocked.h>
#include <math/symmetry/corrtab.h>
#include <chemistry/molecule/findisp.h>

using namespace std;
using namespace sc;

#define DEFAULT_CHECKPOINT  0
#define DEFAULT_RESTART     1

/////////////////////////////////////////////////////////////////
// FinDispMolecularHessian::Impl

EGH::EGH() :
    energy_(0.0), gradient_(0), hessian_(0)
{
}

EGH::EGH(double e, const RefSCVector& g, const RefSymmSCMatrix& h) :
    energy_(e), gradient_(g), hessian_(h)
{
}

EGH::~EGH() {}

template<>
void
sc::FromStateIn<EGH>(EGH& v, StateIn& s, int& count) {
  double e;
  RefSCVector g;
  RefSymmSCMatrix h;

  s.get(e);
  FromStateIn(g,s,count);
  FromStateIn(h,s,count);

  v = EGH(e,g,h);
}

template<>
void
sc::ToStateOut<EGH>(const EGH& v, StateOut& s, int& count) {
  s.put(v.energy());
  ToStateOut<RefSCVector>(v.gradient(),s,count);
  ToStateOut<RefSymmSCMatrix>(v.hessian(),s,count);
}

/////////////////////////////////////////////////////////////////
// FinDispMolecularHessian::Params

ClassDesc FinDispMolecularHessian::Params::class_desc_(
  typeid(FinDispMolecularHessian::Params),"FinDispMolecularHessian::Params",2,"virtual public SavableState",
  create<FinDispMolecularHessian::Params>, create<FinDispMolecularHessian::Params>, create<FinDispMolecularHessian::Params>);

FinDispMolecularHessian::Params::Params()
{
  disp_pg_ = 0;
  disp_ = 1.0e-2;
  use_energies_ = false;
  only_totally_symmetric_ = false;
  // default for eliminate_quadratic_terms will be overridden by FinDispMolecularHessian
  // unless user provided eliminate_quadratic_terms
  user_provided_eliminate_quadratic_terms_ = false;
  eliminate_quadratic_terms_ = true;
  do_null_displacement_ = true;
  double desired_accuracy = 1e-4;
  energy_accuracy_ = desired_accuracy*disp_*disp_;
  gradient_accuracy_ = desired_accuracy*disp_;
  checkpoint_ = DEFAULT_CHECKPOINT;
  checkpoint_file_ = SCFormIO::fileext_to_filename_string(".ckpt.hess");
  restart_ = DEFAULT_RESTART;
  restart_file_ = SCFormIO::fileext_to_filename_string(".ckpt.hess");
  debug_ = 0;
}

FinDispMolecularHessian::Params::Params(const Ref<KeyVal>& keyval)
{
  use_energies_ = keyval->booleanvalue("use_energies", KeyValValueboolean(false));
  debug_ = keyval->intvalue("debug", KeyValValueint(0));
  disp_pg_ << keyval->describedclassvalue("point_group");
  disp_ = keyval->doublevalue("displacement",KeyValValuedouble(1.0e-2));
  only_totally_symmetric_ = keyval->booleanvalue("only_totally_symmetric",
                                                 KeyValValueboolean(false));
  user_provided_eliminate_quadratic_terms_ = keyval->exists("eliminate_quadratic_terms");
  eliminate_quadratic_terms_ = keyval->booleanvalue("eliminate_quadratic_terms",
                                                KeyValValueboolean(true));
  do_null_displacement_ = keyval->booleanvalue("do_null_displacement",
                                               KeyValValueboolean(true));
  const double desired_accuracy = keyval->doublevalue("accuracy", KeyValValuedouble(1e-4));
  energy_accuracy_ = keyval->doublevalue("energy_accuracy",
                                           KeyValValuedouble(desired_accuracy*disp_*disp_));
  gradient_accuracy_ = keyval->doublevalue("gradient_accuracy",
                                           KeyValValuedouble(desired_accuracy*disp_));

  KeyValValueboolean def_checkpoint(DEFAULT_CHECKPOINT);
  checkpoint_ = keyval->booleanvalue("checkpoint", def_checkpoint);
  char *def_file = SCFormIO::fileext_to_filename(".ckpt.hess");
  KeyValValuestring def_restart_file(def_file);
  checkpoint_file_ = keyval->stringvalue("checkpoint_file", def_restart_file);

  KeyValValueboolean def_restart(DEFAULT_RESTART);
  restart_ = keyval->booleanvalue("restart", def_restart);
  restart_file_ = keyval->stringvalue("restart_file", def_restart_file);
}

FinDispMolecularHessian::Params::Params(StateIn& s)
{
  s.get(use_energies_);
  s.get(debug_);
  disp_pg_ << SavableState::restore_state(s);
  s.get(disp_);
  s.get(only_totally_symmetric_);
  s.get(user_provided_eliminate_quadratic_terms_);
  s.get(eliminate_quadratic_terms_);
  s.get(do_null_displacement_);
  s.get(energy_accuracy_);
  s.get(gradient_accuracy_);
  s.get(checkpoint_);
  s.get(checkpoint_file_);
  s.get(restart_);
  s.get(restart_file_);
}

FinDispMolecularHessian::Params::~Params()
{
}

void
FinDispMolecularHessian::Params::save_data_state(StateOut& s)
{
  s.put(use_energies_);
  s.put(debug_);
  SavableState::save_state(disp_pg_.pointer(), s);
  s.put(disp_);
  s.put(only_totally_symmetric_);
  s.put(user_provided_eliminate_quadratic_terms_);
  s.put(eliminate_quadratic_terms_);
  s.put(do_null_displacement_);
  s.put(energy_accuracy_);
  s.put(gradient_accuracy_);
  s.put(checkpoint_);
  s.put(checkpoint_file_);
  s.put(restart_);
  s.put(restart_file_);
}

void
FinDispMolecularHessian::Params::set_desired_accuracy(double acc) {
  energy_accuracy_ = acc*disp_*disp_;
  gradient_accuracy_ = acc*disp_;
}

void
FinDispMolecularHessian::Params::set_disp_size(double d) {
  const double base_acc = gradient_accuracy_ / disp_; // extract the accuracy of the hessian (it's displacement-size independent)
  disp_ = d; // update displacement size
  set_desired_accuracy(base_acc); // recompute desired accuracy of gradients and energies
}

/////////////////////////////////////////////////////////////////
// FinDispMolecularHessian::Impl

ClassDesc FinDispMolecularHessian::Impl::class_desc_(
  typeid(FinDispMolecularHessian::Impl),"FinDispMolecularHessian::Impl",1,"virtual public SavableState",
  0, 0, 0);

FinDispMolecularHessian::Impl::Impl(const Ref<MolecularEnergy>& e,
                                    const Ref<Params>& params) :
  params_(params),
  mole_(e)
{
  init();
}

FinDispMolecularHessian::Impl::Impl(StateIn& s):
  SavableState(s)
{
  restore_displacements(s);
}

void
FinDispMolecularHessian::Impl::save_data_state(StateOut& s)
{
  checkpoint_displacements(s);
}

FinDispMolecularHessian::Impl::~Impl()
{
}

void
FinDispMolecularHessian::Impl::init()
{
  if (mole_.null()) return;

  Ref<Molecule> mol = mole_->molecule();

  if (params_->disp_pg().null()) {
    params_->set_disp_pg(new PointGroup(*mol->point_group().pointer()));
  }

  original_point_group_ = mol->point_group();
  original_geometry_ = matrixkit()->vector(d3natom());

  int i, coor;
  for (i=0, coor=0; i<mol->natom(); i++) {
    for (int j=0; j<3; j++, coor++) {
      original_geometry_(coor) = mol->r(i,j);
    }
  }

  symbasis_ = cartesian_to_symmetry(mol,
                                    params_->disp_pg(),
                                    matrixkit());
}

void
FinDispMolecularHessian::Impl::restart() {
#if 0
  // update geometry to the current displacement; if none had been done re-init
  if (ndisplace()) {
    int irrep;
    std::vector<std::pair<int,double> > disp;
    get_disp(ndisplacements_done(), irrep, disp);
    if (disp.empty() == false && irrep != 0) {
      displace(ndisplacements_done());
      mole_->symmetry_changed();
    }
  }
  else {
    init();
  }
#endif
}

void
FinDispMolecularHessian::Impl::restore_displacements(StateIn& s)
{
  int count;

  params_ = 0;
  params_ << SavableState::restore_state(s);
  original_point_group_ << SavableState::restore_state(s);
  FromStateIn(original_geometry_,s,count);
  FromStateIn(symbasis_,s,count);
  values_.restore(s);
}

void
FinDispMolecularHessian::Impl::checkpoint_displacements(StateOut& s)
{
  int count;

  SavableState::save_state(params_.pointer(),s);
  SavableState::save_state(original_point_group_.pointer(),s);
  ToStateOut<RefSCVector>(original_geometry_,s,count);
  ToStateOut<RefSCMatrix>(symbasis_,s,count);
  values_.save(s);
}

RefSCMatrix
FinDispMolecularHessian::Impl::displacements(int irrep) const
{
  BlockedSCMatrix *bsymbasis = dynamic_cast<BlockedSCMatrix*>(symbasis_.pointer());
  RefSCMatrix block = bsymbasis->block(irrep);
  if (block.null() || (params_->only_totally_symmetric() && irrep > 0)) {
    RefSCDimension zero = new SCDimension(0);
    block = matrixkit()->matrix(zero,zero);
    return block;
    }
  return block.t();
}

void
FinDispMolecularHessian::Impl::displace(const Displacement& disp)
{

  Ref<Molecule> mol = mole_->molecule();
  const unsigned int natom = mol->natom();

  for (int i=0, coor=0; i<natom; i++) {
    for (int j=0; j<3; j++, coor++) {
      mol->r(i,j) = original_geometry_(coor);
    }
  }
  if (disp.empty()) {
    if (mole_) mole_->obsolete();
    return;
  }

  // either one or two coordinates are displaced, in the latter case both must have same symmetry
  const int irrep = this->coor_to_irrep(disp[0].first);

  const size_t nsalc = disp.size();
  for(int s=0; s<nsalc; ++s) { // for every SALC that is part of this displacement
    const double dispsize = disp[s].second;
    const int salc = disp[s].first;
    for (int i=0, coor=0; i<natom; i++) {
      for (int j=0; j<3; j++, coor++) {
        mol->r(i,j) += dispsize * params_->disp_size()
                       * symbasis_.get_element(salc,coor);

      }
    }
  }

  if (irrep == 0) {
    mol->set_point_group(original_point_group_);
  }
  else {
    Ref<PointGroup> oldpg = mol->point_group();
    Ref<PointGroup> newpg = mol->highest_point_group();
    CorrelationTable corrtab;
    if (corrtab.initialize_table(original_point_group_, newpg)) {
      // something went wrong so use c1 symmetry
      newpg = new PointGroup("c1");
    }
    if (!oldpg->equiv(newpg)) {
      mol->set_point_group(newpg);
      mole_->symmetry_changed();
    }
  }

#if 0
  ExEnv::out0() << indent
       << "Displacement point group: " << endl
       << incindent;
  params_->disp_pg()->print();
  ExEnv::out0() << decindent;
  ExEnv::out0() << indent
       << "Displaced molecule: " << endl
       << incindent;
  mol->print();
  ExEnv::out0() << decindent;
#endif

  ExEnv::out0() << indent
       << "Displacement is "
       << params_->disp_pg()->char_table().gamma(irrep).symbol()
       << " in " << params_->disp_pg()->symbol()
       << ".  Using point group "
       << mol->point_group()->symbol()
       << " for displaced molecule."
       << endl;
  if (mole_) mole_->obsolete();
  mol->print();
}

void
FinDispMolecularHessian::Impl::original_geometry()
{
  Ref<Molecule> mol = mole_->molecule();
  const int natom =mol->natom();
  for (int i=0, coor=0; i<natom; i++) {
    for (int j=0; j<3; j++, coor++) {
      mol->r(i,j) = original_geometry_(coor);
    }
  }

  if (!mol->point_group()->equiv(original_point_group_)) {
    mol->set_point_group(original_point_group_);
    mole_->symmetry_changed();
  }

  if (mole_) mole_->obsolete();
}

#if 0
void
FinDispMolecularHessian::Impl::get_disp(size_t disp_index, int &irrep,
                                        std::vector< std::pair<int, double> >& disp)
{
  if (disp_index >= disp.size())
    throw ProgrammingError("bad displacement number",
                           __FILE__, __LINE__, class_desc());

  disp.resize(0);
  typedef Displacements<unsigned int>::citer citer;
  citer i = disps_.find(disp_index);
  disp.resize(i->first.size());
  std::copy(i->first.begin(), i->first.end(), disp.begin());
}
#endif

void
FinDispMolecularHessian::Impl::do_hess_for_irrep(int irrep,
                                                 const RefSymmSCMatrix &dhessian,
                                                 const RefSymmSCMatrix &xhessian)
{
  RefSCMatrix dtrans = displacements(irrep);
  RefSCDimension ddim = dtrans.coldim();
  if (ddim.n() == 0) return;
  if (params_->debug()) {
    std::ostringstream oss; oss << "dhessian for irrep " << irrep;
    dhessian.print(oss.str().c_str());
    dtrans.print("dtrans");
    }
  xhessian.accumulate_transform(dtrans, dhessian);
}

RefSymmSCMatrix
FinDispMolecularHessian::Impl::cartesian_hessian()
{
  Timer tim("hessian");
  RefSymmSCMatrix xhessian = compute_hessian();
  tim.exit("hessian");

  symbasis_ = 0;

  return xhessian;
}

unsigned int
FinDispMolecularHessian::Impl::coor_to_irrep(unsigned int symm_coord) const {
  // row dimension of symbasis_ is blocked according to irreps -> map row to block # -> voila!
  RefSCDimension dsym = symbasis_.rowdim();
  MPQC_ASSERT(symm_coord < dsym.n());
  Ref<SCBlockInfo> bsym = dsym->blocks();
  MPQC_ASSERT(bsym);
  int block, offset;
  bsym->elem_to_block(symm_coord, block, offset);
  return block;
}

/////////////////////////////////////////////////////////////////
// FinDispMolecularHessian::GradientsImpl

ClassDesc FinDispMolecularHessian::GradientsImpl::class_desc_(
  typeid(FinDispMolecularHessian::GradientsImpl),"FinDispMolecularHessian::GradientsImpl",1,"public FinDispMolecularHessian::Impl",
  0, 0, create<FinDispMolecularHessian::GradientsImpl>);

FinDispMolecularHessian::GradientsImpl::GradientsImpl(const Ref<MolecularEnergy>& e,
                                                      const Ref<Params>& params) :
    Impl(e, params)
{
  validate_mole(mole_);
}

FinDispMolecularHessian::GradientsImpl::GradientsImpl(StateIn&s):
  SavableState(s), Impl(s)
{
}

FinDispMolecularHessian::GradientsImpl::~GradientsImpl()
{
}

void
FinDispMolecularHessian::GradientsImpl::save_data_state(StateOut&s)
{
  Impl::save_data_state(s);
}

void
FinDispMolecularHessian::GradientsImpl::validate_mole(const Ref<MolecularEnergy>& e)
{
  if (e.null()) return;
  if (e->gradient_implemented() == 0)
    throw ProgrammingError("FinDispMolecularHessian -- hessian from gradients requested but MolecularEnergy cannot compute gradients", __FILE__, __LINE__);
}

void
FinDispMolecularHessian::GradientsImpl::compute_mole(const Displacement& disp) {

  // check if disp has been computed
  if (values_.has(disp)) {
    if (values_.find(disp).second.gradient())
      return;
  }

  // This produces side-effects in mol and may even change
  // its symmetry.
  ExEnv::out0() << endl << indent
       << "Beginning displacement " << values_.size()+1 << ":" << endl;
  displace(disp);

  // mole_->obsolete(); displace obsoleted mole
  const double original_accuracy = mole_->desired_gradient_accuracy();
  mole_->set_desired_gradient_accuracy(params_->gradient_accuracy());
  RefSCVector gradv = mole_->get_cartesian_gradient();
  const double energy = mole_->energy();
  mole_->set_desired_gradient_accuracy(original_accuracy);
  set_gradient(disp, energy, gradv);

  if (params_->checkpoint()) {
    const char *hessckptfile;
    if (MessageGrp::get_default_messagegrp()->me() == 0) {
      hessckptfile = params_->checkpoint_file().c_str();
    }
    else {
      hessckptfile = "/dev/null";
    }
    StateOutBin so(hessckptfile);
    checkpoint_displacements(so);
  }
}

void
FinDispMolecularHessian::GradientsImpl::set_gradient(const Displacement& disp, double energy, const RefSCVector &grad)
{
  // what's the irrep of the displacement? Only single-coordinate displacements are needed to compute hessian from gradients
  int irrep = 0;
  if (!disp.empty()) {
    MPQC_ASSERT(disp.size() == 1);
    const int int_coor = disp[0].first;
    irrep = coor_to_irrep(int_coor);
  }
  // transform the gradient into symmetrized coordinates
  RefSCVector symm_grad = displacements(irrep).t() * grad;
  values_.push(disp,EGH(energy, symm_grad, 0));
  if (params_->debug()) {
    grad.print("cartesian gradient");
    symm_grad.print("internal gradient");
  }
}

int
FinDispMolecularHessian::GradientsImpl::ndisplace() const {
  int result = 0;
  // start with the totally symmetric displacements
  {
    RefSCMatrix dtrans = displacements(0);
    const int n = dtrans.ncol();
    result += n * (params_->eliminate_quadratic_terms() ? 2 : 1);
    result += (params_->do_null_displacement() && !params_->eliminate_quadratic_terms() )? 1 : 0;
  }
  // process all nonsymmetric displacements
  if (!params_->only_totally_symmetric()) {
    for (int irrep=1; irrep<params_->nirrep(); irrep++) {
      RefSCMatrix dtrans = displacements(irrep);
      const int n = dtrans.ncol();
      result += n;
    }
  }

  return result;
}

RefSymmSCMatrix
FinDispMolecularHessian::GradientsImpl::compute_hessian()
{
  RefSymmSCMatrix dhessian;  // Hessian in the internal coordinates
  RefSymmSCMatrix xhessian = matrixkit()->symmmatrix(d3natom()); // Hessian in the Cartesian coordinates
  xhessian.assign(0.0);

  // need gradient at reference geometry?
  RefSCVector grad_0;
  if (params_->do_null_displacement() && !params_->eliminate_quadratic_terms()) {
    Displacement empty_disp;
    this->compute_mole(empty_disp);
    grad_0 = values_.find(empty_disp).second.gradient();
  }

  // start with the totally symmetric displacements
  RefSCMatrix dtrans = displacements(0);
  RefSCDimension ddim = dtrans.coldim();
  dhessian = matrixkit()->symmmatrix(ddim);
  for (int i=0; i<ddim.n(); i++) {

    RefSCVector grad_i_p, grad_i_m;
    // gradient at +delta along i
    {
      Displacement disp; disp.push_back(make_pair(i,1.0));
      this->compute_mole(disp);
      grad_i_p = values_.find(disp).second.gradient();
    }
    // gradient at -delta along j
    if (params_->eliminate_quadratic_terms()) {
      Displacement disp; disp.push_back(make_pair(i,-1.0));
      this->compute_mole(disp);
      grad_i_m = values_.find(disp).second.gradient();
    }

    for (int j=0; j<=i; j++) {

      RefSCVector grad_j_p, grad_j_m;
      // gradient at +delta along j
      {
        Displacement disp; disp.push_back(make_pair(j,1.0));
        this->compute_mole(disp);
        grad_j_p = values_.find(disp).second.gradient();
      }
      // gradient at -delta along j
      if (params_->eliminate_quadratic_terms()) {
        Displacement disp; disp.push_back(make_pair(j,-1.0));
        this->compute_mole(disp);
        grad_j_m = values_.find(disp).second.gradient();
      }

      double hij = grad_i_p(j) + grad_j_p(i);
      double ncontrib = 2.0;

      // gradient at reference geometry
      if (params_->do_null_displacement() && !params_->eliminate_quadratic_terms()) {
        hij -= grad_0(j) + grad_0(i);
      }

      // gradients at -delta along i and j
      if (params_->eliminate_quadratic_terms()) {
        hij -= grad_i_m(j) + grad_j_m(i);
        ncontrib += 2.0;
      }

      hij /= ncontrib * params_->disp_size();
      dhessian(i,j) = hij;
    }
  }
  do_hess_for_irrep(0, dhessian, xhessian);

  // non-totally symmetric blocks
  if (!params_->only_totally_symmetric()) {
    int coor_offset = ddim.n();
    for (int irrep=1; irrep<params_->nirrep(); irrep++) {
      dtrans = displacements(irrep);
      ddim = dtrans.coldim();
      if (ddim.n() == 0) continue;
      dhessian = matrixkit()->symmmatrix(ddim);
      for (int i=0; i<ddim.n(); i++) {

        RefSCVector grad_i_p;
        // gradient at +delta along i
        {
          Displacement disp; disp.push_back(make_pair(i+coor_offset,1.0));
          this->compute_mole(disp);
          grad_i_p = values_.find(disp).second.gradient();
        }

        for (int j=0; j<=i; j++) {

          RefSCVector grad_j_p;
          // gradient at +delta along j
          {
            Displacement disp; disp.push_back(make_pair(j+coor_offset,1.0));
            this->compute_mole(disp);
            grad_j_p = values_.find(disp).second.gradient();
          }

          dhessian(i,j) = (grad_i_p(j) + grad_j_p(i)) /(2.0 * params_->disp_size());
        }
      }
      do_hess_for_irrep(irrep, dhessian, xhessian);
      coor_offset += ddim.n();
    }
  }

  if (params_->debug()) {
    xhessian.print("xhessian");
  }

  original_geometry();

  return xhessian;
}

/////////////////////////////////////////////////////////////////
// FinDispMolecularHessian::EnergiesImpl

ClassDesc FinDispMolecularHessian::EnergiesImpl::class_desc_(
  typeid(FinDispMolecularHessian::EnergiesImpl),"FinDispMolecularHessian::EnergiesImpl",1,"public FinDispMolecularHessian::Impl",
  0, 0, create<FinDispMolecularHessian::EnergiesImpl>);


FinDispMolecularHessian::EnergiesImpl::EnergiesImpl(const Ref<MolecularEnergy>& e,
                                                    const Ref<Params>& params) :
  Impl(e, params)
{
}

FinDispMolecularHessian::EnergiesImpl::EnergiesImpl(StateIn&s):
  SavableState(s), Impl(s)
{
}

FinDispMolecularHessian::EnergiesImpl::~EnergiesImpl()
{
}

void
FinDispMolecularHessian::EnergiesImpl::save_data_state(StateOut&s)
{
  Impl::save_data_state(s);
}

void
FinDispMolecularHessian::EnergiesImpl::validate_mole(const Ref<MolecularEnergy>& e)
{
  if (e.null()) return;
  if (e->value_implemented() == false)
    throw ProgrammingError("FinDispMolecularHessian -- hessian from energies requested but MolecularEnergy cannot compute energies", __FILE__, __LINE__);
}

int
FinDispMolecularHessian::EnergiesImpl::ndisplace() const {
  int result = 1;  // always need the energy at reference geometry
  // start with the totally symmetric displacements
  {
    RefSCMatrix dtrans = displacements(0);
    const int n = dtrans.ncol();
    // diagonal force constants
    result += (params_->eliminate_quadratic_terms() ? 4 : 2) * n;
    // off-diagonal force constants
    result += (params_->eliminate_quadratic_terms() ? 4 : 2) * n * (n-1) / 2;
  }
  // process all nonsymmetric displacements
  if (!params_->only_totally_symmetric()) {
    for (int irrep=1; irrep<params_->nirrep(); irrep++) {
      RefSCMatrix dtrans = displacements(irrep);
      const int n = dtrans.ncol();
      // diagonal force constants
      result += (params_->eliminate_quadratic_terms() ? 2 : 1) * n;
      // off-diagonal force constants
      result += (params_->eliminate_quadratic_terms() ? 2 : 1) * n * (n-1) / 2;
    }
  }
  return result;
}

RefSymmSCMatrix
FinDispMolecularHessian::EnergiesImpl::compute_hessian()
{
  RefSymmSCMatrix dhessian;  // Hessian in the internal coordinates
  RefSymmSCMatrix xhessian = matrixkit()->symmmatrix(d3natom()); // Hessian in the Cartesian coordinates
  xhessian.assign(0.0);

  //
  // to compute diagonal force constants:
  // Fii = (E(+d_i) + E(-d_i) - 2E(0))/d^2, or extended (5-point) stencil
  // Fii = (- E(+2d_i) - E(-2d_i) + 16 E(+d_i) + 16 E(-d_i) - 30 E(0))/(12 d^2)
  // to compute off-diagonal force constants:
  // Fij = (- E(+d_i-d_j) - E(-d_i+d_j) + E(+d_i) + E(-d_i) + E(+d_j) + E(-d_j) - 2E(0))/ (2 d^2), or extended stencil
  // Fij = (- 30 f[0, 0]
  //        + 16 f[0, -d] + 16 f[0, d] + 16 f[d, 0] + 16 f[-d, 0]
  //        - f[0, -2 d] - f[0, 2 d] - f[-2 d, 0] - f[2 d, 0]
  //        - 16 f[-d, d] - 16 f[d, -d]
  //        + f[-2 d, 2 d] + f[2 d, -2 d] )/(24 d^2)
  //

  // energy at reference geometry
  double energy_0;
  {
    Displacement empty_disp;
    this->compute_mole(empty_disp);
    energy_0 = values_.find(empty_disp).second.energy();
  }

  const double disp_size2 = params_->disp_size() * params_->disp_size();

  // start with the totally symmetric displacements
  RefSCMatrix dtrans = displacements(0);
  RefSCDimension ddim = dtrans.coldim();
  dhessian = matrixkit()->symmmatrix(ddim);
  for (int i=0; i<ddim.n(); i++) {

    Eij e_ii(i, i, *this);

    const double energy_i_p = e_ii(1,0);
    const double energy_i_m = e_ii(-1,0);

    if (!params_->eliminate_quadratic_terms()) {
      dhessian(i,i) = (energy_i_p + energy_i_m - 2.0 * energy_0) / disp_size2;
    }
    else {
      const double energy_i_p2 = e_ii(2,0);
      const double energy_i_m2 = e_ii(-2,0);
      dhessian(i,i) = (- energy_i_p2 - energy_i_m2 + 16.0 * (energy_i_p + energy_i_m) - 30.0 * energy_0) / (12.0 * disp_size2);
    }

    for (int j=0; j<i; j++) {

      Eij e_ij(i, j, *this);

      if (!params_->eliminate_quadratic_terms()) {
        dhessian(i,j) = (- e_ij(+1,-1) - e_ij(-1,+1) + e_ij(+1,0) + e_ij(-1,0) + e_ij(0,+1) + e_ij(0,-1) - 2.0 * energy_0) / (2.0 * disp_size2);
      }
      else {
        dhessian(i,j) = (+ (e_ij(+2,-2) + e_ij(-2,+2))
                         - 16.0 * (e_ij(+1,-1) + e_ij(-1,+1))
                         -        (e_ij(+2,0) + e_ij(-2,0) + e_ij(0,+2) + e_ij(0,-2))
                         + 16.0 * (e_ij(+1,0) + e_ij(-1,0) + e_ij(0,+1) + e_ij(0,-1))
                         - 30.0 * energy_0
                        ) / (24.0 * disp_size2);
      }

    }
  }
  do_hess_for_irrep(0, dhessian, xhessian);

  // non-totally symmetric blocks
  if (!params_->only_totally_symmetric()) {
    int coor_offset = ddim.n();
    for (int irrep=1; irrep<params_->nirrep(); irrep++) {
      dtrans = displacements(irrep);
      ddim = dtrans.coldim();
      if (ddim.n() == 0) continue;
      dhessian = matrixkit()->symmmatrix(ddim);
      for (int i=0; i<ddim.n(); i++) {

        Eij e_ii(i + coor_offset, i + coor_offset, *this);

        const double energy_i_p = e_ii(1,0);
        const double energy_i_m = energy_i_p;

        if (!params_->eliminate_quadratic_terms()) {
          dhessian(i,i) = (energy_i_p + energy_i_m - 2.0 * energy_0) / disp_size2;
        }
        else {
          const double energy_i_p2 = e_ii(2,0);
          const double energy_i_m2 = energy_i_p2;
          dhessian(i,i) = (- energy_i_p2 - energy_i_m2 + 16.0 * (energy_i_p + energy_i_m) - 30.0 * energy_0) / (12.0 * disp_size2);
        }

        for (int j=0; j<i; j++) {

          Eij e_ij(i + coor_offset, j + coor_offset, *this);

          if (!params_->eliminate_quadratic_terms()) {
            dhessian(i,j) = (- e_ij(+1,-1) + e_ij(+1,0) + e_ij(0,+1) - energy_0) / (disp_size2);
          }
          else {
            dhessian(i,j) = (+ (e_ij(+2,-2))
                             - 16.0 * (e_ij(+1,-1))
                             -        (e_ij(+2,0) + e_ij(0,+2))
                             + 16.0 * (e_ij(+1,0) + e_ij(0,+1))
                             - 15.0 * energy_0
                            ) / (12.0 * disp_size2);
          }

        }
      }
      do_hess_for_irrep(irrep, dhessian, xhessian);
      coor_offset += ddim.n();
    }
  }

  if (params_->debug()) {
    xhessian.print("xhessian");
  }

  original_geometry();

  return xhessian;
}

void
FinDispMolecularHessian::EnergiesImpl::compute_mole(const Displacement& disp) {

  // check if disp has been computed
  if (values_.has(disp)) {
    if (values_.find(disp).second.energy() != 0.0)
      ExEnv::out0() << indent << "WARNING: FinDispMolecularHessian(energy) already encountered a dislacement, but energy is 0" << std::endl;
    return;
  }

  // This produces side-effects in mol and may even change
  // its symmetry.
  ExEnv::out0() << endl << indent
       << "Beginning displacement " << values_.size()+1 << ":" << endl;
  displace(disp);

  // mole_->obsolete(); displace obsoleted mole
  const double original_accuracy = mole_->desired_value_accuracy();
  mole_->set_desired_value_accuracy(params_->energy_accuracy());
  const double energy = mole_->energy();
  values_.push(disp,EGH(energy,0,0));
  mole_->set_desired_value_accuracy(original_accuracy);

  if (params_->checkpoint()) {
    const char *hessckptfile;
    if (MessageGrp::get_default_messagegrp()->me() == 0) {
      hessckptfile = params_->checkpoint_file().c_str();
    }
    else {
      hessckptfile = "/dev/null";
    }
    StateOutBin so(hessckptfile);
    checkpoint_displacements(so);
  }
}


/////////////////////////////////////////////////////////////////
// FinDispMolecularHessian

static ClassDesc FinDispMolecularHessian_cd(
  typeid(FinDispMolecularHessian),"FinDispMolecularHessian",1,"public MolecularHessian",
  0, create<FinDispMolecularHessian>, create<FinDispMolecularHessian>);

FinDispMolecularHessian::FinDispMolecularHessian(const Ref<MolecularEnergy> &e)
{
  params_ = new Params;
  //init_pimpl(e);
  mole_init_ = e;
  if (mole_init_) override_default_params();
}

FinDispMolecularHessian::FinDispMolecularHessian(const Ref<KeyVal>&keyval):
  MolecularHessian(keyval)
{
  Ref<MolecularEnergy> e; e << keyval->describedclassvalue("energy");
  params_ = new Params(keyval);
  //init_pimpl(e);
  mole_init_ = e;
  if (mole_init_) override_default_params();
}

FinDispMolecularHessian::FinDispMolecularHessian(StateIn&s):
  SavableState(s),
  MolecularHessian(s)
{
  pimpl_ << SavableState::restore_state(s);
  mole_init_ = 0;
}

FinDispMolecularHessian::~FinDispMolecularHessian()
{
  pimpl_ = 0;
  mole_init_ = 0;
}

void
FinDispMolecularHessian::save_data_state(StateOut&s)
{
  MolecularHessian::save_data_state(s);
  SavableState::save_state(pimpl_.pointer(),s);
}

void
FinDispMolecularHessian::restart()
{
  if (pimpl_.null()) { init_pimpl(mole_init_); mole_init_ = 0; }

  // broadcast contents of restart file
  int statresult, statsize;
  Ref<MessageGrp> grp = MessageGrp::get_default_messagegrp();
  if (grp->me() == 0) {
    struct stat sb;
    statresult = stat(params()->restart_file().c_str(),&sb);
    statsize = (statresult==0) ? sb.st_size : 0;
    }
  grp->bcast(statsize);
  if (statsize) {
    BcastStateInBin si(grp,params()->restart_file().c_str());
    pimpl_->restore_displacements(si);
  }
  else {
    pimpl_->init(); // no restart file? This is not a restart -- initialize pimpl as usual
  }

  params()->set_restart(false);
}

RefSymmSCMatrix
FinDispMolecularHessian::cartesian_hessian()
{
  if (pimpl_.null()) { init_pimpl(mole_init_); mole_init_ = 0; }

  if (params()->restart()) restart();
  else pimpl_->init();  // initialize original_geometry etc.

  Ref<EnergiesImpl> eimpl; eimpl << pimpl_;

  ExEnv::out0() << indent
       << "Computing molecular hessian by finite differences of "
       << (eimpl ? "energies" : "gradients") << " from "
       << pimpl_->ndisplace() << " displacements:" << endl;
  ExEnv::out0() << indent << "Hessian options: " << endl;
  ExEnv::out0() << indent << "  displacement: " << pimpl_->params()->disp_size()
               << " bohr" << endl;
  ExEnv::out0() << indent << "  eliminate_quadratic_terms: "
               << (pimpl_->params()->eliminate_quadratic_terms() ? "yes" : "no") << endl;
  ExEnv::out0() << indent << "  only_totally_symmetric: "
               << (pimpl_->params()->only_totally_symmetric() ? "yes" : "no") << endl;

  return pimpl_->cartesian_hessian();
}

void
FinDispMolecularHessian::set_energy(const Ref<MolecularEnergy>& mole) {
  mole_init_ = mole;
  pimpl_ = 0;
}

void
FinDispMolecularHessian::init_pimpl(const Ref<MolecularEnergy>& mole) {
  if (mole.null()) {
    if (pimpl_) pimpl_->set_mole(mole);
  }
  else {
    pimpl_ = 0;
    if (mole->gradient_implemented() && !params_->use_energies()) {
      pimpl_ = new GradientsImpl(mole, params_);
    }
    else {
      pimpl_ = new EnergiesImpl(mole, params_);
    }
    pimpl_->init();
  }
}

void
FinDispMolecularHessian::set_desired_accuracy(double acc) {
  MolecularHessian::set_desired_accuracy(acc);
  params_->set_desired_accuracy(acc);
}

void
FinDispMolecularHessian::override_default_params()
{
  MPQC_ASSERT(params_);
  MPQC_ASSERT(mole_init_);
  if (params_->user_provided_eliminate_quadratic_terms() == false) {
    if (mole_init_->gradient_implemented() && !params_->use_energies()) {
      // override the default for eliminate_quadratic_terms
      params_->set_eliminate_quadratic_terms(true);
    }
    else {
      // override the default for eliminate_quadratic_terms
      params_->set_eliminate_quadratic_terms(false);
    }
  }
}

/////////////////////////////////////////////////////////////////
// FinDispMolecularGradient

static ClassDesc FinDispMolecularGradient_cd(
  typeid(FinDispMolecularGradient),"FinDispMolecularGradient",1,"public MolecularGradient",
  0, create<FinDispMolecularGradient>, create<FinDispMolecularGradient>);

FinDispMolecularGradient::FinDispMolecularGradient(const Ref<MolecularEnergy> &e):
  mole_(e)
{
  eliminate_quadratic_terms_ = 0;
  disp_ = 1.0e-2;
  debug_ = 0;
  energy_accuracy_ = MolecularGradient::desired_accuracy() * disp_;
  restart_ = DEFAULT_RESTART;
  checkpoint_ = DEFAULT_CHECKPOINT;
  checkpoint_file_ = SCFormIO::fileext_to_filename_string(".ckpt.grad");
  restart_file_ = SCFormIO::fileext_to_filename_string(".ckpt.grad");
}

FinDispMolecularGradient::FinDispMolecularGradient(const Ref<KeyVal>&keyval):
  MolecularGradient(keyval)
{
  mole_ << keyval->describedclassvalue("energy");
  if (mole_.null()) {
    throw InputError("FinDispMolecularGradient KeyVal ctor: did not find a valid value for keyword energy",
                     __FILE__, __LINE__);
  }

  debug_ = keyval->intvalue("debug", KeyValValueint(0));

  // optional
  displacement_point_group_ << keyval->describedclassvalue("point_group");

  disp_ = keyval->doublevalue("displacement",KeyValValuedouble(1.0e-2));

  KeyValValueboolean def_restart(DEFAULT_RESTART);
  KeyValValueboolean def_checkpoint(DEFAULT_CHECKPOINT);
  KeyValValueboolean truevalue(1);
  KeyValValueboolean falsevalue(0);

  restart_ = keyval->booleanvalue("restart", def_restart);
  char *def_file = SCFormIO::fileext_to_filename(".ckpt.grad");
  KeyValValuestring def_restart_file(def_file);
  restart_file_ = keyval->stringvalue("restart_file", def_restart_file);
  checkpoint_ = keyval->booleanvalue("checkpoint", def_checkpoint);
  checkpoint_file_ = keyval->stringvalue("checkpoint_file", def_restart_file);
  eliminate_quadratic_terms_ = keyval->booleanvalue("eliminate_quadratic_terms",
                                                    falsevalue);

  energy_accuracy_ = keyval->doublevalue("energy_accuracy",
                                         KeyValValuedouble(MolecularGradient::desired_accuracy() * disp_));

}

FinDispMolecularGradient::FinDispMolecularGradient(StateIn&s):
  SavableState(s),
  MolecularGradient(s)
{
  mole_ << SavableState::restore_state(s);
  s.get(checkpoint_);
  s.get(debug_);
  s.get(energy_accuracy_);
  s.get(checkpoint_file_);
  s.get(restart_file_);

  restore_displacements(s);
}

FinDispMolecularGradient::~FinDispMolecularGradient()
{
}

void
FinDispMolecularGradient::save_data_state(StateOut&s)
{
  MolecularGradient::save_data_state(s);

  SavableState::save_state(mole_.pointer(),s);
  s.put(checkpoint_);
  s.put(debug_);
  s.put(energy_accuracy_);
  s.put(checkpoint_file_);
  s.put(restart_file_);

  checkpoint_displacements(s);
}

void
FinDispMolecularGradient::init()
{
  if (mole_.null()) return;

  mol_ = mole_->molecule();

  if (displacement_point_group_.null()) {
    displacement_point_group_
      = new PointGroup(*mol_->point_group().pointer());
    }

  const int nirrep = displacement_point_group_->char_table().nirrep();

  original_geometry_ = matrixkit()->vector(d3natom());

  int i, coor;
  for (i=0, coor=0; i<mol_->natom(); i++) {
      for (int j=0; j<3; j++, coor++) {
          original_geometry_(coor) = mol_->r(i,j);
        }
    }

  symbasis_ = MolecularHessian::cartesian_to_symmetry(mol_,
                                                      displacement_point_group_,
                                                      matrixkit());
  energies_.resize(0);
}

void
FinDispMolecularGradient::set_energy(const Ref<MolecularEnergy> &energy)
{
  mole_ = energy;
}

MolecularEnergy*
FinDispMolecularGradient::energy() const
{
  return mole_.pointer();
}

void
FinDispMolecularGradient::restart()
{
  // broadcast restart file from node 0
  int statresult, statsize;
  Ref<MessageGrp> grp = MessageGrp::get_default_messagegrp();
  if (grp->me() == 0) {
    struct stat sb;
    statresult = stat(restart_file_.c_str(),&sb);
    statsize = (statresult==0) ? sb.st_size : 0;
    }
  grp->bcast(statsize);
  if (statsize) {
    BcastStateInBin si(grp,restart_file_.c_str());
    restore_displacements(si);
    mol_ = mole_->molecule();

    if (ndisplacements_done() == ndisplace()) {
        restart_=0;
        return;
      }
    }

  if (ndisplacements_done()) {
    displace(ndisplacements_done());
  }
  else {
    init();
  }
  restart_ = 0;
}

void
FinDispMolecularGradient::restore_displacements(StateIn& s)
{
  int i;
  displacement_point_group_ << SavableState::restore_state(s);
  original_geometry_ = matrixkit()->vector(d3natom());
  original_geometry_.restore(s);

  s.get(disp_);
  s.get(eliminate_quadratic_terms_);
  s.get(energies_);

  if (energies_.size()) {
    RefSCDimension symrow, symcol;
    symrow << SavableState::restore_state(s);
    symcol << SavableState::restore_state(s);
    Ref<SCMatrixKit> symkit = new BlockedSCMatrixKit(matrixkit());
    symbasis_ = symkit->matrix(symrow,symcol);
    symbasis_.restore(s);
    }
}

void
FinDispMolecularGradient::checkpoint_displacements(StateOut& s)
{
  int i;
  SavableState::save_state(displacement_point_group_.pointer(),s);
  original_geometry_.save(s);

  s.put(disp_);
  s.put(eliminate_quadratic_terms_);
  s.put(energies_);

  if (energies_.size()) {
    SavableState::save_state(symbasis_.rowdim().pointer(),s);
    SavableState::save_state(symbasis_.coldim().pointer(),s);
    symbasis_.save(s);
    }
}

RefSCMatrix
FinDispMolecularGradient::displacements(int irrep) const
{
  BlockedSCMatrix *bsymbasis = dynamic_cast<BlockedSCMatrix*>(symbasis_.pointer());
  RefSCMatrix block = bsymbasis->block(irrep);
  if (block.null() || irrep > 0) {  // only totally symmetric displacements are needed
    RefSCDimension zero = new SCDimension(0);
    block = matrixkit()->matrix(zero,zero);
    return block;
    }
  return block.t();
}

void
FinDispMolecularGradient::get_disp(int disp, int &index, double &dispsize)
{
  int disp_offset = 0;
  const int ndisp_per_coord = eliminate_quadratic_terms_ ? 2 : 1; // number of displacements in each direction for each totally-symmetric coordinate
  const int ndisp_per_dir = ndisp_per_coord * displacements(0).ncol(); // total number of totally-symmetric displacements in each direction.

  // check for +ve totally symmetric displacements
  if (disp < disp_offset + ndisp_per_dir) {
    dispsize = (eliminate_quadratic_terms_ && disp%2 == 1) ? 2.0 : 1.0;  // for 4-pt formula odd displacements are + 2 delta, even are + delta
    index = (disp - disp_offset) / ndisp_per_coord;
    return;
    }
  disp_offset += ndisp_per_dir;
  if (disp < disp_offset + ndisp_per_dir) {
    dispsize = (eliminate_quadratic_terms_ && disp%2 == 1) ? -2.0 : -1.0;  // for 4-pt formula odd displacements are + 2 delta, even are + delta
    index = (disp - disp_offset) / ndisp_per_coord;
    return;
  }
  throw ProgrammingError("bad displacement number",
                         __FILE__, __LINE__, class_desc());
}

int
FinDispMolecularGradient::ndisplace() const
{
  const int ndisp = displacements(0).ncol() * (eliminate_quadratic_terms_ ? 4 : 2);
  return ndisp;
}

void
FinDispMolecularGradient::displace(int disp)
{
  const int irrep = 0;
  int index;
  double dispsize;
  get_disp(disp, index, dispsize);

  for (int i=0, coor=0; i<mol_->natom(); i++) {
    for (int j=0; j<3; j++, coor++) {
      if (index >= 0) {
        mol_->r(i,j) = original_geometry_(coor)
                       + dispsize * disp_
                       * displacements(irrep)->get_element(coor,index);

        }
      else {
        mol_->r(i,j) = original_geometry_(coor);
        }
      }
    }

  // symmetry does not change

#if 0
  ExEnv::out0() << indent
       << "Displacement point group: " << endl
       << incindent;
  displacement_point_group_->print();
  ExEnv::out0() << decindent;
  ExEnv::out0() << indent
       << "Displaced molecule: " << endl
       << incindent;
  mol_->print();
  ExEnv::out0() << decindent;
#endif

  ExEnv::out0() << indent
       << "Displacement is "
       << displacement_point_group_->char_table().gamma(irrep).symbol()
       << " in " << displacement_point_group_->symbol()
       << ".  Using point group "
       << mol_->point_group()->symbol()
       << " for displaced molecule."
       << endl;

  if (mole_) mole_->obsolete();
  mol_->print();
}

void
FinDispMolecularGradient::original_geometry()
{
  for (int i=0, coor=0; i<mol_->natom(); i++) {
    for (int j=0; j<3; j++, coor++) {
      mol_->r(i,j) = original_geometry_(coor);
      }
    }
  if (mole_) mole_->obsolete();
}

RefSCVector
FinDispMolecularGradient::compute_gradient()
{
  // start with the totally symmetric displacements
  RefSCMatrix dtrans = displacements(0);
  RefSCDimension ddim = dtrans.coldim();
  RefSCVector igradient = matrixkit()->vector(ddim);
  igradient.assign(0.0);

  for(int d=0; d<ndisplace(); ++d) {
    int coord; double dispsize;
    get_disp(d, coord, dispsize);
    //std::cout << "disp = " << d << "  coord = " << coord << "  dispsize = " << dispsize << "  energy = " << energies_[d] << std::endl;
    double coeff = 0.0;
    if (!eliminate_quadratic_terms_) {
      // 2-pt formula: f' = (f+ - f-)/(2.0 d)
      coeff = (dispsize == 1.0) ? 1.0 : -1.0;
    }
    else {
      // 4-pt formula: f' = (-f2+ + 8 f+ - 8 f- + f2-)/(12.0 d)
      if (dispsize == 1.0) coeff = 8.0;
      else if (dispsize == -1.0) coeff = -8.0;
      else if (dispsize == 2.0) coeff = -1.0;
      else if (dispsize == -2.0) coeff = 1.0;
    }
    igradient.accumulate_element(coord, coeff * energies_[d]);
  }
  if (eliminate_quadratic_terms_)
    igradient.scale( (1.0 / 12.0) / disp_);
  else
    igradient.scale(0.5 / disp_);
  RefSCVector cgradient = dtrans  * igradient;

  if (true || debug_) {
    igradient.print("gradient in internal coordinates");
    cgradient.print("gradient in Cartesian coordinates");
  }

  return cgradient;
}

RefSCVector
FinDispMolecularGradient::cartesian_gradient()
{
  Timer tim("gradient");

  if (restart_) restart();
  else init();

  ExEnv::out0() << indent
       << "Computing molecular gradient by finite differences of energies from "
       << ndisplace() << " displacements:" << endl
       << indent << "Starting at displacement: "
       << ndisplacements_done() << endl;
  ExEnv::out0() << indent << "Gradient options: " << endl;
  ExEnv::out0() << indent << "  displacement: " << disp_
               << " bohr" << endl;
  ExEnv::out0() << indent << "  energy_accuracy: "
               << energy_accuracy_ << " au" << endl;
  ExEnv::out0() << indent << "  eliminate_quadratic_terms: "
               << (eliminate_quadratic_terms_==0?"no":"yes") << endl;

  for (int i=ndisplacements_done(); i<ndisplace(); i++) {
    // This produces side-effects in mol
    ExEnv::out0() << endl << indent
         << "Beginning displacement " << i+1 << ":" << endl;
    displace(i);

    // mole_->obsolete(); displace() obsoleted mole
    double original_accuracy;
    original_accuracy = mole_->desired_value_accuracy();
    mole_->set_desired_value_accuracy(energy_accuracy_);
    const double energy = mole_->energy();
    mole_->set_desired_value_accuracy(original_accuracy);
    energies_.push_back(energy);

    if (checkpoint_) {
      const char *gradckptfile;
      if (MessageGrp::get_default_messagegrp()->me() == 0) {
        gradckptfile = checkpoint_file_.c_str();
        }
      else {
        gradckptfile = "/dev/null";
        }
      StateOutBin so(gradckptfile);
      checkpoint_displacements(so);
      }
    }
  original_geometry();
  RefSCVector gradient = compute_gradient();
  tim.exit("gradient");

  return gradient;
}

void
FinDispMolecularGradient::set_desired_accuracy(double acc) {
  MolecularGradient::set_desired_accuracy(acc);
  energy_accuracy_ = acc * disp_;
}

void
FinDispMolecularGradient::set_disp_size(double d) {
  const double base_acc = energy_accuracy_ / disp_; // extract the accuracy of the gradient (it's displacement-size independent)
  disp_ = d; // update displacement size
  set_desired_accuracy(base_acc); // recompute desired accuracy of energies
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
