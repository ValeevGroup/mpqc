//
// cart.cc
//
// Copyright (C) 2011 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#include <util/class/class.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <chemistry/qc/basis/cart.h>
#include <chemistry/qc/basis/transform.h>

using namespace sc;

CartesianBasisSet::CartesianBasisSet(const Ref<KeyVal> & kv)
{
  integral_ << kv->describedclassvalue("integral").pointer();
  if (integral_ == 0)
    integral_ = Integral::get_default_integral();
  parent_ << kv->describedclassvalue("basis").pointer();
  if (parent_ == 0)
    throw InputError("keyword basis is required", __FILE__, __LINE__, "basis", "null", this->class_desc());

  convert(parent_, integral_);
}

CartesianBasisSet::CartesianBasisSet(StateIn & si) :
    GaussianBasisSet(si)
{
  integral_ << SavableState::restore_state(si);
  parent_ << SavableState::restore_state(si);
}

CartesianBasisSet::~CartesianBasisSet()
{
}

CartesianBasisSet::CartesianBasisSet(const Ref<GaussianBasisSet> & parent,
                                     Ref<Integral> integral) :
                                     parent_(parent),
                                     integral_(integral)
{
  convert(parent_, integral_);
}

void
CartesianBasisSet::save_data_state(StateOut & so)
{
  GaussianBasisSet::save_data_state(so);
  SavableState::save_state(integral_.pointer(), so);
  SavableState::save_state(parent_.pointer(), so);
}

const Ref<GaussianBasisSet> &
CartesianBasisSet::parent() const
{
  return parent_;
}

void
CartesianBasisSet::convert(const Ref<GaussianBasisSet> & parent,
                           const Ref<Integral>& integral)
{
  molecule_ = parent->molecule();

  const int nshell = parent->nshell();

  if (parent->has_pure() == false) {
    static_cast<GaussianBasisSet&>(*this) = *parent;
  }
  else { // something to do
    // create shells
    std::vector<Shell> shells;

    for (int s = 0; s < nshell; ++s) {
      const GaussianShell &shell = parent->shell(s);
      std::vector<bool> puream(shell.ncontraction(), false);
      shells.push_back(Shell(this, parent->shell_to_center(s), GaussianShell(shell.am(), puream,
                                     shell.exponents(), shell.coefficient_unnorm_block(),
                                     GaussianShell::Unnormalized)
                            )
                       );
    }

    std::ostringstream oss;
    oss << "Cartesian(" << parent->name() << ")";
    std::string name = oss.str();
    init(name, name, parent->molecule(), shells);
  }

#if 0
  // compute tform
  RefSCMatrix tform = this->matrixkit()->matrix(this->basisdim(), parent->basisdim());
  tform.assign(0.0);
  for(int s=0; s<nshell; ++s) {
    int bf_this = this->shell_to_function(s);
    int bf_parent = parent->shell_to_function(s);
    const int ncon = this->shell(s).ncontraction();
    const GaussianShell& pshell = parent->shell(s);
    for(int c=0; c<ncon; ++c) {
      if (pshell.is_pure(c)) {
        const int nbf_parent = pshell.nfunction(c);
        const int nbf_this = this->shell(s).nfunction(c);
        const SphericalTransform* st = integral->spherical_transform(pshell.am(c));
        const int ncoefs = st->n();
        for(int i=0; i<ncoefs; ++i) {
          const int cindex = st->cartindex(i);
          const int pindex = st->pureindex(i);
          const double tcoef = st->coef(i);
          tform.set_element(bf_this + cindex, bf_parent + pindex, tcoef);
        }
        bf_this += nbf_this;
        bf_parent += nbf_parent;
      }
      else {
        const int nbf = pshell.nfunction(c);
        for(int bf=0; bf<nbf; ++bf, ++bf_this, ++bf_parent)
          tform.set_element(bf_this, bf_parent, 1.0);
      }
    }
  }

  return tform;
#endif
}






/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
