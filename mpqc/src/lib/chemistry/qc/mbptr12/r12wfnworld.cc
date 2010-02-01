//
// r12wfnworld.cc
//
// Copyright (C) 2003 Edward Valeev
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <stdexcept>
#include <stdlib.h>
#include <util/misc/formio.h>
#include <util/ref/ref.h>
#include <util/misc/string.h>
#include <util/class/scexception.h>
#include <math/scmat/local.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <chemistry/qc/mbptr12/orbitalspace.h>
#include <chemistry/qc/mbptr12/transform_factory.h>
#include <chemistry/qc/mbptr12/registry.timpl.h>
#include <chemistry/qc/mbptr12/ref.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

/*---------------
  R12WavefunctionWorld
 ---------------*/
static ClassDesc R12WavefunctionWorld_cd(
  typeid(R12WavefunctionWorld),"R12WavefunctionWorld",11,"virtual public SavableState",
  0, 0, create<R12WavefunctionWorld>);

R12WavefunctionWorld::R12WavefunctionWorld(
    const Ref<KeyVal>& keyval,
    const Ref<RefWavefunction>& ref) :
    ref_(ref)
{
  // by default use spin-adapted algorithms for closed-shell
  spinadapted_ = keyval->booleanvalue("spinadapted",KeyValValueboolean(true));
  if (!ref->spin_polarized()) spinadapted_ = false;

  bs_aux_ = require_dynamic_cast<GaussianBasisSet*>(
      keyval->describedclassvalue("aux_basis").pointer(),
      "R12Technology::R12Technology\n"
      );
  if (bs_aux_.pointer() == NULL)
      bs_aux_ = ref->basis();

  r12tech_ = new R12Technology(keyval,ref->basis(),ref->uocc(Alpha)->basis(),bs_aux_);
  // Make sure can use the integral factory for R12 calcs
  r12tech_->check_integral_factory(integral());

  initialize();
}

R12WavefunctionWorld::R12WavefunctionWorld(StateIn& si) : SavableState(si)
{
  ref_ << SavableState::restore_state(si);
  bs_aux_ << SavableState::restore_state(si);
  bs_ri_ << SavableState::restore_state(si);

  si.get(spinadapted_);

  abs_space_ << SavableState::restore_state(si);
  ribs_space_ << SavableState::restore_state(si);

  initialize();
}

R12WavefunctionWorld::~R12WavefunctionWorld()
{
}

void R12WavefunctionWorld::save_data_state(StateOut& so)
{
  SavableState::save_state(ref_.pointer(),so);
  SavableState::save_state(bs_aux_.pointer(),so);
  SavableState::save_state(bs_ri_.pointer(),so);
  so.put(spinadapted_);
  SavableState::save_state(abs_space_.pointer(),so);
  SavableState::save_state(ribs_space_.pointer(),so);
}

void
R12WavefunctionWorld::initialize()
{
  // provide hints to the factory about the likely use of transforms
  {
    // 1) if stdapprox is A' or A'' most transforms will never be reused
    // 2) otherwise use persistent transforms
    if (r12tech()->stdapprox() == R12Technology::StdApprox_Ap ||
        r12tech()->stdapprox() == R12Technology::StdApprox_App)
      world()->tfactory()->hints().data_persistent(false);
    else
      world()->tfactory()->hints().data_persistent(true);
  }

  construct_ri_basis_(r12tech()->safety_check());

  obs_eq_ribs_ = basis()->equiv(basis_ri());
  obs_eq_vbs_ = basis()->equiv( ref()->uocc()->basis() );

  {
    // also create AO spaces
    Ref<OrbitalSpaceRegistry> idxreg = world()->tfactory()->orbital_registry();
    Ref<AOSpaceRegistry> aoidxreg = world()->tfactory()->ao_registry();
    Ref<Integral> localints = ref()->integral()->clone();
    if (!obs_eq_ribs()) { // RI-BS
      Ref<OrbitalSpace> mu = new AtomicOrbitalSpace("mu'", "RIBS(AO)", basis_ri(), localints);
      idxreg->add(make_keyspace_pair(mu));
      aoidxreg->add(mu->basis(),mu);
    }
  }
}

void
R12WavefunctionWorld::obsolete() {
  ref_->obsolete(); // this obsolete WavefunctionWorld
  abs_space_ = 0;
  ribs_space_ = 0;
  cabs_space_[Alpha] = 0;
  cabs_space_[Beta] = 0;
  this->initialize();
}

const Ref<OrbitalSpace>&
R12WavefunctionWorld::cabs_space(const SpinCase1& S) const
{
    if (r12tech()->abs_method() == R12Technology::ABS_CABS ||
        r12tech()->abs_method() == R12Technology::ABS_CABSPlus)
	  return cabs_space_[S];
    else
	  throw ProgrammingError("CABS space requested by abs_method set to ABS/ABS+",__FILE__,__LINE__);
}

bool
R12WavefunctionWorld::sdref() const {
  // only references based on OneBodyWavefunction are detected as single-determinant references!
  Ref<SD_RefWavefunction> sd; sd << ref();
  return sd.nonnull();
}

void
R12WavefunctionWorld::print(std::ostream& o) const {

  o << indent << "R12WavefunctionWorld:" << endl;
  o << incindent;

  if (!bs_aux_->equiv(basis())) {
      o << indent << "Auxiliary Basis Set (ABS):" << endl;
      o << incindent; bs_aux_->print(o); o << decindent << endl;
  }
  r12tech()->print(o);
  o << indent << "Spin-adapted algorithm: " << (spinadapted_ ? "true" : "false") << endl;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
