//
// r12wfnworld.h
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
#pragma interface
#endif

#ifndef _chemistry_qc_mbptr12_vxbevalinfo_h
#define _chemistry_qc_mbptr12_vxbevalinfo_h

#include <string>
#include <util/misc/string.h>
#include <util/ref/ref.h>
#include <math/scmat/abstract.h>
#include <util/group/memory.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbptr12/r12technology.h>
#include <chemistry/qc/mbptr12/linearr12.h>
#include <chemistry/qc/mbptr12/ansatz.h>
#include <chemistry/qc/mbptr12/orbitalspace.h>
#include <chemistry/qc/mbptr12/transform_factory.h>
#include <chemistry/qc/mbptr12/ref.h>
#include <chemistry/qc/mbptr12/moints_runtime.h>
#include <chemistry/qc/mbptr12/fockbuild_runtime.h>

namespace sc {

/** Class R12WavefunctionWorld describes the environment of a Wavefunction implementing an R12 method */
class R12WavefunctionWorld : virtual public SavableState {

  // change to 0 to use the old set of OrbitalSpace keys
  static const int USE_NEW_ORBITALSPACE_KEYS = 1;

public:

  /// Describes the method of storing transformed MO integrals.
  typedef MOIntsTransform::StoreMethod StoreMethod;

  R12WavefunctionWorld(StateIn&);
    /** KeyVal constructor uses keywords of R12Technology and the following keywords:

	<dt><tt>spinadapted</tt><dd> This boolean keyword specifies whether spin-adapted algorithms
	should be used, if possible.

        </dl>
    */
  R12WavefunctionWorld(const Ref<KeyVal>& keyval,
                       const Ref<R12RefWavefunction>& ref);
  ~R12WavefunctionWorld();

  void save_data_state(StateOut&);

  /// Return the R12RefWavefunction object
  const Ref<R12RefWavefunction>& ref() const { return ref_; }
  const Ref<WavefunctionWorld>& world() const { return ref()->world(); }
  Wavefunction* wfn() const { return world()->wfn(); }
  Ref<R12Technology> r12tech() const { return r12tech_; }
  /// Returns the orbital basis set (OBS) object
  const Ref<GaussianBasisSet>& basis() const { return ref()->basis(); }
  /// Returns the resolution-of-the-identity basis set (RIBS) object
  const Ref<GaussianBasisSet>& basis_ri() const { return bs_ri_; }
  /// Returns the virtuals basis set (VBS) object
  const Ref<GaussianBasisSet>& basis_vir() const { return ref()->uocc()->basis(); }
  /// Returns true if VBS is equivalent to OBS
  bool obs_eq_vbs() const { return obs_eq_vbs_; }
  /// Returns true if RIBS is equivalent to OBS
  bool obs_eq_ribs() const { return obs_eq_ribs_; }

  bool spinadapted() const { return spinadapted_; }

  OverlapOrthog::OrthogMethod orthog_method() const { return wfn()->orthog_method(); }
  double lindep_tol() const { return wfn()->lindep_tol(); }
  Ref<Integral> integral() const { return world()->integral(); }
  /// is this a single-determinant reference?
  bool sdref() const;

  /// Returns the OrbitalSpace object for ABS
  const Ref<OrbitalSpace>& abs_space() const { return abs_space_; }
  /// Returns the OrbitalSpace object for RI-BS: approximates the identity
  const Ref<OrbitalSpace>& ribs_space() const { return ribs_space_; }
  /// Returns subspace of ribs_space that is the complement to OBS. If abs_method = ABS/ABS+, this will throw
  const Ref<OrbitalSpace>& cabs_space(const SpinCase1& S) const;

  void print(std::ostream& o) const;

private:

  Ref<R12RefWavefunction> ref_;          //!< describes the reference wavefunction
  Ref<R12Technology> r12tech_;    //!< describes the R12 technology
  Ref<GaussianBasisSet> bs_aux_;  //!< the auxiliary basis used for computing the RI basis used in R12
  Ref<GaussianBasisSet> bs_ri_;   //!< the RI basis used in R12
  bool obs_eq_vbs_;
  bool obs_eq_ribs_;

  bool spinadapted_;
  int nlindep_aux_;
  int nlindep_ri_;

  Ref<OrbitalSpace> abs_space_;  // ABS space
  Ref<OrbitalSpace> ribs_space_; // RIBS basis
  Ref<OrbitalSpace> cabs_space_[NSpinCases1]; // CABS basis

  // construct the RI basis based on abs_method
  void construct_ri_basis_(bool safe);
  void construct_ri_basis_ks_(bool safe);
  void construct_ri_basis_ksplus_(bool safe);
  void construct_ri_basis_ev_(bool safe);
  void construct_ri_basis_evplus_(bool safe);
  // Uses ri_basis to construct a basis that spans the orthogonal complement to the OBS
  void construct_ortho_comp_svd_();
  // Returns true if ABS spans OBS
  bool abs_spans_obs_();
  // Construct orthog_aux_
  void construct_orthog_aux_();
  // Construct orthog_ri_
  void construct_orthog_ri_();
  void initialize();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


