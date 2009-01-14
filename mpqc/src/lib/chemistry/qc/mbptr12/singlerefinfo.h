//
// singlerefinfo.h
//
// Copyright (C) 2005 Edward Valeev
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

#ifndef _chemistry_qc_mbptr12_singlerefinfo_h
#define _chemistry_qc_mbptr12_singlerefinfo_h

#include <chemistry/qc/mbptr12/moindexspace.h>
#include <chemistry/qc/mbptr12/spin.h>

namespace sc {

  class SCF;

  /**
     SingleRefInfo maintains orbital information for the single-reference case.
  */
  class SingleRefInfo : virtual public SavableState {
    public:
    SingleRefInfo(StateIn&);
    /// Construct using SCF reference object
    SingleRefInfo(const Ref<SCF>& ref, unsigned int nfzc, unsigned int nfzv, bool delayed_initialization = false);
    ~SingleRefInfo();

    void save_data_state(StateOut&);
    void initialize();

    /// Returns the reference
    const Ref<SCF>& ref() const;
    /// Returns true if alpha and beta densities are not equal. Thus only false for CLSCF references.
    bool spin_polarized() const;
    /// Return number of frozen occupied orbitals
    unsigned int nfzc() const;
    /// Return number of frozen unoccupied orbitals
    unsigned int nfzv() const;

    /// Returns the space of symmetry-blocked orthogonal SOs (spans the entire space of the basis)
    const Ref<MOIndexSpace>& oso_space() const;

    //@{ \group Spin-independent spaces
    /**  Spin-independent variants can throw when used with spin-polarized reference.
    */
    /// Return the space of symmetry-blocked MOs
    const Ref<MOIndexSpace>& orbs_sb() const;
    /// Return the space of energy-sorted MOs
    const Ref<MOIndexSpace>& orbs() const;
    /// Return the space of symmetry-blocked doubly-occupied MOs
    const Ref<MOIndexSpace>& docc_sb() const;
    /// Return the space of doubly-occupied MOs
    const Ref<MOIndexSpace>& docc() const;
    /// Return the space of active doubly-occupied MOs
    const Ref<MOIndexSpace>& docc_act() const;
    /// Return the space of singly-occupied MOs
    const Ref<MOIndexSpace>& socc() const;
    /// Return the space of symmetry-blocked unoccupied (virtual) MOs
    const Ref<MOIndexSpace>& uocc_sb() const;
    /// Return the space of unoccupied (virtual) MOs
    const Ref<MOIndexSpace>& uocc() const;
    /// Return the space of active unoccupied (virtual) MOs
    const Ref<MOIndexSpace>& uocc_act() const;
    //@}

    /// Return the space of symmetry-blocked MOs of the given spin
    const Ref<MOIndexSpace>& orbs_sb(SpinCase1 spin) const;
    /// Return the space of energy-sorted MOs of the given spin
    const Ref<MOIndexSpace>& orbs(SpinCase1 spin) const;
    /// Return the space of symmery-blocked occupied MOs of the given spin
    const Ref<MOIndexSpace>& occ_sb(SpinCase1 spin) const;
    /// Return the space of occupied MOs of the given spin
    const Ref<MOIndexSpace>& occ(SpinCase1 spin) const;
    /// Return the space of occupied MOs of the given spin
    const Ref<MOIndexSpace>& occ_act(SpinCase1 spin) const;
    /// Return the space of symmetry-blocked unoccupied (virtual) MOs of the given spin
    const Ref<MOIndexSpace>& uocc_sb(SpinCase1 spin) const;
    /// Return the space of unoccupied (virtual) MOs of the given spin
    const Ref<MOIndexSpace>& uocc(SpinCase1 spin) const;
    /// Return the space of unoccupied (virtual) MOs of the given spin
    const Ref<MOIndexSpace>& uocc_act(SpinCase1 spin) const;

    private:
    /// initialized?
    bool initialized_;
    /// The reference function
    Ref<SCF> ref_;
    /// Number of occupied orbitals not used in correlated method
    unsigned int nfzc_;
    /// Number of unoccupied orbitals not used in correlated method
    unsigned int nfzv_;

    /// Following data structure is defined for each spin case
    typedef struct {
      Ref<MOIndexSpace> orbs_sb_;
      Ref<MOIndexSpace> orbs_;
      Ref<MOIndexSpace> occ_sb_;
      Ref<MOIndexSpace> occ_;
      Ref<MOIndexSpace> occ_act_;
      Ref<MOIndexSpace> uocc_sb_;
      Ref<MOIndexSpace> uocc_;
      Ref<MOIndexSpace> uocc_act_;
      /// "constructor"
      void init(SpinCase1 spin, const Ref<GaussianBasisSet>& bs,
                const Ref<Integral>& integral,
                const RefDiagSCMatrix& evals, const RefSCMatrix& evecs,
                const std::vector<double>& occs, unsigned int nfzc, unsigned int nfzv);
    } SpinSpaces;

    //@{
    /** see corresponding public member function
    */
    Ref<MOIndexSpace> orbs_sb_;
    Ref<MOIndexSpace> orbs_;
    Ref<MOIndexSpace> docc_sb_;
    Ref<MOIndexSpace> docc_;
    Ref<MOIndexSpace> docc_act_;
    Ref<MOIndexSpace> socc_;
    Ref<MOIndexSpace> uocc_sb_;
    Ref<MOIndexSpace> uocc_;
    Ref<MOIndexSpace> uocc_act_;
    SpinSpaces spinspaces_[NSpinCases1];
    //}@

    /// initialize spin-specific spaces
    void init_spinspecific_spaces();
    /// initialize spin-independent spaces
    void init_spinindependent_spaces();

    /// throws if ref_ is spin-polarized
    void throw_if_spin_polarized() const;

  };

};

#endif

