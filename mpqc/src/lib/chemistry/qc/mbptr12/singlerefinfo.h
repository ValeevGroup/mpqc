//
// singlerefinfo.h
//
// Copyright (C) 2005 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#include <chemistry/qc/mbptr12/moindexspace.h>

#ifndef _chemistry_qc_mbptr12_singlerefinfo_h
#define _chemistry_qc_mbptr12_singlerefinfo_h


namespace sc {
  
  class SCF;
  
  /**
     SingleRefInfo maintains orbital information for the single-reference case.
  */
  class SingleRefInfo : virtual public SavableState {
    public:
    /// Spin cases
    typedef enum {AlphaSpin=0, BetaSpin=1} SpinCase;
    
    SingleRefInfo(StateIn&);
    /// Construct using SCF reference object
    SingleRefInfo(const Ref<SCF>& ref);
    ~SingleRefInfo();
    
    /// Returns the reference
    const Ref<SCF>& ref() const;
    
    /// Returns the space of symmetry-blocked orthogonal SOs (spans the entire space of the basis)
    const Ref<MOIndexSpace>& oso_space() const;
    
    //@{ \group Spin-independent spaces
    /**  Spin-independent variants can throw when used with spin-polarized reference.
    */
    /// Return the space of symmetry-blocked MOs
    const Ref<MOIndexSpace>& symblk_mo() const;
    /// Return the space of energy-sorted MOs
    const Ref<MOIndexSpace>& energy_mo() const;
    /// Return the space of doubly-occupied MOs
    const Ref<MOIndexSpace>& docc() const;
    /// Return the space of singly-occupied MOs
    const Ref<MOIndexSpace>& socc() const;
    /// Return the space of unoccupied (virtual) MOs
    const Ref<MOIndexSpace>& uocc() const;
    //@}
    
    /// Return the space of symmetry-blocked MOs of the given spin
    const Ref<MOIndexSpace>& symblk_mo(SpinCase spin) const;
    /// Return the space of energy-sorted MOs of the given spin
    const Ref<MOIndexSpace>& energy_mo(SpinCase spin) const;
    /// Return the space of occupied MOs of the given spin
    const Ref<MOIndexSpace>& occ(SpinCase spin) const;
    /// Return the space of unoccupied (virtual) MOs of the given spin
    const Ref<MOIndexSpace>& uocc(SpinCase spin) const;
    
    private:
    /// The reference function
    Ref<SCF> ref_;
    
    /// Following data structure is defined for each spin case
    typedef struct {
      Ref<MOIndexSpace> symblk_mo_;
      Ref<MOIndexSpace> energy_mo_;
      Ref<MOIndexSpace> occ_;
      Ref<MOIndexSpace> uocc_;
      /// "constructor"
      void init(const std::string& prefix, const Ref<GaussianBasisSet>& bs,
                const RefDiagSCMatrix& evals, const RefSCMatrix& evecs,
                const std::vector<double>& occs);
    } SpinSpaces;
    
    //@{
    /** see corresponding public member function
    */
    Ref<MOIndexSpace> symblk_mo_;
    Ref<MOIndexSpace> energy_mo_;
    Ref<MOIndexSpace> docc_;
    Ref<MOIndexSpace> socc_;
    Ref<MOIndexSpace> uocc_;
    SpinSpaces spinspaces_[2];
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

