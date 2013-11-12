//
// nbwfn.h
//
// Copyright (C) 2012 Edward Valeev
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
//#pragma interface // LLVM, Intel also define __GNUG__
#endif

#ifndef _mpqc_src_lib_chemistry_qc_nbody_nbwfn_h
#define _mpqc_src_lib_chemistry_qc_nbody_nbwfn_h

#include <chemistry/qc/nbody/ref.h>

namespace sc {

  /// @addtogroup ChemistryElectronicStructureNBody
  /// @{

  /**
   * ManyBodyWavefunction is a Wavefunction obtained from
   * a reference OneBodyWavefunction (its orbitals or more).
   */
  class ManyBodyWavefunction : public Wavefunction {
    public:
      /** A KeyVal constructor is used to generate a ManyBodyWavefunction
          object from a KeyVal object. This constructor accepts all keywords
          of the KeyVal constructor of the Wavefunction class, plus the additional
          keywords listed below.

          <table border="1">

          <tr><td><b>%Keyword</b><td><b>Type</b><td><b>Default</b><td><b>Description</b>

          <tr><td><tt>reference</tt><td>RefWavefunction<td>none<td>
	  specifies the reference wavefunction.
	  The most common choice is an SD_RefWavefunction.

          <tr><td><tt>world</tt><td>WavefunctionWorld<td>see the notes<td>
	  the WavefunctionWorld object that this Wavefunction belongs to.
	  If not given, this object will live in its own WavefunctionWorld.
          Ordinarily one does need to specify this.

          </table>
       */
      ManyBodyWavefunction(const Ref<KeyVal>& kv);
      ManyBodyWavefunction(StateIn&);
      virtual ~ManyBodyWavefunction();
      void save_data_state(StateOut&);

      const Ref<WavefunctionWorld>& world() const { return world_; }
      const Ref<RefWavefunction>& refwfn() const { return refwfn_; }

      double ref_energy() { return refwfn_->energy(); }
      double corr_energy() {
        return energy() - ref_energy();
      }

      void print(std::ostream& o=ExEnv::out0()) const;

      /// overloads MolecularEnergy::purge()
      void purge();
      void obsolete();

      void symmetry_changed();
      void set_desired_value_accuracy(double acc);

      /// specifies the ratio of the desired accuracy of RefWavefunction to the
      /// desired accuracy of this object
      static double ref_to_corr_acc() {
        return 0.01; // reference need to compute 100 times more precisely than this object. totally empirical.
      }

    private:
      static ClassDesc class_desc_;

      Ref<WavefunctionWorld> world_;
      Ref<RefWavefunction> refwfn_;
  };

  /// @}
  // end of addtogroup ChemistryElectronicStructureNBody

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
