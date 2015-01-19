//
// moinfo.h
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

#ifndef _mpqc_src_lib_extern_moinfo_moinfo_h
#define _mpqc_src_lib_extern_moinfo_moinfo_h

#include <vector>
#include <string>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/wfn/rdm.h>
#include <math/scmat/abstract.h>

namespace sc {

  /// Reads MO information from a text file
  /// Note that the MO ordering in the external file may not be the same as in MPQC
  /// For example, irreducible representations may be ordered differently in different programs
  /// Thus MOs will be reordered to be consistent with MPQC rules, and a map from the native
  /// to MPQC representation will be provided so that other files produced by the external program
  /// can be interpreted
  class ExternMOInfo : public DescribedClass {
    public:
      ExternMOInfo(std::string filename,
                   Ref<Integral> integral = Integral::get_default_integral()->clone(),
                   std::string basisname = std::string());
      ~ExternMOInfo() {}

      typedef OrderedOrbitalSpace<SymmetryMOOrder> SymmOrbitalSpace;
      typedef OrderedOrbitalSpace<EnergyMOOrder<std::less<double> > > CorrOrbitalSpace;
      /// This object reports orbitals in symmetry-blocked order
      const Ref<SymmOrbitalSpace>& orbs() const { return orbs_sb_; }
      /** maps the MO indices assumed by the contents of the data file
          to that of orbs() */
      const std::vector<unsigned int>& indexmap() const;
      /** same as @c occindexmap_sb(), except it maps the occupied subset of orbs only */
      const std::vector<unsigned int>& occindexmap_occ() const;
      /** same as @c actindexmap_sb(), except it maps the active subset of orbs only */
      const std::vector<unsigned int>& actindexmap_occ() const;

      const std::vector<unsigned int>& fzcpi() const;
      const std::vector<unsigned int>& fzvpi() const;
      const std::vector<unsigned int>& inactpi() const;
      const std::vector<unsigned int>& actpi() const;
      const std::vector<unsigned int>& corrpi() const;
      const std::vector<unsigned int>& mopi() const;

    private:
      static ClassDesc class_desc_;

      std::vector<unsigned int> indexmap_; //< all file order -> all mpqc order
      std::vector<unsigned int> occindexmap_occ_; //< occ file order -> occ mpqc order
      std::vector<unsigned int> actindexmap_occ_; //< act file order -> act mpqc order

      Ref<SymmOrbitalSpace> orbs_sb_;
      Ref<CorrOrbitalSpace> orbs_;
      std::vector<unsigned int> mopi_;
      std::vector<unsigned int> fzcpi_;
      std::vector<unsigned int> fzvpi_;
      std::vector<unsigned int> inactpi_;
      std::vector<unsigned int> actpi_;
      std::vector<unsigned int> corrpi_; // explicitly correlated orbitals
  };

  /// Reads 1-RDM from a text file
  class ExternSpinFreeRDMOne : public SpinFreeRDM<One>{
    public:
      /// reads 1-rdm from filename
      /// assumes that 1-rdm is expressed in orbs_sb

      /**
       * reads 1-rdm from filename, assumes that 1-rdm is expressed in orbs_sb,
       * indices in file are mapped to orbs_sb via indexmap
       *
       * @param filename file that contains an ASCII text specification of 1-RDM.
       *                 the file contains 3 columns: row index, column index, value.
       * @param indexmap maps the indices assumed in the file to orbs_sb
       * @param orbs_sb
       */
      ExternSpinFreeRDMOne(const std::string& filename,
                           const std::vector<unsigned int>& indexmap,
                           const Ref<OrbitalSpace>& orbs);
      /// receives 1-rdm as a constructor argument
      /// assumes that 1-rdm is expressed in orbs_sb
      ExternSpinFreeRDMOne(const RefSymmSCMatrix& rdm,
                           const Ref<OrbitalSpace>& orbs);
      virtual ~ExternSpinFreeRDMOne();

      /// cannot be obsoleted
      void obsolete() {}
      /// already computed
      void compute() {}
      /// the orbital space in which the density is reported
      Ref<OrbitalSpace> orbs() const { return orbs_; }
      /// density matrix
      RefSymmSCMatrix scmat() const { return scmat_; }

    private:
      static ClassDesc class_desc_;

      Ref<OrbitalSpace> orbs_;
  };

  /// Reads 2-RDM from a text file
  class ExternSpinFreeRDMTwo : public SpinFreeRDM<Two>{
    public:
      /// this assumes that the file reports 2-rdm in the entire occupied space.
      /// act_occ_indexmap maps the indices assumed in the file to occupied MPQC orbs
      ExternSpinFreeRDMTwo(const std::string& filename,
                           const std::vector<unsigned int>& indexmap,
                           const Ref<OrbitalSpace>& occ_orbs);
      virtual ~ExternSpinFreeRDMTwo();

      /// cannot be obsoleted
      void obsolete() {}
      /// already computed
      void compute() {}
      /// the orbital space of spincase s in which the density is reported
      Ref<OrbitalSpace> orbs() const { return orbs_; }
      /// density matrix
      RefSymmSCMatrix scmat() const { return scmat_; }
      /// density matrix in DistArray4 format
      const Ref<DistArray4>& da4() const;
      Ref< SpinFreeRDM<One> > rdm_m_1() const;

    private:
      static ClassDesc class_desc_;

    private:
      Ref<OrbitalSpace> orbs_;
      std::string filename_; // filename from which this was constructed -- may be useful to find rdm1 file

      void init_from_rdm2_occspace(const std::vector<unsigned int>& indexmap,
                                   const Ref<OrbitalSpace> & occ_orbs);
      void init_from_rdm2_actspace(const std::vector<unsigned int>& indexmap,
                                   const Ref<OrbitalSpace> & occ_orbs);

      mutable Ref<ExternSpinFreeRDMOne> rdm1_; // cached value to avoid recomputation
  };

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
