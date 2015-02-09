//
// molcas_pt2r12.h
//
// Copyright (C) 2014 Chong Peng
//
// Authors: Chong Peng
// Maintainer: Chong Peng and Edward Valeev
//
// This file is part of the MPQC Toolkit.
//
// The MPQC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The MPQC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the MPQC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//


#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_molcas_pt2r12_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_molcas_pt2r12_h

#include <chemistry/qc/mbptr12/extern_pt2r12.h>

namespace sc{
  class MolcasPT2R12 : public MolecularEnergy {
    public:

      MolcasPT2R12(const Ref<KeyVal>& kv);

      void compute();

      /// initialize values needed for ExternPT2R12
      void initialize();

      // run molcas input file
      void run_molcas();

      /// read molcas log file and get the energy needed
      void read_energy();

      // purge
      void purge();
      // set value implemented
      int value_implemented() const { return 1; }

    private:
      /// prefix of molcas input file
      std::string prefix_;
      std::string xyz_file_;
      /// location of molcas program
      std::string molcas_;
      /// name of molcas input
      std::string molcas_input_;
      /// molcas command line options
      std::string molcas_options_;

      double rasscf_energy_;
      double caspt2_energy_;

      Ref<ExternPT2R12> extern_pt2r12_;
      Ref<AssignedKeyVal> extern_pt2r12_akv_;

      static ClassDesc class_desc_;
    };
} // end of namespace sc

#endif //end of header guard
