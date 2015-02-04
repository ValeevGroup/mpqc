//
// extern_pt2r12.h
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

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_externpt2r12_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_externpt2r12_h

#include <chemistry/qc/mbptr12/pt2r12.h>
#include <extern/moinfo/moinfo.h>

namespace sc {

  /// ExternPT2R12 is a PT2R12 wave function computed from external MO info and 2-RDM
  class ExternPT2R12 : public Wavefunction {
    public:
      /** A KeyVal constructor is used to generate a ExternPT2R12
          object from the input. In addition to all keywords of Wavefunction,
          the following additional keywords will be queried:

          <table border="1">

          <tr><td><b>%Keyword</b><td><b>Type</b><td><b>Default</b><td><b>Description</b>

          <tr><td><tt>world</tt><td>WavefunctionWorld<td>none<td>this wavefunction will own this world (@sa WavefunctionWorld::set_wfn )

          <tr><td><tt>rdm2</tt><td>SpinFreeRDM<Two><td>none<td>the 2-particle reduced density matrix will be provided by this object

          <tr><td><tt>orbs_info</tt><td>ExternMOInfo<td>none<td>the information about the orbitals used by this object

          <tr><td><tt>cabs</tt><td>string<td>none<td>specifies the name of the RI basis from which CABS will be constructed (@sa GaussianBasisSet::GaussianBasisSet)

          <tr><td><tt>f12exp</tt><td>string<td>none<td>specifies the exponent of the Slater geminal

          </table>
       */
      ExternPT2R12(const Ref<KeyVal>& kv);
      ~ExternPT2R12();

      void compute();
      int nelectron();
      RefSymmSCMatrix density();
      double magnetic_moment() const;
      int value_implemented() const { return 1; }
      void set_desired_value_accuracy(double acc);
      void print(std::ostream& os=ExEnv::out0());
      void initialize();

    private:
      static ClassDesc class_desc_;
      static const unsigned int debug_print_ = 0; // set to 1 to print out some debugging info

      // need to initialize with initialize()
      Ref<PT2R12> pt2r12_;


    protected:
      // provided by the user
      Ref<WavefunctionWorld> world_;
      Ref<ExternMOInfo> orbs_info_;
      Ref<SpinFreeRDM<Two> > rdm2_;

      std::string cabs_name_;
      std::string obs_name_;
      std::string dfbs_name_;
      std::string f12exp_str_;
      std::string r12_str_;

      #if defined(HAVE_MPQC3_RUNTIME)
          std::string singles_str_;
          std::string partition_str_;
          std::string cabs_singles_name_;
      #endif

      bool cabs_contraction_;

  };


} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
