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
  // //////////////////////////////////////////////////////////////////////

/** The MolcasPT2R12 class is a interface between Molcas and MPQC to perform
* CASPT2F12 calculations. The Molcas program should provide prefix.pt2r12.dat and
* prefix.pt2r12.rdm2.dat files for input to MPQC. Those two inputfiles are used
* to create the object ExternPT2R12 in MolcasPT2R12.
* Here is an example of Molcas input file. It should perform a CASPT2 calculation.
* Please put the option and value in one single line
<pre>
  *-------------------------------------------------------------------------------
  * Molecule: HF
  * Basis: cc-pVDZ-F12
  * Symmetry:
  *-------------------------------------------------------------------------------

  &GATEWAY
    coord=prefix.xyz
    basis=cc-pVDZ-F12
    Group=C1
  *-------------------------------------------------------------------------------
  &SEWARD
  *-------------------------------------------------------------------------------
  &SCF &END
  End of input
  *-------------------------------------------------------------------------------
  &RASSCF &END
    nActEl= 6 0 0
    Symmetry= 1
    Spin= 1
    Inactive= 2 0 0 0
    Ras2= 2 1 0 1
    LumOrb
  End of input
  *-------------------------------------------------------------------------------
  &CASPT2 &End
    Convergence = 1.0D-8
  End of Input
  *-------------------------------------------------------------------------------
</pre>
 *
 * the must have option in the input file is
 *
<pre>
  &GATEWAY
    coord=
    basis=
    Group=
  &SEWARD
  &SCF
  &RASSCF
    nActEl=
    Inactive=
    Ras2=
  &CASPT2
</pre>

*/


  class MolcasPT2R12 : public MolecularEnergy {
    public:

    /** The KeyVal constructor reads the following keywords:

         <dl>

          All the KeyVal keywords needed for MolcasPT2R12 class

          <table border="1">

          <tr><td><b>%Keyword</b><td><b>Type</b><td><b>Default</b><td><b>Description</b>

          <tr><td><tt>prefix</tt><td> string <td>none<td> mandatory filename prefix,
          will look for files: prefix.pt2r12.dat and prefix.pt2r12.rdm2.dat\n

          <tr><td><tt>molcas</tt><td> string <td>molcas<td> molcas program name,
          the default is molcas

          <tr><td><tt>molcas_input</tt><td> string <td>none<td> molcas input file,
          it should perform CASPT2 calculation and provide prefix.pt2r12.dat and
          prefix.pt2r12.rdm2.dat files after molcas calculation

          <tr><td><tt>molcas_options</tt><td> string <td> -f <td> molcas command
          line options, default is -f, which will create prefix.log file

          <tr><td><tt>xyz_file</tt><td> string <td>none<td> geometry xyz file, should
          not define geometry inside Molcas or MPQC input file because Molcas and
          MPQC will rely on the same geometry xyz file.

          <tr><td><tt>obs</tt><td> string <td>none<td>  name for the orbital basis
          set; optional, if given will be used to set the defaults for cabs, dfbs,
          and f12exp

          <tr><td><tt>cabs</tt><td>string<td>based on obs <td>specifies
          the name of the RI basis from which CABS will be constructed. If not given
          it will use obs to choose the default.

          <tr><td><tt>cabs_contraction</tt><td>string<td> true <td> if use contracted cabs
          basis sets

          <tr><td><tt>dfbs</tt><td>string<td>based on obs<td>  name for DFBS; default: no
          density fitting; use "none" to override the default for the obs

          <tr><td><tt>f12exp</tt><td>string<td> based on obs <td> specifies the
          exponent of the Slater geminal, default constructed based on obs, if
          not, the default is 1.0

          <tr><td><tt>r12</tt><td>string<td> true <td> if compute [2]_R12 correction

          </table>

          The following keywords are MPQC_NEW_FEATHRES only, which include keywords
          for [2]_S calculation:

          <table border="1">

          <tr><td><b>%Keyword</b><td><b>Type</b><td><b>Default</b><td><b>Description</b>

          <tr><td><tt>cabs_singles </tt><td>string<td> true <td> if compute [2]_S correction

          <tr><td><tt>singles_basis </tt><td>string<td> none <td> if use a different
          CABS basis for [2]_S calculation, the default is to use the same CABS basis
          sets as [2]_R12

          <tr><td><tt>partitionH </tt><td>string<td> fock <td> the way to How to
          partition Hamiltonian in [2]_S: fock, dyall_1, dyall_2; default: fock"

          </table>

         </dl> */
      MolcasPT2R12(const Ref<KeyVal>& kv);
      ~MolcasPT2R12();

      void compute();

      /// initialize values needed for ExternPT2R12
      void initialize();

      // run molcas input file
      void run_molcas();

      /// read molcas log file and get the energy needed
      void read_energy();

      /// convert input to c1 symmetry when the symmetry from MPQC
      /// is different than original symmetry
      void convert_c1_symmetry();

      /// restore the original input file after change of symmetry
      void restore_molcas_input();

      /// purge
      void purge();

      /// obsolete
      void obsolete();

    /// set value implemented
      int value_implemented() const { return 1; }

    /// print function
      void print(std::ostream& os=ExEnv::out0()) const;

    private:
      /// prefix of molcas input file
      std::string prefix_;
      /// xyz file name
      std::string xyz_file_;
      /// location of molcas program
      std::string molcas_;
      /// name of molcas input file
      std::string molcas_input_;
      /// molcas command line options
      std::string molcas_options_;
      /// inactive orbital input
      std::vector<int> inactive_;
      /// active orbital input
      std::vector<int> active_;
      /// number of active electron
      std::vector<int> active_electron_;
      /// sysmetry
      std::string symmetry_;

      double rasscf_energy_;
      double caspt2_energy_;

      /// ref to ExternPT2R12 object
      Ref<ExternPT2R12> extern_pt2r12_;
      /// keyval to store key needed for ExternPT2R12
      Ref<AssignedKeyVal> extern_pt2r12_akv_;

      static ClassDesc class_desc_;
    };
} // end of namespace sc

#endif //end of header guard
