//
// etrain.h
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

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _mpqc_src_lib_chemistry_qc_etrain_etrain_h
#define _mpqc_src_lib_chemistry_qc_etrain_etrain_h

#include <util/class/scexception.h>
#include <util/misc/runnable.h>
#include <math/mmisc/grid.h>
#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/scf/clscf.h>
#include <chemistry/qc/basis/orthog.h>
#include <chemistry/qc/wfn/eht.h>
#include <chemistry/qc/wfn/orbitalspace.h>

namespace sc {

///
/// Class ETraIn evaluates transfer and overlap matrix in the basis of monomer SCF wave functions
///
class ETraIn: public Function, public Runnable {

  public:
  /** A KeyVal constructor is used to generate a name
      object from the input. The full list of keywords
      that are accepted is below.

      <table border="1">

      <tr><td><b>%Keyword</b><td><b>Type</b><td><b>Default</b><td><b>Description</b>

      <tr><td><tt>wfn1</tt><td>OneBodyWavefunction<td>none<td>Specifies how to compute fragment 1 orbitals
      <tr><td><tt>wfn2</tt><td>OneBodyWavefunction<td>none<td>Specifies how to compute fragment 2 orbitals
      <tr><td><tt>wfn12</tt><td>OneBodyWavefunction<td>none<td>Specifies how to compute the Fock operator
                                (typically it's a wavefunction for the combination of fragment 1 and 2)

      <tr><td><tt>nocc</tt><td>int<td>-1<td>Number of occupied orbitals from each fragment to consider. The default, -1, means to include all.
      <tr><td><tt>nuocc</tt><td>int<td>-1<td>Number of unoccupied orbitals from each fragment to consider. The default, -1, means to include all.

      <tr><td><tt>grid</tt><td>Grid<td>null<td>If specified, ETraIn will evaluate fragment 1 and 2 orbitals on this Grid.

      </table>
   */
  ETraIn(const Ref<KeyVal>&);
  // Only KeyVal constructor is provided. SavableState functionality will not be used

  // Implementation of Runnable::run()
  void run();
  // Implementation of Compute::compute()
  void compute(void);
  // Overload of Function::obsolete()
  void obsolete(void);

  private:

  // n-mer wave function
  Ref<OneBodyWavefunction> obwfn12_;
  // Monomer wave functions
  Ref<OneBodyWavefunction> obwfn1_;
  Ref<OneBodyWavefunction> obwfn2_;
  // because the computational frames for monomers may not coincide with that of the n-mer
  // compute the atom maps from monomers to the n-mer
  std::vector<unsigned int> atom_map1_;
  std::vector<unsigned int> atom_map2_;

  // Number of monomer HOMOs and LUMOs specifies the basis in which to compute transfer matrices
  int nocc_, nuocc_;
  // debug level
  unsigned int debug_;
  // grid
  Ref<Grid> grid_;

  // this function computes and prints out transfer integral (and overlap) matrix
  void compute_train();

};

}; // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
