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
  // Only KeyVal constructor is provided. SavableState functionality will not be used
  ETraIn(const Ref<KeyVal>&);

  // Implementation of Runnable::run()
  void run();
  // Implementation of Compute::compute()
  void compute(void);
  // Overload of Function::obsolete()
  void obsolete(void);

  private:

  // n-mer wave function
  Ref<OneBodyWavefunction> scf12_;
  // Monomer wave functions
  Ref<OneBodyWavefunction> scf1_;
  Ref<OneBodyWavefunction> scf2_;
  // because the computational frames for monomers may not coincide with that of the n-mer
  // compute the atom maps from monomers to the n-mer
  std::vector<unsigned int> atom_map1_;
  std::vector<unsigned int> atom_map2_;

  // Number of monomer HOMOs and LUMOs specifies the basis in which to compute transfer matrices
  int nhomo_, nlumo_;
  // debug level
  unsigned int debug_;

  // this function computes and prints out transfer integral (and overlap) matrix
  void compute_train();

};

}; // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
