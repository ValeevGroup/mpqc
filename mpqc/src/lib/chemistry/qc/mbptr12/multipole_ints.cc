//
// multipole_ints.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdexcept>
#include <math/scmat/blocked.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbptr12/transform_factory.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>

using namespace std;
using namespace sc;

void
R12IntEvalInfo::compute_multipole_ints(const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2,
                                       RefSCMatrix& MX, RefSCMatrix& MY, RefSCMatrix& MZ,
				       RefSCMatrix& MXX, RefSCMatrix& MYY, RefSCMatrix& MZZ)
{
  if (!space1->integral()->equiv(space2->integral()))
    throw ProgrammingError("two OrbitalSpaces use incompatible Integral factories");
  const Ref<GaussianBasisSet> bs1 = space1->basis();
  const Ref<GaussianBasisSet> bs2 = space2->basis();
  const bool bs1_eq_bs2 = (bs1 == bs2);
  int nshell1 = bs1->nshell();
  int nshell2 = bs2->nshell();

  RefSCMatrix vec1t = space1->coefs().t();
  RefSCMatrix vec2 = space2->coefs();

  Ref<Integral> localints = space1->integral()->clone();
  localints->set_basis(bs1,bs2);

  Ref<OneBodyInt> m1_ints = localints->dipole(0);
  Ref<OneBodyInt> m2_ints = localints->quadrupole(0);

  // form AO moment matrices
  RefSCDimension aodim1 = vec1t.coldim();
  RefSCDimension aodim2 = vec2.rowdim();
  Ref<SCMatrixKit> aokit = bs1->so_matrixkit();
  RefSCMatrix mx(aodim1, aodim2, aokit);
  RefSCMatrix my(aodim1, aodim2, aokit);
  RefSCMatrix mz(aodim1, aodim2, aokit);
  RefSCMatrix mxx(aodim1, aodim2, aokit);
  RefSCMatrix myy(aodim1, aodim2, aokit);
  RefSCMatrix mzz(aodim1, aodim2, aokit);
  mx.assign(0.0);
  my.assign(0.0);
  mz.assign(0.0);
  mxx.assign(0.0);
  myy.assign(0.0);
  mzz.assign(0.0);
  
  for(int sh1=0; sh1<nshell1; sh1++) {
    int bf1_offset = bs1->shell_to_function(sh1);
    int nbf1 = bs1->shell(sh1).nfunction();

    int sh2max;
    if (bs1_eq_bs2)
      sh2max = sh1;
    else
      sh2max = nshell2-1;
    
    for(int sh2=0; sh2<=sh2max; sh2++) {
      int bf2_offset = bs2->shell_to_function(sh2);
      int nbf2 = bs2->shell(sh2).nfunction();
      
      m1_ints->compute_shell(sh1,sh2);
      const double *m1intsptr = m1_ints->buffer();

      m2_ints->compute_shell(sh1,sh2);
      const double *m2intsptr = m2_ints->buffer();

      int bf1_index = bf1_offset;
      for(int bf1=0; bf1<nbf1; bf1++, bf1_index++, m1intsptr+=3*nbf2, m2intsptr+=6*nbf2) {
	int bf2_index = bf2_offset;
	const double *ptr1 = m1intsptr;
        const double *ptr2 = m2intsptr;
	int bf2max;
        if (bs1_eq_bs2 && sh1 == sh2)
          bf2max = bf1;
        else
	  bf2max = nbf2-1;
	for(int bf2=0; bf2<=bf2max; bf2++, bf2_index++) {

	  mx.set_element(bf1_index, bf2_index, *(ptr1++));
	  my.set_element(bf1_index, bf2_index, *(ptr1++));
	  mz.set_element(bf1_index, bf2_index, *(ptr1++));

          mxx.set_element(bf1_index, bf2_index, *(ptr2++));
          ptr2 += 2;
          myy.set_element(bf1_index, bf2_index, *(ptr2++));
          ptr2++;
          mzz.set_element(bf1_index, bf2_index, *(ptr2++));

        }
      }
    }
  }

  // and clean up a bit
  m1_ints = 0;
  m2_ints = 0;

  // Symmetrize matrices, if necessary
  if (bs1_eq_bs2) {

    const int nbasis = bs1->nbasis();
    
    for(int bf1=0; bf1<nbasis; bf1++)
      for(int bf2=0; bf2<=bf1; bf2++) {
        mx(bf2,bf1) = mx(bf1,bf2);
        my(bf2,bf1) = my(bf1,bf2);
        mz(bf2,bf1) = mz(bf1,bf2);
        mxx(bf2,bf1) = mxx(bf1,bf2);
        myy(bf2,bf1) = myy(bf1,bf2);
        mzz(bf2,bf1) = mzz(bf1,bf2);
      }

  }
      

  // finally, transform
  MX = vec1t * mx * vec2;
  MY = vec1t * my * vec2;
  MZ = vec1t * mz * vec2;
  MXX = vec1t * mxx * vec2;
  MYY = vec1t * myy * vec2;
  MZZ = vec1t * mzz * vec2;
  
  // and clean up a bit
  mx = 0;
  my = 0;
  mz = 0;
  mxx = 0;
  myy = 0;
  mzz = 0;

  //if (debug_ > 1) {
  //  MX.print("mu(X)");
  //  MY.print("mu(Y)");
  //  MZ.print("mu(Z)");
  //  MXX.print("mu(XX)");
  //  MYY.print("mu(YY)");
  //  MZZ.print("mu(ZZ)");
  //}
}


void
R12IntEvalInfo::compute_overlap_ints(const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2,
                                     RefSCMatrix& S)
{
  if (!space1->integral()->equiv(space2->integral()))
    throw ProgrammingError("two OrbitalSpaces use incompatible Integral factories");
  const Ref<GaussianBasisSet> bs1 = space1->basis();
  const Ref<GaussianBasisSet> bs2 = space2->basis();
  const bool bs1_eq_bs2 = (bs1 == bs2);
  int nshell1 = bs1->nshell();
  int nshell2 = bs2->nshell();

  RefSCMatrix vec1t = space1->coefs().t();
  RefSCMatrix vec2 = space2->coefs();

  Ref<Integral> localints = space1->integral()->clone();
  localints->set_basis(bs1,bs2);

  Ref<OneBodyInt> ov_ints = localints->overlap();

  // form AO moment matrices
  RefSCDimension aodim1 = vec1t.coldim();
  RefSCDimension aodim2 = vec2.rowdim();
  Ref<SCMatrixKit> aokit = bs1->so_matrixkit();
  RefSCMatrix s(aodim1, aodim2, aokit);
  s.assign(0.0);

  for(int sh1=0; sh1<nshell1; sh1++) {
    int bf1_offset = bs1->shell_to_function(sh1);
    int nbf1 = bs1->shell(sh1).nfunction();

    int sh2max;
    if (bs1_eq_bs2)
      sh2max = sh1;
    else
      sh2max = nshell2-1;

    for(int sh2=0; sh2<=sh2max; sh2++) {
      int bf2_offset = bs2->shell_to_function(sh2);
      int nbf2 = bs2->shell(sh2).nfunction();

      ov_ints->compute_shell(sh1,sh2);
      const double *ovintsptr = ov_ints->buffer();

      int bf1_index = bf1_offset;
      for(int bf1=0; bf1<nbf1; bf1++, bf1_index++, ovintsptr+=nbf2) {
        int bf2_index = bf2_offset;
        const double *ptr = ovintsptr;
        int bf2max;
        if (bs1_eq_bs2 && sh1 == sh2)
          bf2max = bf1;
        else
          bf2max = nbf2-1;
        for(int bf2=0; bf2<=bf2max; bf2++, bf2_index++) {

          s.set_element(bf1_index, bf2_index, *(ptr++));

        }
      }
    }
  }

  // and clean up a bit
  ov_ints = 0;

  // Symmetrize matrices, if necessary
  if (bs1_eq_bs2) {

    const int nbasis = bs1->nbasis();

    for(int bf1=0; bf1<nbasis; bf1++)
      for(int bf2=0; bf2<=bf1; bf2++) {
        s(bf2,bf1) = s(bf1,bf2);
      }

  }


    // finally, transform
    S = vec1t * s * vec2;

    // and clean up a bit
    s = 0;

    //if (debug_ > 1) {
    //  S.print("Overlap");
    //}
}

///////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
