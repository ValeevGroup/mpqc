//
// multipole_ints.cc
//
// Copyright (C) 2003 Edward Valeev
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdexcept>
#include <math/scmat/blocked.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>

using namespace std;
using namespace sc;

void
R12IntEvalInfo::compute_multipole_ints(RefSymmSCMatrix& MX, RefSymmSCMatrix& MY, RefSymmSCMatrix& MZ,
				       RefSymmSCMatrix& MXX, RefSymmSCMatrix& MYY, RefSymmSCMatrix& MZZ)
{
  Ref<PetiteList> pl = ref_->integral()->petite_list();
  int nshell = bs_->nshell();

  // Have to convert scf_vec_ to unblocked form, get_subblock, and convert back to blocked form
  // because get_subblock doesn't like blocked matrices
  RefSCDimension aodim = pl->AO_basisdim();
  RefSCDimension modim;
  {
    int nmo = scf_vec_.rowdim().n();
    int blksize[1];  blksize[0] = nmo;
    modim = new SCDimension(nmo,1,blksize,"Sorted MO dimension");
  }
  RefSCMatrix unblvec(modim,aodim,bs_->matrixkit());
  unblvec->convert(scf_vec_.pointer());
  RefSCMatrix subblock = unblvec.get_subblock(nfzc_,nocc_-1,0,bs_->nbasis()-1);
  RefSCDimension occactdim;
  {
    int blksize[1]; blksize[0] = nocc_ - nfzc_;
    occactdim = new SCDimension(nocc_-nfzc_,1,blksize,"Active occupied MO dimension");
    occactdim->blocks()->set_subdim(0,new SCDimension(nocc_-nfzc_,"occact MO dim"));
  }
  RefSCMatrix OccAct_Vec(occactdim,aodim,bs_->so_matrixkit());
  OccAct_Vec->convert(subblock.pointer());
  subblock = 0;
  unblvec=0;

  Ref<OneBodyInt> m1_ints = integral_->dipole(0);

  // form AO dipole moment matrices
  RefSymmSCMatrix mx(pl->AO_basisdim(), bs_->so_matrixkit());
  RefSymmSCMatrix my(pl->AO_basisdim(), bs_->so_matrixkit());
  RefSymmSCMatrix mz(pl->AO_basisdim(), bs_->so_matrixkit());
  mx.assign(0.0);
  my.assign(0.0);
  mz.assign(0.0);
    
  for(int sh1=0; sh1<nshell; sh1++) {
    int bf1_offset = bs_->shell_to_function(sh1);
    int nbf1 = bs_->shell(sh1).nfunction();
    for(int sh2=0; sh2<=sh1; sh2++) {
      int bf2_offset = bs_->shell_to_function(sh2);
      int nbf2 = bs_->shell(sh2).nfunction();
      
      m1_ints->compute_shell(sh1,sh2);
      const double *m1intsptr = m1_ints->buffer();
      int bf1_index = bf1_offset;
      for(int bf1=0; bf1<nbf1; bf1++, bf1_index++, m1intsptr+=3*nbf2) {
	int bf2_index = bf2_offset;
	const double *ptr1 = m1intsptr;
	int bf2max;
	if (sh1 != sh2)
	  bf2max = nbf2-1;
	else
	  bf2max = bf1;
	for(int bf2=0; bf2<=bf2max; bf2++, bf2_index++) {
	  mx.set_element(bf1_index, bf2_index, *(ptr1++));
	  my.set_element(bf1_index, bf2_index, *(ptr1++));
	  mz.set_element(bf1_index, bf2_index, *(ptr1++));
	}
      }
    }
  }

  // and clean up a bit
  m1_ints = 0;

  // finally, transform
  MX = bs_->so_matrixkit()->symmmatrix(occactdim);
  MY = bs_->so_matrixkit()->symmmatrix(occactdim);
  MZ = bs_->so_matrixkit()->symmmatrix(occactdim);
  MX.assign(0.0);
  MY.assign(0.0);
  MZ.assign(0.0);
  MX.accumulate_transform(OccAct_Vec,mx);
  mx = 0;
  MY.accumulate_transform(OccAct_Vec,my);
  my = 0;
  MZ.accumulate_transform(OccAct_Vec,mz);
  mz = 0;

  // same for quadrupole integrals
  Ref<OneBodyInt> m2_ints = integral_->quadrupole(0);
  RefSymmSCMatrix mxx(pl->AO_basisdim(), bs_->so_matrixkit());
  RefSymmSCMatrix myy(pl->AO_basisdim(), bs_->so_matrixkit());
  RefSymmSCMatrix mzz(pl->AO_basisdim(), bs_->so_matrixkit());
  mxx.assign(0.0);
  myy.assign(0.0);
  mzz.assign(0.0);
    
  for(int sh1=0; sh1<nshell; sh1++) {
    int bf1_offset = bs_->shell_to_function(sh1);
    int nbf1 = bs_->shell(sh1).nfunction();
    for(int sh2=0; sh2<=sh1; sh2++) {
      int bf2_offset = bs_->shell_to_function(sh2);
      int nbf2 = bs_->shell(sh2).nfunction();
      
      m2_ints->compute_shell(sh1,sh2);
      const double *m2intsptr = m2_ints->buffer();
      int bf1_index = bf1_offset;
      for(int bf1=0; bf1<nbf1; bf1++, bf1_index++, m2intsptr+=6*nbf2) {
	int bf2_index = bf2_offset;
	const double *ptr2 = m2intsptr;
	int bf2max;
	if (sh1 != sh2)
	  bf2max = nbf2-1;
	else
	  bf2max = bf1;
	for(int bf2=0; bf2<=bf2max; bf2++, bf2_index++) {
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
  m2_ints = 0;
  pl = 0;

  // transform
  MXX = bs_->so_matrixkit()->symmmatrix(occactdim);
  MYY = bs_->so_matrixkit()->symmmatrix(occactdim);
  MZZ = bs_->so_matrixkit()->symmmatrix(occactdim);
  MXX.assign(0.0);
  MYY.assign(0.0);
  MZZ.assign(0.0);
  MXX.accumulate_transform(OccAct_Vec,mxx);
  mxx = 0;
  MYY.accumulate_transform(OccAct_Vec,myy);
  myy = 0;
  MZZ.accumulate_transform(OccAct_Vec,mzz);
  mzz = 0;

  if (debug_ > 1) {
    MX.print("mu(X) in active occupied MO basis");
    MY.print("mu(Y) in active occupied MO basis");
    MZ.print("mu(Z) in active occupied MO basis");
    MXX.print("mu(XX) in active occupied MO basis");
    MYY.print("mu(YY) in active occupied MO basis");
    MZZ.print("mu(ZZ) in active occupied MO basis");
  }

  OccAct_Vec = 0;

  return;
}

///////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
