//
// fock.cc
//
// Copyright (C) 2004 Edward Valeev
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
#include <chemistry/qc/mbptr12/r12int_eval.h>

using namespace std;
using namespace sc;

RefSCMatrix
R12IntEval::fock_(const Ref<MOIndexSpace>& occ_space, const Ref<MOIndexSpace>& bra_space,
                  const Ref<MOIndexSpace>& ket_space, double scale_J, double scale_K)
{
  const Ref<GaussianBasisSet> bs1 = bra_space->basis();
  const Ref<GaussianBasisSet> bs2 = ket_space->basis();
  const bool bs1_eq_bs2 = (bs1 == bs2);
  int nshell1 = bs1->nshell();
  int nshell2 = bs2->nshell();

  RefSCMatrix vec1t = bra_space->coefs().t();
  RefSCMatrix vec2 = ket_space->coefs();

  Ref<Integral> localints = r12info_->integral()->clone();
  localints->set_basis(bs1,bs2);

  Ref<OneBodyInt> h_ints = localints->hcore();

  // form AO moment matrices
  RefSCDimension aodim1 = vec1t.coldim();
  RefSCDimension aodim2 = vec2.rowdim();
  Ref<SCMatrixKit> aokit = bs1->so_matrixkit();
  RefSCMatrix h(aodim1, aodim2, aokit);
  h.assign(0.0);
  
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
      
      h_ints->compute_shell(sh1,sh2);
      const double *hintsptr = h_ints->buffer();

      int bf1_index = bf1_offset;
      for(int bf1=0; bf1<nbf1; bf1++, bf1_index++, hintsptr+=nbf2) {
	int bf2_index = bf2_offset;
	const double *ptr = hintsptr;
	int bf2max;
        if (bs1_eq_bs2 && sh1 == sh2)
          bf2max = bf1;
        else
	  bf2max = nbf2-1;
	for(int bf2=0; bf2<=bf2max; bf2++, bf2_index++, ptr++) {

	  h.set_element(bf1_index, bf2_index, *ptr);

        }
      }
    }
  }

  // Symmetrize matrices, if necessary
  if (bs1_eq_bs2) {
    const int nbasis = bs1->nbasis();
    for(int bf1=0; bf1<nbasis; bf1++)
      for(int bf2=0; bf2<=bf1; bf2++) {
        h(bf2,bf1) = h(bf1,bf2);
      }
  }

  // finally, transform
  RefSCMatrix F = vec1t * h * vec2;

  // add coulomb and exchange parts
  if (scale_J != 0.0) {
    RefSCMatrix J = coulomb_(occ_space,bra_space,ket_space);
    J.scale(2.0*scale_J); F.accumulate(J); J = 0;
  }
  if (scale_K != 0.0) {
    RefSCMatrix K = exchange_(occ_space,bra_space,ket_space);
    K.scale(-1.0*scale_K); F.accumulate(K); K = 0;
  }
  
  // and clean up a bit
  h_ints = 0;
  h = 0;

  if (debug_ > 1) {
    F.print("Fock matrix");
  }
  
  return F;
}

void
R12IntEval::compute_norms_(const RefSCMatrix& A, const std::string& label, std::ostream& os)
{
  Ref<SCElementMaxAbs> maxabs_op(new SCElementMaxAbs);
  A.element_op(maxabs_op);
  const double maxabs = maxabs_op->result();

  Ref<SCElementKNorm> onenorm_op(new SCElementKNorm(1.0));
  A.element_op(onenorm_op);
  const double onenorm = onenorm_op->result();

  Ref<SCElementKNorm> twonorm_op(new SCElementKNorm(2.0));
  A.element_op(twonorm_op);
  const double twonorm = twonorm_op->result();

  os << indent << "Norms of " << label << endl;
  os << indent << "------------------------" << endl;
  os << indent << "||A||_{\\infty} = " << scprintf("%10.5lf",maxabs) << endl;
  os << indent << "||A||_1        = " << scprintf("%10.5lf",onenorm) << endl;
  os << indent << "||A||_2        = " << scprintf("%10.5lf",twonorm) << endl << endl;
}

///////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
