//
// svd.cc
//
// Copyright (C) 2004 Edward Valeev
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
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/basis/petite.h>

using namespace std;
using namespace sc;

void
SCF::svd_product_basis()
{
  Ref<GaussianBasisSet> bs = basis();
  int nao = bs->nbasis();
  Ref<PetiteList> pl = integral()->petite_list(bs);
  Ref<SCMatrixKit> ao_mkit = bs->matrixkit();
  Ref<TwoBodyInt> grt_eval = integral()->grt<4>();
  const double* ints = grt_eval->buffer(TwoBodyOper::eri);
  int* blocksizes = new int[1];
  blocksizes[0] = nao*(nao+1)/2;
  RefSCDimension ao2_dim = new SCDimension(blocksizes[0],1,blocksizes);

  RefSCMatrix G(ao2_dim,ao2_dim,ao_mkit);
  RefSCMatrix U(ao2_dim,ao2_dim,ao_mkit);
  RefSCMatrix V(ao2_dim,ao2_dim,ao_mkit);
  RefDiagSCMatrix Sigma(ao2_dim,ao_mkit);
  int nshell = bs->nshell();
  for(int si=0; si<nshell; si++) {
    int ni = bs->shell(si).nfunction();
    for(int sj=0; sj<=si; sj++) {
      int nj = bs->shell(sj).nfunction();

      for(int sk=0; sk<nshell; sk++) {
        int nk = bs->shell(sk).nfunction();
        for(int sl=0; sl<=sk; sl++) {
          int nl = bs->shell(sl).nfunction();

          grt_eval->compute_shell(si,sj,sk,sl);

          int ii = bs->shell_to_function(si);
          int jj = bs->shell_to_function(sj);
          int kk = bs->shell_to_function(sk);
          int ll = bs->shell_to_function(sl);
          for(int i=0; i<ni; i++) {
            int jmax = (si == sj) ? i : nj-1;
            for(int j=0; j<=jmax; j++) {
              int ij = (ii+i)*(ii+i+1)/2 + (jj+j);
              for(int k=0; k<nk; k++) {
                int lmax = (sk == sl) ? k : nl-1;
                for(int l=0; l<=lmax; l++) {
                  int kl = (kk+k)*(kk+k+1)/2 + (ll+l);

                  int ijkl = ((i*nj+j)*nk+k)*nl+l;
                  double value = ints[ijkl];
                  G.set_element(ij,kl,value);
                }
              }
            }
          }
        }
      }
    }
  }

  G.svd(U,Sigma,V);
  Sigma.print("Sigmas");

  throw std::runtime_error("RI-CLHF not yet implemented");
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:

