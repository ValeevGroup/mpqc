//
// dftest.cc
//
// Copyright (C) 2001 Edward Valeev
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

#include <stdlib.h>
#include <string.h>

#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <util/keyval/keyval.h>
#include <util/group/message.h>
#include <util/group/pregtime.h>
#include <chemistry/qc/intv3/intv3.h>
#include <chemistry/qc/df/df.h>

using namespace std;
using namespace sc;

RefSCMatrix twobody_matrix(const Ref<TwoBodyInt>& tbint, TwoBodyInt::tbint_type tbtype = TwoBodyInt::eri);

int main(int argc, char **argv) {

  Ref<MessageGrp> msg = MessageGrp::initial_messagegrp(argc, argv);
  if (msg.null())
    msg = new ProcMessageGrp();
  MessageGrp::set_default_messagegrp(msg);

  Ref<RegionTimer> tim = new ParallelRegionTimer(msg, "dftest", 1, 1);

  char *infile = new char[strlen(SRCDIR) + strlen("/dftest.in") + 1];
  sprintf(infile,SRCDIR "/dftest.in");
  if (argc == 2) {
    infile = argv[1];
  }

  Ref<KeyVal> pkv(new ParsedKeyVal(infile));
  Ref<KeyVal> tkeyval(new PrefixKeyVal(pkv, ":test"));

  Ref<GaussianBasisSet>
  bs =
    require_dynamic_cast<GaussianBasisSet*> (
        tkeyval->describedclassvalue("basis").pointer(), "main\n");

  Ref<GaussianBasisSet>
  fbs =
    require_dynamic_cast<GaussianBasisSet*> (
        tkeyval->describedclassvalue("fitbasis").pointer(), "main\n");

  int tproc = tkeyval->intvalue("test_processor");
  if (tproc >= msg->n())
    tproc = 0;
  int me = msg->me();

  if (me == tproc)
    cout << "testing on processor " << tproc << endl;

  tim->enter("DensityFitting");
  Ref<Integral> integral = new IntegralV3(bs);

  // construct the Coulomb matrix
  // and reconstruct it using DensityFitting
  integral->set_basis(bs, bs, bs, bs);
  RefSCMatrix C = twobody_matrix(integral->electron_repulsion());
  C.print("Coulomb operator");
  Ref<DensityFitting> df = new DensityFitting(integral, bs, bs, fbs);
  df->compute();
  RefSCMatrix C_df = df->C().t() * df->conjugateC();
  C_df.print("Reconstructed Coulomb operator");
  (C-C_df).print("Reconstruction error");

  tim->exit();

  tim->print();
  return 0;
}


RefSCMatrix twobody_matrix(const Ref<TwoBodyInt>& tbint, TwoBodyInt::tbint_type tbtype)
{

  Ref<GaussianBasisSet> b1 = tbint->basis1();
  Ref<GaussianBasisSet> b2 = tbint->basis2();
  Ref<GaussianBasisSet> b3 = tbint->basis3();
  Ref<GaussianBasisSet> b4 = tbint->basis4();

  RefSCDimension dim12 = new SCDimension(b1->nbasis() * b2->nbasis());
  RefSCDimension dim34 = new SCDimension(b3->nbasis() * b4->nbasis());

  const int nbasis2 = b2->nbasis();
  const int nbasis4 = b4->nbasis();

  Ref<SCMatrixKit> kit = SCMatrixKit::default_matrixkit();
  RefSCMatrix result = kit->matrix(dim12,dim34);

  const double* buffer = tbint->buffer(tbtype);
  for (int s1 = 0; s1 < b1->nshell(); ++s1) {
    const int s1offset = b1->shell_to_function(s1);
    const int nf1 = b1->shell(s1).nfunction();
    for (int s2 = 0; s2 < b2->nshell(); ++s2) {
      const int s2offset = b2->shell_to_function(s2);
      const int nf2 = b2->shell(s2).nfunction();

      for (int s3 = 0; s3 < b3->nshell(); ++s3) {
        const int s3offset = b3->shell_to_function(s3);
        const int nf3 = b3->shell(s3).nfunction();
        for (int s4 = 0; s4 < b4->nshell(); ++s4) {
          const int s4offset = b4->shell_to_function(s4);
          const int nf4 = b4->shell(s4).nfunction();

          // compute shell quartet
          tbint->compute_shell(s1, s2, s3, s4);

          // copy buffer into result
          const double* bufptr = buffer;
          int f12 = 0;
          for(int f1=0; f1<nf1; ++f1) {
            const int s12offset = (s1offset+f1) * nbasis2 + s2offset;
            for(int f2=0; f2<nf2; ++f2, ++f12) {

              int f34 = 0;
              for(int f3=0; f3<nf3; ++f3) {
                const int s34offset = (s3offset+f3) * nbasis4 + s4offset;
                for(int f4=0; f4<nf4; ++f4, ++f34, ++bufptr) {
                  result.set_element(f2+s12offset, f4+s34offset, *bufptr);
                }
              }
            }
          }

        }
      }
    }
  }

  return result;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
