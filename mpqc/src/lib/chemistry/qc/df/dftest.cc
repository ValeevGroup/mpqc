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
#include <math/scmat/local.h>
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

  int me = msg->me();
  int nproc = msg->n();
  cout << "testing on " << nproc << " processors" << endl;

  // use test::DensityFitting to test first
  {
    using sc::test::DensityFitting;

    tim->enter("test::DensityFitting");
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
    (C - C_df).print("Reconstruction error");

    tim->exit();
  }

  // use DensityFitting to test first
  {
    using sc::DensityFitting;

    tim->enter("DensityFitting");
    Ref<Integral> integral = new IntegralV3(bs);

    // construct the Coulomb matrix
    // and reconstruct it using DensityFitting
    integral->set_basis(bs, bs, bs, bs);
    RefSCMatrix C = twobody_matrix(integral->electron_repulsion());
    C.print("Coulomb operator");

    Ref<OrbitalSpace> bs_space = new AtomicOrbitalSpace("mu", "AO(OBS)", bs, integral);
    Ref<OrbitalSpace> fbs_space = new AtomicOrbitalSpace("Mu", "AO(FBS)", fbs, integral);
    OrbitalSpaceRegistry::instance()->add(make_keyspace_pair(bs_space));
    OrbitalSpaceRegistry::instance()->add(make_keyspace_pair(fbs_space));
    AOSpaceRegistry::instance()->add(bs, bs_space);
    AOSpaceRegistry::instance()->add(fbs, fbs_space);
    Ref<MOIntsTransformFactory> factory = new MOIntsTransformFactory(integral);
    factory->set_ints_method(MOIntsTransformFactory::StoreMethod::posix);
    Ref<DensityFitting> df = new DensityFitting(factory, "1/r_{12}", bs_space, bs_space, fbs);
    df->compute();

    RefSCMatrix C_df = C.clone(); C_df.assign(0.0);
    // assemble reconstructed matrix
    {
      Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
      Ref<R12IntsAcc> C = df->C();
      C->activate();
      Ref<R12IntsAcc> cC = df->conjugateC();
      C->activate();
      const int nbs = bs->nbasis();
      const int nfbs = fbs->nbasis();
      RefSCDimension bsdim = new SCDimension(nbs);
      RefSCDimension fbsdim = new SCDimension(nfbs);
      RefSCMatrix C_jR = localkit->matrix(bsdim, fbsdim);
      RefSCMatrix cC_jR = localkit->matrix(bsdim, fbsdim);
      RefSCMatrix cCC = localkit->matrix(bsdim, bsdim);
      RefSCMatrix C_df_localcopy = localkit->matrix(C_df.rowdim(),
                                                    C_df.coldim());
      C_df_localcopy.assign(0.0);
      if (me == 0) {
        for (int i1 = 0; i1 < nbs; ++i1) {
          const double* C_jR_buf = C->retrieve_pair_block(0, i1,
                                                          TwoBodyInt::eri);
          C_jR.assign(C_jR_buf);
          C->release_pair_block(0, i1, TwoBodyInt::eri);
          RefSCMatrix C_jR_t = C_jR.t();
          for (int i2 = 0; i2 < nbs; ++i2) {
            const double* cC_jR_buf = cC->retrieve_pair_block(0, i2,
                                                              TwoBodyInt::eri);
            cC_jR.assign(cC_jR_buf);
            cC->release_pair_block(0, i2, TwoBodyInt::eri);

            cCC.assign(0.0);
            cCC.accumulate_product(cC_jR, C_jR_t);
            C_df_localcopy.assign_subblock(cCC, i2 * nbs, (i2 + 1) * nbs - 1,
                                           i1 * nbs, (i1 + 1) * nbs - 1);

          }
        }
      }
      C_df->convert(C_df_localcopy);
    }

    C_df.print("Reconstructed Coulomb operator");
    (C - C_df).print("Reconstruction error");

    tim->exit();
  }

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
