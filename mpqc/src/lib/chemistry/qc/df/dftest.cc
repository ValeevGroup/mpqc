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
#include <chemistry/qc/mbptr12/transform_ixjy_df.h>
#include <chemistry/qc/mbptr12/fockbuilder.h>

using namespace std;
using namespace sc;

RefSCMatrix twobody_matrix(const Ref<TwoBodyInt>& tbint,
                           TwoBodyOper::type tbtype = TwoBodyOper::eri);
RefSCMatrix twobody_matrix(const Ref<TwoBodyThreeCenterInt>& tbint,
                           TwoBodyOper::type tbtype = TwoBodyOper::eri);
RefSCMatrix twobody_matrix(const Ref<TwoBodyTwoCenterInt>& tbint,
                           TwoBodyOper::type tbtype = TwoBodyOper::eri);

int main(int argc, char **argv) {

  Ref<MessageGrp> msg = MessageGrp::initial_messagegrp(argc, argv);
  if (msg.null())
    msg = new ProcMessageGrp();
  MessageGrp::set_default_messagegrp(msg);

  Ref<Integral> integral = Integral::initial_integral(argc, argv);
  if (integral.null())
    integral = new IntegralV3();
  Integral::set_default_integral(integral);

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
    Ref<Integral> integral = Integral::get_default_integral();

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

#if 0
  // now, do the same by using DensityFitting directly
  {
    using sc::DensityFitting;

    tim->enter("DensityFitting");
    Ref<Integral> integral = Integral::get_default_integral();

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
                                                          TwoBodyOper::eri);
          C_jR.assign(C_jR_buf);
          C->release_pair_block(0, i1, TwoBodyOper::eri);
          RefSCMatrix C_jR_t = C_jR.t();
          for (int i2 = 0; i2 < nbs; ++i2) {
            const double* cC_jR_buf = cC->retrieve_pair_block(0, i2,
                                                              TwoBodyOper::eri);
            cC_jR.assign(cC_jR_buf);
            cC->release_pair_block(0, i2, TwoBodyOper::eri);

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
#endif

#if 0
  // now, do the same by using TwoBodyMOIntsTransform_ixjy_df
  {
    tim->enter("DensityFitting ixjy transform");

    Ref<Integral> integral = Integral::get_default_integral();
    // construct the Coulomb matrix
    // and reconstruct it using DensityFitting
    integral->set_basis(bs, bs, bs, bs);
    RefSCMatrix C = twobody_matrix(integral->electron_repulsion());
    C.print("Coulomb operator");

    Ref<OrbitalSpace> bs_space = new AtomicOrbitalSpace("mu", "AO(OBS)", bs, integral);
    OrbitalSpaceRegistry::instance()->add(make_keyspace_pair(bs_space));
    AOSpaceRegistry::instance()->add(bs, bs_space);
    Ref<MOIntsTransformFactory> factory = new MOIntsTransformFactory(integral);
    factory->set_ints_method(MOIntsTransformFactory::StoreMethod::posix);
    Ref<TwoBodyMOIntsTransform> ixjy_tform =
      new TwoBodyMOIntsTransform_ixjy_df("test", factory,
                                         new TwoBodyIntDescrERI(integral),
                                         bs_space, bs_space,
                                         bs_space, bs_space,
                                         fbs, fbs);

    ixjy_tform->compute();
    Ref<R12IntsAcc> ixjy_acc = ixjy_tform->ints_acc();
    ixjy_acc->activate();

    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    const int nbs = bs->nbasis();
    RefSCDimension bsdim = new SCDimension(nbs);
    RefSCMatrix C_df_ij = localkit->matrix(bsdim, bsdim);
    RefSCMatrix C_df = C.clone(); C_df.assign(0.0);
    RefSCMatrix C_df_localcopy = localkit->matrix(C_df.rowdim(),
                                                  C_df.coldim());
    C_df_localcopy.assign(0.0);
    if (me == 0) {
      for (int i1 = 0; i1 < nbs; ++i1) {
        for (int i2 = 0; i2 < nbs; ++i2) {
          const double* buf = ixjy_acc->retrieve_pair_block(i1, i2,
                                                            TwoBodyOper::eri);
          C_df_ij.assign(buf);
          ixjy_acc->release_pair_block(i1, i2, TwoBodyOper::eri);

          C_df_localcopy.assign_subblock(C_df_ij,
                                         i1 * nbs, (i1 + 1) * nbs - 1,
                                         i2 * nbs, (i2 + 1) * nbs - 1);

        }
      }
    }
    ixjy_acc->deactivate();

    C_df->convert(C_df_localcopy);
    C_df.print("Reconstructed Coulomb operator");
    (C - C_df).print("Reconstruction error");

    tim->exit();
  }
#endif

  typedef IntParamsG12::PrimitiveGeminal PrimitiveGeminal;
  std::vector<PrimitiveGeminal> geminal(1, std::make_pair(2.0, 1.0));
  Ref<IntParamsG12> g12params = new IntParamsG12(geminal);

  // use test::DensityFitting to test Gaussian geminal integrals
  {
    using sc::test::DensityFitting;

    tim->enter("test::DensityFitting G12");
    Ref<Integral> integral = Integral::get_default_integral();

    // construct the G12 matrix
    // and reconstruct it using DensityFitting
    integral->set_basis(bs, bs, bs, bs);
    RefSCMatrix C = twobody_matrix(integral->g12nc<4>(g12params), TwoBodyOper::r12_0_g12);
    C.print("G12 operator");
    Ref<DensityFitting> df = new DensityFitting(integral, bs, bs, fbs);
    df->compute();
    integral->set_basis(fbs,fbs);
    RefSCMatrix kernel = twobody_matrix(integral->g12nc<2>(g12params),
                                        TwoBodyOper::r12_0_g12);

    {
      RefSCMatrix C_df = df->C().t() * kernel * df->C();
      C_df.print("Reconstructed G12 operator (non-robust formula)");
      (C - C_df).print("Reconstruction error (non-robust formula)");
    }

    {
      integral->set_basis(bs, bs, fbs);
      RefSCMatrix G = twobody_matrix(integral->g12nc<3>(g12params),
                                     TwoBodyOper::r12_0_g12);

      RefSCMatrix C_df = G * df->C() + df->C().t() * G.t()
                         - df->C().t() * kernel * df->C();
      C_df.print("Reconstructed G12 operator (robust formula)");
      (C - C_df).print("Reconstruction error (robust formula)");
    }

    tim->exit();
  }

  // now, use TwoBodyMOIntsTransform_ixjy_df to reconstruct a Gaussian geminal matrix
  {
    tim->enter("DensityFitting ixjy transform G12NC");

    Ref<TwoBodyIntDescr> descr = IntDescrFactory::make<4>(integral,
        TwoBodyOperSet::G12NC,
        g12params);
    const unsigned int num_te_types = descr->num_sets();

    Ref<Integral> integral = Integral::get_default_integral();

    Ref<OrbitalSpace> bs_space = new AtomicOrbitalSpace("mu", "AO(OBS)", bs, integral);
    OrbitalSpaceRegistry::instance()->add(make_keyspace_pair(bs_space));
    AOSpaceRegistry::instance()->add(bs, bs_space);
    Ref<MOIntsTransformFactory> factory = new MOIntsTransformFactory(integral);
    factory->set_ints_method(MOIntsTransform::StoreMethod::posix);
    Ref<TwoBodyMOIntsTransform> ixjy_tform =
      new TwoBodyMOIntsTransform_ixjy_df("test", factory,
                                         descr,
                                         bs_space, bs_space,
                                         bs_space, bs_space,
                                         fbs, fbs);

    ixjy_tform->compute();
    Ref<R12IntsAcc> ixjy_acc = ixjy_tform->ints_acc();
    ixjy_acc->activate();

    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    const int nbs = bs->nbasis();
    RefSCDimension bsdim = new SCDimension(nbs);

    // try density fitting for each operator type
    for(int te_type = 0; te_type<num_te_types; ++te_type) {
      const TwoBodyOper::type opertype = descr->intset(te_type);
      integral->set_basis(bs, bs, bs, bs);
      RefSCMatrix C = twobody_matrix(integral->g12nc<4>(g12params), opertype);
      std::ostringstream oss; oss << "G12NC operator te_type " << te_type;
      const std::string operlabel = oss.str();
      C.print(operlabel.c_str());

      RefSCMatrix C_df_ij = localkit->matrix(bsdim, bsdim);
      RefSCMatrix C_df = C.clone(); C_df.assign(0.0);
      RefSCMatrix C_df_localcopy = localkit->matrix(C_df.rowdim(),
                                                    C_df.coldim());
      C_df_localcopy.assign(0.0);
      if (me == 0) {
        for (int i1 = 0; i1 < nbs; ++i1) {
          for (int i2 = 0; i2 < nbs; ++i2) {
            const double* buf = ixjy_acc->retrieve_pair_block(i1, i2, te_type);
            C_df_ij.assign(buf);
            ixjy_acc->release_pair_block(i1, i2, te_type);

            C_df_localcopy.assign_subblock(C_df_ij,
                                           i1 * nbs, (i1 + 1) * nbs - 1,
                                           i2 * nbs, (i2 + 1) * nbs - 1);

          }
        }
      }

      C_df->convert(C_df_localcopy);
      C_df.print((std::string("Reconstructed ") + operlabel).c_str());
      (C - C_df).print("Reconstruction error");

    }

    ixjy_acc->deactivate();

    tim->exit();
  }

#define TEST_TWOBODYTWOCENTERINTOP 0
#if TEST_TWOBODYTWOCENTERINTOP
  {
    integral->set_basis(fbs,fbs);
    RefSCMatrix kernel_ref = twobody_matrix(integral->g12nc<2>(g12params),
                                            TwoBodyOper::eri);
    RefSymmSCMatrix kernel_so = detail::twobody_twocenter_coulomb(fbs,integral);
    integral->set_basis(fbs);
    Ref<PetiteList> pl = integral->petite_list();
    RefSymmSCMatrix kernel_ao = pl->to_AO_basis(kernel_so);
    Ref<SCMatrixKit> kit = SCMatrixKit::default_matrixkit();
    RefSymmSCMatrix kernel = kit->symmmatrix(new SCDimension(fbs->nbasis()));
    kernel->convert(kernel_ao);
    kernel_ref.print("testing TwoBodyTwoCenterIntOp: eri kernel(ref)");
    kernel.print("testing TwoBodyTwoCenterIntOp: eri kernel(intop)");
  }
#endif


  tim->print();
  return 0;
}


RefSCMatrix twobody_matrix(const Ref<TwoBodyInt>& tbint, TwoBodyOper::type tbtype)
{

  Ref<GaussianBasisSet> b1 = tbint->basis1();
  Ref<GaussianBasisSet> b2 = tbint->basis2();
  Ref<GaussianBasisSet> b3 = tbint->basis3();
  Ref<GaussianBasisSet> b4 = tbint->basis4();

  const int nbasis2 = b2->nbasis();
  const int nbasis4 = b4->nbasis();

  RefSCDimension dim12 = new SCDimension(b1->nbasis() * b2->nbasis());
  RefSCDimension dim34 = new SCDimension(b3->nbasis() * b4->nbasis());
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

RefSCMatrix twobody_matrix(const Ref<TwoBodyThreeCenterInt>& tbint, TwoBodyOper::type tbtype)
{

  const Ref<GaussianBasisSet>& b1 = tbint->basis1();
  const Ref<GaussianBasisSet>& b2 = tbint->basis2();
  const Ref<GaussianBasisSet>& b3 = tbint->basis3();

  RefSCDimension rowdim = new SCDimension(b1->nbasis() * b2->nbasis(), "");
  RefSCDimension coldim = new SCDimension(b3->nbasis(), "");
  Ref<SCMatrixKit> kit = SCMatrixKit::default_matrixkit();
  RefSCMatrix result = kit->matrix(rowdim,coldim);

  const int nbasis2 = b2->nbasis();
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

        // compute shell triplet
        tbint->compute_shell(s1, s2, s3);

        // copy buffer into kernel_
        const double* bufptr = buffer;
        for(int f1=0; f1<nf1; ++f1) {
          const int s12offset = (s1offset+f1) * nbasis2 + s2offset;
          for(int f2=0; f2<nf2; ++f2) {
            for(int f3=0; f3<nf3; ++f3, ++bufptr) {
              result.set_element(f2+s12offset, f3+s3offset, *bufptr);
            }
          }
        }

      }
    }
  }

  return result;
}

RefSCMatrix twobody_matrix(const Ref<TwoBodyTwoCenterInt>& tbint, TwoBodyOper::type tbtype)
{
  Ref<GaussianBasisSet> b1 = tbint->basis1();
  Ref<GaussianBasisSet> b2 = tbint->basis2();

  RefSCDimension rowdim = new SCDimension(b1->nbasis(), "");
  RefSCDimension coldim = new SCDimension(b2->nbasis(), "");
  Ref<SCMatrixKit> kit = SCMatrixKit::default_matrixkit();
  RefSCMatrix result = kit->matrix(rowdim,coldim);

  const double* buffer = tbint->buffer(tbtype);
  for (int s1 = 0; s1 < b1->nshell(); ++s1) {
    const int s1offset = b1->shell_to_function(s1);
    const int nf1 = b1->shell(s1).nfunction();
    for (int s2 = 0; s2 < b2->nshell(); ++s2) {
      const int s2offset = b2->shell_to_function(s2);
      const int nf2 = b2->shell(s2).nfunction();

      // compute shell doublet
      tbint->compute_shell(s1, s2);

      // copy buffer into kernel_
      const double* bufptr = buffer;
      for(int f1=0; f1<nf1; ++f1) {
        for(int f2=0; f2<nf2; ++f2, ++bufptr) {
          result.set_element(f1+s1offset, f2+s2offset, *bufptr);
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
