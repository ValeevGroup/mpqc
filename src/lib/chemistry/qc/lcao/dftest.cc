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

#include <mpqc_config.h>
#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <util/misc/consumableresources.h>
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>
#include <util/group/message.h>
#include <util/group/pregtime.h>
#include <chemistry/qc/intv3/intv3.h>
#include <chemistry/qc/lcao/df.h>
#include <chemistry/qc/lcao/transform_ixjy_df.h>
#include <chemistry/qc/lcao/fockbuilder.h>

#include <chemistry/qc/basis/linkage.h>
#ifdef HAVE_LIBINT2
#include <chemistry/qc/libint2/linkage.h>
#endif

using namespace std;
using namespace sc;

RefSCMatrix twobody_matrix(const Ref<TwoBodyInt>& tbint,
                           TwoBodyOper::type tbtype = TwoBodyOper::eri);
RefSCMatrix twobody_matrix(const Ref<TwoBodyThreeCenterInt>& tbint,
                           TwoBodyOper::type tbtype = TwoBodyOper::eri);
RefSCMatrix twobody_matrix(const Ref<TwoBodyTwoCenterInt>& tbint,
                           TwoBodyOper::type tbtype = TwoBodyOper::eri);

const int debug = 1;
void compare_2body_oper(const Ref<GaussianBasisSet>& bs,
                        const RefSCMatrix& C,
                        const RefSCMatrix& Capprox,
                        const std::string& descr,
                        int lmax = -1);

const char* kernel_key = "exp(20)";
const char* dkernel_key = "exp(20)";
const char* pkernel_key = kernel_key;

int main(int argc, char **argv) {

  Ref<MessageGrp> msg = MessageGrp::initial_messagegrp(argc, argv);
  if (msg.null())
    msg = new ProcMessageGrp();
  MessageGrp::set_default_messagegrp(msg);

  Ref<Integral> integral = Integral::initial_integral(argc, argv);
  if (integral.null())
    integral = new IntegralV3();
  Integral::set_default_integral(integral);

  Ref<ConsumableResources> res = ConsumableResources::initial_instance(argc, argv);
  if (res)
    ConsumableResources::set_default_instance(res);
  res->print();

  Ref<RegionTimer> tim = new ParallelRegionTimer(msg, "dftest", 1, 1);

  char *infile = new char[strlen(SRCDIR) + strlen("/dftest.in") + 1];
  sprintf(infile,SRCDIR "/dftest.in");
  if (argc == 2) {
    delete[] infile;
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


    using sc::DensityFitting;
    tim->enter("DensityFitting");

    // construct the Coulomb matrix
    // and reconstruct it using DensityFitting
    integral->set_basis(bs, bs, bs, bs);
    RefSCMatrix C = twobody_matrix(integral->electron_repulsion());
    RefSCMatrix Cref = C;
    if (debug > 1) C.print("Coulomb operator");

    Ref<MOIntsTransformFactory> factory = new MOIntsTransformFactory(integral);
    factory->set_ints_method(MOIntsTransform::StoreMethod::posix);
    Ref<DensityFittingParams> dfparams = new DensityFittingParams(fbs,
                                                                  kernel_key,
                                                                  "bunchkaufman_refine");
    Ref<DensityFitting::MOIntsRuntime> rtime_u = new DensityFitting::MOIntsRuntime(factory);

    Ref<DensityFittingRuntime> dfruntime = new DensityFittingRuntime(new DensityFitting::MOIntsRuntime(factory),
                                                                     dfparams.pointer());
    const DensityFittingInfo* df_info = new DensityFittingInfo(dfparams, dfruntime);
    factory->df_info(df_info);

    Ref<OrbitalSpace> bs_space = new AtomicOrbitalSpace("mu", "AO(OBS)", bs, integral);
    Ref<OrbitalSpace> fbs_space = new AtomicOrbitalSpace("Mu", "AO(FBS)", fbs, integral);
    Ref<OrbitalSpaceRegistry> oreg = factory->orbital_registry();
    Ref<AOSpaceRegistry> aoreg = factory->ao_registry();
    oreg->add(make_keyspace_pair(bs_space));
    oreg->add(make_keyspace_pair(fbs_space));
    aoreg->add(bs, bs_space);
    aoreg->add(fbs, fbs_space);

    Ref<DensityFitting> df = new DensityFitting(rtime_u, dfparams->kernel_key(), dfparams->solver(),
                                                bs_space, bs_space, fbs);
    df->compute();

    RefSCMatrix C_df = C.clone(); C_df.assign(0.0);
    // assemble reconstructed matrix
    {
      Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
      Ref<DistArray4> C = df->C();
      C->activate();

      std::string cC_key = ParsedTwoBodyThreeCenterIntKey::key(bs_space->id(),
                                                               aoreg->value(fbs)->id(),
                                                               bs_space->id(),
                                                               "ERI", "");
      Ref<TwoBodyThreeCenterMOIntsTransform> tform3 = dfruntime->moints_runtime()->runtime_3c()->get(cC_key);
      tform3->compute();
      Ref<DistArray4> cC = tform3->ints_acc();
      cC->activate();

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
          const double* C_jR_buf = C->retrieve_pair_block(0, i1, 0);
          C_jR.assign(C_jR_buf);
          C->release_pair_block(0, i1, TwoBodyOper::eri);
          RefSCMatrix C_jR_t = C_jR.t();
          for (int i2 = 0; i2 < nbs; ++i2) {
            const double* cC_jR_buf = cC->retrieve_pair_block(0, i2, 0);
            cC_jR.assign(cC_jR_buf);
            cC->release_pair_block(0, i2, 0);

            cCC.assign(0.0);
            cCC.accumulate_product(cC_jR, C_jR_t);
            C_df_localcopy.assign_subblock(cCC, i2 * nbs, (i2 + 1) * nbs - 1,
                                           i1 * nbs, (i1 + 1) * nbs - 1);

          }
        }
      }
      C_df->convert(C_df_localcopy);
    }

    std::ostringstream oss; oss << "kernel=" << kernel_key << ",std";
    compare_2body_oper(bs, C, C_df, oss.str());

    tim->exit();

  // now, do the same by using TwoBodyMOIntsTransform_ixjy_df
  if (0) {
    tim->enter("DensityFitting ixjy transform");

    Ref<Integral> integral = Integral::get_default_integral();
    // construct the Coulomb matrix
    // and reconstruct it using DensityFitting
    integral->set_basis(bs, bs, bs, bs);
    RefSCMatrix C = twobody_matrix(integral->electron_repulsion());

    Ref<TwoBodyMOIntsTransform> ixjy_tform =
      new TwoBodyMOIntsTransform_ixjy_df("test", df_info,
                                         new TwoBodyIntDescrERI(integral),
                                         bs_space, bs_space,
                                         bs_space, bs_space);

    ixjy_tform->compute();
    Ref<DistArray4> ixjy_acc = ixjy_tform->ints_distarray4();
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
          const double* buf = ixjy_acc->retrieve_pair_block(i1, i2, 0);
          C_df_ij.assign(buf);
          ixjy_acc->release_pair_block(i1, i2, 0);

          C_df_localcopy.assign_subblock(C_df_ij,
                                         i1 * nbs, (i1 + 1) * nbs - 1,
                                         i2 * nbs, (i2 + 1) * nbs - 1);

        }
      }
    }
    ixjy_acc->deactivate();

    C_df->convert(C_df_localcopy);

    std::ostringstream oss; oss << "kernel=" << kernel_key << ",robust";
    compare_2body_oper(bs, C, C_df, oss.str());
    tim->exit();
  }

  ///////////////////////////
  // experiments with fitting
  ///////////////////////////
  if (1) {

    //
    // check how close to exact are the charges (maybe other multipoles also?) -> this essentially checks the long range of the potential
    //
    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    Ref<DistArray4> C = df->C();
    C->activate();

    // need to use the unit basis .. check if it is already in the registry
    Ref<OrbitalSpace> unit_space;
    if (not oreg->key_exists("1")) {
      unit_space = new AtomicOrbitalSpace("1", "AO(1)", GaussianBasisSet::unit(), integral);
      oreg->add(make_keyspace_pair(unit_space));
      aoreg->add(unit_space->basis(), unit_space);
    }

    const int nbs = bs->nbasis();
    const int nfbs = fbs->nbasis();
    RefSCDimension bsdim = new SCDimension(nbs);
    RefSCDimension fbsdim = new SCDimension(nfbs);
    RefSCMatrix S = sc::overlap(*bs_space, *bs_space, localkit);
    RefSCMatrix S_df = sc::overlap(*fbs_space, *unit_space, localkit);
    RefSCMatrix C_jX = localkit->matrix(bsdim, fbs_space->dim());
    RefSCMatrix N_ij = localkit->matrix(bs_space->dim(), bs_space->dim());
    N_ij.assign(0.0);

    if (me == 0) {
      for (int i1 = 0; i1 < nbs; ++i1) {
        const double* C_jX_buf = C->retrieve_pair_block(0, i1, 0);
        C_jX.assign(C_jX_buf);
        C->release_pair_block(0, i1, 0);

        RefSCMatrix CS = C_jX * S_df;

        for (int i2 = 0; i2 < nbs; ++i2) {

          N_ij.set_element(i1, i2, CS.get_element(i2,0));
        }
      }
    }

    if (debug > 1) N_ij.print("charges of the fitted densities");
    //S_df.print("overlap in DF basis");
    //SCFormIO::setverbose(ExEnv::out0(),1);
    //fbs->print();

    // charges of the exact densities are the 2-center overlaps
    if (debug > 1) S.print("charges of the exact product densities");

    if (debug > 0) (N_ij - S).print("fitting errors in the \"charges\" of the product densities");

    {
      RefSCMatrix error = N_ij - S;

      Ref<SCElementMaxAbs> maxabs_op = new SCElementMaxAbs;
      error->element_op(maxabs_op);
      ExEnv::out0() << indent << "Charge error maxabs = " << maxabs_op->result() << std::endl;

      Ref<SCElementKNorm> norm2_op = new SCElementKNorm(2);
      error->element_op(norm2_op);
      ExEnv::out0() << indent << "Charge error norm2 = " << norm2_op->result() << std::endl;
    }

    //
    // check how close to exact (in overlap sense) are the densities
    //

    // dRho_mu,nu = ||Rho_mu,nu - ~Rho_mu,nu||_2 = ||Rho_mu,nu||_2 + ~Rho_mu,nu||_2 - 2 <Rho_mu,nu | ~Rho_mu,nu>
    //            = <mu,nu|mu,nu> + C^X_mu,nu <X|Y> C^Y_mu,nu - 2 <mu,nu|X> C^X_mu,nu
    {

      integral->set_basis(bs, bs, bs, bs);
      RefSCMatrix S_oooo = twobody_matrix(integral->delta_function<4>(), TwoBodyOper::delta);
      RefSCMatrix S_dd = sc::overlap(*fbs_space, *fbs_space, localkit);

      std::string S_ood_key = ParsedTwoBodyThreeCenterIntKey::key(bs_space->id(),
                                                                  fbs_space->id(),
                                                                  bs_space->id(),
                                                                  "Delta", "");
      Ref<TwoBodyThreeCenterMOIntsTransform> tform3 = dfruntime->moints_runtime()->runtime_3c()->get(S_ood_key);
      tform3->compute();
      Ref<DistArray4> S_ood = tform3->ints_acc();
      S_ood->activate();

      RefSCMatrix C_jX = localkit->matrix(bsdim, fbs_space->dim());
      RefSCMatrix S_ood_jX = localkit->matrix(bsdim, fbs_space->dim());
      RefSCMatrix S_ij = localkit->matrix(bs_space->dim(), bs_space->dim());
      S_ij.assign(0.0);

      if (me == 0) {
        for (int i1 = 0; i1 < nbs; ++i1) {
          const double* C_jX_buf = C->retrieve_pair_block(0, i1, 0);
          C_jX.assign(C_jX_buf);
          C->release_pair_block(0, i1, 0);
          const double* S_ood_jX_buf = S_ood->retrieve_pair_block(0, i1, 0);
          S_ood_jX.assign(S_ood_jX_buf);
          S_ood->release_pair_block(0, i1, 0);

          RefSCMatrix CSC = C_jX * S_dd * C_jX.t();
          RefSCMatrix CS = C_jX * S_ood_jX.t();

          for (int i2 = 0; i2 < nbs; ++i2) {
            const size_t i1i2 = i1 * nbs + i2;

            S_ij.set_element(i1, i2, sqrt(S_oooo.get_element(i1i2, i1i2) +  CSC.get_element(i2,i2) - 2 * CS.get_element(i2,i2)));
            //S_ij.set_element(i1, i2, S_oooo.get_element(i1i2, i1i2) - CS.get_element(i2,i2));
          }
        }
      }

      if (debug > 0) S_ij.print("fitting errors in the \"norms\" of the product densities");
      //S_oooo.print("4-center overlap ... the diagonal elements are the norms of the product densities");
      //RefSCMatrix S_oo = sc::overlap(*bs_space, *bs_space, localkit);
      //S_oo.print("2-center overlap");
      //SCFormIO::setverbose(ExEnv::out0(),1);
      //bs->print();

      if (0) {
        std::string S_ooo_key = ParsedTwoBodyThreeCenterIntKey::key(bs_space->id(),
                                                                    bs_space->id(),
                                                                    bs_space->id(),
                                                                    "Delta", "");
        Ref<TwoBodyThreeCenterMOIntsTransform> tform3 = dfruntime->moints_runtime()->runtime_3c()->get(S_ooo_key);
        tform3->compute();
        Ref<DistArray4> S_ooo = tform3->ints_acc();
        S_ooo->activate();
        RefSCMatrix S_ooo_jX = localkit->matrix(bsdim, bsdim);
        if (me == 0) {
          ExEnv::out0() << "3-center overlap" << std::endl;
          for (int i1 = 0; i1 < nbs; ++i1) {
            const double* S_ooo_jX_buf = S_ooo->retrieve_pair_block(0, i1, 0);
            S_ooo_jX.assign(S_ooo_jX_buf);
            S_ooo->release_pair_block(0, i1, 0);

            ExEnv::out0() << "block # " << i1 << std::endl;
            S_ooo_jX.print("");
          }
        }

      }

      //
      // combine potential fitting and density fitting
      //
      if (1) {
        Ref<DensityFitting> df = new DensityFitting(rtime_u, dkernel_key, DensityFitting::SolveMethod_RefinedBunchKaufman,
                                                    bs_space, bs_space, fbs);
        df->compute();
        Ref<DistArray4> P = df->C(); P->activate();
        Ref<DensityFitting> pf = new DensityFitting(rtime_u, pkernel_key, DensityFitting::SolveMethod_RefinedBunchKaufman,
                                                    bs_space, bs_space, fbs);
        pf->compute();
        Ref<DistArray4> VP = pf->C(); VP->activate();

        std::string V_ff_key = ParsedTwoBodyTwoCenterIntKey::key(fbs_space->id(), fbs_space->id(),
                                                                 "ERI", "");
        RefSCMatrix V_ff = dfruntime->moints_runtime()->runtime_2c()->get(V_ff_key);
        RefSCMatrix V_ff_local = localkit->matrix(fbsdim, fbsdim);
        //V_ff.print("V_ff");
        V_ff_local->convert(V_ff);
        //V_ff_local.print("V_ff_local");

        // assemble reconstructed coulomb matrix
        C_df = Cref.clone(); C_df.assign(0.0);
        {

          RefSCMatrix df_jR = localkit->matrix(bsdim, fbsdim);
          RefSCMatrix pf_jR = localkit->matrix(bsdim, fbsdim);
          if (me == 0) {
            for (int i1 = 0; i1 < nbs; ++i1) {
              const double* df_jR_buf = P->retrieve_pair_block(0, i1, 0);
              df_jR.assign(df_jR_buf);
              P->release_pair_block(0, i1, 0);

              for (int i2 = 0; i2 < nbs; ++i2) {
                const double* pf_jR_buf = VP->retrieve_pair_block(0, i2, 0);
                pf_jR.assign(pf_jR_buf);
                VP->release_pair_block(0, i2, 0);

                RefSCMatrix C_j = df_jR * V_ff_local * pf_jR.t();
                for (int j1 = 0; j1 < nbs; ++j1) {
                  const int i1j1 = i1 * nbs + j1;
                  for (int j2 = 0; j2 < nbs; ++j2) {
                    const int i2j2 = i2 * nbs + j2;

                    C_df.set_element(i1j1, i2j2, C_j.get_element(j1, j2));
                  }
                }
              }

            }
          }
        }

        // symmetrize
        RefSCMatrix C_df_symm = C_df.copy();
        C_df_symm.accumulate(C_df.t());
        C_df_symm.scale(0.5);
        C_df = C_df_symm;

        std::ostringstream oss; oss << "kernel(pot)=" << pkernel_key << ",kernel(dens)=" << dkernel_key << ",std";
        compare_2body_oper(bs, Cref, C_df, oss.str());

      }

      //
      // check the errors in intermolecular integrals assuming purely local fitting
      //
      if (1) {

        // params
        const char* bsname = "Def2-SVP/JK";
        const char* fbsname = "Def2-QZVPP/JK";
        //const char* fbsname = "cc-pVDZ-RI";
        const double xdisp = 0.0;
        const double zdisp = 2.0;
        const double R = 1.0;

        Ref<Molecule> mol[3]; // fragments 0 and 1, and their union
        {
          const double Ro2sqrt2 = R / (2 * sqrt(2.0));
          for(int i=0; i<3; ++i) mol[i] = new Molecule;
          mol[0]->add_atom(static_cast<int>(1.0), 0.0, 0.0,  R/2);
          mol[0]->add_atom(static_cast<int>(1.0), 0.0, 0.0, -R/2);
          mol[1]->add_atom(static_cast<int>(1.0), xdisp + Ro2sqrt2, 0.0,  zdisp + Ro2sqrt2);
          mol[1]->add_atom(static_cast<int>(1.0), xdisp - Ro2sqrt2, 0.0,  zdisp - Ro2sqrt2);
          //mol[1] = mol[0];
          for(int a=0; a<2; ++a) mol[2]->add_atom(mol[0]->Z(0), mol[0]->r(a, 0), mol[0]->r(a, 1), mol[0]->r(a, 2));
          for(int a=0; a<2; ++a) mol[2]->add_atom(mol[1]->Z(0), mol[1]->r(a, 0), mol[1]->r(a, 1), mol[1]->r(a, 2));
        }
        Ref<GaussianBasisSet> bs[2], fbs[3];
        Ref<AtomicOrbitalSpace> bs_space[2], fbs_space[3];
        for(int i=0; i<2; ++i) {
          Ref<AssignedKeyVal> akv = new AssignedKeyVal;
          akv->assign("molecule", mol[i].pointer());
          akv->assign("name", bsname);
          bs[i] = new GaussianBasisSet(Ref<KeyVal>(akv));

          std::ostringstream oss; oss << i;
          std::string bs_id = "mu" + oss.str();
          std::string bs_descr = "AO(OBS" + oss.str() + ")";
          bs_space[i] = new AtomicOrbitalSpace(bs_id, bs_descr, bs[i], integral);
          oreg->add(make_keyspace_pair(bs_space[i]));
          aoreg->add(bs[i], bs_space[i]);
        }
        for(int i=0; i<3; ++i) {
          Ref<AssignedKeyVal> akv = new AssignedKeyVal;
          akv->assign("molecule", mol[i].pointer());
          akv->assign("name", fbsname);
          fbs[i] = new GaussianBasisSet(Ref<KeyVal>(akv));

          std::ostringstream oss; oss << i;
          std::string fbs_id = "Mu" + oss.str();
          std::string fbs_descr = "AO(DFBS" + oss.str() + ")";
          fbs_space[i] = new AtomicOrbitalSpace(fbs_id, fbs_descr, fbs[i], integral);
          oreg->add(make_keyspace_pair(fbs_space[i]));
          aoreg->add(fbs[i], fbs_space[i]);
        }

        for(int fitmethod=0; fitmethod<2; ++fitmethod) { // do 2 kinds of fits -- local and nonlocal

          ExEnv::out0() << "Dimer fitting " << (fitmethod == 0 ? "locally" : "globally") << std::endl;

          // redirect to the proper fitting basis
          Ref<AtomicOrbitalSpace> fitspace[2];
          for(int i=0; i<2; ++i)
            fitspace[i] = fitmethod == 0 ? fbs_space[i] : fbs_space[2];

          Ref<DensityFitting> df[2], pf[2];
          for(int i=0; i<2; ++i) {
            df[i] = new DensityFitting(rtime_u, dkernel_key, DensityFitting::SolveMethod_RefinedBunchKaufman,
                                       bs_space[i], bs_space[i], fitspace[i]->basis());
            df[i]->compute();
            pf[i] = new DensityFitting(rtime_u, pkernel_key, DensityFitting::SolveMethod_RefinedBunchKaufman,
                                       bs_space[i], bs_space[i], fitspace[i]->basis());
            pf[i]->compute();
          }

          std::string V_1_2_ff_key = ParsedTwoBodyTwoCenterIntKey::key(fitspace[0]->id(), fitspace[1]->id(),
                                                                       "ERI", "");
          RefSCMatrix V_1_2_ff = dfruntime->moints_runtime()->runtime_2c()->get(V_1_2_ff_key);
          RefSCMatrix V_1_2_ff_local = localkit->matrix(fitspace[0]->dim(), fitspace[1]->dim());
          V_1_2_ff_local->convert(V_1_2_ff);

          Ref<DistArray4> P0 = df[0]->C(); P0->activate();
          Ref<DistArray4> VP0 = pf[0]->C(); VP0->activate();
          Ref<DistArray4> P1 = df[1]->C(); P1->activate();
          Ref<DistArray4> VP1 = pf[1]->C(); VP1->activate();

          // compute the exact coulomb matrix
          integral->set_basis(bs[0], bs[0], bs[1], bs[1]);
          RefSCMatrix C0011 = twobody_matrix(integral->electron_repulsion());

          // assemble reconstructed coulomb matrix
          RefSCMatrix C0011_df = C0011.clone(); C0011_df.assign(0.0);
          {
            enum type {Density, Potential};
            type bra[] = {Density, Potential};
            type ket[] = {Potential, Density};

            for(int k=0; k<2; ++k) {

              const int nbra = bs_space[0]->rank();
              const int nket = bs_space[1]->rank();

              const Ref<DistArray4>& bra_data = bra[k] == Density ? P0 : VP0;
              const Ref<DistArray4>& ket_data = ket[k] == Density ? P1 : VP1;

              RefSCMatrix bra_jR = localkit->matrix(bs_space[0]->dim(), fitspace[0]->dim());
              RefSCMatrix ket_jR = localkit->matrix(bs_space[1]->dim(), fitspace[1]->dim());

              if (me == 0) {
                for (int i1 = 0; i1 < nbra; ++i1) {
                  const double* bra_jR_buf = bra_data->retrieve_pair_block(0, i1, 0);
                  bra_jR.assign(bra_jR_buf);
                  bra_data->release_pair_block(0, i1, 0);

                  for (int i2 = 0; i2 < nket; ++i2) {
                    const double* ket_jR_buf = ket_data->retrieve_pair_block(0, i2, 0);
                    ket_jR.assign(ket_jR_buf);
                    ket_data->release_pair_block(0, i2, 0);

                    RefSCMatrix C_j = bra_jR * V_1_2_ff_local * ket_jR.t();
                    for (int j1 = 0; j1 < nbra; ++j1) {
                      const int i1j1 = i1 * nbra + j1;
                      for (int j2 = 0; j2 < nket; ++j2) {
                        const int i2j2 = i2 * nket + j2;

                        C0011_df.accumulate_element(i1j1, i2j2, C_j.get_element(j1, j2));
                      }
                    }
                  }
                }
              }
            }

            C0011_df.scale(0.5);

          }

          std::ostringstream oss; oss << "kernel(pot)=" << pkernel_key << ",kernel(dens)=" << dkernel_key << ",std";
          compare_2body_oper(bs[0], C0011, C0011_df, oss.str());
        }
      }

    }

  }


  typedef IntParamsG12::PrimitiveGeminal PrimitiveGeminal;
  std::vector<PrimitiveGeminal> geminal(1, std::make_pair(2.0, 1.0));
  Ref<IntParamsG12> g12params = new IntParamsG12(geminal);


  // now, use TwoBodyMOIntsTransform_ixjy_df to reconstruct a Gaussian geminal matrix
  try { if (0) {
    tim->enter("DensityFitting ixjy transform G12NC");

    Ref<TwoBodyIntDescr> descr = IntDescrFactory::make<4>(integral,
        TwoBodyOperSet::G12NC,
        g12params);
    const unsigned int num_te_types = descr->num_sets();

    Ref<Integral> integral = Integral::get_default_integral();

    Ref<MOIntsTransformFactory> factory = new MOIntsTransformFactory(integral);
    factory->set_ints_method(MOIntsTransform::StoreMethod::posix);
    Ref<OrbitalSpace> bs_space = new AtomicOrbitalSpace("mu", "AO(OBS)", bs, integral);
    factory->orbital_registry()->add(make_keyspace_pair(bs_space));
    factory->ao_registry()->add(bs, bs_space);
    Ref<TwoBodyMOIntsTransform> ixjy_tform =
      new TwoBodyMOIntsTransform_ixjy_df("test", df_info,
                                         descr,
                                         bs_space, bs_space,
                                         bs_space, bs_space);

    ixjy_tform->compute();
    Ref<DistArray4> ixjy_acc = ixjy_tform->ints_distarray4();
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
  }}
  catch (...) {
    ExEnv::out0() << "Failed to test GTG fitting" << std::endl;
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
    RefSymmSCMatrix kernel_key = kit->symmmatrix(new SCDimension(fbs->nbasis()));
    kernel_key->convert(kernel_ao);

    kernel_ref.print("testing TwoBodyTwoCenterIntOp: eri kernel(ref)");
    kernel_key.print("testing TwoBodyTwoCenterIntOp: eri kernel(intop)");

    // try the new runtime
    Ref<MOIntsTransformFactory> factory = new MOIntsTransformFactory(integral);
    factory->set_ints_method(MOIntsTransform::StoreMethod::posix);
    Ref<TwoBodyTwoCenterMOIntsRuntime> rtime2 = new TwoBodyTwoCenterMOIntsRuntime(factory);
    const std::string key = ParsedTwoBodyTwoCenterIntKey::key("Mu","Mu","ERI","");
    RefSCMatrix kernel_from_rtime = rtime2->get(key);
    kernel_from_rtime.print("testing TwoBodyTwoCenterIntOp: eri kernel(factory)");
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

void compare_2body_oper(const Ref<GaussianBasisSet>& bs,
                        const RefSCMatrix& C,
                        const RefSCMatrix& Capprox,
                        const std::string& descr,
                        const int lmax) {

  std::string label = std::string("Coulomb operator: (") + descr + ")";
  if (debug > 1) Capprox.print(label.c_str());

  RefSCMatrix error = C - Capprox;
  label = std::string("Coulomb operator error: (") + descr + ")";
  if (debug > 1) error.print(label.c_str());

  Ref<SCElementMaxAbs> maxabs_op = new SCElementMaxAbs;
  error->element_op(maxabs_op);
  ExEnv::out0() << indent << label << ": maxabs = " << maxabs_op->result() << std::endl;

  Ref<SCElementKNorm> norm2_op = new SCElementKNorm(2);
  error->element_op(norm2_op);
#if 0
  // correct for some elements in C AND in Capprox zero
  const int nrow = C.nrow();
  const int ncol = C.ncol();
  size_t num_nonzero_elem = 0.0;
  for(int r=0; r<nrow; ++r)
    for(int c=0; c<ncol; ++c)
      if (C(r,c) > 1e-9)
        ++num_nonzero_elem;
#endif
  ExEnv::out0() << indent << label << ": norm2 = " << norm2_op->result() << std::endl;

  // 2-norm of the geometric mean error
  RefSCMatrix gerror = error.clone();
  const int nrow = gerror.nrow();
  const int ncol = gerror.ncol();
  for(int r=0; r<nrow; ++r) {
    for(int c=0; c<ncol; ++c) {
      gerror.set_element(r, c, sqrt(fabs(error(r,c) * C(r,c))));
    }
  }
  Ref<SCElementKNorm> gnorm2_op = new SCElementKNorm(2);
  gerror->element_op(gnorm2_op);
  ExEnv::out0() << indent << label << ": geo-norm2 = " << gnorm2_op->result() << std::endl;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
