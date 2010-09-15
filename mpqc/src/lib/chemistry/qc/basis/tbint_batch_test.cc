//
// tbint_batch_test.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Ida Nielsen <ida@kemi.aau.dk>
// Maintainer: LPS
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

#ifdef HAVE_CONFIG_H
#include <scconfig.h>
#endif

#include <string.h>

#include <sys/stat.h>
#include <unistd.h>
#include <new>

#include <util/class/scexception.h>
#include <util/keyval/keyval.h>
#include <util/group/message.h>
#include <util/group/pregtime.h>
#include <util/misc/bug.h>
#include <util/misc/formio.h>
#include <util/misc/autovec.h>
#include <util/state/stateio.h>
#include <util/state/state_bin.h>

#include <math/scmat/repl.h>

#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/energy.h>

#include <chemistry/qc/scf/clhf.h>
#include <chemistry/qc/scf/hsoshf.h>

using namespace std;
using namespace sc;

class MP2: public Wavefunction {
    Ref<OneBodyWavefunction> ref_mp2_wfn_;
    double compute_mp2_energy();

  public:
    MP2(const Ref<KeyVal> &);
    MP2(StateIn &);
    void save_data_state(StateOut &);
    void compute(void);
    void obsolete(void);
    int nelectron(void);
    RefSymmSCMatrix density(void);
    int spin_polarized(void);
    int value_implemented(void) const;
};

static ClassDesc MP2_cd(typeid(MP2), "MP2", 1, "public Wavefunction", 0,
                        create<MP2> , create<MP2> );

MP2::MP2(const Ref<KeyVal> &keyval) :
  Wavefunction(keyval) {
  ref_mp2_wfn_ << keyval->describedclassvalue("reference");
  if (ref_mp2_wfn_.null()) {
    throw InputError("require a OneBodyWavefunction object", __FILE__,
                     __LINE__, "reference", 0, class_desc());
  }
}

MP2::MP2(StateIn &statein) :
  Wavefunction(statein) {
  ref_mp2_wfn_ << SavableState::restore_state(statein);
}

void MP2::save_data_state(StateOut &stateout) {
  Wavefunction::save_data_state(stateout);

  SavableState::save_state(ref_mp2_wfn_.pointer(), stateout);
}

void MP2::compute(void) {
  if (gradient_needed()) {
    throw FeatureNotImplemented("no gradients yet", __FILE__, __LINE__,
                                class_desc());
  }

  double extra_hf_acc = 10.;
  ref_mp2_wfn_->set_desired_value_accuracy(desired_value_accuracy()
      / extra_hf_acc);
  double refenergy = ref_mp2_wfn_->energy();
  double mp2energy = compute_mp2_energy();

  ExEnv::out0() << indent << "MP2 Energy = " << mp2energy << endl;

  set_value(refenergy + mp2energy);
  set_actual_value_accuracy(ref_mp2_wfn_->actual_value_accuracy()
      * extra_hf_acc);
}

void MP2::obsolete(void) {
  Wavefunction::obsolete();
  ref_mp2_wfn_->obsolete();
}

int MP2::nelectron(void) {
  return ref_mp2_wfn_->nelectron();
}

RefSymmSCMatrix MP2::density(void) {
  throw FeatureNotImplemented("no density yet", __FILE__, __LINE__,
                              class_desc());
  return 0;
}

int MP2::spin_polarized(void) {
  return 0;
}

int MP2::value_implemented(void) const {
  return 1;
}

double MP2::compute_mp2_energy() {
  if (molecule()->point_group()->char_table().order() != 1) {
    throw FeatureNotImplemented("C1 symmetry only", __FILE__, __LINE__,
                                class_desc());
  }

  RefSCMatrix vec = ref_mp2_wfn_->eigenvectors();

  int nao = vec.nrow();
  int nmo = vec.ncol();
  int nocc = ref_mp2_wfn_->nelectron() / 2;
  int nvir = nmo - nocc;

  auto_vec<double> cvec_av(new double[vec.nrow() * vec.ncol()]);
  double *cvec = cvec_av.get();
  vec->convert(cvec);

  auto_vec<double> pqrs_av(new double[nao * nao * nao * nao]);
  double *pqrs = pqrs_av.get();
  for (int n = 0; n < nao * nao * nao * nao; n++)
    pqrs[n] = 0.0;

  Ref<TwoBodyInt> twoint = integral()->electron_repulsion();
  const double *buffer = twoint->buffer();

  Ref<GaussianBasisSet> basis = this->basis();

  int nshell = basis->nshell();
  for (int P = 0; P < nshell; P++) {
    int nump = basis->shell(P).nfunction();

    for (int Q = 0; Q < nshell; Q++) {
      int numq = basis->shell(Q).nfunction();

      for (int R = 0; R < nshell; R++) {
        int numr = basis->shell(R).nfunction();

        for (int S = 0; S < nshell; S++) {
          int nums = basis->shell(S).nfunction();

          twoint->compute_shell(P, Q, R, S);

          int index = 0;
          for (int p = 0; p < nump; p++) {
            int op = basis->shell_to_function(P) + p;

            for (int q = 0; q < numq; q++) {
              int oq = basis->shell_to_function(Q) + q;

              for (int r = 0; r < numr; r++) {
                int oor = basis->shell_to_function(R) + r;

                for (int s = 0; s < nums; s++, index++) {
                  int os = basis->shell_to_function(S) + s;

                  int ipqrs = (((op * nao + oq) * nao + oor) * nao + os);

                  pqrs[ipqrs] = buffer[index];

                }
              }
            }
          }

        }
      }
    }
  }

  twoint = 0;

  auto_vec<double> ijkl_av(new double[nmo * nmo * nmo * nmo]);
  double *ijkl = ijkl_av.get();

  int idx = 0;
  for (int i = 0; i < nmo; i++) {
    for (int j = 0; j < nmo; j++) {
      for (int k = 0; k < nmo; k++) {
        for (int l = 0; l < nmo; l++, idx++) {

          ijkl[idx] = 0.0;

          int index = 0;
          for (int p = 0; p < nao; p++) {
            for (int q = 0; q < nao; q++) {
              for (int r = 0; r < nao; r++) {
                for (int s = 0; s < nao; s++, index++) {

                  ijkl[idx] += cvec[p * nmo + i] * cvec[q * nmo + j] * cvec[r
                      * nmo + k] * cvec[s * nmo + l] * pqrs[index];

                }
              }
            }
          }

        }
      }
    }
  }

  pqrs_av.release();
  pqrs = 0;
  cvec_av.release();
  cvec = 0;

  auto_vec<double> evals_av(new double[nmo]);
  double *evals = evals_av.get();
  ref_mp2_wfn_->eigenvalues()->convert(evals);

  double energy = 0.0;
  for (int i = 0; i < nocc; i++) {
    for (int j = 0; j < nocc; j++) {
      for (int a = nocc; a < nmo; a++) {
        for (int b = nocc; b < nmo; b++) {

          int iajb = (((i * nmo + a) * nmo + j) * nmo + b);
          int ibja = (((i * nmo + b) * nmo + j) * nmo + a);

          energy += (2 * ijkl[iajb] - ijkl[ibja]) * ijkl[iajb] / (evals[i]
              + evals[j] - evals[a] - evals[b]);

        }
      }
    }
  }

  ijkl_av.release();
  ijkl = 0;
  evals_av.release();
  evals = 0;

  return energy;
}

// Force linkages:
static ForceLink<CLHF> fl0a;
static ForceLink<MP2> fl0e;
static ForceLink<ReplSCMatrixKit> fl6;
static ForceLink<ProcMessageGrp> fl9;

Ref<MessageGrp> grp;

static Ref<MessageGrp> init_mp(const Ref<KeyVal>& keyval, int &argc,
                               char **&argv) {
  grp << keyval->describedclassvalue("message");

  if (grp.null())
    grp = MessageGrp::initial_messagegrp(argc, argv);

  if (grp.null()) {
    grp << keyval->describedclassvalue("messagegrp");
  }

  if (grp.null())
    grp = MessageGrp::get_default_messagegrp();

  if (grp.null()) {
    std::cerr << indent << "Couldn't initialize MessageGrp\n";
    abort();
  }

  MessageGrp::set_default_messagegrp(grp);

  Ref<Debugger> debugger;
  debugger << keyval->describedclassvalue(":debug");
  // Let the debugger know the name of the executable and the node
  if (debugger.nonnull()) {
    debugger->set_exec("test");
    debugger->set_prefix(grp->me());
    debugger->debug("curt is a hog");
  }

  RegionTimer::set_default_regiontimer(new ParallelRegionTimer(grp, "test",
                                                               1, 0));
  
  SCFormIO::set_printnode(0);
  SCFormIO::init_mp(grp->me());
  //SCFormIO::set_debug(1);

  SCFormIO::setindent(ExEnv::outn(), 2);
  SCFormIO::setindent(cerr, 2);
  
  return grp;
}

main(int argc, char**argv)
{
  const char *input = (argc > 1)? argv[1] : SRCDIR "/tbint_batch_test.in";
  const char *keyword = (argc > 2)? argv[2] : "mole";
  const char *optkeyword = (argc > 3)? argv[3] : "opt";

  // open keyval input
  Ref<KeyVal> rpkv(new ParsedKeyVal(input));

  init_mp(rpkv, argc, argv);

  Timer tim;
  tim.enter("input");

  int do_gradient = rpkv->booleanvalue("gradient");

  if (rpkv->exists("matrixkit")) {
    Ref<SCMatrixKit> kit; kit << rpkv->describedclassvalue("matrixkit");
    SCMatrixKit::set_default_matrixkit(kit);
  }

  struct stat sb;
  Ref<MolecularEnergy> mole;
  mole << rpkv->describedclassvalue("mole");

  tim.exit("input");

  if (mole.nonnull()) {
    ExEnv::out0() << indent << "energy: " << mole->energy() << endl;
    if (mole->value_implemented()) {
      ExEnv::out0() << indent
      << scprintf("value of mole is %20.15f\n\n", mole->energy());
    }

    mole->print(ExEnv::out0());
  }

  StateOutBin so("tbint_batch_test.wfn");
  SavableState::save_state(mole.pointer(),so);

  tim.print(ExEnv::out0());

  grp = 0;
  RegionTimer::set_default_regiontimer(0);
  MessageGrp::set_default_messagegrp(0);
  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
