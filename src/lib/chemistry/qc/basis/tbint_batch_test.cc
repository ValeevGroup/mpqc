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
#include <mpqc_config.h>
#endif

#include <string.h>

#include <sys/stat.h>
#include <unistd.h>
#include <new>

#include <util/misc/scexception.h>
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

#include <chemistry/qc/basis/tbint_batch.h>

using namespace std;
using namespace sc;

// jtf: Example class to deliver integrals in batches

template <unsigned int N, typename T> struct tuple {
    T data[N];
    T operator[](size_t i) const {
      return data[i];
    }
    T& operator[](size_t i) {
      return data[i];
    }
};

typedef tuple<4, unsigned int> ShellDesc;

class Int_Batch {

  public:

    // Initialize with a buffer size s and a basis set b on which to compute ints
    Int_Batch(int s, const Ref<TwoBodyInt>& tbint);
    ~Int_Batch();
    bool next();

    const std::vector<ShellDesc>& shells_in_batch() {
      return shells_in_batch_;
    }
    const std::vector<ShellDesc>& pqrs_start() {
      return pqrs_start_;
    }
    const std::vector<ShellDesc>& pqrs_len() {
      return pqrs_len_;
    }

    const double* buffer() const {
      return &(buffer_[0]);
    }


  private:
    Ref<TwoBodyInt> twoint_;
    std::vector<double> buffer_;
    ShellDesc s_; // Shell Quartet state - persists between calls to next()
    std::vector<ShellDesc> shells_in_batch_;
    std::vector<ShellDesc> pqrs_start_;
    std::vector<ShellDesc> pqrs_len_;


};

Int_Batch::Int_Batch(int s, const Ref<TwoBodyInt>& t) :
  twoint_(t) {
  buffer_.reserve(s);
  s_[0] = 0;
  s_[1] = 0;
  s_[2] = 0;
  s_[3] = 0;
}



Int_Batch::~Int_Batch() {

}

bool Int_Batch::next() {
  const double *in_buf = twoint_->buffer();
  const Ref<GaussianBasisSet> basis = twoint_->basis();

  buffer_.clear();
  shells_in_batch_.clear();
  pqrs_start_.clear();
  pqrs_len_.clear();

  const int nshell = basis->nshell();

  ShellDesc s = s_;
  ShellDesc p_i;
  ShellDesc p_l;
  for (; s[0] < nshell; ++s[0], s[1]=0) {
    p_l[0] = basis->shell(s[0]).nfunction();
    p_i[0] = basis->shell_to_function(s[0]);
    for (; s[1] < nshell; ++s[1], s[2]=0) {
      p_l[1] = basis->shell(s[1]).nfunction();
      p_i[1] = basis->shell_to_function(s[1]);
      for (; s[2] < nshell; ++s[2], s[3]=0) {
        p_l[2] = basis->shell(s[2]).nfunction();
        p_i[2] = basis->shell_to_function(s[2]);
        for (; s[3] < nshell; ++s[3]) {
          p_l[3] = basis->shell(s[3]).nfunction();
          p_i[3] = basis->shell_to_function(s[3]);

          const int n = p_l[0]*p_l[1]*p_l[2]*p_l[3];
          const size_t new_size = buffer_.size() + n;
          s_ = s;
          if (new_size >= buffer_.capacity()) {
            return true; //return 1 if buffer is full and there are more to come
          }

          shells_in_batch_.push_back(s_);
          pqrs_start_.push_back(p_i);
          pqrs_len_.push_back(p_l);
          printf("computing shell %d %d %d %d\n", s_[0], s_[1], s_[2], s_[3]);
          twoint_->compute_shell(s_[0], s_[1], s_[2], s_[3]);
          std::vector<double>::iterator end = buffer_.end();
          buffer_.resize(new_size);
          std::copy(in_buf, in_buf + n, end);

        }
      }
    }
  }

  return false; // if finish these loops, all integrals are done and returned
}
// jtf: End of Int_Batch class declaration

class MP2:public Wavefunction {
    Ref<OneBodyWavefunction> ref_mp2_wfn_;
    double compute_mp2_energy();

  public:
    // constructors to initialize new from input stream or old from savable state
    MP2(const Ref<KeyVal> &);
    MP2(StateIn &);

    // checkpoint the object
    void save_data_state(StateOut &);

    // compute the energy
    void compute(void);

    // deem computed quantites obsolete, require re-computation if needed again
    void obsolete(void);

    // obtain number of electrons - why is this a function, and not a stored property of the molecule?
    int nelectron(void);

    // unimplemented computation of density matrix
    RefSymmSCMatrix density(void);

    // unimplemented(?) check (always returns 0) on whether to use a spin polarized wfn
    int spin_polarized(void);

    // check to determine if this class has implemented a value member function
    int value_implemented(void) const;
};

static ClassDesc MP2_cd(typeid(MP2), "MP2", 1, "public Wavefunction", 0,
    create<MP2> , create<MP2> );

// constructor from input stream
MP2::MP2(const Ref<KeyVal> &keyval) :
  Wavefunction(keyval) {
  ref_mp2_wfn_ << keyval->describedclassvalue("reference");
  if (ref_mp2_wfn_ == 0) {
    throw InputError("require a OneBodyWavefunction object", __FILE__,
        __LINE__, "reference", 0, class_desc());
  }
}

// constructor from checkpointed state
MP2::MP2(StateIn &statein) :
  Wavefunction(statein) {
  ref_mp2_wfn_ << SavableState::restore_state(statein);
}
void MP2::compute(void) {
  // exit if gradient requested
  if (gradient_needed()) {
    throw FeatureNotImplemented("no gradients yet", __FILE__, __LINE__,
        class_desc());
  }

  double extra_hf_acc = 10.;
  ref_mp2_wfn_->set_desired_value_accuracy(desired_value_accuracy()
      / extra_hf_acc);

  // first obtain the reference wavefunction and energy
  double refenergy = ref_mp2_wfn_->energy();

  // then compute the MP2 energy from MO coefficients and newly computed (direct) integrals
  double mp2energy = compute_mp2_energy();

  // dump to output
  ExEnv::out0() << indent << "MP2 Energy = " << mp2energy << endl;

  // change the function.value() to total energy
  set_value(refenergy + mp2energy);
  set_actual_value_accuracy(ref_mp2_wfn_->actual_value_accuracy()
      * extra_hf_acc);
}

// function to checkpoint the current state
void MP2::save_data_state(StateOut &stateout) {
  Wavefunction::save_data_state(stateout);

  SavableState::save_state(ref_mp2_wfn_.pointer(), stateout);
}

// principal function to obtain MP2 energy from reference wavefunction
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
  typedef detail::tuple<4, unsigned int> IntTuple;

  RefSCMatrix vec = ref_mp2_wfn_->eigenvectors();
  //Int_Batch batch(1024, integral()->electron_repulsion());

  //TwoBodyIntBatchGeneric<4> *batch = new TwoBodyIntBatchGeneric<4>( integral()->electron_repulsion() );

  TwoBodyIntBatchGeneric<4> batch(integral()->electron_repulsion());

  int nao = vec.nrow();
  int nmo = vec.ncol();
  int nocc = ref_mp2_wfn_->nelectron() / 2;
  int nvir = nmo - nocc;
  int i;
  auto_vec<double> cvec_av(new double[vec.nrow() * vec.ncol()]);
  double *cvec = cvec_av.get();
  vec->convert(cvec);

#if 0
  auto_vec<double> pqrs_av(new double[nao * nao * nao * nao]);
  double *pqrs = pqrs_av.get();
  for (int n = 0; n < nao * nao * nao * nao; n++)
  pqrs[n] = 0.0;

  Ref<TwoBodyInt> twoint_ = integral()->electron_repulsion();
  const double *buffer = twoint_->buffer();

  TwoBodyIntBatch<4>* tmp = 0;

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

          twoint_->compute_shell(P, Q, R, S);
          // Close loops, accumulate several shell quartet worth of ints into buffer
          // Open next loops, use integrals - not into pqrs[], but into mp2 energy calculation

        }
      }
    }
  }

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

  twoint_ = 0;

#endif

  int nmo4 = nmo * nmo * nmo * nmo;
  auto_vec<double> ijkl_av(new double[nmo4]);
  double *ijkl = ijkl_av.get();
  // need to zero the ijkl[] array?
  bzero(ijkl, nmo4 * sizeof(double));
  int bn = 0;


  int ismore = 1;

  while (ismore) {
    ismore = batch.next();
    printf("batch number = %d\n", ++bn);

    const double *pqrs = batch.buffer();
    int index = 0;

    // jf OK
    std::vector<IntTuple> sib = batch.current_batch(); // shells in this batch
    std::vector<IntTuple>::iterator it;

    // jf OK
    for (it = sib.begin(), i = 0; it < sib.end(); it++, i++) {

      IntTuple start = batch.start()[i];
      IntTuple fence = batch.fence()[i];
      TensorIndexRangeIterator<4> function_range(start, fence);

      // do 4-index transformation
      for (function_range.init();
          function_range.in_range();
          function_range.next(), ++index) {
        const IntTuple& current = function_range.current();
        // at this point, user has information and access to:
        // p, q, r, s indices as p[]
        // integral value at pqrs[idx] note lack of permutational symmetry, etc.
        // offsets po, etc., for locating position in mo coefficient matrix cvec.
        int mo_4idx = 0;
        for (int i = 0; i < nmo; i++) {
          for (int j = 0; j < nmo; j++) {
            for (int k = 0; k < nmo; k++) {
              for (int l = 0; l < nmo; l++, mo_4idx++) {

                // place each ao integral in the batch into the mo integral ijkl[idx]
                ijkl[mo_4idx] += cvec[current[0]*nmo + i] * cvec[current[1]*nmo + j] * cvec[current[2]*nmo + k]
                    * cvec[current[3]*nmo + l] * pqrs[index];

              }
            }
          }
        }
      }
    }
  }

#if 0
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

#endif

  // don't need pqrs autovec any more; using storage within Int_Batch batch
  //  pqrs_av.release();
  //  pqrs = 0;
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

  if (grp == 0) grp = MessageGrp::initial_messagegrp(argc, argv);

  if (grp == 0) {
    grp << keyval->describedclassvalue("messagegrp");
  }

  if (grp == 0) grp = MessageGrp::get_default_messagegrp();

  if (grp == 0) {
    std::cerr << indent << "Couldn't initialize MessageGrp\n";
    abort();
  }

  MessageGrp::set_default_messagegrp(grp);

  Ref<Debugger> debugger;
  debugger << keyval->describedclassvalue(":debug");
  // Let the debugger know the name of the executable and the node
  if (debugger) {
    debugger->set_exec("test");
    debugger->set_prefix(grp->me());
    debugger->debug("curt is a hog");
  }

  RegionTimer::set_default_regiontimer(new ParallelRegionTimer(grp, "test", 1,
      0));

  SCFormIO::set_printnode(0);
  SCFormIO::init_mp(grp->me());
  //SCFormIO::set_debug(1);

  SCFormIO::setindent(ExEnv::outn(), 2);
  SCFormIO::setindent(cerr, 2);

  return grp;
}

int main(int argc, char**argv) {
  const char *input = (argc > 1)? argv[1] : SRCDIR"/tbint_batch_test.in";
  const char *keyword = (argc > 2) ? argv[2] : "mole";
  const char *optkeyword = (argc > 3) ? argv[3] : "opt";
  double val;

  // open keyval input
  Ref<KeyVal> rpkv(new ParsedKeyVal(input));

  init_mp(rpkv, argc, argv);

  Timer tim;
  tim.enter("input");

  // int do_gradient = rpkv->booleanvalue("gradient");

  if (rpkv->exists("matrixkit")) {
    Ref<SCMatrixKit> kit;
    kit << rpkv->describedclassvalue("matrixkit");
    SCMatrixKit::set_default_matrixkit(kit);
  }

  struct stat sb;
  Ref<MolecularEnergy> mole;
  mole << rpkv->describedclassvalue("mole");

  tim.exit("input");

  if (mole) {
    // this line performs the entire computation, as well as printing the energy
    // by accessing mole->energy(), I think
    val = mole->value();

    ExEnv::out0() << indent << "energy: " << mole->energy() << endl;

    if (mole->value_implemented()) {
      ExEnv::out0() << indent << scprintf("value of mole is %20.15f\n\n",
          mole->energy());
    }

    mole->print(ExEnv::out0());
  }

  StateOutBin so("tbint_batch_test.wfn");
  SavableState::save_state(mole.pointer(), so);

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
