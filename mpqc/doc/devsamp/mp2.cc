
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/scf/clhf.h>

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

static ClassDesc MP2_cd(typeid(MP2), "MP2", 1, "public Wavefunction",
                        0, create<MP2>, create<MP2>);

MP2::MP2(const Ref<KeyVal> &keyval):Wavefunction(keyval) {
  ref_mp2_wfn_ << keyval->describedclassvalue("reference");
  if(ref_mp2_wfn_.null()) {
    ExEnv::out0() << "reference is null" << endl;
    abort();
  }
}

MP2::MP2(StateIn &statein):Wavefunction(statein)
{
  ref_mp2_wfn_ << SavableState::restore_state(statein);
}

void
MP2::save_data_state(StateOut &stateout) {
  Wavefunction::save_data_state(stateout);

  SavableState::save_state(ref_mp2_wfn_.pointer(),stateout);
}

void
MP2::compute(void)
{
  if(gradient_needed()) {
    ExEnv::out0() << "No gradients yet" << endl;
    abort();
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

void
MP2::obsolete(void) {
  Wavefunction::obsolete();
  ref_mp2_wfn_->obsolete();
}

int
MP2::nelectron(void) {
  return ref_mp2_wfn_->nelectron();
}

RefSymmSCMatrix
MP2::density(void) {
  ExEnv::out0() << "No density yet" << endl;
  abort();
  return 0;
}

int
MP2::spin_polarized(void) {
  return 0;
}

int
MP2::value_implemented(void) const {
  return 1;
}

double
MP2::compute_mp2_energy()
{
  if(molecule()->point_group()->char_table().order() != 1) {
    ExEnv::out0() << "C1 symmetry only" << endl;
    abort();
  }

  RefSCMatrix vec = ref_mp2_wfn_->eigenvectors();

  int nao = vec.nrow();
  int nmo = vec.ncol();
  int nocc = ref_mp2_wfn_->nelectron()/2;
  int nvir = nmo - nocc;

  double *cvec = new double [vec.nrow() * vec.ncol()];
  vec->convert(cvec);

  double *pqrs = new double [nao * nao * nao * nao];
  for(int n = 0; n < nao*nao*nao*nao; n++) pqrs[n] = 0.0;

  Ref<TwoBodyInt> twoint = integral()->electron_repulsion();
  const double *buffer = twoint->buffer();

  Ref<GaussianBasisSet> basis = this->basis();

  int nshell = basis->nshell();
  for(int P = 0; P < nshell; P++) {
    int nump = basis->shell(P).nfunction();

    for(int Q = 0; Q < nshell; Q++) {
      int numq = basis->shell(Q).nfunction();

      for(int R = 0; R < nshell; R++) {
	int numr = basis->shell(R).nfunction();

	for(int S = 0; S < nshell; S++) {
	  int nums = basis->shell(S).nfunction();

	  twoint->compute_shell(P,Q,R,S);

	  int index = 0;
	  for(int p=0; p < nump; p++) {
	    int op = basis->shell_to_function(P)+p;

	    for(int q = 0; q < numq; q++) {
	      int oq = basis->shell_to_function(Q)+q;

	      for(int r = 0; r < numr; r++) {
		int oor = basis->shell_to_function(R)+r;

		for(int s = 0; s < nums; s++,index++) {
		  int os = basis->shell_to_function(S)+s;

		  int ipqrs = (((op*nao+oq)*nao+oor)*nao+os);

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

  double *ijkl = new double [nmo * nmo * nmo * nmo];

  int idx = 0;
  for(int i = 0; i < nmo; i++) {
    for(int j = 0; j < nmo; j++) {
      for(int k = 0; k < nmo; k++) {
	for(int l = 0; l < nmo; l++, idx++) {

	  ijkl[idx] = 0.0;

	  int index = 0;
	  for(int p = 0; p < nao; p++) {
	    for(int q = 0; q < nao; q++) {
	      for(int r = 0; r < nao; r++) {
		for(int s = 0; s < nao; s++,index++) {

		  ijkl[idx] += cvec[p*nmo + i] * cvec[q*nmo +j]
		    * cvec[r*nmo + k] * cvec[s*nmo + l]
		    * pqrs[index];

		}
	      }
	    }
	  }

	}
      }
    }
  }

  delete [] pqrs;
  delete [] cvec;

  double *evals = new double [nmo];
  ref_mp2_wfn_->eigenvalues()->convert(evals);

  double energy = 0.0;
  for(int i=0; i < nocc; i++) {
    for(int j=0; j < nocc; j++) {
      for(int a=nocc; a < nmo; a++) {
	for(int b=nocc; b < nmo; b++) {

	  int iajb = (((i*nmo+a)*nmo+j)*nmo+b);
	  int ibja = (((i*nmo+b)*nmo+j)*nmo+a);

	  energy += (2 * ijkl[iajb] - ijkl[ibja]) * ijkl[iajb]/
	    (evals[i] + evals[j] - evals[a] - evals[b]);

	}
      }
    }
  }

  delete [] ijkl;
  delete [] evals;

  return energy;
}
