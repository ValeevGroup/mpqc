
#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>

#include <chemistry/qc/wfn/hcore.h>
#include <chemistry/qc/scf/clscf.h>

///////////////////////////////////////////////////////////////////////////
// CLSCF

#define CLASSNAME CLSCF
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public SCF
#include <util/class/classi.h>
void *
CLSCF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCF::_castdown(cd);
  return do_castdowns(casts,cd);
}

void
CLSCF::init()
{
}

CLSCF::CLSCF(StateIn& s) :
  SCF(s)
{
  s.get(ndocc_);
}

CLSCF::CLSCF(const RefKeyVal& keyval) :
  SCF(keyval)
{
  if (keyval->exists("ndocc")) {
    ndocc_ = keyval->intvalue("ndocc");
  } else {
    int Z=0;
    PointBag_double *z = _mol->charges();
  
    for (Pix i=z->first(); i; z->next(i)) Z += (int) z->get(i);

    ndocc_ = Z/2;
    if (Z%2) {
      fprintf(stderr,"CLSCF::init: Warning, there's a leftover electron.\n");
      fprintf(stderr,"  total nuclear charge = %d, %d closed shells\n",
              Z, ndocc_);
      fprintf(stderr,"  total charge = %d\n\n",Z-2*ndocc_);
    }
  }

  // check to see if this was done in SCF(keyval)
  if (!keyval->exists("maxiter"))
    maxiter_ = 40;
}

CLSCF::~CLSCF()
{
}

void
CLSCF::save_data_state(StateOut& s)
{
  s.put(ndocc_);
}

double
CLSCF::occupation(int i)
{
  if (i < ndocc_) return 2.0;
  return 0.0;
}

int
CLSCF::value_implemented()
{
  return 1;
}

int
CLSCF::gradient_implemented()
{
  return 0;
}

int
CLSCF::hessian_implemented()
{
  return 0;
}

void
CLSCF::print(ostream&o)
{
  OneBodyWavefunction::print(o);
}

/////////////////////////////////////////////////////////////////////////////

void
CLSCF::init_vector()
{
  // calculate the core Hamiltonian
  RefAccumHCore hc = new AccumHCore();
  hc->init(basis(), integral());

  cl_hcore_ = matrixkit()->symmmatrix(basis_dimension());
  hc->accum(cl_hcore_);
  hc=0;
  
  // allocate storage for other temp matrices
  cl_dens_ = cl_hcore_.clone();
  cl_dens_.assign(0.0);
  
  cl_dens_diff_ = cl_hcore_.clone();
  cl_dens_diff_.assign(0.0);

  cl_gmat_ = cl_hcore_.clone();
  cl_gmat_.assign(0.0);

  cl_fock_ = cl_hcore_.clone();
  cl_fock_.assign(0.0);
  
  // test to see if we need a guess vector
  if (_eigenvectors.result_noupdate().null()) {
    RefSCMatrix vec(basis_dimension(), basis_dimension(), matrixkit());
    RefDiagSCMatrix val(basis_dimension(), matrixkit());
    cl_hcore_.diagonalize(val,vec);

    // transform to S**(-1/2) basis
    vec = ao_to_orthog_ao() * vec;
    
    vec->schmidt_orthog(overlap().pointer(), ndocc_);
    
    _eigenvectors = vec;
  }

  scf_vector_ = _eigenvectors.result_noupdate();
}

void
CLSCF::done_vector()
{
  cl_hcore_ = 0;
  cl_fock_ = 0;
  cl_gmat_ = 0;
  cl_dens_ = 0;
  cl_dens_diff_ = 0;

  scf_vector_ = 0;
}

void
CLSCF::reset_density()
{
  cl_gmat_.assign(0.0);
  cl_dens_diff_.assign(cl_dens_);
}

double
CLSCF::new_density()
{
  int nbasis = basis()->nbasis();

  double delta=0;
  
  // find out what type of matrices we're dealing with
  int ij=0;
  for (int i=0; i < nbasis; i++) {
    for (int j=0; j <= i; j++,ij++) {
      double pt=0;
      int k;
      for (k=0; k < ndocc_; k++)
        pt += 2.0*scf_vector_->get_element(i,k)*scf_vector_->get_element(j,k);

      double dlt = pt-cl_dens_->get_element(i,j);
      delta += dlt*dlt;
      
      cl_dens_diff_->set_element(i,j,dlt);
      cl_dens_->set_element(i,j,pt);
    }
  }

  delta = sqrt(delta/ij);
  return delta;
}

double
CLSCF::scf_energy()
{
  RefSymmSCMatrix t = cl_fock_.copy();
  t.accumulate(cl_hcore_);

  double eelec=0;
  for (int i=0; i < t->n(); i++) {
    for (int j=0; j < i; j++) {
      eelec += cl_dens_.get_element(i,j)*t.get_element(i,j);
    }
    eelec += 0.5*cl_dens_.get_element(i,i)*t.get_element(i,i);
  }

  return eelec;
}

void
CLSCF::ao_fock()
{
  // calculate G
  cl_dens_diff_->scale(2.0);
  cl_dens_diff_->scale_diagonal(0.5);
  ao_gmat();
  cl_dens_diff_->scale(0.5);
  cl_dens_diff_->scale_diagonal(2.0);
  
  cl_fock_.assign(cl_hcore_);
  cl_fock_.accumulate(cl_gmat_);
}

void
CLSCF::make_contribution(int i, int j, int k, int l, double val, int type)
{
  SymmSCMatrix& gmat = *cl_gmat_.pointer();
  SymmSCMatrix& pmat = *cl_dens_diff_.pointer();
  
  switch(type) {
  case 1:
    gmat.accumulate_element(i, j, val*pmat.get_element(k,l));
    gmat.accumulate_element(k, l, val*pmat.get_element(i,j));
    break;
    
  case 2:
    gmat.accumulate_element(i, j, -0.25*val*pmat.get_element(k,l));
    gmat.accumulate_element(k, l, -0.25*val*pmat.get_element(i,j));
    break;
    
  case 3:
    gmat.accumulate_element(i, j, -0.5*val*pmat.get_element(k,l));
    gmat.accumulate_element(k, l, -0.5*val*pmat.get_element(i,j));
    break;
    
  case 4:
    gmat.accumulate_element(i, j, 0.75*val*pmat.get_element(k,l));
    gmat.accumulate_element(k, l, 0.75*val*pmat.get_element(i,j));
    break;
    
  case 5:
    gmat.accumulate_element(i, j, 0.5*val*pmat.get_element(k,l));
    gmat.accumulate_element(k, l, 0.5*val*pmat.get_element(i,j));
    break;
    
  default:
    fprintf(stderr,"CLSCF::make_contribution: invalid type %d\n",type);
    abort();
  }

}

RefSCExtrapError
CLSCF::extrap_error()
{
  RefSymmSCMatrix moerror = effective_fock();

  for (int i=0; i < moerror->n(); i++) {
    double occi = occupation(i);

    for (int j=0; j <= i; j++) {
      double occj = occupation(j);
      if (occi==occj)
        moerror.set_element(i,j,0.0);
    }
  }

  RefSymmSCMatrix aoerror = moerror.clone();
  aoerror.assign(0.0);
  aoerror.accumulate_transform(scf_vector_,moerror);
  moerror=0;

  RefSCExtrapError error = new SymmSCMatrixSCExtrapError(aoerror);
  aoerror=0;

  return error;
}

RefSCExtrapData
CLSCF::extrap_data()
{
  RefSCExtrapData data = new SymmSCMatrixSCExtrapData(cl_fock_);
  return data;
}

RefSymmSCMatrix
CLSCF::effective_fock()
{
  // just return MO fock matrix
  RefSymmSCMatrix mofock = cl_fock_.clone();
  mofock.assign(0.0);
  mofock.accumulate_transform(scf_vector_.t(), cl_fock_);
  return mofock;
}


/////////////////////////////////////////////////////////////////////////////

void
CLSCF::init_gradient()
{
}

void
CLSCF::done_gradient()
{
}

void
CLSCF::init_hessian()
{
}

void
CLSCF::done_hessian()
{
}

