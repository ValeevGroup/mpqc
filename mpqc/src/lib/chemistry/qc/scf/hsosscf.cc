
#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>

#include <chemistry/qc/wfn/hcore.h>
#include <chemistry/qc/scf/hsosscf.h>

///////////////////////////////////////////////////////////////////////////
// HSOSSCF

#define CLASSNAME HSOSSCF
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public SCF
#include <util/class/classi.h>
void *
HSOSSCF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCF::_castdown(cd);
  return do_castdowns(casts,cd);
}

void
HSOSSCF::init()
{
}

HSOSSCF::HSOSSCF(StateIn& s) :
  SCF(s)
{
  s.get(ndocc_);
  s.get(nsocc_);
}

HSOSSCF::HSOSSCF(const RefKeyVal& keyval) :
  SCF(keyval)
{
  // check to see if occupation numbers have been given
  if (keyval->exists("ndocc") && keyval->exists("nsocc")) {
    ndocc_ = keyval->intvalue("ndocc");
    nsocc_ = keyval->intvalue("nsocc");

    // if only the number of singly occupied orbitals is given, then try
    // to calculate the number of doubly occupied
  } else if (keyval->exists("nsocc")) {

    // if we only have the number of doubly occupied orbitals, then try
    // to figure out the number of singly occupied.
  } else if (keyval->exists("nsocc")) {

    // we have no occupation information.  
  } else {
  }
}

HSOSSCF::~HSOSSCF()
{
}

void
HSOSSCF::save_data_state(StateOut& s)
{
  SCF::save_data_state(s);
  s.put(ndocc_);
  s.put(nsocc_);
}

double
HSOSSCF::occupation(int ir, int i)
{
  if (i < ndocc_) return 2.0;
  if (i < ndocc_ + nsocc_) return 1.0;
  return 0.0;
}

int
HSOSSCF::value_implemented()
{
  return 1;
}

int
HSOSSCF::gradient_implemented()
{
  return 0;
}

int
HSOSSCF::hessian_implemented()
{
  return 0;
}

void
HSOSSCF::print(ostream&o)
{
  SCF::print(o);
}

//////////////////////////////////////////////////////////////////////////////

void
HSOSSCF::init_vector()
{
  // calculate the core Hamiltonian
  cl_hcore_ = core_hamiltonian();
  
  // allocate storage for other temp matrices
  cl_dens_ = cl_hcore_.clone();
  cl_dens_.assign(0.0);
  
  cl_dens_diff_ = cl_hcore_.clone();
  cl_dens_diff_.assign(0.0);

  cl_gmat_ = cl_hcore_.clone();
  cl_gmat_.assign(0.0);

  cl_fock_ = cl_hcore_.clone();
  cl_fock_.assign(0.0);

  op_dens_ = cl_hcore_.clone();
  op_dens_.assign(0.0);
  
  op_dens_diff_ = cl_hcore_.clone();
  op_dens_diff_.assign(0.0);

  op_gmat_ = cl_hcore_.clone();
  op_gmat_.assign(0.0);

  op_fock_ = cl_hcore_.clone();
  op_fock_.assign(0.0);
  
  // test to see if we need a guess vector
  if (eigenvectors_.result_noupdate().null())
    eigenvectors_ = hcore_guess();

  scf_vector_ = eigenvectors_.result_noupdate();
}

void
HSOSSCF::done_vector()
{
  // save these if we're doing the gradient or hessian
  if (!gradient_needed() && !hessian_needed()) {
    cl_fock_ = 0;
    op_fock_ = 0;
  }

  cl_hcore_ = 0;
  cl_gmat_ = 0;
  cl_dens_ = 0;
  cl_dens_diff_ = 0;
  op_gmat_ = 0;
  op_dens_ = 0;
  op_dens_diff_ = 0;

  scf_vector_ = 0;
}

void
HSOSSCF::reset_density()
{
  cl_gmat_.assign(0.0);
  cl_dens_diff_.assign(cl_dens_);
  op_gmat_.assign(0.0);
  op_dens_diff_.assign(op_dens_);
}

double
HSOSSCF::new_density()
{
  int nbasis = basis()->nbasis();

  double delta=0;
  
  int ij=0;
  for (int i=0; i < nbasis; i++) {
    for (int j=0; j <= i; j++,ij++) {
      double pt=0, po=0;
      int k;
      for (k=0; k < ndocc_; k++)
        pt += 2.0*scf_vector_->get_element(i,k)*scf_vector_->get_element(j,k);
      for (k=ndocc_; k < ndocc_+nsocc_; k++)
        po += scf_vector_->get_element(i,k)*scf_vector_->get_element(j,k);

      double dlt = pt+po-cl_dens_->get_element(i,j);
      double dlto = po-op_dens_->get_element(i,j);
      delta += dlt*dlt;
      
      cl_dens_diff_->set_element(i,j,dlt);
      cl_dens_->set_element(i,j,pt);
      op_dens_diff_->set_element(i,j,dlto);
      op_dens_->set_element(i,j,po);
    }
  }

  delta = sqrt(delta/ij);
  return delta;
}

double
HSOSSCF::scf_energy()
{
  RefSymmSCMatrix t = cl_fock_.copy();
  t.accumulate(cl_hcore_);

  double eelec=0;
  for (int i=0; i < t->n(); i++) {
    for (int j=0; j < i; j++) {
      eelec += cl_dens_.get_element(i,j)*t.get_element(i,j)
               - op_dens_.get_element(i,j)*op_fock_.get_element(i,j);
    }
    eelec += 0.5*(cl_dens_.get_element(i,i)*t.get_element(i,i)
               - op_dens_.get_element(i,i)*op_fock_.get_element(i,i));
  }
  return eelec;
}

void
HSOSSCF::ao_fock()
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
HSOSSCF::make_contribution(int i, int j, int k, int l, double val, int type)
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
    fprintf(stderr,"HSOSSCF::make_contribution: invalid type %d\n",type);
    abort();
  }

}

RefSCExtrapError
HSOSSCF::extrap_error()
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
HSOSSCF::extrap_data()
{
  RefSCExtrapData data = new SymmSCMatrixSCExtrapData(cl_fock_);
  return data;
}

RefSymmSCMatrix
HSOSSCF::effective_fock()
{
  // just return MO fock matrix
  RefSymmSCMatrix mofock = cl_fock_.clone();
  mofock.assign(0.0);
  mofock.accumulate_transform(scf_vector_.t(), cl_fock_);
  return mofock;
}


/////////////////////////////////////////////////////////////////////////////

void
HSOSSCF::init_gradient()
{
  // presumably the eigenvectors have already been computed by the time
  // we get here
  scf_vector_ = eigenvectors_.result_noupdate();
}

void
HSOSSCF::done_gradient()
{
  // save these if we're doing the hessian
  if (!hessian_needed()) {
    cl_fock_=0;
  }

  scf_vector_ = 0;
}

RefSymmSCMatrix
HSOSSCF::lagrangian()
{
  RefSymmSCMatrix ewdens(basis_dimension(), basis_matrixkit());
  RefSymmSCMatrix evals = effective_fock();
  
  for (int i=0; i < scf_vector_->nrow(); i++) {
    for (int j=0; j <= i; j++) {
      double pt=0;
      for (int k=0; k < ndocc_; k++)
        pt += scf_vector_->get_element(i,k)*scf_vector_->get_element(j,k)*
          evals->get_element(k,k);
      
      ewdens->set_element(i,j,pt);
    }
  }
  ewdens->scale(-2.0);

  return ewdens;
}

RefSymmSCMatrix
HSOSSCF::gradient_density()
{
  cl_dens_ = basis_matrixkit()->symmmatrix(basis_dimension());
  
  for (int i=0; i < scf_vector_->nrow(); i++) {
    for (int j=0; j <= i; j++) {
      double pt=0;
      for (int k=0; k < ndocc_; k++)
        pt += scf_vector_->get_element(i,k)*scf_vector_->get_element(j,k);
      
      cl_dens_->set_element(i,j,2.0*pt);
    }
  }

  return cl_dens_;
}

void
HSOSSCF::make_gradient_contribution()
{
}

/////////////////////////////////////////////////////////////////////////////

void
HSOSSCF::init_hessian()
{
}

void
HSOSSCF::done_hessian()
{
}

