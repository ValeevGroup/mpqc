
#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <math/scmat/blocked.h>
#include <math/scmat/local.h>
#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>

#include <chemistry/qc/basis/petite.h>
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

void
CLSCF::set_occupations(const RefDiagSCMatrix& ev)
{
  if (user_occupations_)
    return;
  
  int i,j;
  
  RefDiagSCMatrix evals;
  
  if (ev.null())
    evals = core_hamiltonian().eigvals();
  else
    evals = ev;

  // first convert evals to something we can deal with easily
  LocalDiagSCMatrix *lvals = new LocalDiagSCMatrix(basis()->basisdim(),
                                                   new LocalSCMatrixKit());
  lvals->convert(evals);

  RefPetiteList pl = integral()->petite_list(basis());
  
  double **vals = new double*[nirrep_];
  for (i=j=0; i < nirrep_; i++) {
    int nf=pl->nfunction(i);
    if (nf) {
      vals[i] = new double[nf];
      for (int k=0; k < nf; k++,j++)
        vals[i][k] = lvals->get_element(j);
    } else {
      vals[i] = 0;
    }
  }

  delete lvals;
  
  // now loop to find the tndocc_ lowest eigenvalues and populate those
  // MO's
  int *newocc = new int[nirrep_];
  memset(newocc,0,sizeof(int)*nirrep_);

  for (i=0; i < tndocc_; i++) {
    // find lowest eigenvalue
    int lir,ln;
    double lowest=999999999;

    for (int ir=0; ir < nirrep_; ir++) {
      int nf=pl->nfunction(ir);
      if (!nf)
        continue;
      for (j=0; j < nf; j++) {
        if (vals[ir][j] < lowest) {
          lowest=vals[ir][j];
          lir=ir;
          ln=j;
        }
      }
    }
    vals[lir][ln]=999999999;
    newocc[lir]++;
  }

  // get rid of vals
  for (i=0; i < nirrep_; i++)
    if (vals[i])
      delete[] vals[i];
  delete[] vals;

  if (!ndocc_) {
    ndocc_=newocc;
  } else {
    // test to see if newocc is different from ndocc_
    for (i=0; i < nirrep_; i++) {
      if (ndocc_[i] != newocc[i]) {
        fprintf(stderr,"  CLSCF::set_occupations:  WARNING!!!!\n");
        fprintf(stderr,"    occupations for irrep %d have changed\n",i+1);
        fprintf(stderr,"    ndocc was %d, changed to %d\n",
                ndocc_[i],newocc[i]);
      }
    }

    memcpy(ndocc_,newocc,sizeof(int)*nirrep_);
    delete[] newocc;
  }
}

CLSCF::CLSCF(StateIn& s) :
  SCF(s)
{
  s.get(user_occupations_);
  s.get(tndocc_);
  s.get(nirrep_);
  s.get(ndocc_);
}

CLSCF::CLSCF(const RefKeyVal& keyval) :
  SCF(keyval)
{
  // calculate the total nuclear charge
  int Znuc=0;
  PointBag_double *z = molecule()->charges();
  
  for (Pix i=z->first(); i; z->next(i)) Znuc += (int) z->get(i);

  // check to see if this is to be a charged molecule
  int charge = keyval->intvalue("total_charge");
  
  // now see if ndocc was specified
  if (keyval->exists("ndocc")) {
    tndocc_ = keyval->intvalue("ndocc");
  } else {
    tndocc_ = (Znuc-charge)/2;
    if ((Znuc-charge)%2) {
      fprintf(stderr,
              "\n  CLSCF::init: Warning, there's a leftover electron.\n");
      fprintf(stderr,"    total_charge = %d\n",charge);
      fprintf(stderr,"    total nuclear charge = %d\n", Znuc);
      fprintf(stderr,"    ndocc_ = %d\n", tndocc_);
    }
  }

  printf("\n  CLSCF::init: total charge = %d\n\n", Znuc-2*tndocc_);

  nirrep_ = molecule()->point_group().char_table().ncomp();

  if (keyval->exists("docc")) {
    ndocc_ = new int[nirrep_];
    user_occupations_=1;
    for (int i=0; i < nirrep_; i++) {
      ndocc_[i]=0;

      if (keyval->exists("docc",i))
        ndocc_[i] = keyval->intvalue("docc",i);
    }
  } else {
    ndocc_=0;
    user_occupations_=0;
    set_occupations(0);
  }

  printf("  docc = [");
  for (int i=0; i < nirrep_; i++)
    printf(" %d",ndocc_[i]);
  printf(" ]\n");

  // check to see if this was done in SCF(keyval)
  if (!keyval->exists("maxiter"))
    maxiter_ = 40;
}

CLSCF::~CLSCF()
{
  if (ndocc_) {
    delete[] ndocc_;
    ndocc_=0;
  }
}

void
CLSCF::save_data_state(StateOut& s)
{
  SCF::save_data_state(s);
  s.put(user_occupations_);
  s.put(tndocc_);
  s.put(nirrep_);
  s.put(ndocc_,nirrep_);
}

double
CLSCF::occupation(int ir, int i)
{
  if (i < ndocc_[ir]) return 2.0;
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
  SCF::print(o);
}

/////////////////////////////////////////////////////////////////////////////

void
CLSCF::init_vector()
{
  // calculate the core Hamiltonian
  cl_hcore_ = core_hamiltonian();
  //cl_hcore_.print("core hamiltonian");
  
  // allocate storage for other temp matrices
  cl_dens_ = cl_hcore_.clone();
  cl_dens_.assign(0.0);
  
  cl_dens_diff_ = cl_hcore_.clone();
  cl_dens_diff_.assign(0.0);

  cl_fock_ = cl_hcore_.clone();
  cl_fock_.assign(0.0);
  
  // gmat is in AO basis
  cl_gmat_ = basis()->matrixkit()->symmmatrix(basis()->basisdim());
  cl_gmat_.assign(0.0);

  // test to see if we need a guess vector
  if (eigenvectors_.result_noupdate().null())
    eigenvectors_ = hcore_guess();

  scf_vector_ = eigenvectors_.result_noupdate();
  //scf_vector_.print("guess vector");
}

void
CLSCF::done_vector()
{
  // save these if we're doing the gradient or hessian
  if (!gradient_needed() && !hessian_needed())
    cl_fock_ = 0;

  cl_hcore_ = 0;
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
  BlockedSCMatrix *vecp = BlockedSCMatrix::require_castdown(
    scf_vector_, "CLSCF::new_density: scf_vector");

  BlockedSymmSCMatrix *densp = BlockedSymmSCMatrix::require_castdown(
    cl_dens_, "CLSCF::new_density: density");

  BlockedSymmSCMatrix *ddensp = BlockedSymmSCMatrix::require_castdown(
    cl_dens_diff_, "CLSCF::new_density: density difference");
  
  RefPetiteList pl = integral()->petite_list(basis());
  
  int ij=0;
  double delta=0;

  for (int ir=0; ir < vecp->nblocks(); ir++) {
    int nbasis = pl->nfunction(ir);

    RefSCMatrix vir = vecp->block(ir);
    RefSymmSCMatrix dir = densp->block(ir);
    RefSymmSCMatrix ddir = ddensp->block(ir);
  
    for (int i=0; i < nbasis; i++) {
      for (int j=0; j <= i; j++,ij++) {
        double pt=0;
        int k;
        for (k=0; k < ndocc_[ir]; k++)
          pt += 2.0*vir->get_element(i,k)*vir->get_element(j,k);

        double dlt = pt - dir->get_element(i,j);
        delta += dlt*dlt;
      
        ddir->set_element(i,j,dlt);
        dir->set_element(i,j,pt);
      }
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

  BlockedSymmSCMatrix *densp = BlockedSymmSCMatrix::require_castdown(
    cl_dens_, "CLSCF::scf_energy: density");

  BlockedSymmSCMatrix *fockp = BlockedSymmSCMatrix::require_castdown(
    t, "CLSCF::new_density: H+F");

  double eelec=0;
  for (int ir=0; ir < fockp->nblocks(); ir++) {
    RefSymmSCMatrix dir = densp->block(ir);
    RefSymmSCMatrix fhir = fockp->block(ir);
    
    for (int i=0; i < fhir.n(); i++) {
      for (int j=0; j < i; j++) {
        eelec += dir.get_element(i,j)*fhir.get_element(i,j);
      }
      eelec += 0.5*dir.get_element(i,i)*fhir.get_element(i,i);
    }
  }

  return eelec;
}

void
CLSCF::ao_fock()
{
  RefPetiteList pl = integral()->petite_list(basis());
  
  // calculate G
  RefSymmSCMatrix dd = cl_dens_diff_;
  cl_dens_diff_ = pl->to_AO_basis(dd);
  cl_dens_diff_->scale(2.0);
  cl_dens_diff_->scale_diagonal(0.5);
  //cl_dens_diff_.print("dens diff");
  ao_gmat();

  //cl_gmat_.print("ao gmat");
  cl_dens_diff_ = dd;

  RefSymmSCMatrix foo = cl_gmat_.copy();
  foo.scale(1.0/(double)pl->order());
  //foo.print("foo");
  
  pl->symmetrize(foo,dd);
  //dd = pl->to_SO_basis(cl_gmat_);
  
  //dd.print("SO gmat");
  
  cl_fock_.assign(cl_hcore_);
  cl_fock_.accumulate(dd);

  //cl_fock_.print("SO fock");
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
    fprintf(stderr,"  CLSCF::make_contribution: invalid type %d\n",type);
    abort();
  }

}

RefSCExtrapError
CLSCF::extrap_error()
{
  RefSymmSCMatrix mofock = effective_fock();
  
  BlockedSymmSCMatrix *moerror = BlockedSymmSCMatrix::require_castdown(
    mofock,"CLSCF::extrap_error: moerror");

  for (int ir=0; ir < moerror->nblocks(); ir++) {
    RefSymmSCMatrix moeir = moerror->block(ir);
    
    for (int i=0; i < moeir.n(); i++) {
      double occi = occupation(ir,i);

      for (int j=0; j <= i; j++) {
        double occj = occupation(ir,j);
        if (occi==occj)
          moeir.set_element(i,j,0.0);
      }
    }
  }

  RefSymmSCMatrix aoerror = cl_fock_.clone();
  aoerror.assign(0.0);
  aoerror.accumulate_transform(scf_vector_,mofock);
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
  //mofock.print("MO fock");
  return mofock;
}


/////////////////////////////////////////////////////////////////////////////

void
CLSCF::init_gradient()
{
  // presumably the eigenvectors have already been computed by the time
  // we get here
  scf_vector_ = eigenvectors_.result_noupdate();
}

void
CLSCF::done_gradient()
{
  // save these if we're doing the hessian
  if (!hessian_needed()) {
    cl_fock_=0;
  }

  scf_vector_ = 0;
}

RefSymmSCMatrix
CLSCF::lagrangian()
{
  RefSymmSCMatrix ewdens(basis_dimension(), basis_matrixkit());
  RefSymmSCMatrix evals = effective_fock();
  
#if 0
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
#endif

  return ewdens;
}

RefSymmSCMatrix
CLSCF::gradient_density()
{
  cl_dens_ = basis_matrixkit()->symmmatrix(basis_dimension());
  
#if 0
  for (int i=0; i < scf_vector_->nrow(); i++) {
    for (int j=0; j <= i; j++) {
      double pt=0;
      for (int k=0; k < ndocc_; k++)
        pt += scf_vector_->get_element(i,k)*scf_vector_->get_element(j,k);
      
      cl_dens_->set_element(i,j,2.0*pt);
    }
  }
#endif
  
  return cl_dens_;
}

void
CLSCF::make_gradient_contribution()
{
}

/////////////////////////////////////////////////////////////////////////////

void
CLSCF::init_hessian()
{
}

void
CLSCF::done_hessian()
{
}

