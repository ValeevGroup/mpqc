
#ifdef __GNUC__
#pragma implementation
#pragma implementation "tccont.h"
#endif

#include <iostream.h>
#include <iomanip.h>

#include <math.h>

#include <util/misc/timer.h>
#include <util/misc/formio.h>

#include <math/scmat/block.h>
#include <math/scmat/blocked.h>
#include <math/scmat/local.h>

#include <math/optimize/scextrapmat.h>

#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/scf/scflocal.h>
#include <chemistry/qc/scf/scfden.h>
#include <chemistry/qc/scf/effh.h>

#include <chemistry/qc/scf/tcscf.h>
#include <chemistry/qc/scf/tccont.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/ltbgrad.h>

///////////////////////////////////////////////////////////////////////////

#ifdef __GNUC__
template class GBuild<LocalTCContribution>;
template class LocalGBuild<LocalTCContribution>;

template class TBGrad<LocalTCGradContribution>;
template class LocalTBGrad<LocalTCGradContribution>;
#endif

///////////////////////////////////////////////////////////////////////////
// TCSCF

#define CLASSNAME TCSCF
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public SCF
#include <util/class/classi.h>
void *
TCSCF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCF::_castdown(cd);
  return do_castdowns(casts,cd);
}

TCSCF::TCSCF(StateIn& s) :
  SCF(s),
  focka_(this),
  fockb_(this),
  ka_(this),
  kb_(this)
  maybe_SavableState(s)
{
  focka_.result_noupdate() = basis_matrixkit()->symmmatrix(basis_dimension());
  focka_.restore_state(s);
  focka_.result_noupdate().restore(s);
  
  fockb_.result_noupdate() = basis_matrixkit()->symmmatrix(basis_dimension());
  fockb_.restore_state(s);
  fockb_.result_noupdate().restore(s);
  
  ka_.result_noupdate() = basis_matrixkit()->symmmatrix(basis_dimension());
  ka_.restore_state(s);
  ka_.result_noupdate().restore(s);
  
  kb_.result_noupdate() = basis_matrixkit()->symmmatrix(basis_dimension());
  kb_.restore_state(s);
  kb_.result_noupdate().restore(s);
  
  s.get(user_occupations_);
  s.get(tndocc_);
  s.get(nirrep_);
  s.get(ndocc_);
  s.get(osa_);
  s.get(osb_);
  s.get(occa_);
  s.get(occb_);
  s.get(ci1_);
  s.get(ci2_);

  // now take care of memory stuff
  init_mem(8);
}

TCSCF::TCSCF(const RefKeyVal& keyval) :
  SCF(keyval),
  focka_(this),
  fockb_(this),
  ka_(this),
  kb_(this)
{
  int me = scf_grp_->me();
  
  focka_.compute()=0;
  focka_.computed()=0;
  
  fockb_.compute()=0;
  fockb_.computed()=0;
  
  ka_.compute()=0;
  ka_.computed()=0;
  
  kb_.compute()=0;
  kb_.computed()=0;
  
  // calculate the total nuclear charge
  int Znuc=0;
  PointBag_double *z = molecule()->charges();
  
  for (Pix i=z->first(); i; z->next(i)) Znuc += (int) z->get(i);

  // check to see if this is to be a charged molecule
  int charge = keyval->intvalue("total_charge");
  int nelectrons = Znuc-charge;

  // figure out how many doubly occupied shells there are
  if (keyval->exists("ndocc")) {
    tndocc_ = keyval->intvalue("ndocc");
  } else {
    tndocc_ = (nelectrons-2)/2;
    if ((nelectrons-2)%2) {
      if (me==0) {
        cerr << endl << indent
             << "TCSCF::init: Warning, there's a leftover electron.\n";
        cerr << incindent;
        cerr << indent << "total_charge = " << charge << endl;
        cerr << indent << "total nuclear charge = " << Znuc << endl;
        cerr << indent << "ndocc_ = " << tndocc_ << endl << decindent;
      }
    }
  }

  if (me==0)
    cout << endl << indent << "TCSCF::init: total charge = "
         << Znuc-2*tndocc_-2 << endl << endl;

  nirrep_ = molecule()->point_group().char_table().ncomp();

  if (nirrep_==1) {
    cerr << indent << "TCSCF::init: cannot do C1 symmetry\n";
    abort();
  }

  occa_=occb_=1.0;
  ci1_=ci2_ = 0.5*sqrt(2.0);
  
  if (keyval->exists("ci1")) {
    ci1_ = keyval->doublevalue("ci1");
    ci2_ = sqrt(1.0 - ci1_*ci1_);
    occa_ = 2.0*ci1_*ci1_;
    occb_ = 2.0*ci2_*ci2_;
  }

  if (keyval->exists("occa")) {
    occa_ = keyval->doublevalue("occa");
    ci1_ = sqrt(occa_/2.0);
    ci2_ = sqrt(1.0 - ci1_*ci1_);
    occb_ = 2.0*ci2_*ci2_;
  }

  osa_=-1;
  osb_=-1;

  if (keyval->exists("docc") && keyval->exists("socc")) {
    ndocc_ = new int[nirrep_];
    user_occupations_=1;
    for (int i=0; i < nirrep_; i++) {
      ndocc_[i] = keyval->intvalue("docc",i);
      int nsi = keyval->intvalue("socc",i);
      if (nsi && osa_<0)
        osa_==i;
      else if (nsi && osb_<0)
        osb_==i;
      else if (nsi) {
        cerr << indent << "TCSCF::init: too many open shells\n";
        abort();
      }
    }
  } else {
    ndocc_=0;
    user_occupations_=0;
    set_occupations(0);
  }

  if (me==0) {
    cout << indent << "docc = [";
    for (int i=0; i < nirrep_; i++)
      cout << " " << ndocc_[i];
    cout << " ]\n";

    cout << indent << "socc = [";
    for (int i=0; i < nirrep_; i++)
      cout << " " << (i==osa_ || i==osb_) ? 1 : 0;
    cout << " ]\n";
  }

  // check to see if this was done in SCF(keyval)
  if (!keyval->exists("maxiter"))
    maxiter_ = 200;

  if (!keyval->exists("level_shift"))
    level_shift_ = 1.0;

  // now take care of memory stuff
  init_mem(8);
}

TCSCF::~TCSCF()
{
}

void
TCSCF::save_data_state(StateOut& s)
{
  SCF::save_data_state(s);
  
  focka_.save_data_state(s);
  focka_.result_noupdate().save(s);
  
  fockb_.save_data_state(s);
  fockb_.result_noupdate().save(s);
  
  ka_.save_data_state(s);
  ka_.result_noupdate().save(s);
  
  kb_.save_data_state(s);
  kb_.result_noupdate().save(s);
  
  s.put(user_occupations_);
  s.put(tndocc_);
  s.put(nirrep_);
  s.put(ndocc_,nirrep_);
  s.put(osa_);
  s.put(osb_);
  s.put(occa_);
  s.put(occb_);
  s.put(ci1_);
  s.put(ci2_);
}

double
TCSCF::occupation(int ir, int i)
{
  if (i < ndocc_[ir])
    return 2.0;
  else if (ir==osa_ && i==ndocc_[ir])
    return occa_;
  else if (ir==osb_ && i==ndocc_[ir])
    return occb_;
  else
    return 0.0;
}

int
TCSCF::n_fock_matrices() const
{
  return 2;
}

RefSymmSCMatrix
TCSCF::fock(int n)
{
  if (n > 4) {
    cerr << indent << "TCSCF::fock: there are only four fock matrices, "
         << "but fock(" << n << ") was requested" << endl;
    abort();
  }

  if (n==0)
    return focka_.result();
  else if (n==1)
    return fockb_.result();
  else if (n==2)
    return ka_.result();
  else
    return kb_.result();
}

int
TCSCF::value_implemented()
{
  return 1;
}

int
TCSCF::gradient_implemented()
{
  return 1;
}

int
TCSCF::hessian_implemented()
{
  return 0;
}

void
TCSCF::print(ostream&o)
{
  SCF::print(o);
  if (scf_grp_->me()==0) {
    o << indent << "TCSCF Parameters:\n" << incindent;
    o << indent << "ndocc = " << tndocc_ << endl;
    o << indent << "occa = " << occa_ << endl;
    o << indent << "occb = " << occb_ << endl;
    o << indent << "ci1 = " << ci1_ << endl;
    o << indent << "ci2 = " << ci2_ << endl;
    o << indent << "docc = [";
    for (int i=0; i < nirrep_; i++)
      o << " " << ndocc_[i];
    o << " ]" << endl;
    o << indent << "socc = [";
    for (int i=0; i < nirrep_; i++)
      o << " " << (i==osa_ || i==osb_) ? 1 : 0;
    o << " ]" << endl << decindent << endl;
  }
}

//////////////////////////////////////////////////////////////////////////////

void
TCSCF::set_occupations(const RefDiagSCMatrix& ev)
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
  BlockedDiagSCMatrix *evalsb = BlockedDiagSCMatrix::require_castdown(evals,
                                                 "TCSCF::set_occupations");
  
  RefPetiteList pl = integral()->petite_list(basis());
  
  double **vals = new double*[nirrep_];
  for (i=0; i < nirrep_; i++) {
    int nf=pl->nfunction(i);
    if (nf) {
      vals[i] = new double[nf];
      evalsb->block(i)->convert(vals[i]);
    } else {
      vals[i] = 0;
    }
  }

  // now loop to find the tndocc_ lowest eigenvalues and populate those
  // MO's
  int *newdocc = new int[nirrep_];
  memset(newdocc,0,sizeof(int)*nirrep_);

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
    newdocc[lir]++;
  }

  int osa=-1, osb=-1;
  
  for (i=0; i < 2; i++) {
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

    if (!i) {
      osa=lir;
    } else {
      if (lir==osa) {
        i--;
        continue;
      }
      osb=lir;
    }
  }

   if (osa > osb) {
     int tmp=osa;
     osa=osb;
     osb=tmp;
   }
  
  // get rid of vals
  for (i=0; i < nirrep_; i++)
    if (vals[i])
      delete[] vals[i];
  delete[] vals;

  if (!ndocc_) {
    ndocc_=newdocc;
    osa_=osa;
    osb_=osb;
  } else {
    // test to see if newocc is different from ndocc_
    for (i=0; i < nirrep_; i++) {
      if (ndocc_[i] != newdocc[i] && scf_grp_->me()==0) {
        cerr << indent << "TCSCF::set_occupations:  WARNING!!!!\n";
        cerr << incindent << indent <<
          "occupations for irrep " << i+1 << " have changed\n";
        cerr << indent << "ndocc was " << ndocc_[i] << ", changed to "
             << newdocc[i] << endl << decindent;
      }
      if (((osa != osa_ && osa != osb_) || (osb != osb_ && osb != osa_))
          && scf_grp_->me()==0) {
        cerr << indent << "TCSCF::set_occupations:  WARNING!!!!\n";
        cerr << incindent << indent
             << "open shell occupations have changed"
             << endl << decindent;
        osa_=osa;
        osb_=osb;
        reset_density();
      }
    }

    memcpy(ndocc_,newdocc,sizeof(int)*nirrep_);
    
    delete[] newdocc;
  }
}

//////////////////////////////////////////////////////////////////////////////
//
// scf things
//

void
TCSCF::init_vector()
{
  // initialize the two electron integral classes
  tbi_ = integral()->electron_repulsion();

  // calculate the core Hamiltonian
  cl_hcore_ = core_hamiltonian();
  
  // allocate storage for other temp matrices
  cl_dens_ = cl_hcore_.clone();
  cl_dens_.assign(0.0);
  
  cl_dens_diff_ = cl_hcore_.clone();
  cl_dens_diff_.assign(0.0);

  op_densa_ = cl_hcore_.clone();
  op_densa_.assign(0.0);
  
  op_densa_diff_ = cl_hcore_.clone();
  op_densa_diff_.assign(0.0);

  op_densb_ = cl_hcore_.clone();
  op_densb_.assign(0.0);
  
  op_densb_diff_ = cl_hcore_.clone();
  op_densb_diff_.assign(0.0);

  // gmat is in AO basis
  ao_gmata_ = basis()->matrixkit()->symmmatrix(basis()->basisdim());
  ao_gmata_.assign(0.0);

  ao_gmatb_ = ao_gmata_.clone();
  ao_gmatb_.assign(0.0);

  ao_ka_ = ao_gmata_.clone();
  ao_ka_.assign(0.0);

  ao_kb_ = ao_gmata_.clone();
  ao_kb_.assign(0.0);

  // test to see if we need a guess vector
  if (eigenvectors_.result_noupdate().null()) {
    initial_vector();

    focka_ = cl_hcore_.clone();
    focka_.result_noupdate().assign(0.0);
    fockb_ = cl_hcore_.clone();
    fockb_.result_noupdate().assign(0.0);
    ka_ = cl_hcore_.clone();
    ka_.result_noupdate().assign(0.0);
    kb_ = cl_hcore_.clone();
    kb_.result_noupdate().assign(0.0);
  }

  scf_vector_ = eigenvectors_.result_noupdate();
}

void
TCSCF::done_vector()
{
  tbi_=0;
  
  cl_hcore_ = 0;
  cl_dens_ = 0;
  cl_dens_diff_ = 0;
  op_densa_ = 0;
  op_densa_diff_ = 0;
  op_densb_ = 0;
  op_densb_diff_ = 0;

  ao_gmata_ = 0;
  ao_gmatb_ = 0;
  ao_ka_ = 0;
  ao_kb_ = 0;

  scf_vector_ = 0;
}

void
TCSCF::reset_density()
{
  cl_dens_diff_.assign(cl_dens_);
  
  ao_gmata_.assign(0.0);
  op_densa_diff_.assign(op_densa_);

  ao_gmatb_.assign(0.0);
  op_densb_diff_.assign(op_densb_);

  ao_ka_.assign(0.0);
  ao_kb_.assign(0.0);
}

double
TCSCF::new_density()
{
  // copy current density into density diff and scale by -1.  later we'll
  // add the new density to this to get the density difference.
  cl_dens_diff_.assign(cl_dens_);
  cl_dens_diff_.scale(-1.0);

  op_densa_diff_.assign(op_densa_);
  op_densa_diff_.scale(-1.0);

  op_densb_diff_.assign(op_densb_);
  op_densb_diff_.scale(-1.0);

  cl_dens_.assign(0.0);
  RefSCElementOp op = new SCFDensity(this, scf_vector_, 2.0);
  cl_dens_.element_op(op);
  cl_dens_.scale(2.0);

  op_densa_.assign(0.0);
  op = new SCFDensity(this, scf_vector_, occa_);
  op_densa_.element_op(op);
  BlockedSymmSCMatrix::castdown(op_densa_)->block(osb_)->assign(0.0);
  op_densa_.scale(2.0);

  op_densb_.assign(0.0);
  op = new SCFDensity(this, scf_vector_, occb_);
  op_densb_.element_op(op);
  BlockedSymmSCMatrix::castdown(op_densb_)->block(osa_)->assign(0.0);
  op_densb_.scale(2.0);

  cl_dens_diff_.accumulate(cl_dens_);
  op_densa_diff_.accumulate(op_densa_);
  op_densb_diff_.accumulate(op_densb_);

  RefSymmSCMatrix del = cl_dens_diff_.copy();
  del.accumulate(op_densa_diff_);
  del.accumulate(op_densb_diff_);
  
  RefSCElementScalarProduct sp(new SCElementScalarProduct);
  del.element_op(sp, del);
  
  double delta = sp->result();
  delta = sqrt(delta/i_offset(cl_dens_diff_.n()));

  return delta;
}

double
TCSCF::scf_energy()
{
  // first calculate the elements of the CI matrix
  SCFEnergy *eop = new SCFEnergy;
  eop->reference();
  RefSCElementOp2 op = eop;

  RefSymmSCMatrix t = focka_.result_noupdate().copy();
  t.accumulate(cl_hcore_);

  RefSymmSCMatrix d = cl_dens_.copy();
  d.accumulate(op_densa_);

  t.element_op(op, d);
  double h11 = eop->result();

  t.assign(fockb_.result_noupdate().copy());
  t.accumulate(cl_hcore_);

  d.assign(cl_dens_);
  d.accumulate(op_densb_);

  eop->reset();
  t.element_op(op, d);
  double h22 = eop->result();

  t = ka_.result_noupdate();
  eop->reset();
  t.element_op(op, op_densb_);
  double h21 = eop->result();

  t = kb_.result_noupdate();
  eop->reset();
  t.element_op(op, op_densa_);
  double h12 = eop->result();
  
  op=0;
  eop->dereference();
  delete eop;

  // now diagonalize the CI matrix to get the coefficients
  RefSCDimension l2 = new SCDimension(2);
  RefSCMatrixKit lkit = new LocalSCMatrixKit;
  RefSymmSCMatrix h = lkit->symmmatrix(l2);
  RefSCMatrix hv = lkit->matrix(l2,l2);
  RefDiagSCMatrix hl = lkit->diagmatrix(l2);
  
  h.set_element(0,0,h11);
  h.set_element(1,1,h22);
  h.set_element(1,0,h12);
  h.diagonalize(hl,hv);

  ci1_ = hv.get_element(0,0);
  ci2_ = hv.get_element(1,0);
  double c1c2 = ci1_*ci2_;

  if (scf_grp_->me()==0) {
    cout.setf(ios::fixed);
    cout << indent
         << "ci1 = " << setprecision(7) << ci1_ << " "
         << "ci2 = " << setprecision(7) << ci2_ << endl;
  }
  
  occa_ = 2*ci1_*ci1_;
  occb_ = 2*ci2_*ci2_;
  
  double eelec = 0.5*occa_*h11 + 0.5*occb_*h22 + 2.0*c1c2*h12;
  
  return eelec;
}

RefSCExtrapData
TCSCF::extrap_data()
{
  RefSymmSCMatrix *m = new RefSymmSCMatrix[4];
  m[0] = focka_.result_noupdate();
  m[1] = fockb_.result_noupdate();
  m[2] = ka_.result_noupdate();
  m[3] = kb_.result_noupdate();
  
  RefSCExtrapData data = new SymmSCMatrixNSCExtrapData(4, m);
  delete[] m;
  
  return data;
}

RefSymmSCMatrix
TCSCF::effective_fock()
{
  // use fock() instead of cl_fock_ just in case this is called from
  // someplace outside SCF::compute_vector()
  RefSymmSCMatrix mofocka = fock(0).clone();
  mofocka.assign(0.0);
  mofocka.accumulate_transform(scf_vector_.t(), fock(0));
  mofocka.scale(ci1_*ci1_);

  RefSymmSCMatrix mofockb = mofocka.clone();
  mofockb.assign(0.0);
  mofockb.accumulate_transform(scf_vector_.t(), fock(1));
  mofockb.scale(ci2_*ci2_);

  RefSymmSCMatrix moka = mofocka.clone();
  moka.assign(0.0);
  moka.accumulate_transform(scf_vector_.t(), fock(3));
  moka.scale(ci1_*ci2_);

  RefSymmSCMatrix mokb = mofocka.clone();
  mokb.assign(0.0);
  mokb.accumulate_transform(scf_vector_.t(), fock(4));
  mokb.scale(ci1_*ci2_);

  RefSymmSCMatrix mofock = mofocka.copy();
  mofock.accumulate(mofockb);

  BlockedSymmSCMatrix *F = BlockedSymmSCMatrix::castdown(mofock);
  BlockedSymmSCMatrix *Fa = BlockedSymmSCMatrix::castdown(mofocka);
  BlockedSymmSCMatrix *Fb = BlockedSymmSCMatrix::castdown(mofockb);
  BlockedSymmSCMatrix *Ka = BlockedSymmSCMatrix::castdown(moka);
  BlockedSymmSCMatrix *Kb = BlockedSymmSCMatrix::castdown(mokb);
  
  double scalea = (fabs(ci1_) < fabs(ci2_)) ? 1.0/(ci1_*ci1_ + 0.05) : 1.0;
  double scaleb = (fabs(ci2_) < fabs(ci1_)) ? 1.0/(ci2_*ci2_ + 0.05) : 1.0;

  for (int b=0; b < Fa->nblocks(); b++) {
    if (b==osa_) {
      RefSymmSCMatrix f = F->block(b);
      RefSymmSCMatrix fa = Fa->block(b);
      RefSymmSCMatrix fb = Fb->block(b);
      RefSymmSCMatrix kb = Kb->block(b);

      int i,j;

      i=ndocc_[b];
      for (j=0; j < ndocc_[b]; j++) 
        f->set_element(i,j,
                       scaleb*(fb->get_element(i,j)-kb->get_element(i,j)));

      j=ndocc_[b];
      for (i=ndocc_[b]+1; i < f->n(); i++)
        f->set_element(i,j,
                       scalea*(fa->get_element(i,j)+kb->get_element(i,j)));
      
    } else if (b==osb_) {
      RefSymmSCMatrix f = F->block(b);
      RefSymmSCMatrix fa = Fa->block(b);
      RefSymmSCMatrix fb = Fb->block(b);
      RefSymmSCMatrix ka = Ka->block(b);

      int i,j;

      double scale=1.0/(ci2_*ci2_ + 0.05);

      i=ndocc_[b];
      for (j=0; j < ndocc_[b]; j++) 
        f->set_element(i,j,
                       scalea*(fa->get_element(i,j)-ka->get_element(i,j)));

      j=ndocc_[b];
      for (i=ndocc_[b]+1; i < f->n(); i++)
        f->set_element(i,j,
                       scaleb*(fb->get_element(i,j)+ka->get_element(i,j)));
    }
  }

  return mofock;
}

//////////////////////////////////////////////////////////////////////////////

void
TCSCF::ao_fock()
{
  RefPetiteList pl = integral()->petite_list(basis());
  
  // calculate G.  First transform cl_dens_diff_ to the AO basis, then
  // scale the off-diagonal elements by 2.0
  RefSymmSCMatrix da = pl->to_AO_basis(cl_dens_diff_);
  RefSymmSCMatrix db = da.copy();
  RefSymmSCMatrix oda = pl->to_AO_basis(op_densa_diff_);
  RefSymmSCMatrix odb = pl->to_AO_basis(op_densb_diff_);
  da.accumulate(oda);
  db.accumulate(odb);

  da->scale(2.0);
  da->scale_diagonal(0.5);
  
  db->scale(2.0);
  db->scale_diagonal(0.5);
  
  oda->scale(2.0);
  oda->scale_diagonal(0.5);
  
  odb->scale(2.0);
  odb->scale_diagonal(0.5);
  
  // now try to figure out the matrix specialization we're dealing with
  // if we're using Local matrices, then there's just one subblock, or
  // see if we can convert G and P to local matrices
  if (local_ || local_dens_) {

    // grab the data pointers from the G and P matrices
    double *gmata, *gmatb, *kmata, *kmatb, *pmata, *pmatb, *opmata, *opmatb;
    RefSymmSCMatrix gatmp = get_local_data(ao_gmata_, gmata, SCF::Accum);
    RefSymmSCMatrix patmp = get_local_data(da, pmata, SCF::Read);
    RefSymmSCMatrix gbtmp = get_local_data(ao_gmatb_, gmatb, SCF::Accum);
    RefSymmSCMatrix pbtmp = get_local_data(db, pmatb, SCF::Read);
    RefSymmSCMatrix katmp = get_local_data(ao_ka_, kmata, SCF::Accum);
    RefSymmSCMatrix opatmp = get_local_data(oda, opmata, SCF::Read);
    RefSymmSCMatrix kbtmp = get_local_data(ao_kb_, kmatb, SCF::Accum);
    RefSymmSCMatrix opbtmp = get_local_data(odb, opmatb, SCF::Read);
    
    char * pmax = init_pmax(pmata);
    char * pmaxb = init_pmax(pmatb);
  
    for (int i=0; i < i_offset(basis()->nshell()); i++) {
      if (pmaxb[i] > pmax[i])
        pmax[i]=pmaxb[i];
    }
    
    delete[] pmaxb;
    
    LocalTCContribution lclc(gmata, pmata, gmatb, pmatb,
                             kmata, opmata, kmatb, opmatb);
    LocalGBuild<LocalTCContribution>
      gb(lclc, tbi_, integral(), basis(), scf_grp_, pmax);
    gb.build_gmat(desired_value_accuracy()/100.0);

    delete[] pmax;

    // if we're running on multiple processors, then sum the G matrices
    if (scf_grp_->n() > 1) {
      scf_grp_->sum(gmata, i_offset(basis()->nbasis()));
      scf_grp_->sum(gmatb, i_offset(basis()->nbasis()));
      scf_grp_->sum(kmata, i_offset(basis()->nbasis()));
      scf_grp_->sum(kmatb, i_offset(basis()->nbasis()));
    }
    
    // if we're running on multiple processors, or we don't have local
    // matrices, then accumulate gtmp back into G
    if (!local_ || scf_grp_->n() > 1) {
      ao_gmata_->convert_accumulate(gatmp);
      ao_gmatb_->convert_accumulate(gbtmp);
      ao_ka_->convert_accumulate(katmp);
      ao_kb_->convert_accumulate(kbtmp);
    }
  }

  // for now quit
  else {
    cerr << indent << "Cannot yet use anything but Local matrices\n";
    abort();
  }
  
  da=0;
  db=0;
  oda=0;
  odb=0;

  // now symmetrize the skeleton G matrix, placing the result in dd
  RefSymmSCMatrix skel_gmat = ao_gmata_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  pl->symmetrize(skel_gmat,focka_.result_noupdate());
  
  skel_gmat = ao_gmatb_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  pl->symmetrize(skel_gmat,fockb_.result_noupdate());
  
  skel_gmat = ao_ka_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  pl->symmetrize(skel_gmat,ka_.result_noupdate());
  
  skel_gmat = ao_kb_.copy();
  skel_gmat.scale(1.0/(double)pl->order());
  pl->symmetrize(skel_gmat,kb_.result_noupdate());
  
  // Fa = H+Ga
  focka_.result_noupdate().accumulate(cl_hcore_);

  // Fb = H+Gb
  fockb_.result_noupdate().accumulate(cl_hcore_);

  focka_.computed()=1;
  fockb_.computed()=1;
  ka_.computed()=1;
  kb_.computed()=1;
}

/////////////////////////////////////////////////////////////////////////////

void
TCSCF::init_gradient()
{
  // presumably the eigenvectors have already been computed by the time
  // we get here
  scf_vector_ = eigenvectors_.result_noupdate();
}

void
TCSCF::done_gradient()
{
  cl_dens_=0;
  op_densa_=0;
  op_densb_=0;
  scf_vector_ = 0;
}

/////////////////////////////////////////////////////////////////////////////

// MO lagrangian
//       c    o   v
//  c  |2*FC|2*FC|0|
//     -------------
//  o  |2*FC| FO |0|
//     -------------
//  v  | 0  |  0 |0|
//
RefSymmSCMatrix
TCSCF::lagrangian()
{
  RefSymmSCMatrix mofocka = focka_.result_noupdate().clone();
  mofocka.assign(0.0);
  mofocka.accumulate_transform(scf_vector_.t(), focka_.result_noupdate());
  mofocka.scale(ci1_*ci1_);

  RefSymmSCMatrix mofockb = mofocka.clone();
  mofockb.assign(0.0);
  mofockb.accumulate_transform(scf_vector_.t(), fockb_.result_noupdate());
  mofockb.scale(ci2_*ci2_);

  // FOa = c1^2*Fa + c1c2*Kb
  RefSymmSCMatrix moka = mofocka.clone();
  moka.assign(0.0);
  moka.accumulate_transform(scf_vector_.t(), kb_.result_noupdate());
  moka.scale(ci1_*ci2_);
  moka.accumulate(mofocka);

  // FOb = c1^2*Fb + c1c2*Ka
  RefSymmSCMatrix mokb = mofocka.clone();
  mokb.assign(0.0);
  mokb.accumulate_transform(scf_vector_.t(), ka_.result_noupdate());
  mokb.scale(ci1_*ci2_);
  mokb.accumulate(mofockb);

  BlockedSymmSCMatrix::castdown(moka)->block(osb_)->assign(0.0);
  BlockedSymmSCMatrix::castdown(mokb)->block(osa_)->assign(0.0);
  
  moka.accumulate(mokb);
  mokb=0;

  // FC = c1^2*Fa + c2^2*Fb
  mofocka.accumulate(mofockb);
  mofockb=0;
  
  RefSCElementOp2 op = new MOLagrangian(this);
  mofocka.element_op(op, moka);
  moka=0;
  mofocka.scale(2.0);

  // transform MO lagrangian to SO basis
  RefSymmSCMatrix so_lag(basis_dimension(), basis_matrixkit());
  so_lag.assign(0.0);
  so_lag.accumulate_transform(scf_vector_, mofocka);
  
  // and then from SO to AO
  RefPetiteList pl = integral()->petite_list();
  RefSymmSCMatrix ao_lag = pl->to_AO_basis(so_lag);

  ao_lag.scale(-1.0);

  return ao_lag;
}

RefSymmSCMatrix
TCSCF::gradient_density()
{
  cl_dens_ = basis_matrixkit()->symmmatrix(basis_dimension());
  op_densa_ = cl_dens_.clone();
  op_densb_ = cl_dens_.clone();
  
  cl_dens_.assign(0.0);
  op_densa_.assign(0.0);
  op_densb_.assign(0.0);
  
  RefSCElementOp op = new SCFDensity(this, scf_vector_, 2.0);
  cl_dens_.element_op(op);
  cl_dens_.scale(2.0);
  
  op = new SCFDensity(this, scf_vector_, occa_);
  op_densa_.element_op(op);
  op_densa_.scale(occa_);
  
  op = new SCFDensity(this, scf_vector_, occb_);
  op_densb_.element_op(op);
  op_densb_.scale(occb_);
  
  BlockedSymmSCMatrix::castdown(op_densa_)->block(osb_)->assign(0.0);
  BlockedSymmSCMatrix::castdown(op_densb_)->block(osa_)->assign(0.0);
  
  RefPetiteList pl = integral()->petite_list(basis());
  
  cl_dens_ = pl->to_AO_basis(cl_dens_);
  op_densa_ = pl->to_AO_basis(op_densa_);
  op_densb_ = pl->to_AO_basis(op_densb_);

  RefSymmSCMatrix tdens = cl_dens_.copy();
  tdens.accumulate(op_densa_);
  tdens.accumulate(op_densb_);

  op_densa_.scale(2.0/occa_);
  op_densb_.scale(2.0/occb_);
  
  return tdens;
}

/////////////////////////////////////////////////////////////////////////////

void
TCSCF::two_body_deriv(double * tbgrad)
{
  RefSCElementMaxAbs m = new SCElementMaxAbs;
  cl_dens_.element_op(m);
  double pmax = m->result();
  m=0;

  // now try to figure out the matrix specialization we're dealing with.
  // if we're using Local matrices, then there's just one subblock, or
  // see if we can convert P to local matrices

  if (local_ || local_dens_) {
    // grab the data pointers from the P matrices
    double *pmat, *pmata, *pmatb;
    RefSymmSCMatrix ptmp = get_local_data(cl_dens_, pmat, SCF::Read);
    RefSymmSCMatrix patmp = get_local_data(op_densa_, pmata, SCF::Read);
    RefSymmSCMatrix pbtmp = get_local_data(op_densb_, pmatb, SCF::Read);
  
    LocalTCGradContribution l(pmat,pmata,pmatb,ci1_,ci2_);
    LocalTBGrad<LocalTCGradContribution> tb(l, integral(), basis(), scf_grp_);
    tb.build_tbgrad(tbgrad, pmax, desired_gradient_accuracy());
  }

  // for now quit
  else {
    cerr << indent << "TCSCF::two_body_deriv: can't do gradient yet\n";
    abort();
  }
}

/////////////////////////////////////////////////////////////////////////////

void
TCSCF::init_hessian()
{
}

void
TCSCF::done_hessian()
{
}
