//
// uscf.cc --- implementation of the UnrestrictedSCF abstract base class
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>
#include <limits.h>
#include <sys/stat.h>

#include <util/misc/timer.h>
#include <util/misc/formio.h>

#include <math/scmat/local.h>
#include <math/scmat/repl.h>
#include <math/scmat/offset.h>
#include <math/scmat/block.h>
#include <math/scmat/blocked.h>
#include <math/scmat/blkiter.h>
#include <math/scmat/local.h>

#include <math/optimize/scextrapmat.h>
#include <math/optimize/diis.h>

#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/scf/scfops.h>
#include <chemistry/qc/scf/scflocal.h>
#include <chemistry/qc/scf/uscf.h>

///////////////////////////////////////////////////////////////////////////
// UnrestrictedSCF

#define CLASSNAME UnrestrictedSCF
#define PARENTS public SCF
#include <util/class/classia.h>
void *
UnrestrictedSCF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCF::_castdown(cd);
  return do_castdowns(casts,cd);
}

UnrestrictedSCF::UnrestrictedSCF(StateIn& s) :
  SCF(s),
  cb_(this),
  eb_(this),
  focka_(this),
  fockb_(this)
  maybe_SavableState(s)
{
  need_vec_ = 1;
  compute_guess_ = 0;

  cb_.result_noupdate() =
    basis_matrixkit()->matrix(basis_dimension(), basis_dimension());
  cb_.restore_state(s);
  cb_.result_noupdate().restore(s);

  eb_.result_noupdate() =
    basis_matrixkit()->diagmatrix(basis_dimension());
  eb_.restore_state(s);
  eb_.result_noupdate().restore(s);

  focka_.result_noupdate() =
    basis_matrixkit()->symmmatrix(basis_dimension());
  focka_.restore_state(s);
  focka_.result_noupdate().restore(s);

  fockb_.result_noupdate() =
    basis_matrixkit()->symmmatrix(basis_dimension());
  fockb_.restore_state(s);
  fockb_.result_noupdate().restore(s);

  s.get(user_occupations_);
  s.get(tnalpha_);
  s.get(tnbeta_);
  s.get(nirrep_);
  s.get(nalpha_);
  s.get(nbeta_);

  init_mem(4);
}

UnrestrictedSCF::UnrestrictedSCF(const RefKeyVal& keyval) :
  SCF(keyval),
  cb_(this),
  eb_(this),
  focka_(this),
  fockb_(this)
{
  int i;
  
  double acc = eigenvectors_.desired_accuracy();
  cb_.set_desired_accuracy(acc);
  eb_.set_desired_accuracy(acc);

  if (cb_.desired_accuracy() < DBL_EPSILON) {
    cb_.set_desired_accuracy(DBL_EPSILON);
    eb_.set_desired_accuracy(DBL_EPSILON);
  }

  focka_.compute()=0;
  focka_.computed()=0;
  fockb_.compute()=0;
  fockb_.computed()=0;

  // calculate the total nuclear charge
  int Znuc=0;
  PointBag_double *z = molecule()->charges();
  
  for (Pix p=z->first(); p; z->next(p)) Znuc += (int) z->get(p);
  delete z;

  // check to see if this is to be a charged molecule
  int charge = keyval->intvalue("total_charge");
  int nelectrons = Znuc-charge;

  // first let's try to figure out how many open shells there are
  if (keyval->exists("multiplicity")) {
    int mult = keyval->intvalue("multiplicity");
    if (mult < 1) {
      cerr << node0 << endl << indent
           << "USCF::init: bad value for multiplicity: " << mult << endl
           << indent << "assuming singlet" << endl;
      mult=1;
    }
    
    // for singlet, triplet, etc. we need an even number of electrons
    // for doublet, quartet, etc. we need an odd number of electrons
    if ((mult%2 && nelectrons%2) || (!(mult%2) && !(nelectrons%2))) {
      cerr << node0 << endl << indent
           << "USCF::init: Warning, there's a leftover electron..."
           << " I'm going to get rid of it" << endl
           << incindent << indent << "total_charge = " << charge << endl
           << indent << "total nuclear charge = " << Znuc << endl
           << indent << "multiplicity = " << mult << endl << decindent;
      nelectrons--;
    }
    if (mult%2)
      tnalpha_ = nelectrons/2 + (mult-1)/2;
    else
      tnalpha_ = nelectrons/2 + mult/2;

  } else {
    // if there's an odd number of electrons, then do a doublet, otherwise
    // do a triplet
    tnalpha_=nelectrons/2+1;
  }

  tnbeta_ = nelectrons-tnalpha_;
  
  cout << node0 << endl << indent
       << "USCF::init: total charge = " << Znuc-tnalpha_-tnbeta_
       << endl << endl;

  nirrep_ = molecule()->point_group().char_table().ncomp();

  if (keyval->exists("alpha") && keyval->exists("beta")) {
    nalpha_ = new int[nirrep_];
    nbeta_ = new int[nirrep_];
    tnalpha_ = 0;
    tnbeta_ = 0;
    user_occupations_=1;
    for (i=0; i < nirrep_; i++) {
      nalpha_[i] = keyval->intvalue("alpha",i);
      nbeta_[i] = keyval->intvalue("beta",i);
      tnalpha_ += nalpha_[i];
      tnbeta_ += nbeta_[i];
    }
  } else {
    nalpha_=0;
    nbeta_=0;
    user_occupations_=0;
    set_occupations(0,0);
  }

  cout << node0 << indent << "alpha = [";
  for (i=0; i < nirrep_; i++)
    cout << node0 << " " << nalpha_[i];
  cout << node0 << " ]\n";

  cout << node0 << indent << "beta  = [";
  for (i=0; i < nirrep_; i++)
    cout << node0 << " " << nbeta_[i];
  cout << node0 << " ]\n";

  // check to see if this was done in SCF(keyval)
  if (!keyval->exists("maxiter"))
    maxiter_ = 100;

  if (!keyval->exists("level_shift"))
    level_shift_ = 1.0;

  // now take care of memory stuff
  init_mem(4);
}

UnrestrictedSCF::~UnrestrictedSCF()
{
  if (nalpha_) {
    delete[] nalpha_;
    nalpha_=0;
  }
  if (nbeta_) {
    delete[] nbeta_;
    nbeta_=0;
  }
}

void
UnrestrictedSCF::save_data_state(StateOut& s)
{
  SCF::save_data_state(s);
  cb_.save_data_state(s);
  cb_.result_noupdate().save(s);
  eb_.save_data_state(s);
  eb_.result_noupdate().save(s);
  focka_.save_data_state(s);
  focka_.result_noupdate().save(s);
  fockb_.save_data_state(s);
  fockb_.result_noupdate().save(s);

  s.put(user_occupations_);
  s.put(tnalpha_);
  s.put(tnbeta_);
  s.put(nirrep_);
  s.put(nalpha_, nirrep_);
  s.put(nbeta_, nirrep_);
}

double
UnrestrictedSCF::occupation(int ir, int i)
{
  abort();
}

double
UnrestrictedSCF::alpha_occupation(int ir, int i)
{
  if (i < nalpha_[ir]) return 1.0;
  return 0.0;
}

double
UnrestrictedSCF::beta_occupation(int ir, int i)
{
  if (i < nbeta_[ir]) return 1.0;
  return 0.0;
}

RefSCMatrix
UnrestrictedSCF::eigenvectors()
{
  abort();
}

RefDiagSCMatrix
UnrestrictedSCF::eigenvalues()
{
  abort();
}

RefSCMatrix
UnrestrictedSCF::alpha_eigenvectors()
{
  return eigenvectors_.result();
}

RefDiagSCMatrix
UnrestrictedSCF::alpha_eigenvalues()
{
  return eigenvalues_.result();
}

RefSCMatrix
UnrestrictedSCF::beta_eigenvectors()
{
  return cb_.result();
}

RefDiagSCMatrix
UnrestrictedSCF::beta_eigenvalues()
{
  return eb_.result();
}

int
UnrestrictedSCF::spin_polarized()
{
  return 1;
}

int
UnrestrictedSCF::spin_unrestricted()
{
  return 1;
}

int
UnrestrictedSCF::n_fock_matrices() const
{
  return 2;
}

RefSymmSCMatrix
UnrestrictedSCF::fock(int n)
{
  if (n > 1) {
    cerr << node0 << indent
         << "USCF::fock: there are only two fock matrices, "
         << scprintf("but fock(%d) was requested\n",n);
    abort();
  }

  if (n==0)
    return focka_.result();
  else
    return fockb_.result();
}

void
UnrestrictedSCF::print(ostream&o)
{
  int i;

  SCF::print(o);
  o << node0 << indent << "UnrestrictedSCF Parameters:\n" << incindent
    << indent << "charge = " << molecule()->nuclear_charge()
                                - tnalpha_ - tnbeta_ << endl
    << indent << "nalpha = " << tnalpha_ << endl
    << indent << "nbeta = " << tnbeta_ << endl
    << indent << "alpha = [";

  for (i=0; i < nirrep_; i++)
    o << node0 << " " << nalpha_[i];
  o << node0 << " ]" << endl;

  o << node0 << indent << "beta  = [";
  for (i=0; i < nirrep_; i++)
    o << node0 << " " << nbeta_[i];
  o << node0 << " ]" << endl << decindent << endl;
}

//////////////////////////////////////////////////////////////////////////////

void
UnrestrictedSCF::initial_vector(int needv)
{
  if (need_vec_) {
    if (eigenvectors_.result_noupdate().null()) {
      // if guess_wfn_ is non-null then try to get a guess vector from it.
      // First check that the same basis is used...if not, then project the
      // guess vector into the present basis.
      // right now the check is crude...there should be an equiv member in
      // GaussianBasisSet
      if (guess_wfn_.nonnull()) {
        if (guess_wfn_->basis()->nbasis() == basis()->nbasis()) {
          cout << node0 << indent
               << "Using guess wavefunction as starting vector" << endl;

          // indent output of eigenvectors() call if there is any
          cout << incindent << incindent;
          UnrestrictedSCF *ug =
            UnrestrictedSCF::castdown(guess_wfn_.pointer());
          if (!ug || compute_guess_) {
            eigenvectors_ = guess_wfn_->alpha_eigenvectors().copy();
            eigenvalues_ = guess_wfn_->alpha_eigenvalues().copy();
            cb_ = guess_wfn_->beta_eigenvectors().copy();
            eb_ = guess_wfn_->beta_eigenvalues().copy();
          } else if (ug) {
            eigenvectors_ = ug->eigenvectors_.result_noupdate().copy();
            eigenvalues_ = ug->eigenvalues_.result_noupdate().copy();
            cb_ = ug->cb_.result_noupdate().copy();
            eb_ = ug->eb_.result_noupdate().copy();
            eigenvectors_.result_noupdate()->schmidt_orthog(
              overlap().pointer(), basis()->nbasis());
            cb_.result_noupdate()->schmidt_orthog(
              overlap().pointer(), basis()->nbasis());
          }
          cout << decindent << decindent;
        } else {
          cout << node0 << indent
               << "Projecting guess wavefunction into the present basis set"
               << endl;

          // indent output of projected_eigenvectors() call if there is any
          cout << incindent << incindent;
          eigenvectors_ = projected_eigenvectors(guess_wfn_, 1);
          eigenvalues_ = projected_eigenvalues(guess_wfn_, 1);
          cb_ = projected_eigenvectors(guess_wfn_, 0);
          eb_ = projected_eigenvalues(guess_wfn_, 0);
          cout << decindent << decindent;
        }

        // we should only have to do this once, so free up memory used
        // for the old wavefunction
        guess_wfn_=0;

        cout << node0 << endl;
      
      } else {
        cout << node0 << indent << "Starting from core Hamiltonian guess\n"
             << endl;
        eigenvectors_ = hcore_guess();
        eigenvalues_ = core_hamiltonian().eigvals();
        cb_ = eigenvectors_.result_noupdate().copy();
        eb_ = eigenvalues_.result_noupdate().copy();
      }
    } else {
      // this is just an old vector, so orthogonalize it
      eigenvectors_.result_noupdate()->schmidt_orthog(overlap().pointer(),
                                                      basis()->nbasis());
      cb_.result_noupdate()->schmidt_orthog(overlap().pointer(),
                                                      basis()->nbasis());
    }
  }

  need_vec_=needv;
}

//////////////////////////////////////////////////////////////////////////////

void
UnrestrictedSCF::set_occupations(const RefDiagSCMatrix& ev)
{
  abort();
}

void
UnrestrictedSCF::set_occupations(const RefDiagSCMatrix& eva,
                                 const RefDiagSCMatrix& evb)
{
  if (user_occupations_)
    return;
  
  if (nirrep_==1) {
    if (!nalpha_) {
      nalpha_=new int[1];
      nalpha_[0] = tnalpha_;
    }
    if (!nbeta_) {
      nbeta_=new int[1];
      nbeta_[0] = tnbeta_;
    }
    return;
  }
  
  int i,j;
  
  RefDiagSCMatrix evalsa, evalsb;
  
  if (eva.null()) {
    initial_vector(0);
    evalsa = eigenvalues_.result_noupdate();
    evalsb = eb_.result_noupdate();
  }
  else {
    evalsa = eva;
    evalsb = evb;
  }

  // first convert evals to something we can deal with easily
  BlockedDiagSCMatrix *bevalsa = BlockedDiagSCMatrix::require_castdown(evalsa,
                                "UnrestrictedSCF::set_occupations");
  BlockedDiagSCMatrix *bevalsb = BlockedDiagSCMatrix::require_castdown(evalsb,
                                "UnrestrictedSCF::set_occupations");
  
  RefPetiteList pl = integral()->petite_list(basis());
  
  double **valsa = new double*[nirrep_];
  double **valsb = new double*[nirrep_];
  for (i=0; i < nirrep_; i++) {
    int nf=pl->nfunction(i);
    if (nf) {
      valsa[i] = new double[nf];
      valsb[i] = new double[nf];
      bevalsa->block(i)->convert(valsa[i]);
      bevalsb->block(i)->convert(valsb[i]);
    } else {
      valsa[i] = 0;
      valsb[i] = 0;
    }
  }

  // now loop to find the tnalpha_ lowest eigenvalues and populate those
  // MO's
  int *newalpha = new int[nirrep_];
  memset(newalpha,0,sizeof(int)*nirrep_);

  for (i=0; i < tnalpha_; i++) {
    // find lowest eigenvalue
    int lir,ln;
    double lowest=999999999;

    for (int ir=0; ir < nirrep_; ir++) {
      int nf=pl->nfunction(ir);
      if (!nf)
        continue;
      for (j=0; j < nf; j++) {
        if (valsa[ir][j] < lowest) {
          lowest=valsa[ir][j];
          lir=ir;
          ln=j;
        }
      }
    }
    valsa[lir][ln]=999999999;
    newalpha[lir]++;
  }

  int *newbeta = new int[nirrep_];
  memset(newbeta,0,sizeof(int)*nirrep_);

  for (i=0; i < tnbeta_; i++) {
    // find lowest eigenvalue
    int lir,ln;
    double lowest=999999999;

    for (int ir=0; ir < nirrep_; ir++) {
      int nf=pl->nfunction(ir);
      if (!nf)
        continue;
      for (j=0; j < nf; j++) {
        if (valsb[ir][j] < lowest) {
          lowest=valsb[ir][j];
          lir=ir;
          ln=j;
        }
      }
    }
    valsb[lir][ln]=999999999;
    newbeta[lir]++;
  }

  // get rid of vals
  for (i=0; i < nirrep_; i++) {
    if (valsa[i])
      delete[] valsa[i];
    if (valsb[i])
      delete[] valsb[i];
  }
  delete[] valsa;
  delete[] valsb;

  if (!nalpha_) {
    nalpha_=newalpha;
    nbeta_=newbeta;
  } else {
    // test to see if newocc is different from nalpha_
    for (i=0; i < nirrep_; i++) {
      if (nalpha_[i] != newalpha[i]) {
        cerr << node0 << indent << "UnrestrictedSCF::set_occupations:  WARNING!!!!\n"
             << incindent << indent
             << scprintf("occupations for irrep %d have changed\n",i+1)
             << indent
             << scprintf("nalpha was %d, changed to %d", nalpha_[i], newalpha[i])
             << endl << decindent;
      }
      if (nbeta_[i] != newbeta[i]) {
        cerr << node0 << indent << "UnrestrictedSCF::set_occupations:  WARNING!!!!\n"
             << incindent << indent
             << scprintf("occupations for irrep %d have changed\n",i+1)
             << indent
             << scprintf("nbeta was %d, changed to %d", nbeta_[i], newbeta[i])
             << endl << decindent;
      }
    }

    memcpy(nalpha_,newalpha,sizeof(int)*nirrep_);
    memcpy(nbeta_,newbeta,sizeof(int)*nirrep_);
    delete[] newalpha;
    delete[] newbeta;
  }
}

//////////////////////////////////////////////////////////////////////////////
//
// scf things
//

void
UnrestrictedSCF::init_vector()
{
  // initialize the two electron integral classes
  tbi_ = integral()->electron_repulsion();
  tbi_->set_integral_storage(integral()->storage_unused());

  // calculate the core Hamiltonian
  hcore_ = core_hamiltonian();
  
  // allocate storage for other temp matrices
  densa_ = hcore_.clone();
  densa_.assign(0.0);
  
  diff_densa_ = hcore_.clone();
  diff_densa_.assign(0.0);

  densb_ = hcore_.clone();
  densb_.assign(0.0);
  
  diff_densb_ = hcore_.clone();
  diff_densb_.assign(0.0);

  // gmat is in AO basis
  gmata_ = basis()->matrixkit()->symmmatrix(basis()->basisdim());
  gmata_.assign(0.0);

  gmatb_ = gmata_.clone();
  gmatb_.assign(0.0);

  if (focka_.result_noupdate().null()) {
    focka_ = hcore_.clone();
    focka_.result_noupdate().assign(0.0);
    fockb_ = hcore_.clone();
    fockb_.result_noupdate().assign(0.0);
  }

  // set up trial vector
  initial_vector(1);
    
  scf_vector_ = eigenvectors_.result_noupdate();
  scf_vectorb_ = cb_.result_noupdate();
}

void
UnrestrictedSCF::done_vector()
{
  tbi_=0;
  
  hcore_ = 0;
  gmata_ = 0;
  densa_ = 0;
  diff_densa_ = 0;
  gmatb_ = 0;
  densb_ = 0;
  diff_densb_ = 0;

  scf_vector_ = 0;
  scf_vectorb_ = 0;
}

RefSymmSCMatrix
UnrestrictedSCF::alpha_density()
{
  RefSymmSCMatrix dens(basis_dimension(), basis_matrixkit());
  so_density(dens, 1.0, 1);
  return dens;
}

RefSymmSCMatrix
UnrestrictedSCF::beta_density()
{
  RefSymmSCMatrix dens(basis_dimension(), basis_matrixkit());
  so_density(dens, 1.0, 0);
  return dens;
}

void
UnrestrictedSCF::reset_density()
{
  gmata_.assign(0.0);
  diff_densa_.assign(densa_);

  gmatb_.assign(0.0);
  diff_densb_.assign(densb_);
}

double
UnrestrictedSCF::new_density()
{
  // copy current density into density diff and scale by -1.  later we'll
  // add the new density to this to get the density difference.
  diff_densa_.assign(densa_);
  diff_densa_.scale(-1.0);

  diff_densb_.assign(densb_);
  diff_densb_.scale(-1.0);

  so_density(densa_, 1.0, 1);
  so_density(densb_, 1.0, 0);

  diff_densa_.accumulate(densa_);
  diff_densb_.accumulate(densb_);

  RefSymmSCMatrix d = diff_densa_ + diff_densb_;

  RefSCElementScalarProduct sp(new SCElementScalarProduct);
  d.element_op(sp, d);
  d=0;
  
  double delta = sp->result();
  delta = sqrt(delta/i_offset(diff_densa_.n()));

  return delta;
}

RefSymmSCMatrix
UnrestrictedSCF::density()
{
  if (!density_.computed()) {
    RefSymmSCMatrix densa(basis_dimension(), basis_matrixkit());
    RefSymmSCMatrix densb(basis_dimension(), basis_matrixkit());
    so_density(densa, 1.0, 1);
    so_density(densb, 1.0, 0);
    densa.accumulate(densb);
    densb=0;
    
    density_ = densa;
    density_.computed() = 1;
  }

  return density_.result_noupdate();
}

double
UnrestrictedSCF::scf_energy()
{
  SCFEnergy *eop = new SCFEnergy;
  eop->reference();
  RefSCElementOp2 op = eop;
  focka_.result_noupdate().element_op(op, densa_);
  double ea = eop->result();
  
  eop->reset();
  fockb_.result_noupdate().element_op(op, densb_);
  double eb = eop->result();

  RefSymmSCMatrix denst = densa_+densb_;
  eop->reset();
  hcore_.element_op(op, denst);
  double ec = eop->result();
  denst=0;
  
  op=0;
  eop->dereference();
  delete eop;

  return ec+ea+eb;
}

RefSymmSCMatrix
UnrestrictedSCF::effective_fock()
{
  abort();
}

////////////////////////////////////////////////////////////////////////////

class UAExtrapErrorOp : public BlockedSCElementOp {
  private:
    UnrestrictedSCF *scf_;

  public:
    UAExtrapErrorOp(UnrestrictedSCF *s) : scf_(s) {}
    ~UAExtrapErrorOp() {}

    int has_side_effects() { return 1; }

    void process(SCMatrixBlockIter& bi) {
      int ir=current_block();
      
      for (bi.reset(); bi; bi++) {
        int i=bi.i();
        int j=bi.j();
        if (scf_->alpha_occupation(ir,i) == scf_->alpha_occupation(ir,j))
          bi.set(0.0);
      }
    }
};

class UBExtrapErrorOp : public BlockedSCElementOp {
  private:
    UnrestrictedSCF *scf_;

  public:
    UBExtrapErrorOp(UnrestrictedSCF *s) : scf_(s) {}
    ~UBExtrapErrorOp() {}

    int has_side_effects() { return 1; }

    void process(SCMatrixBlockIter& bi) {
      int ir=current_block();
      
      for (bi.reset(); bi; bi++) {
        int i=bi.i();
        int j=bi.j();
        if (scf_->beta_occupation(ir,i) == scf_->beta_occupation(ir,j))
          bi.set(0.0);
      }
    }
};

RefSCExtrapData
UnrestrictedSCF::extrap_data()
{
  RefSCExtrapData data =
    new SymmSCMatrix2SCExtrapData(focka_.result_noupdate(),
                                  fockb_.result_noupdate());
  return data;
}

RefSCExtrapError
UnrestrictedSCF::extrap_error()
{
  // form Error_a
  RefSymmSCMatrix moa = hcore_.clone();
  moa.assign(0.0);
  moa.accumulate_transform(scf_vector_, focka_.result_noupdate(),
                              SCMatrix::TransposeTransform);
  
  RefSCElementOp op = new UAExtrapErrorOp(this);
  moa.element_op(op);
  
  // form Error_b
  RefSymmSCMatrix mob = hcore_.clone();
  mob.assign(0.0);
  mob.accumulate_transform(scf_vectorb_, fockb_.result_noupdate(),
                              SCMatrix::TransposeTransform);
  
  op = new UBExtrapErrorOp(this);
  mob.element_op(op);

  RefSymmSCMatrix aoa = moa.clone();
  aoa.assign(0.0);
  aoa.accumulate_transform(scf_vector_, moa);
  
  RefSymmSCMatrix aob = moa;
  moa=0;
  aob.assign(0.0);
  aob.accumulate_transform(scf_vectorb_,mob);
  mob=0;

  aoa.accumulate(aob);
  aob=0;

  RefSCExtrapError error = new SymmSCMatrixSCExtrapError(aoa);
  aoa=0;

  return error;
}

///////////////////////////////////////////////////////////////////////////

void
UnrestrictedSCF::compute_vector(double& eelec)
{
  tim_enter("vector");
  int i;

  // reinitialize the extrapolation object
  extrap_->reinitialize();
  
  // create level shifter
  ALevelShift *alevel_shift = new ALevelShift(this);
  alevel_shift->reference();
  BLevelShift *blevel_shift = new BLevelShift(this);
  blevel_shift->reference();
  
  // set up subclass for vector calculation
  init_vector();
  
  // calculate the nuclear repulsion energy
  double nucrep = molecule()->nuclear_repulsion_energy();
  cout << node0 << indent
       << scprintf("nuclear repulsion energy = %20.15f", nucrep)
       << endl << endl;

  RefDiagSCMatrix evalsa(basis_dimension(), basis_matrixkit());
  RefDiagSCMatrix evalsb(basis_dimension(), basis_matrixkit());

  int iter;
  for (iter=0; iter < maxiter_; iter++) {
    // form the density from the current vector 
    tim_enter("density");
    double delta = new_density();
    tim_exit("density");
    
    // check convergence
    if (delta < desired_value_accuracy())
      break;

    // reset the density from time to time
    if (iter && !(iter%dens_reset_freq_))
      reset_density();
      
    // form the AO basis fock matrix
    tim_enter("fock");
    ao_fock();
    tim_exit("fock");

    // calculate the electronic energy
    eelec = scf_energy();
    cout << node0 << indent
         << scprintf("iter %5d energy = %20.15f delta = %10.5e",
                     iter+1, eelec+nucrep, delta)
         << endl;

    // now extrapolate the fock matrix
    tim_enter("extrap");
    RefSCExtrapData data = extrap_data();
    RefSCExtrapError error = extrap_error();
    extrap_->extrapolate(data,error);
    data=0;
    error=0;
    tim_exit("extrap");

    // diagonalize effective MO fock to get MO vector
    tim_enter("evals");

    RefSymmSCMatrix moa = hcore_.clone();
    moa.assign(0.0);
    moa.accumulate_transform(scf_vector_, focka_.result_noupdate(),
                             SCMatrix::TransposeTransform);
    
    RefSymmSCMatrix mob = hcore_.clone();
    mob.assign(0.0);
    mob.accumulate_transform(scf_vectorb_, fockb_.result_noupdate(),
                             SCMatrix::TransposeTransform);
    
    RefSCMatrix nvectora = scf_vector_.clone();
    RefSCMatrix nvectorb = scf_vector_.clone();
  
    // level shift effective fock
    alevel_shift->set_shift(level_shift_);
    moa.element_op(alevel_shift);
    blevel_shift->set_shift(level_shift_);
    mob.element_op(blevel_shift);
    
    moa.diagonalize(evalsa,nvectora);
    mob.diagonalize(evalsb,nvectorb);
    moa=0;
    mob=0;
    tim_exit("evals");

    // now un-level shift eigenvalues
    alevel_shift->set_shift(-level_shift_);
    evalsa.element_op(alevel_shift);
    blevel_shift->set_shift(-level_shift_);
    evalsb.element_op(blevel_shift);
    
    if (reset_occ_)
      set_occupations(evalsa, evalsb);

    // transform MO vector to AO basis
    scf_vector_ = scf_vector_ * nvectora;
    scf_vectorb_ = scf_vectorb_ * nvectorb;
    nvectora=0;
    nvectorb=0;
    
    // and orthogonalize vector
    tim_enter("schmidt");
    scf_vector_->schmidt_orthog(overlap().pointer(),basis()->nbasis());
    scf_vectorb_->schmidt_orthog(overlap().pointer(),basis()->nbasis());
    tim_exit("schmidt");
  }
      
  eigenvalues_ = evalsa;
  eigenvalues_.computed() = 1;
  evalsa = 0;
  
  eigenvectors_ = scf_vector_;
  eigenvectors_.computed() = 1;
  
  cb_ = scf_vectorb_;
  cb_.computed() = 1;

  eb_ = evalsb;
  eb_.computed() = 1;
  evalsb = 0;
  
  {
    // compute spin contamination
    RefSCMatrix Sab = scf_vector_.t() * overlap() * scf_vectorb_;
    //Sab.print("Sab");
    BlockedSCMatrix *pSab = BlockedSCMatrix::castdown(Sab);
    double s2=0;
    for (int ir=0; ir < nirrep_; ir++) {
      RefSCMatrix Sab_ir=pSab->block(0);
      if (Sab_ir.nonnull()) {
        for (i=0; i < nalpha_[ir]; i++)
          for (int j=0; j < nbeta_[ir]; j++)
            s2 += Sab_ir.get_element(i,j)*Sab_ir.get_element(i,j);
      }
    }

    double S2real = (double)(tnalpha_-tnbeta_)/2.;
    S2real = S2real*(S2real+1);
    double S2 = S2real + tnbeta_ - s2;

    cout << node0 << endl
         << indent << scprintf("<S^2>exact = %f", S2real) << endl
         << indent << scprintf("<S^2>      = %f", S2) << endl;
  }
  
  // now clean up
  done_vector();

  alevel_shift->dereference();
  delete alevel_shift;
  blevel_shift->dereference();
  delete blevel_shift;

  tim_exit("vector");
  //tim_print(0);
}

////////////////////////////////////////////////////////////////////////////

void
UnrestrictedSCF::init_gradient()
{
  // presumably the eigenvectors have already been computed by the time
  // we get here
  scf_vector_ = eigenvectors_.result_noupdate();
  scf_vectorb_ = cb_.result_noupdate();
}

void
UnrestrictedSCF::done_gradient()
{
  densa_=0;
  densb_=0;
  scf_vector_ = 0;
  scf_vectorb_ = 0;
}

/////////////////////////////////////////////////////////////////////////////

RefSymmSCMatrix
UnrestrictedSCF::lagrangian()
{
  RefDiagSCMatrix ea = alpha_eigenvalues().copy();
  RefDiagSCMatrix eb = beta_eigenvalues().copy();
  
  BlockedDiagSCMatrix *eab = BlockedDiagSCMatrix::castdown(ea);
  BlockedDiagSCMatrix *ebb = BlockedDiagSCMatrix::castdown(eb);

  RefPetiteList pl = integral()->petite_list(basis());

  for (int ir=0; ir < nirrep_; ir++) {
    RefDiagSCMatrix eair = eab->block(ir);
    RefDiagSCMatrix ebir = ebb->block(ir);

    if (eair.null())
      continue;

    int i;
    for (i=nalpha_[ir]; i < pl->nfunction(ir); i++)
      eair.set_element(i,0.0);
    for (i=nbeta_[ir]; i < pl->nfunction(ir); i++)
      ebir.set_element(i,0.0);
  }
  
  RefSymmSCMatrix la = basis_matrixkit()->symmmatrix(basis_dimension());
  la.assign(0.0);
  la.accumulate_transform(scf_vector_, ea);

  RefSymmSCMatrix lb = la.clone();
  lb.assign(0.0);
  lb.accumulate_transform(scf_vectorb_, eb);

  la.accumulate(lb);
  
  la = pl->to_AO_basis(la);
  la->scale(-1.0);

  return la;
}

RefSymmSCMatrix
UnrestrictedSCF::gradient_density()
{
  densa_ = basis_matrixkit()->symmmatrix(basis_dimension());
  densb_ = densa_.clone();
  
  so_density(densa_, 1.0, 1);
  so_density(densb_, 1.0, 0);
  
  RefPetiteList pl = integral()->petite_list(basis());
  
  densa_ = pl->to_AO_basis(densa_);
  densb_ = pl->to_AO_basis(densb_);

  RefSymmSCMatrix tdens = densa_.copy();
  tdens.accumulate(densb_);
  return tdens;
}

//////////////////////////////////////////////////////////////////////////////

void
UnrestrictedSCF::init_hessian()
{
}

void
UnrestrictedSCF::done_hessian()
{
}

//////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
