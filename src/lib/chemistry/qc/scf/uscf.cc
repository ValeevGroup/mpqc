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

#include <numeric>
#include <math.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <util/state/stateio.h>

#include <util/misc/regtime.h>
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
#include <chemistry/qc/scf/ltbgrad.h>
#include <chemistry/qc/scf/uhftmpl.h>
#include <chemistry/qc/wfn/femo.h>

using namespace std;
using namespace sc;

namespace sc {

///////////////////////////////////////////////////////////////////////////
// UnrestrictedSCF

static ClassDesc UnrestrictedSCF_cd(
  typeid(UnrestrictedSCF),"UnrestrictedSCF",2,"public SCF",
  0, 0, 0);

UnrestrictedSCF::UnrestrictedSCF(StateIn& s) :
  SavableState(s),
  SCF(s),
  oso_eigenvectors_beta_(this),
  eigenvalues_beta_(this),
  focka_(this),
  fockb_(this)
{
  compute_guess_ = 0;

  oso_eigenvectors_beta_.result_noupdate() =
    basis_matrixkit()->matrix(so_dimension(), oso_dimension());
  oso_eigenvectors_beta_.restore_state(s);
  oso_eigenvectors_beta_.result_noupdate().restore(s);

  eigenvalues_beta_.result_noupdate() =
    basis_matrixkit()->diagmatrix(oso_dimension());
  eigenvalues_beta_.restore_state(s);
  eigenvalues_beta_.result_noupdate().restore(s);

  focka_.result_noupdate() =
    basis_matrixkit()->symmmatrix(so_dimension());
  focka_.restore_state(s);
  focka_.result_noupdate().restore(s);

  fockb_.result_noupdate() =
    basis_matrixkit()->symmmatrix(so_dimension());
  fockb_.restore_state(s);
  fockb_.result_noupdate().restore(s);

  s.get(user_occupations_);
  s.get(tnalpha_);
  s.get(tnbeta_);
  s.get(nirrep_);
  s.get(nalpha_);
  s.get(nbeta_);

  if (s.version(::class_desc<UnrestrictedSCF>()) >= 2) {
    s.get(initial_nalpha_);
    s.get(initial_nbeta_);
    most_recent_pg_ << SavableState::restore_state(s);
  } else {
    initial_nalpha_ = new int[nirrep_];
    memcpy(initial_nalpha_, nalpha_, sizeof(int)*nirrep_);
    initial_nbeta_ = new int[nirrep_];
    memcpy(initial_nbeta_, nbeta_, sizeof(int)*nirrep_);
  }

  init_mem(4);
}

UnrestrictedSCF::UnrestrictedSCF(const Ref<KeyVal>& keyval) :
  SCF(keyval),
  oso_eigenvectors_beta_(this),
  eigenvalues_beta_(this),
  focka_(this),
  fockb_(this)
{
  int i;

  double acc = oso_eigenvectors_.desired_accuracy();
  oso_eigenvectors_beta_.set_desired_accuracy(acc);
  eigenvalues_beta_.set_desired_accuracy(acc);

  if (oso_eigenvectors_beta_.desired_accuracy() < DBL_EPSILON) {
    oso_eigenvectors_beta_.set_desired_accuracy(DBL_EPSILON);
    eigenvalues_beta_.set_desired_accuracy(DBL_EPSILON);
  }

  focka_.compute()=0;
  focka_.computed()=0;
  fockb_.compute()=0;
  fockb_.computed()=0;

  // calculate the total nuclear charge
  double Znuc=molecule()->total_charge();

  // check to see if this is to be a charged molecule
  double charge = keyval->doublevalue("total_charge");
  int nelectrons = (int)(Znuc-charge+1.0e-4);

  bool multiplicity_given = false;
  // first let's try to figure out how many open shells there are
  if (keyval->exists("multiplicity")) {
    int mult = keyval->intvalue("multiplicity");
    if (mult < 1) {
      ExEnv::err0() << endl << indent
           << "USCF::init: bad value for multiplicity: " << mult << endl
           << indent << "will guess automatically" << endl;
    }
    else {
      multiplicity_given = true;
    }

    // for singlet, triplet, etc. we need an even number of electrons
    // for doublet, quartet, etc. we need an odd number of electrons
    if ((mult%2 && nelectrons%2) || (!(mult%2) && !(nelectrons%2))) {
      ExEnv::err0() << endl << indent
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

  ExEnv::out0() << endl << indent
       << "USCF::init: total charge = " << Znuc-tnalpha_-tnbeta_
       << endl << endl;

  nirrep_ = molecule()->point_group()->char_table().ncomp();

  nalpha_ = read_occ(keyval, "alpha", nirrep_);
  nbeta_ = read_occ(keyval, "beta", nirrep_);
  if (nalpha_ && nbeta_) {
    tnalpha_ = 0;
    tnbeta_ = 0;
    user_occupations_=1;
    for (i=0; i < nirrep_; i++) {
      tnalpha_ += nalpha_[i];
      tnbeta_ += nbeta_[i];
    }
    initial_nalpha_ = new int[nirrep_];
    memcpy(initial_nalpha_, nalpha_, sizeof(int)*nirrep_);
    initial_nbeta_ = new int[nirrep_];
    memcpy(initial_nbeta_, nbeta_, sizeof(int)*nirrep_);
  }
  else if ((nalpha_ && (!nbeta_)) || ((!nalpha_) && nbeta_)) {
    ExEnv::out0() << "ERROR: USCF: only one of alpha and beta specified: "
                 << "give both or none" << endl;
    abort();
  }
  else {
    initial_nalpha_=0;
    initial_nbeta_=0;
    nalpha_=0;
    nbeta_=0;
    user_occupations_=0;
    // if multiplicity was not given, guess
    set_occupations(0,0, multiplicity_given ? false : true);
  }

  ExEnv::out0() << indent << "alpha = [";
  for (i=0; i < nirrep_; i++)
    ExEnv::out0() << " " << nalpha_[i];
  ExEnv::out0() << " ]\n";

  ExEnv::out0() << indent << "beta  = [";
  for (i=0; i < nirrep_; i++)
    ExEnv::out0() << " " << nbeta_[i];
  ExEnv::out0() << " ]\n";

  // check to see if this was done in SCF(keyval)
  if (!keyval->exists("maxiter"))
    maxiter_ = 100;

  if (!keyval->exists("level_shift"))
    level_shift_ = 0.25;

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
  delete[] initial_nalpha_;
  delete[] initial_nbeta_;
}

void
UnrestrictedSCF::save_data_state(StateOut& s)
{
  SCF::save_data_state(s);
  oso_eigenvectors_beta_.save_data_state(s);
  oso_eigenvectors_beta_.result_noupdate().save(s);
  eigenvalues_beta_.save_data_state(s);
  eigenvalues_beta_.result_noupdate().save(s);
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

  s.put(initial_nalpha_,initial_pg_->char_table().ncomp());
  s.put(initial_nbeta_,initial_pg_->char_table().ncomp());
  SavableState::save_state(most_recent_pg_.pointer(),s);
}

double
UnrestrictedSCF::occupation(int ir, int i)
{
  abort();
  return 0;
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
  return 0;
}

RefDiagSCMatrix
UnrestrictedSCF::eigenvalues()
{
  abort();
  return 0;
}

RefSCMatrix
UnrestrictedSCF::oso_alpha_eigenvectors()
{
  return oso_eigenvectors_.result();
}

RefSCMatrix
UnrestrictedSCF::alpha_eigenvectors()
{
  return so_to_orthog_so().t() * oso_eigenvectors_.result();
}

RefDiagSCMatrix
UnrestrictedSCF::alpha_eigenvalues()
{
  return eigenvalues_.result();
}

RefSCMatrix
UnrestrictedSCF::oso_beta_eigenvectors()
{
  return oso_eigenvectors_beta_.result();
}

RefSCMatrix
UnrestrictedSCF::beta_eigenvectors()
{
  return so_to_orthog_so().t() * oso_eigenvectors_beta_.result();
}

RefDiagSCMatrix
UnrestrictedSCF::beta_eigenvalues()
{
  return eigenvalues_beta_.result();
}

double
UnrestrictedSCF::magnetic_moment() const
{
  const int mm = std::accumulate(nalpha_, nalpha_+nirrep_, 0) -
                 std::accumulate(nbeta_, nbeta_+nirrep_, 0);
  return static_cast<double>(mm);
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
    ExEnv::err0() << indent
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
UnrestrictedSCF::print(ostream&o) const
{
  int i;

  SCF::print(o);
  o << indent << "UnrestrictedSCF Parameters:\n" << incindent
    << indent << "charge = " << molecule()->total_charge()
                                - tnalpha_ - tnbeta_ << endl
    << indent << "nalpha = " << tnalpha_ << endl
    << indent << "nbeta = " << tnbeta_ << endl
    << indent << "alpha = [";

  for (i=0; i < nirrep_; i++)
    o << " " << nalpha_[i];
  o << " ]" << endl;

  o << indent << "beta  = [";
  for (i=0; i < nirrep_; i++)
    o << " " << nbeta_[i];
  o << " ]" << endl << decindent << endl;
}

//////////////////////////////////////////////////////////////////////////////

void
UnrestrictedSCF::initial_vector()
{
  const bool vector_is_null = oso_eigenvectors_.result_noupdate() == 0;
  if (vector_is_null) {
    // if guess_wfn_ is non-null then try to get a guess vector from it.
    // First check that the same basis is used...if not, then project the
    // guess vector into the present basis.
    // right now the check is crude...there should be an equiv member in
    // GaussianBasisSet
    if (guess_wfn_) {
      if (basis()->equiv(guess_wfn_->basis())
          &&orthog_method() == guess_wfn_->orthog_method()
          &&oso_dimension()->equiv(guess_wfn_->oso_dimension().pointer())) {
        ExEnv::out0() << indent
                      << "Using " << guess_wfn_->class_name()
                      << " guess wavefunction as starting vector" << endl;

        // indent output of eigenvectors() call if there is any
        ExEnv::out0() << incindent << incindent;
        UnrestrictedSCF *ug =
          dynamic_cast<UnrestrictedSCF*>(guess_wfn_.pointer());
        if (!ug || compute_guess_) {
          oso_eigenvectors_ = guess_wfn_->oso_alpha_eigenvectors().copy();
          eigenvalues_ = guess_wfn_->alpha_eigenvalues().copy();
          oso_eigenvectors_beta_ = guess_wfn_->oso_beta_eigenvectors().copy();
          eigenvalues_beta_ = guess_wfn_->beta_eigenvalues().copy();
        } else if (ug) {
          oso_eigenvectors_ = ug->oso_eigenvectors_.result_noupdate().copy();
          eigenvalues_ = ug->eigenvalues_.result_noupdate().copy();
          oso_eigenvectors_beta_ = ug->oso_eigenvectors_beta_.result_noupdate().copy();
          eigenvalues_beta_ = ug->eigenvalues_beta_.result_noupdate().copy();
        }
        ExEnv::out0() << decindent << decindent;
      } else {
        ExEnv::out0() << indent
                      << "Projecting guess wavefunction ("
                      << guess_wfn_->class_name()
                      << ") into the present basis set"
                      << endl;

        // indent output of projected_eigenvectors() call if there is any
        ExEnv::out0() << incindent << incindent;
        oso_eigenvectors_ = projected_eigenvectors(guess_wfn_, 1);
        eigenvalues_ = projected_eigenvalues(guess_wfn_, 1);
        oso_eigenvectors_beta_ = projected_eigenvectors(guess_wfn_, 0);
        eigenvalues_beta_ = projected_eigenvalues(guess_wfn_, 0);
        ExEnv::out0() << decindent << decindent;
      }

      // we should only have to do this once, so free up memory used
      // for the old wavefunction, unless told otherwise
      if (!keep_guess_wfn_) guess_wfn_=0;

      ExEnv::out0() << endl;

    } else {
      ExEnv::out0() << indent << "Starting from core Hamiltonian guess\n"
                    << endl;
      oso_eigenvectors_ = hcore_guess(eigenvalues_.result_noupdate());
      oso_eigenvectors_beta_ = oso_eigenvectors_.result_noupdate().copy();
      eigenvalues_beta_ = eigenvalues_.result_noupdate().copy();
    }
  } else {
    // if vector exists do nothing
  }
} // end of initial_vector()

//////////////////////////////////////////////////////////////////////////////

void
UnrestrictedSCF::set_occupations(const RefDiagSCMatrix& ev)
{
  abort();
}

void
UnrestrictedSCF::set_occupations(const RefDiagSCMatrix& eva,
                                 const RefDiagSCMatrix& evb,
                                 bool can_change_multiplicity)
{
  if (user_occupations_ || (initial_nalpha_ && eva == 0)) {
    if (form_occupations(nalpha_, initial_nalpha_)) {
      form_occupations(nbeta_, initial_nbeta_);
      most_recent_pg_ = new PointGroup(molecule()->point_group());
      return;
    }
    ExEnv::out0() << indent
         << "UnrestrictedSCF: WARNING: reforming occupation vector from scratch" << endl;
  }

  int i,j;

  RefDiagSCMatrix evalsa, evalsb;

  if (eva == 0) {
    initial_vector();
    evalsa = eigenvalues_.result_noupdate();
    evalsb = eigenvalues_beta_.result_noupdate();
  }
  else {
    evalsa = eva;
    evalsb = evb;
  }

  // if can change the multiplicity, use HundsFEMOSeeker to get the FEMO with maximum multiplicity
  Ref<FEMO> femo;
  if (can_change_multiplicity) {
    //                            # electron,       Etol, allow closed-shell?
    HundsFEMOSeeker femoseeker(tnalpha_ + tnbeta_, HundsFEMOSeeker::tolerance, true,
                               evalsa, evalsb);
    femo = femoseeker.result();
  }
  else {
    femo = new FEMO(tnalpha_, tnbeta_,evalsa, evalsb);
  }
  // copy into local arrays
  int *newalpha = new int[nirrep_];
  int *newbeta = new int[nirrep_];
  tnalpha_ = tnbeta_ = 0;
  for(int g=0; g<nirrep_; ++g) {
    const int na = femo->nalpha(g);
    const int nb = femo->nbeta(g);
    newalpha[g] = na;
    newbeta[g] = nb;
    tnalpha_ += na;
    tnbeta_ += nb;
  }

  if (!nalpha_) {
    nalpha_=newalpha;
    nbeta_=newbeta;
  } else if (most_recent_pg_
             && most_recent_pg_->equiv(molecule()->point_group())) {
    // test to see if newocc is different from nalpha_
    for (i=0; i < nirrep_; i++) {
      if (nalpha_[i] != newalpha[i]) {
        ExEnv::err0() << indent << "UnrestrictedSCF::set_occupations:  WARNING!!!!\n"
             << incindent << indent
             << scprintf("occupations for irrep %d have changed\n",i+1)
             << indent
             << scprintf("nalpha was %d, changed to %d", nalpha_[i], newalpha[i])
             << endl << decindent;
      }
      if (nbeta_[i] != newbeta[i]) {
        ExEnv::err0() << indent << "UnrestrictedSCF::set_occupations:  WARNING!!!!\n"
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

  if (initial_pg_->equiv(molecule()->point_group())) {
    delete[] initial_nalpha_;
    initial_nalpha_ = new int[nirrep_];
    memcpy(initial_nalpha_,nalpha_,sizeof(int)*nirrep_);
  }

  if (initial_pg_->equiv(molecule()->point_group())) {
    delete[] initial_nbeta_;
    initial_nbeta_ = new int[nirrep_];
    memcpy(initial_nbeta_,nbeta_,sizeof(int)*nirrep_);
  }

  most_recent_pg_ = new PointGroup(molecule()->point_group());
}

void
UnrestrictedSCF::symmetry_changed()
{
  SCF::symmetry_changed();
  nirrep_ = molecule()->point_group()->char_table().ncomp();
  oso_eigenvectors_beta_.result_noupdate() = 0;
  eigenvalues_beta_.result_noupdate() = 0;
  focka_.result_noupdate() = 0;
  fockb_.result_noupdate() = 0;
  // do not change multiplicity
  set_occupations(0,0,false);
}

//////////////////////////////////////////////////////////////////////////////
//
// scf things
//

void
UnrestrictedSCF::init_vector()
{
  init_threads();

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

  if (focka_.result_noupdate() == 0) {
    focka_ = hcore_.clone();
    focka_.result_noupdate().assign(0.0);
    fockb_ = hcore_.clone();
    fockb_.result_noupdate().assign(0.0);
  }

  // make sure trial vector is set up
  initial_vector();

  oso_scf_vector_ = oso_eigenvectors_.result_noupdate();
  oso_scf_vector_beta_ = oso_eigenvectors_beta_.result_noupdate();
}

void
UnrestrictedSCF::done_vector()
{
  done_threads();

  hcore_ = 0;
  gmata_ = 0;
  densa_ = 0;
  diff_densa_ = 0;
  gmatb_ = 0;
  densb_ = 0;
  diff_densb_ = 0;

  oso_scf_vector_ = 0;
  oso_scf_vector_beta_ = 0;
}

RefSymmSCMatrix
UnrestrictedSCF::alpha_density()
{
  RefSymmSCMatrix dens(so_dimension(), basis_matrixkit());
  so_density(dens, 1.0, 1);
  return dens;
}

RefSymmSCMatrix
UnrestrictedSCF::beta_density()
{
  RefSymmSCMatrix dens(so_dimension(), basis_matrixkit());
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

  Ref<SCElementScalarProduct> sp(new SCElementScalarProduct);
  d.element_op(sp.pointer(), d);
  d=0;

  double delta = sp->result();
  delta = sqrt(delta/i_offset(diff_densa_.n()));

  return delta;
}

RefSymmSCMatrix
UnrestrictedSCF::density()
{
  if (!density_.computed()) {
    RefSymmSCMatrix densa(so_dimension(), basis_matrixkit());
    RefSymmSCMatrix densb(so_dimension(), basis_matrixkit());
    so_density(densa, 1.0, 1);
    so_density(densb, 1.0, 0);
    densa.accumulate(densb);
    densb=0;

    density_ = densa;
    // only flag the density as computed if the calc is converged
    if (!value_needed()) density_.computed() = 1;
  }

  return density_.result_noupdate();
}

double
UnrestrictedSCF::scf_energy()
{
  SCFEnergy *eop = new SCFEnergy;
  eop->reference();
  Ref<SCElementOp2> op = eop;
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
  return 0;
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

Ref<SCExtrapData>
UnrestrictedSCF::extrap_data()
{
  Ref<SCExtrapData> data =
    new SymmSCMatrix2SCExtrapData(focka_.result_noupdate(),
                                  fockb_.result_noupdate());
  return data;
}

Ref<SCExtrapError>
UnrestrictedSCF::extrap_error()
{
  RefSCMatrix so_to_ortho_so_tr = so_to_orthog_so().t();

  // form Error_a
  RefSymmSCMatrix moa(oso_dimension(), basis_matrixkit());
  moa.assign(0.0);
  moa.accumulate_transform(so_to_ortho_so_tr * oso_scf_vector_,
                           focka_.result_noupdate(),
                           SCMatrix::TransposeTransform);

  Ref<SCElementOp> op = new UAExtrapErrorOp(this);
  moa.element_op(op.pointer());

  // form Error_b
  RefSymmSCMatrix mob(oso_dimension(), basis_matrixkit());
  mob.assign(0.0);
  mob.accumulate_transform(so_to_ortho_so_tr * oso_scf_vector_beta_,
                           fockb_.result_noupdate(),
                           SCMatrix::TransposeTransform);

  op = new UBExtrapErrorOp(this);
  mob.element_op(op);

  RefSymmSCMatrix aoa(so_dimension(), basis_matrixkit());
  aoa.assign(0.0);
  aoa.accumulate_transform(so_to_ortho_so_tr * oso_scf_vector_, moa);
  moa = 0;

  RefSymmSCMatrix aob(so_dimension(), basis_matrixkit());
  aob.assign(0.0);
  aob.accumulate_transform(so_to_ortho_so_tr * oso_scf_vector_beta_,mob);
  mob=0;

  aoa.accumulate(aob);
  aob=0;

  Ref<SCExtrapError> error = new SymmSCMatrixSCExtrapError(aoa);
  aoa=0;

  return error;
}

///////////////////////////////////////////////////////////////////////////

double
UnrestrictedSCF::compute_vector(double& eelec, double nucrep)
{
  Timer tim("vector");
  int i;

    // reinitialize the extrapolation object
  extrap_->reinitialize();

  // create level shifter
  ALevelShift *alevel_shift = new ALevelShift(this);
  alevel_shift->reference();
  BLevelShift *blevel_shift = new BLevelShift(this);
  blevel_shift->reference();

  // calculate the core Hamiltonian
  hcore_ = core_hamiltonian();

  // add density independant contributions to Hcore
  accumdih_->accum(hcore_);

  // set up subclass for vector calculation
  init_vector();

  RefDiagSCMatrix evalsa(oso_dimension(), basis_matrixkit());
  RefDiagSCMatrix evalsb(oso_dimension(), basis_matrixkit());

  double delta = 1.0;
  int iter, iter_since_reset = 0;
  double accuracy = 1.0;

  ExEnv::out0() << indent
                << "Beginning iterations.  Basis is "
                << basis()->label() << '.' << std::endl;
  for (iter=0; iter < maxiter_; iter++, iter_since_reset++) {
    // form the density from the current vector
    tim.enter("density");
    delta = new_density();
    tim.exit("density");

    // reset the density from time to time
    if (iter_since_reset && !(iter_since_reset%dens_reset_freq_)) {
      reset_density();
      iter_since_reset = 0;
    }

    // form the AO basis fock matrix
    tim.enter("fock");
    double base_accuracy = delta;
    if (base_accuracy < desired_value_accuracy())
      base_accuracy = desired_value_accuracy();
    double new_accuracy = 0.01 * base_accuracy;
    if (new_accuracy > 0.001) new_accuracy = 0.001;
    if (iter == 0) accuracy = new_accuracy;
    else if (new_accuracy < accuracy) {
      accuracy = new_accuracy/10.0;
      if (iter_since_reset > 0) {
        reset_density();
        iter_since_reset = 0;
      }
    }
    ao_fock(accuracy);
    tim.exit("fock");

    // calculate the electronic energy
    eelec = scf_energy();
    ExEnv::out0() << indent
         << scprintf("iter %5d energy = %15.10f delta = %10.5e",
                     iter+1, eelec+nucrep, delta)
         << endl;

    // check convergence
    if (delta < desired_value_accuracy()
        && iter+1 >= miniter_)
      break;

    // now extrapolate the fock matrix
    tim.enter("extrap");
    Ref<SCExtrapData> data = extrap_data();
    Ref<SCExtrapError> error = extrap_error();
    extrap_->extrapolate(data,error);
    data=0;
    error=0;
    tim.exit("extrap");

    // diagonalize effective MO fock to get MO vector
    tim.enter("evals");

    RefSCMatrix so_to_oso_tr = so_to_orthog_so().t();

    RefSymmSCMatrix moa(oso_dimension(), basis_matrixkit());
    moa.assign(0.0);
    moa.accumulate_transform(so_to_oso_tr * oso_scf_vector_,
                             focka_.result_noupdate(),
                             SCMatrix::TransposeTransform);

    RefSymmSCMatrix mob(oso_dimension(), basis_matrixkit());
    mob.assign(0.0);
    mob.accumulate_transform(so_to_oso_tr * oso_scf_vector_beta_,
                             fockb_.result_noupdate(),
                             SCMatrix::TransposeTransform);

    RefSCMatrix nvectora(oso_dimension(), oso_dimension(), basis_matrixkit());
    RefSCMatrix nvectorb(oso_dimension(), oso_dimension(), basis_matrixkit());

    // level shift effective fock in the mo basis
    alevel_shift->set_shift(level_shift_);
    moa.element_op(alevel_shift);
    blevel_shift->set_shift(level_shift_);
    mob.element_op(blevel_shift);

    // transform back to the oso basis to do the diagonalization
    RefSymmSCMatrix osoa(oso_dimension(), basis_matrixkit());
    osoa.assign(0.0);
    osoa.accumulate_transform(oso_scf_vector_,moa);
    moa = 0;
    osoa.diagonalize(evalsa,oso_scf_vector_);
    osoa = 0;

    RefSymmSCMatrix osob(oso_dimension(), basis_matrixkit());
    osob.assign(0.0);
    osob.accumulate_transform(oso_scf_vector_beta_,mob);
    mob = 0;
    osob.diagonalize(evalsb,oso_scf_vector_beta_);
    osob = 0;

    tim.exit("evals");

    // now un-level shift eigenvalues
    alevel_shift->set_shift(-level_shift_);
    evalsa.element_op(alevel_shift);
    blevel_shift->set_shift(-level_shift_);
    evalsb.element_op(blevel_shift);

    if (reset_occ_)
      // Maintain multiplicity
      set_occupations(evalsa, evalsb, false);

    savestate_iter(iter);
    }

  eigenvalues_ = evalsa;
  eigenvalues_.computed() = 1;
  eigenvalues_.set_actual_accuracy(delta);
  evalsa = 0;

  oso_eigenvectors_ = oso_scf_vector_;
  oso_eigenvectors_.computed() = 1;
  oso_eigenvectors_.set_actual_accuracy(delta);

  oso_eigenvectors_beta_ = oso_scf_vector_beta_;
  oso_eigenvectors_beta_.computed() = 1;
  oso_eigenvectors_beta_.set_actual_accuracy(delta);

  eigenvalues_beta_ = evalsb;
  eigenvalues_beta_.computed() = 1;
  eigenvalues_beta_.set_actual_accuracy(delta);
  evalsb = 0;

  {
    // compute spin contamination
    RefSCMatrix so_to_oso_tr = so_to_orthog_so().t();
    RefSCMatrix Sab
      = (so_to_oso_tr * oso_scf_vector_).t()
      * overlap()
      * (so_to_oso_tr * oso_scf_vector_beta_);
    //Sab.print("Sab");
    BlockedSCMatrix *pSab = dynamic_cast<BlockedSCMatrix*>(Sab.pointer());
    double s2=0;
    for (int ir=0; ir < nirrep_; ir++) {
      RefSCMatrix Sab_ir=pSab->block(0);
      if (Sab_ir) {
        for (i=0; i < nalpha_[ir]; i++)
          for (int j=0; j < nbeta_[ir]; j++)
            s2 += Sab_ir.get_element(i,j)*Sab_ir.get_element(i,j);
      }
    }

    double S2real = (double)(abs(tnalpha_-tnbeta_))/2.;
    S2real = S2real*(S2real+1);
    double S2 = S2real + min(tnalpha_,tnbeta_) - s2;

    ExEnv::out0() << endl
         << indent << scprintf("<S^2>exact = %f", S2real) << endl
         << indent << scprintf("<S^2>      = %f", S2) << endl;
  }

  // now clean up
  done_vector();

  alevel_shift->dereference();
  delete alevel_shift;
  blevel_shift->dereference();
  delete blevel_shift;

  tim.exit("vector");
  //tim.print();

  return delta;
}

////////////////////////////////////////////////////////////////////////////

void
UnrestrictedSCF::init_gradient()
{
  // presumably the eigenvectors have already been computed by the time
  // we get here
  oso_scf_vector_ = oso_eigenvectors_.result_noupdate();
  oso_scf_vector_beta_ = oso_eigenvectors_beta_.result_noupdate();
}

void
UnrestrictedSCF::done_gradient()
{
  densa_=0;
  densb_=0;
  oso_scf_vector_ = 0;
  oso_scf_vector_beta_ = 0;
}

/////////////////////////////////////////////////////////////////////////////

RefSymmSCMatrix
UnrestrictedSCF::lagrangian()
{
  RefSCMatrix so_to_oso_tr = so_to_orthog_so().t();

  RefDiagSCMatrix ea = eigenvalues_.result_noupdate().copy();
  RefDiagSCMatrix eb = eigenvalues_beta_.result_noupdate().copy();

  BlockedDiagSCMatrix *eab = dynamic_cast<BlockedDiagSCMatrix*>(ea.pointer());
  BlockedDiagSCMatrix *ebb = dynamic_cast<BlockedDiagSCMatrix*>(eb.pointer());

  Ref<PetiteList> pl = integral()->petite_list(basis());

  for (int ir=0; ir < nirrep_; ir++) {
    RefDiagSCMatrix eair = eab->block(ir);
    RefDiagSCMatrix ebir = ebb->block(ir);

    if (eair == 0)
      continue;

    int i;
    for (i=nalpha_[ir]; i < eair.dim().n(); i++)
      eair.set_element(i,0.0);
    for (i=nbeta_[ir]; i < ebir.dim().n(); i++)
      ebir.set_element(i,0.0);
  }

  RefSymmSCMatrix la = basis_matrixkit()->symmmatrix(so_dimension());
  la.assign(0.0);
  la.accumulate_transform(so_to_oso_tr * oso_scf_vector_, ea);

  RefSymmSCMatrix lb = la.clone();
  lb.assign(0.0);
  lb.accumulate_transform(so_to_oso_tr * oso_scf_vector_beta_, eb);

  la.accumulate(lb);

  la = pl->to_AO_basis(la);
  la->scale(-1.0);

  return la;
}

RefSymmSCMatrix
UnrestrictedSCF::gradient_density()
{
  densa_ = basis_matrixkit()->symmmatrix(so_dimension());
  densb_ = densa_.clone();

  so_density(densa_, 1.0, 1);
  so_density(densb_, 1.0, 0);

  Ref<PetiteList> pl = integral()->petite_list(basis());

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

void
UnrestrictedSCF::two_body_deriv_hf(double * tbgrad, double exchange_fraction)
{
  Ref<SCElementMaxAbs> m = new SCElementMaxAbs;
  densa_.element_op(m.pointer());
  double pmax = m->result();
  m=0;

  // now try to figure out the matrix specialization we're dealing with.
  // if we're using Local matrices, then there's just one subblock, or
  // see if we can convert P to local matrices

  if (local_ || local_dens_) {
    // grab the data pointers from the P matrices
    double *pmata, *pmatb;
    RefSymmSCMatrix ptmpa = get_local_data(densa_, pmata, SCF::Read);
    RefSymmSCMatrix ptmpb = get_local_data(densb_, pmatb, SCF::Read);

    Ref<PetiteList> pl = integral()->petite_list();
    LocalUHFGradContribution l(pmata,pmatb);

    int i;
    int na3 = molecule()->natom()*3;
    int nthread = threadgrp_->nthread();
    double **grads = new double*[nthread];
    Ref<TwoBodyDerivInt> *tbis = new Ref<TwoBodyDerivInt>[nthread];
    for (i=0; i < nthread; i++) {
      tbis[i] = integral()->electron_repulsion_deriv();
      grads[i] = new double[na3];
      memset(grads[i], 0, sizeof(double)*na3);
    }

    LocalTBGrad<LocalUHFGradContribution> **tblds =
      new LocalTBGrad<LocalUHFGradContribution>*[nthread];

    for (i=0; i < nthread; i++) {
      tblds[i] = new LocalTBGrad<LocalUHFGradContribution>(
        l, tbis[i], pl, basis(), scf_grp_, grads[i], pmax,
        desired_gradient_accuracy(), nthread, i, exchange_fraction);
      threadgrp_->add_thread(i, tblds[i]);
    }

    if (threadgrp_->start_threads() < 0
        ||threadgrp_->wait_threads() < 0) {
      ExEnv::err0() << indent
           << "USCF: error running threads" << endl;
      abort();
    }

    for (i=0; i < nthread; i++) {
      for (int j=0; j < na3; j++)
        tbgrad[j] += grads[i][j];

      delete[] grads[i];
      delete tblds[i];
      tbis[i] = 0;
    }
    delete[] tbis;

    scf_grp_->sum(tbgrad,3 * basis()->molecule()->natom());
  }

  // for now quit
  else {
    ExEnv::err0() << indent
         << "USCF::two_body_deriv_hf: can't do gradient yet\n";
    abort();
  }
}

void
UnrestrictedSCF::set_desired_value_accuracy(double eps)
{
  OneBodyWavefunction::set_desired_value_accuracy(eps);
  oso_eigenvectors_beta_.set_desired_accuracy(eps);
  eigenvalues_beta_.set_desired_accuracy(eps);
}

void
UnrestrictedSCF::obsolete_vector() {
  SCF::obsolete_vector();
  oso_eigenvectors_beta_ = RefSCMatrix(0);
}

//////////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
