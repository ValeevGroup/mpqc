//
// molfreq.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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

#include <util/misc/math.h>
#include <util/misc/scexception.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <util/group/message.h>
#include <math/symmetry/corrtab.h>
#include <math/scmat/local.h>
#include <math/scmat/blocked.h>
#include <chemistry/molecule/molfreq.h>
#include <chemistry/molecule/molrender.h>

using namespace std;
using namespace sc;

#undef DEBUG

static ClassDesc MolecularFrequencies_cd(
  typeid(MolecularFrequencies),"MolecularFrequencies",3,"public SavableState",
  0, create<MolecularFrequencies>, create<MolecularFrequencies>);

MolecularFrequencies::MolecularFrequencies(const Ref<KeyVal>& keyval)
{
  mol_ << keyval->describedclassvalue("molecule");
  if (mol_ == 0) {
      throw InputError("missing required input of type Molecule",
                       __FILE__, __LINE__, "molecule", 0,
                       class_desc());
    }
  KeyValValueRefDescribedClass def_pg(mol_->point_group().pointer());
  pg_ << keyval->describedclassvalue("point_group", def_pg);
  nirrep_ = pg_->char_table().nirrep();
  debug_ = keyval->booleanvalue("debug");
  nfreq_ = 0;
  freq_ = 0;
}

MolecularFrequencies::MolecularFrequencies(const Ref<Molecule>& mol) :
    mol_(mol),
    pg_(mol->point_group()),
    debug_(0),
    nirrep_(mol->point_group()->char_table().nirrep()),
    nfreq_(0), freq_(0)
{
}

MolecularFrequencies::~MolecularFrequencies()
{
  delete[] nfreq_;
  if (freq_) {
      for (int i=0; i<nirrep_; i++) {
          delete[] freq_[i];
        }
      delete[] freq_;
    }
}

MolecularFrequencies::MolecularFrequencies(StateIn& si):
  SavableState(si)
{
  int i;

  if (si.version(::class_desc<MolecularFrequencies>()) < 3) {
      throw FileOperationFailed("cannot restore from old version",
                                __FILE__, __LINE__, 0,
                                FileOperationFailed::Corrupt,
                                class_desc());
    }

  mol_ << SavableState::restore_state(si);
  pg_ << SavableState::restore_state(si);

  si.get(nirrep_);
  si.get(nfreq_);
  for (i=0; i<nirrep_; i++) si.get(freq_[i]);
}

void
MolecularFrequencies::save_data_state(StateOut& so)
{
  int i;

  SavableState::save_state(mol_.pointer(),so);
  SavableState::save_state(pg_.pointer(),so);

  so.put(nirrep_);
  so.put(nfreq_,nirrep_);
  for (i=0; i<nirrep_; i++) so.put(freq_[i],nfreq_[i]);
}

void
MolecularFrequencies::compute_frequencies(const RefSymmSCMatrix &xhessian)
{
  int i, coor;

  RefSCMatrix symmbasis
      = MolecularHessian::cartesian_to_symmetry(mol_,pg_);
  BlockedSCMatrix *bsymmbasis = dynamic_cast<BlockedSCMatrix*>(symmbasis.pointer());

  kit_ = xhessian->kit();
  d3natom_ = xhessian->dim();
  symkit_ = symmbasis->kit();
  bd3natom_ = symmbasis->coldim();
  disym_ = symmbasis->rowdim();

  ExEnv::out0() << endl
               << indent << "Frequencies (cm-1; negative is imaginary):"
               << endl;

  // initialize the frequency tables
  if (nfreq_) delete[] nfreq_;
  nfreq_ = new int[nirrep_];
  if (freq_) delete[] freq_;
  freq_ = new double*[nirrep_];

  // initialize normal coordinate matrix
  normco_ = symmatrixkit()->matrix(bd3natom_, disym_);

  // find the inverse sqrt mass matrix
  RefDiagSCMatrix m(d3natom_, matrixkit());
  for (i=0,coor=0; i<mol_->natom(); i++) {
      for (int j=0; j<3; j++, coor++) {
          m(coor) = 1.0/sqrt(mol_->mass(i)*(1.0/5.48579903e-4));
        }
    }

  RefSymmSCMatrix dhessian;

  for (int irrep=0; irrep<nirrep_; irrep++) {
      RefSCMatrix dtranst = bsymmbasis->block(irrep);
      RefSCDimension ddim = dtranst.rowdim();
      nfreq_[irrep] = ddim.n();
      freq_[irrep] = new double[nfreq_[irrep]];
      if (ddim.n() == 0) continue;
      dhessian = matrixkit()->symmmatrix(ddim);
      dhessian.assign(0.0);
      dhessian.accumulate_transform(dtranst,xhessian);
      do_freq_for_irrep(irrep, m, dhessian, dtranst);
    }
}

void
MolecularFrequencies::do_freq_for_irrep(
    int irrep,
    const RefDiagSCMatrix &m,
    const RefSymmSCMatrix &dhessian,
    const RefSCMatrix &dtranst)
{
  int i;
  RefSCMatrix dtrans = dtranst.t();
  RefSCDimension ddim = dtrans.coldim();
  if (ddim.n() == 0) return;
  if (debug_) {
      dhessian.print("dhessian");
      dtrans.print("dtrans");
    }
  // find the basis for the normal coordinates
  RefSCMatrix ncbasis = m * dtrans;
  // use the SVD to orthogonalize and check this basis
  RefSCMatrix basU(d3natom_, d3natom_, matrixkit());
  RefSCMatrix basV(ddim, ddim, matrixkit());
  RefDiagSCMatrix bassigma(ddim, matrixkit());
  ncbasis.svd(basU, bassigma, basV);
  for (i=0; i<ddim.n(); i++) {
      if (bassigma(i) < 1.e-3) {
          throw ToleranceExceeded("singular value too small: "
                                  "displacements don't span coordinates",
                                  __FILE__, __LINE__, 1.e-3, bassigma(i),
                                  class_desc());
        }
    }
  ncbasis.assign_subblock(basU, 0, d3natom_.n()-1, 0, ddim.n()-1, 0, 0);
  // a transform from disp to x to q (mass weighted x) to disp
  RefSCMatrix dxqd = ncbasis.t() * m * dtrans;
  // transform the dhessian to the mass weighted dhessian
  RefSymmSCMatrix mdhessian = matrixkit()->symmmatrix(dxqd.rowdim());
  mdhessian.assign(0.0);
  mdhessian.accumulate_transform(dxqd, dhessian);
  if (debug_) {
      mdhessian.print("mass weighted dhessian");
    }
  // diagonalize the hessian
  RefDiagSCMatrix freqs(ddim,matrixkit());
  RefSCMatrix eigvecs(ddim,ddim,matrixkit());
  mdhessian.diagonalize(freqs,eigvecs);
  // convert the eigvals to frequencies in wavenumbers
  for (i=0; i<freqs.n(); i++) {
      if (freqs(i) >=0.0) freqs(i) = sqrt(freqs(i));
      else freqs(i) = -sqrt(-freqs(i));
      freq_[irrep][i] = freqs(i);
      freqs(i) = freqs->get_element(i) * 219474.63;
    }

  ExEnv::out0() << indent
               << pg_->char_table().gamma(irrep).symbol() << endl;
  int ifreqoff = 1;
  for (i=0; i<irrep; i++) ifreqoff += nfreq_[i];
  for (i=0; i<freqs.n(); i++) {
      double freq = freqs(freqs.n()-i-1);
      ExEnv::out0() << indent
                   << scprintf("%4d % 8.2f",i+ifreqoff,freq)
                   << endl;
    }
  ExEnv::out0() << endl;

  if (debug_) {
      eigvecs.print("eigenvectors");
      ncbasis.print("ncbasis");
      (ncbasis*eigvecs).print("ncbasis*eigvecs");
    }
  dynamic_cast<BlockedSCMatrix*>(
      normco_.pointer())->block(irrep).assign(ncbasis*eigvecs);
}

void
MolecularFrequencies::thermochemistry(int degeneracy, double T, double P)
{
  int i;
  double tmpvar;

  if (!nfreq_) return;

  // default values for temperature T and pressure P are
  // 298.15 K and 1 atm (=101325.0 Pa), respectively 

  // 1986 CODATA
  const double NA = 6.0221367e23;  // Avogadro's number
  const double k  = 1.380658e-23;  // Boltzmann's constant (J/K)
  const double h  = 6.6260755e-34; // Planck's constant (J*s)
  const double R  = 8.314510;      // gas constant (J/(mol*K))  (R=k*NA)
  const double pi = M_PI;
  const double hartree_to_hertz = 6.5796838e15; // (hertz/hartree)

  const double hartree_to_joule = 4.3597482e-18; // (J/hartree)
  const double hartree_to_joule_per_mol = hartree_to_joule*NA;
                                            // (J/(mol*hartree))
  const double amu_to_kg = 1.6605402e-27; // (kg/amu)
  const double angstrom_to_meter = 1.0e-10;
  const double atm_to_Pa = 101325.0; // (Pa/atm)


  ////////////////////////////////////////////////////////////////////////
  // compute the molar entropy using formulas for ideal polyatomic gasses
  // from McQuarrie, Statistical Mechanics, 1976, Ch. 8; [use (8-27) for
  // linear and (8-33) for non-linear molecules]
  // S = S_trans + S_rot + S_vib + S_el
  ////////////////////////////////////////////////////////////////////////

  // compute the mass of the molecule (in kg)
  double mass = 0.0;
  for (i=0; i<mol_->natom(); i++) {
      mass += mol_->mass(i);
      }
  mass *= amu_to_kg;

  // compute principal moments of inertia (pmi) in amu*angstrom^2
  double pmi[3];
  mol_->principal_moments_of_inertia(pmi);

  // find out if molecule is linear (if smallest pmi < 1.0e-5 amu angstrom^2)
  // (elements of pmi are sorted in order smallest to largest)
  int linear = 0;
  if (pmi[0] < 1.0e-5) linear = 1;

  // compute the symmetry number sigma;
  // for linear molecules: sigma = 2 (D_inf_h), sigma = 1 (C_inf_v)
  // for non-linear molecules: sigma = # of rot. in pt. grp, including E
  int sigma;
  CharacterTable ct = pg_->char_table();
  if (linear) {
      //if (D_inf_h) sigma = 2;
      if (ct.symbol()[0] == 'D' ||
          ct.symbol()[0] == 'd') sigma = 2;
      else if (ct.symbol()[0] == 'C' ||
               ct.symbol()[0] == 'c') sigma = 1;
      else {
          throw InputError("for linear molecules "
                           " the specified point group must be Cnv or Dnh",
                           __FILE__, __LINE__, 0, 0, class_desc());
          }
      }
  else if ((ct.symbol()[0] == 'C' ||
            ct.symbol()[0] == 'c') &&
           (ct.symbol()[1] >= '1'  &&
            ct.symbol()[1] <= '8') &&
            ct.symbol()[2] == '\0') {
      sigma = ct.order();  // group is a valid CN
      }
  else if ((ct.symbol()[0] == 'D' ||
            ct.symbol()[0] == 'd') &&
           (ct.symbol()[1] >= '2'  &&
            ct.symbol()[1] <= '6') &&
            ct.symbol()[2] == '\0') {
      sigma = ct.order();  // group is a valid DN
      }
  else if ((ct.symbol()[0] == 'T' ||
            ct.symbol()[0] == 't') &&
            ct.symbol()[1] == '\0') {
      sigma = ct.order();  // group is T
      }
  else sigma = (int)(0.5*ct.order()); // group is not pure rot. group (CN, DN, or T)

  // compute S_trans
  double S_trans;
  tmpvar = pow(2*pi*mass*k*T/(h*h),1.5);
  S_trans = R*(log(tmpvar*R*T/(P*atm_to_Pa)) + 2.5 - log(NA));

  // compute S_rot
  double S_rot;
  double theta[3]; // rotational temperatures (K)
  if (linear) {
      theta[1] = h*h/(8*pi*pi*pmi[1]*amu_to_kg*pow(angstrom_to_meter,2.0)*k);
      S_rot = log(T/(sigma*theta[1])) + 1.0;
      }
  else {
      theta[0] = h*h/(8*pi*pi*pmi[0]*amu_to_kg*pow(angstrom_to_meter,2.0)*k);
      theta[1] = h*h/(8*pi*pi*pmi[1]*amu_to_kg*pow(angstrom_to_meter,2.0)*k);
      theta[2] = h*h/(8*pi*pi*pmi[2]*amu_to_kg*pow(angstrom_to_meter,2.0)*k);
      tmpvar = theta[0]*theta[1]*theta[2];
      S_rot = log(pow(pi*T*T*T/tmpvar,0.5)/sigma) + 1.5;
      }
  S_rot *= R;

  // compute S_vib
  double S_vib = 0.0;
  for (i=0; i<nirrep_; i++) {
      for (int j=0; j<nfreq_[i]; j++) {
          if (freq_[i][j] > 0.0) {
              tmpvar = hartree_to_hertz*h*freq_[i][j]/(k*T);
              double expval = exp(-tmpvar);
              S_vib += tmpvar*expval/(1.0-expval) - log(1.0-expval);
              }
          }
      }
   S_vib *= R;

  // compute S_el
  double S_el;
  S_el = R*log(double(degeneracy));

  // compute total molar entropy S (in J/(mol*K))
  double S;
  S = S_trans + S_rot + S_vib + S_el;

  
  //////////////////////////////////////////////
  // compute the molar enthalpy (nonelectronic) 
  //////////////////////////////////////////////

  int n_zero_or_imaginary = 0;
  double E0vib = 0.0;
  for (i=0; i<nirrep_; i++) {
      for (int j=0; j<nfreq_[i]; j++) {
          if (freq_[i][j] > 0.0) E0vib += freq_[i][j] * hartree_to_joule_per_mol;
          else n_zero_or_imaginary++;
        }
    }
  E0vib *= 0.5;

  double EvibT = 0.0;
  for (i=0; i<nirrep_; i++) {
      for (int j=0; j<nfreq_[i]; j++) {
          if (freq_[i][j] > 0.0) {
              double expval = exp(-freq_[i][j]*hartree_to_joule/(k*T));
              EvibT += freq_[i][j] * hartree_to_joule_per_mol
                     * expval/(1.0-expval);
            }
        }
    }

  double EPV = NA*k*T;

  int nexternal = 6;
  if (mol_->natom() == 1) nexternal = 3;
  else if (mol_->is_linear()) nexternal = 5;

  double Erot;
  if (nexternal == 3) {
      // atom
      Erot = 0.0;
    }
  else if (nexternal == 5) {
      // linear
      Erot = EPV;
    }
  else if (nexternal == 6) {
      // nonlinear
      Erot = 1.5 * EPV;
    }
  else {
      ExEnv::errn() << "Strange number of external coordinates: " << nexternal
           << ".  Setting Erot to 0.0" << endl;
      Erot = 0.0;
    }

  double Etrans = 1.5 * EPV;

  ////////////////////////////////////////////////
  // Print out results of thermodynamic analysis
  ////////////////////////////////////////////////

  ExEnv::out0() << indent << "THERMODYNAMIC ANALYSIS:" << endl << endl
       << indent << scprintf("Contributions to the nonelectronic enthalpy at %.2lf K:\n",T)
       << indent << "                   kJ/mol       kcal/mol"<< endl
       << indent << scprintf("  E0vib        = %9.4lf    %9.4lf\n",
          E0vib/1000, E0vib/(4.184*1000))
       << indent << scprintf("  Evib(T)      = %9.4lf    %9.4lf\n",
          EvibT/1000, EvibT/(4.184*1000))
       << indent << scprintf("  Erot(T)      = %9.4lf    %9.4lf\n",
          Erot/1000, Erot/(4.184*1000))
       << indent << scprintf("  Etrans(T)    = %9.4lf    %9.4lf\n",
          Etrans/1000, Etrans/(4.184*1000))
       << indent << scprintf("  PV(T)        = %9.4lf    %9.4lf\n",
          EPV/1000, EPV/(4.184*1000))
       << indent << scprintf("  Total nonelectronic enthalpy:\n")
       << indent << scprintf("  H_nonel(T)   = %9.4lf    %9.4lf\n",
         (E0vib+EvibT+Erot+Etrans+EPV)/1000,
         (E0vib+EvibT+Erot+Etrans+EPV)/(4.184*1000))
       << endl
       << indent
       << scprintf("Contributions to the entropy at %.2lf K and %.1lf atm:\n",
                   T, P)
       << indent << "                   J/(mol*K)    cal/(mol*K)"<< endl
       << indent
       << scprintf("  S_trans(T,P) = %9.4lf    %9.4lf\n",
                   S_trans, S_trans/4.184)
       << indent
       << scprintf("  S_rot(T)     = %9.4lf    %9.4lf\n", S_rot,S_rot/4.184)
       << indent
       << scprintf("  S_vib(T)     = %9.4lf    %9.4lf\n", S_vib,S_vib/4.184)
       << indent
       << scprintf("  S_el         = %9.4lf    %9.4lf\n", S_el,S_el/4.184)
       << indent << scprintf("  Total entropy:\n")
       << indent << scprintf("  S_total(T,P) = %9.4lf    %9.4lf\n", S, S/4.184)
       << indent << endl

       << indent << "Various data used for thermodynamic analysis:" << endl
       << indent << endl;

  if (linear) ExEnv::out0() << indent << "Linear molecule" << endl;
  else ExEnv::out0() << indent << "Nonlinear molecule" << endl;

  ExEnv::out0() << indent
       << scprintf("Principal moments of inertia (amu*angstrom^2):"
          " %.5lf, %.5lf, %.5lf\n", pmi[0], pmi[1], pmi[2])
       << indent << "Point group: " << ct.symbol()
       << endl
       << indent << "Order of point group: " << ct.order() << endl
       << indent << "Rotational symmetry number: " << sigma << endl;

  if (linear) {
      ExEnv::out0() << indent
           << scprintf("Rotational temperature (K): %.4lf\n", theta[1]);
    }
  else {
      ExEnv::out0() << indent
           << scprintf("Rotational temperatures (K): %.4lf, %.4lf, %.4lf\n",
                       theta[0], theta[1], theta[2]);
    }

  ExEnv::out0() << indent << "Electronic degeneracy: " << degeneracy
       << endl << endl;
}

void
MolecularFrequencies::animate(const Ref<Render>& render,
                              const Ref<MolFreqAnimate>& anim)
{
  int i,j, symoff = 0;
  for (i=0; i<nirrep_; i++) {
      int nfreq = disym_->blocks()->size(i);
      for (j=0; j<nfreq; j++) {
          char name[128];
          sprintf(name,"%02d.%s",
                  nfreq-j+symoff, pg_->char_table().gamma(i).symbol_ns());
          anim->set_name(name);
          anim->set_mode(i,j);
          render->animate(anim.pointer());
        }
      symoff += nfreq;
    }
}

/////////////////////////////////////////////////////////////////////////////
// MolFreqAnimate

static ClassDesc MolFreqAnimate_cd(
  typeid(MolFreqAnimate),"MolFreqAnimate",1,"public AnimatedObject",
  0, create<MolFreqAnimate>, 0);

MolFreqAnimate::MolFreqAnimate(const Ref<KeyVal> &keyval):
  AnimatedObject(keyval)
{
  renmol_ << keyval->describedclassvalue("rendered");
  molfreq_ << keyval->describedclassvalue("freq");
  dependent_mole_ << keyval->describedclassvalue("dependent_mole");
  irrep_ = keyval->intvalue("irrep");
  mode_ = keyval->intvalue("mode");
  KeyValValueint default_nframe(10);
  nframe_ = keyval->intvalue("nframe",default_nframe);
  KeyValValuedouble default_disp(0.2);
  disp_ = keyval->doublevalue("displacement", default_disp);
}

MolFreqAnimate::~MolFreqAnimate()
{
}

int
MolFreqAnimate::nobject()
{
  return nframe_;
}

Ref<RenderedObject>
MolFreqAnimate::object(int iobject)
{
  BlockedSCMatrix *normco
      = dynamic_cast<BlockedSCMatrix*>(molfreq_->normal_coordinates().pointer());
  Ref<Molecule> mol = renmol_->molecule();
  Ref<Molecule> molcopy = new Molecule(*mol.pointer());

  double scale = disp_ * cos(M_PI*(iobject+0.5)/(double)nframe_);

  RefSCMatrix irrepblock = normco->block(irrep_);
  int ixyz, iatom, icoor=0;
  for (iatom=0; iatom<mol->natom(); iatom++) {
      for (ixyz=0; ixyz<3; ixyz++, icoor++) {
          mol->r(iatom,ixyz) += scale
                                   * irrepblock->get_element(icoor,mode_);
        }
    }

  if (dependent_mole_) dependent_mole_->obsolete();
  renmol_->init();

  char name[64];
  sprintf(name,"%02d",iobject);
  renmol_->set_name(name);

  // restore the original molecule
  mol->operator = (*molcopy.pointer());
  if (dependent_mole_) dependent_mole_->obsolete();

  return renmol_.pointer();
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
