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

#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>
#include <util/misc/formio.h>
#include <math/scmat/local.h>
#include <math/scmat/blocked.h>
#include <chemistry/molecule/molfreq.h>
#include <chemistry/molecule/molrender.h>

#define CLASSNAME MolecularFrequencies
#define VERSION 2
#define PARENTS public SavableState
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
MolecularFrequencies::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

MolecularFrequencies::MolecularFrequencies(const RefKeyVal& keyval)
{
  mole_ = keyval->describedclassvalue("energy");

  debug_ = keyval->booleanvalue("debug");

  if (mole_.null()) {
      mol_ = keyval->describedclassvalue("molecule");
      kit_ = SCMatrixKit::default_matrixkit();
    }
  else {
      mol_ = mole_->molecule();
      kit_ = mole_->matrixkit();
      d3natom_ = mole_->moldim();
    }
  if (d3natom_.null()) d3natom_ = new SCDimension(3*mol_->natom());
  bd3natom_ = new SCDimension(3*mol_->natom());
  bd3natom_->blocks()->set_subdim(0,d3natom_);
  symkit_ = new BlockedSCMatrixKit(kit_);

  displacement_point_group_ = keyval->describedclassvalue("point_group");
  if (displacement_point_group_.null()) {
      displacement_point_group_
          = new PointGroup(*mol_->point_group().pointer());
    }

  nirrep_ = displacement_point_group_->char_table().nirrep();
  displacements_ = new RefSCMatrix[nirrep_];

  disp_ = keyval->doublevalue("displacement");
  if (keyval->error() != KeyVal::OK) disp_ = 0.001;

  gradients_ = 0;

  nfreq_ = 0;
  freq_ = 0;
  ndisp_ = 0;
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
  delete[] displacements_;
  delete[] gradients_;
}

MolecularFrequencies::MolecularFrequencies(StateIn& si):
  SavableState(si)
{
  int i;

  if (si.version(static_class_desc()) < 2) {
      cerr << "MolecularFrequencies: cannot restore from old version" << endl;
      abort();
    }

  original_point_group_.restore_state(si);

  if (si.version(static_class_desc()) >= 2) {
      displacement_point_group_.restore_state(si);
    }

  mol_.restore_state(si);
  mole_.restore_state(si);

  if (mole_.null()) {
      kit_ = SCMatrixKit::default_matrixkit();
    }
  else {
      kit_ = mole_->matrixkit();
      d3natom_ = mole_->moldim();
    }
  if (d3natom_.null()) d3natom_ = new SCDimension(3*mol_->natom());
  bd3natom_ = new SCDimension(3*mol_->natom());
  bd3natom_->blocks()->set_subdim(0,d3natom_);
  symkit_ = new BlockedSCMatrixKit(kit_);

  si.get(disp_);
  si.get(ndisp_);
  si.get(nirrep_);
  si.get(nexternal_);
  displacements_ = new RefSCMatrix[nirrep_];

  for (i=0; i < nirrep_; i++) {
      int ndisp;
      si.get(ndisp);
      RefSCDimension ddisp = new SCDimension(ndisp);
      displacements_[i] = matrixkit()->matrix(d3natom_,ddisp);
      displacements_[i].restore(si);
    }

  freq_ = 0;
  nfreq_ = 0;
  debug_ = 0;

  gradients_ = new RefSCVector[ndisplace()];
  for (i=0; i < ndisp_; i++) {
      int ndisp;
      si.get(ndisp);
      RefSCDimension ddisp = new SCDimension(ndisp);
      gradients_[i] = matrixkit()->vector(ddisp);
      gradients_[i].restore(si);
    }

  original_geometry_ = matrixkit()->vector(d3natom_);
  original_geometry_.restore(si);
  disym_.restore_state(si);
}

void
MolecularFrequencies::save_data_state(StateOut& so)
{
  int i;
  original_point_group_.save_state(so);

  displacement_point_group_.save_state(so);
  mol_.save_state(so);
  mole_.save_state(so);
  so.put(disp_);
  so.put(ndisp_);
  so.put(nirrep_);
  so.put(nexternal_);
  for (i=0; i < nirrep_; i++) {
      so.put(displacements_[i].ncol());
      displacements_[i].save(so);
    }

  for (i=0; i < ndisp_; i++) {
      so.put(gradients_[i].n());
      gradients_[i].save(so);
    }

  original_geometry_.save(so);
  disym_.save_state(so);
}

void
MolecularFrequencies::compute_displacements()
{
  // create the character table for the point group
  CharacterTable ct = displacement_point_group_->char_table();

  int ng = ct.order();
  int natom = mol_->natom();

  original_point_group_ = mol_->point_group();
  original_geometry_ = matrixkit()->vector(d3natom_);

  int i, coor;
  for (i=0, coor=0; i<mol_->natom(); i++) {
      for (int j=0; j<3; j++, coor++) {
          original_geometry_(coor) = mol_->r(i,j);
        }
    }

  // Form the matrix of displacements in cartesian coordinates
  RefSCMatrix cartdisp(d3natom_,d3natom_,matrixkit());
  cartdisp.assign(0.0);
  for (i=0; i<3*natom; i++) {
      cartdisp(i,i) = 1.0;
    }

  // Project out translations and rotations
  RefSCDimension dext(new SCDimension(6));
  // form a basis for the translation and rotation coordinates
  RefSCMatrix externaldisp(d3natom_,dext,matrixkit());
  externaldisp.assign(0.0);
  for (i=0; i<natom; i++) {
      SCVector3 atom(mol_->r(i));
      for (int j=0; j<3; j++) {
          externaldisp(i*3 + j,j) = 1.0;
        }
      externaldisp(i*3 + 1, 3 + 0) =  atom[2];
      externaldisp(i*3 + 2, 3 + 0) = -atom[1];
      externaldisp(i*3 + 0, 3 + 1) =  atom[2];
      externaldisp(i*3 + 2, 3 + 1) = -atom[0];
      externaldisp(i*3 + 0, 3 + 2) =  atom[1];
      externaldisp(i*3 + 1, 3 + 2) = -atom[0];
    }
  // do an SVD on the external displacements
  RefSCMatrix Uext(d3natom_,d3natom_,matrixkit());
  RefSCMatrix Vext(dext,dext,matrixkit());
  RefSCDimension min;
  nexternal_ = dext.n();
  if (d3natom_.n()<dext.n()) min = d3natom_;
  else min = dext;
  int nmin = min.n();
  RefDiagSCMatrix sigmaext(min,matrixkit());
  externaldisp.svd(Uext,sigmaext,Vext);
  // find the epsilon rank
  const double epsilonext = 1.0e-4;
  int rankext = 0;
  for (i=0; i<nmin; i++) {
      if (sigmaext(i) > epsilonext) rankext++;
    }
  cout << node0 << indent << "The external rank is " << rankext << endl;
  // find the projection onto the externaldisp perp space
  if (rankext) {
      RefSCDimension drankext_tilde = new SCDimension(d3natom_.n() - rankext);
      RefSCMatrix Uextr_tilde(d3natom_,drankext_tilde,matrixkit());
      Uextr_tilde.assign_subblock(Uext,
                                  0, d3natom_.n()-1,
                                  0, drankext_tilde.n()-1,
                                  0, rankext);
      RefSymmSCMatrix projext_perp(d3natom_, matrixkit());
      projext_perp.assign(0.0);
      projext_perp.accumulate_symmetric_product(Uextr_tilde);
      cartdisp = projext_perp * cartdisp;
    }

  // Form the mapping of atom numbers to transformed atom number
  int **atom_map = new int*[natom];
  for (i=0; i < natom; i++) atom_map[i] = new int[ng];
  // loop over all centers
  for (i=0; i < natom; i++) {
      SCVector3 ac(mol_->r(i));
      // then for each symop in the pointgroup, transform the coordinates of
      // center "i" and see which atom it maps into
      for (int g=0; g < ng; g++) {
          double np[3];
          SymmetryOperation so = ct.symm_operation(g);
          for (int ii=0; ii < 3; ii++) {
              np[ii] = 0;
              for (int jj=0; jj < 3; jj++) np[ii] += so(ii,jj) * ac[jj];
            }
          atom_map[i][g] = mol_->atom_at_position(np, 0.05);
          if (atom_map[i][g] < 0) {
              cerr << node0 << indent
                   << "MolecularFrequencies: atom mapping bad" << endl;
              abort();
            }
        }
    }

  int *dims = new int[nirrep_];

  // Project the cartesian displacements into each irrep
  SymmetryOperation so;
  for (i=0; i<nirrep_; i++) {
      IrreducibleRepresentation irrep = ct.gamma(i);
      RefSCMatrix *components = new RefSCMatrix[irrep.degeneracy()];
      // loop over the components of this irrep
      int j;
      for (j=0; j<irrep.degeneracy(); j++) {
          // form the projection matrix for this component of this irrep
          RefSCMatrix projmat(d3natom_,d3natom_,matrixkit());
          projmat.assign(0.0);
          // form the projection matrix for irrep i component j
          // loop over the symmetry operators
          for (int g=0; g < ng; g++) {
              double coef = ((double)irrep.character(g)*irrep.degeneracy())/ng;
              so = ct.symm_operation(g);
              for (int atom=0; atom<natom; atom++) {
                  for (int ii=0; ii < 3; ii++) {
                      for (int jj=0; jj < 3; jj++) {
                          projmat.accumulate_element(atom_map[atom][g]*3+ii,
                                                     atom*3 + jj,
                                                     coef * so(ii,jj));
                        }
                    }
                }
            }
          // projection matrix for irrep i, component j is formed
          RefSCMatrix cartdisp_ij = projmat * cartdisp;
          RefSCMatrix U(d3natom_, d3natom_, matrixkit());
          RefSCMatrix V(d3natom_, d3natom_, matrixkit());
          RefDiagSCMatrix sigma(d3natom_, matrixkit());
          cartdisp_ij.svd(U, sigma, V);
          // Compute the epsilon rank of cartdisp ij
          const double epsilon = 1.0e-3;
          int k, rank = 0;
          for (k=0; k<3*natom; k++) {
              if (sigma(k) > epsilon) rank++;
            }
          if (!rank) continue;
          // Find an orthogonal matrix that spans the range of cartdisp ij
          RefSCDimension drank = new SCDimension(rank);
          RefSCMatrix Ur(d3natom_,drank,matrixkit());
          Ur.assign_subblock(U,0, d3natom_.n()-1, 0, drank.n()-1, 0, 0);
          // Reassign cartdisp_ij to the orthonormal displacement
          cartdisp_ij = Ur;
          if (debug_) {
              cout << node0 << indent
                   << "Irrep " << irrep.symbol() << " component " << j << endl;
              cartdisp_ij.print("cartdisp:",cout);
            }
          components[j] = cartdisp_ij;
        }
      int ndisp = 0;
      for (j=0; j<irrep.degeneracy(); j++) ndisp += components[j].ncol();

      dims[i] = ndisp;
      RefSCDimension ddisp = new SCDimension(ndisp);
      displacements_[i] = matrixkit()->matrix(d3natom_,ddisp);
      int offset = 0;
      for (j=0; j<irrep.degeneracy(); j++) {
          displacements_[i]->assign_subblock(
              components[j],
              0, d3natom_.n()-1,
              offset, offset+components[j].ncol()-1,
              0, 0);
          offset += components[j].ncol();
        }
      delete[] components;
    }

  int total = 0;
  for (i=0; i<nirrep_; i++) {
      total += dims[i];
    }
  RefSCBlockInfo bi = new SCBlockInfo(total, nirrep_, dims);
  delete[] dims;
  for (i=0; i<nirrep_; i++) {
      bi->set_subdim(i, displacements_[i]->coldim());
    }
  disym_ = new SCDimension(bi);

  for (i=0; i<natom; i++) delete[] atom_map[i];
  delete[] atom_map;

  gradients_ = new RefSCVector[ndisplace()];
}

void
MolecularFrequencies::get_disp(int disp, int &irrep, int &index, double &coef)
{
  int disp_offset = 0;

  // check for +ve totally symmetric displacements
  if (disp < disp_offset + displacements_[0].ncol()) {
      irrep = 0;
      coef = 1.0;
      index = disp - disp_offset;
      return;
    }
  disp_offset += displacements_[0].ncol();
  // check for -ve totally symmetric displacements
  if (disp < disp_offset + displacements_[0].ncol()) {
      irrep = 0;
      coef = -1.0;
      index = disp - disp_offset;
      return;
    }
  disp_offset += displacements_[0].ncol();
  for (int i=1; i<nirrep_; i++) {
      if (disp < disp_offset + displacements_[i].ncol()) {
          irrep = i;
          coef = 1.0;
          index = disp - disp_offset;
          return;
        }
      disp_offset += displacements_[i].ncol();
    }
  cerr << node0 << indent
       << "MolecularFrequencies::get_disp: bad disp number" << endl;
  abort();
}

int
MolecularFrequencies::ndisplace() const
{
  int ndisp = 2 * displacements_[0].ncol();
  for (int i=1; i<nirrep_; i++) {
      ndisp += displacements_[i].ncol();
    }
  return ndisp;
}

void
MolecularFrequencies::displace(int disp)
{
  int irrep, index;
  double coef;
  get_disp(disp, irrep, index, coef);

  if (mole_.nonnull()) mole_->obsolete();

  for (int i=0, coor=0; i<mol_->natom(); i++) {
      for (int j=0; j<3; j++, coor++) {
          mol_->r(i,j) = original_geometry_(coor)
                           + coef * disp_
                            * displacements_[irrep]->get_element(coor,index);
        }
    }

  if (irrep == 0) {
      mol_->set_point_group(original_point_group_);
    }
  else {
      // Future work: doesn't need to be reduced to c1 symmetry here
      mol_->set_point_group(new PointGroup("c1"));
    }
}

void
MolecularFrequencies::original_geometry()
{
  if (mole_.nonnull()) mole_->obsolete();

  for (int i=0, coor=0; i<mol_->natom(); i++) {
      for (int j=0; j<3; j++, coor++) {
          mol_->r(i,j) = original_geometry_(coor);
        }
    }

  mol_->set_point_group(original_point_group_);
}

void
MolecularFrequencies::set_energy(const RefMolecularEnergy &mole)
{
  mole_ = mole;
  mol_ = mole_->molecule();
}

void
MolecularFrequencies::set_gradient(int disp, const RefSCVector &grad)
{
  int irrep, index;
  double coef;
  get_disp(disp, irrep, index, coef);

  // transform the gradient into symmetrized coordinates
  gradients_[disp] = displacements_[irrep].t() * grad;
  if (debug_) {
      grad.print("cartesian gradient");
      gradients_[disp].print("internal gradient");
    }

  ndisp_++;
}

void
MolecularFrequencies::compute_frequencies_from_gradients()
{
  int i, coor;

  cout << node0 << endl
       << indent << "Frequencies (cm-1; negative is imaginary):";

  // initialize the frequency tables
  nfreq_ = new int[nirrep_];
  freq_ = new double*[nirrep_];

  // initialize normal cooridinate matrix
  normco_ = symmatrixkit()->matrix(bd3natom_, disym_);

  // find the inverse sqrt mass matrix
  RefDiagSCMatrix m(d3natom_, matrixkit());
  for (i=0,coor=0; i<mol_->natom(); i++) {
      for (int j=0; j<3; j++, coor++) {
          m(coor) = 1.0/sqrt(mol_->mass(i)*(1.0/5.48579903e-4));
        }
    }

  RefSymmSCMatrix dhessian;

  RefSymmSCMatrix xhessian;
  if (debug_) {
      xhessian = matrixkit()->symmmatrix(d3natom_);
      xhessian.assign(0.0);
    }

  // start with the totally symmetry frequencies
  RefSCMatrix dtrans = displacements_[0];
  RefSCDimension ddim = dtrans.coldim();
  nfreq_[0] = ddim.n();
  freq_[0] = new double[nfreq_[0]];
  dhessian = matrixkit()->symmmatrix(ddim);
  for (i=0; i<ddim.n(); i++) {
      for (int j=0; j<=i; j++) {
          dhessian(i,j) = (gradients_[i](j) - gradients_[i+ddim.n()](j)
                         + gradients_[j](i) - gradients_[j+ddim.n()](i))
                         /(4.0*disp_);
        }
    }
  do_freq_for_irrep(0, m, dhessian, xhessian);

  int offset = 2*ddim.n();
  for (int irrep=1; irrep<nirrep_; irrep++) {
      dtrans = displacements_[irrep];
      ddim = dtrans.coldim();
      nfreq_[irrep] = ddim.n();
      freq_[irrep] = new double[nfreq_[irrep]];
      if (ddim.n() == 0) continue;
      dhessian = matrixkit()->symmmatrix(ddim);
      for (i=0; i<ddim.n(); i++) {
          for (int j=0; j<=i; j++) {
              dhessian(i,j) = (gradients_[i+offset](j)
                             + gradients_[j+offset](i))
                             /(2.0*disp_);
            }
        }
      do_freq_for_irrep(irrep, m, dhessian, xhessian);
      offset += ddim.n();
    }

  if (debug_) {
      xhessian.print("xhessian");

      RefSCMatrix mrect(d3natom_,d3natom_,matrixkit());
      mrect.assign(0.0);
      mrect->accumulate(m.pointer());

      RefSymmSCMatrix mxhessian(d3natom_,matrixkit());
      mxhessian.assign(0.0);
      mxhessian.accumulate_transform(mrect,xhessian);
      mxhessian.print("mass weighted cartesian hessian");
      RefDiagSCMatrix freqs(d3natom_, matrixkit());
      RefSCMatrix eigvec(d3natom_, d3natom_, matrixkit());
      mxhessian.diagonalize(freqs,eigvec);
      // convert the eigvals to frequencies in wavenumbers
      for (i=0; i<freqs.n(); i++) {
          if (freqs(i) >=0.0) freqs(i) = sqrt(freqs(i));
          else freqs(i) = -sqrt(-freqs(i));
          freqs(i) = freqs->get_element(i) * 219474.63;
        }
      freqs.print("Frequencies from cartesian hessian");

      eigvec.print("Mass weighted cartesian hessian eigenvectors.");

      normco_.print("Normal Coordinates");
    }

}

void
MolecularFrequencies::do_freq_for_irrep(int irrep,
                                        const RefDiagSCMatrix &m,
                                        const RefSymmSCMatrix &dhessian,
                                        const RefSymmSCMatrix &xhessian)
{
  int i;
  RefSCMatrix dtrans = displacements_[irrep];
  RefSCDimension ddim = dtrans.coldim();
  if (ddim.n() == 0) return;
  if (debug_) {
      dhessian.print("dhessian");
      dtrans.print("dtrans");
      xhessian.accumulate_transform(dtrans, dhessian);
    }
  // find the basis for the normal coordinates
  RefSCMatrix ncbasis = m * dtrans;
  // use the SVD to orthogonalize and check this basis
  RefSCMatrix basU(d3natom_, d3natom_, matrixkit());
  RefSCMatrix basV(ddim, ddim, matrixkit());
  RefDiagSCMatrix bassigma(ddim, matrixkit());
  ncbasis.svd(basU, bassigma, basV);
  for (i=0; i<ddim.n(); i++) {
      if (bassigma(i) < 1.e0-3) {
          cerr << node0 << indent
               << "MolecularFrequencies: displacements don't span"
               << " normal coordinates"
               << endl;
          abort();
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
  freqs.print(displacement_point_group_->char_table().gamma(irrep).symbol());
  if (debug_) {
      eigvecs.print("eigenvectors");
      ncbasis.print("ncbasis");
      (ncbasis*eigvecs).print("ncbasis*eigvecs");
    }
  BlockedSCMatrix::castdown(
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
  const double pi = 3.14159265358979323846;
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
  CharacterTable ct = displacement_point_group_->char_table();
  if (linear) {
      //if (D_inf_h) sigma = 2;
      if (displacement_point_group_->symbol()[0] == 'D' ||
          displacement_point_group_->symbol()[0] == 'd') sigma = 2;
      else if (displacement_point_group_->symbol()[0] == 'C' ||
               displacement_point_group_->symbol()[0] == 'c') sigma = 1;
      else {
          cerr << "MolecularFrequencies: For linear molecules"
               << " the specified point group must be Cnv or Dnh"
               << endl;
          abort();
          }
      }
  else if ((displacement_point_group_->symbol()[0] == 'C' ||
            displacement_point_group_->symbol()[0] == 'c') &&
           (displacement_point_group_->symbol()[1] >= '1'  &&
            displacement_point_group_->symbol()[1] <= '8') &&
            displacement_point_group_->symbol()[2] == '\0') {
      sigma = ct.order();  // group is a valid CN
      }
  else if ((displacement_point_group_->symbol()[0] == 'D' ||
            displacement_point_group_->symbol()[0] == 'd') &&
           (displacement_point_group_->symbol()[1] >= '2'  &&
            displacement_point_group_->symbol()[1] <= '6') &&
            displacement_point_group_->symbol()[2] == '\0') {
      sigma = ct.order();  // group is a valid DN
      }
  else if ((displacement_point_group_->symbol()[0] == 'T' ||
            displacement_point_group_->symbol()[0] == 't') &&
            displacement_point_group_->symbol()[1] == '\0') {
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
              S_vib += tmpvar/(exp(tmpvar)-1) - log(1-exp(-tmpvar));
              }
          }
      }
   S_vib *= R;

  // compute S_el
  double S_el;
  S_el = R*log(degeneracy);

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
          if (freq_[i][j] > 0.0)
              EvibT += freq_[i][j] * hartree_to_joule_per_mol
                      /(exp(freq_[i][j]*hartree_to_joule/k*T-1));
        }
    }

  double EPV = NA*k*T;

  double Erot;
  if (nexternal_ == 3) {
      // atom
      Erot = 0.0;
    }
  else if (nexternal_ == 5) {
      // linear
      Erot = EPV;
    }
  else if (nexternal_ == 6) {
      // nonlinear
      Erot = 1.5 * EPV;
    }
  else {
      cerr << "Strange number of external coordinates: " << nexternal_
           << ".  Setting Erot to 0.0" << endl;
      Erot = 0.0;
    }

  double Etrans = 1.5 * EPV;

  ////////////////////////////////////////////////
  // Print out results of thermodynamic analysis
  ////////////////////////////////////////////////

  cout << node0 << "THERMODYNAMIC ANALYSIS:" << endl;
  cout << node0 << endl;
  cout << node0 << scprintf("Contributions to the nonelectronic enthalpy at %.2lf K:\n",T);
  cout << node0 << "                   kJ/mol       kcal/mol"<< endl;
  cout << node0 << scprintf("  E0vib        = %9.4lf    %9.4lf\n",
          E0vib/1000, E0vib/(4.184*1000));
  cout << node0 << scprintf("  Evib(T)      = %9.4lf    %9.4lf\n",
          EvibT/1000, EvibT/(4.184*1000));
  cout << node0 << scprintf("  Erot(T)      = %9.4lf    %9.4lf\n",
          Erot/1000, Erot/(4.184*1000));
  cout << node0 << scprintf("  Etrans(T)    = %9.4lf    %9.4lf\n",
          Etrans/1000, Etrans/(4.184*1000));
  cout << node0 << scprintf("  PV(T)        = %9.4lf    %9.4lf\n",
          EPV/1000, EPV/(4.184*1000));
  cout << node0 << scprintf("  Total nonelectronic enthalpy:\n");
  cout << node0 << scprintf("  H_nonel(T)   = %9.4lf    %9.4lf\n",
         (E0vib+EvibT+Erot+Etrans+EPV)/1000,
         (E0vib+EvibT+Erot+Etrans+EPV)/(4.184*1000));
  cout << node0 << endl;

  cout << node0
       << scprintf("Contributions to the entropy at %.2lf K and %.1lf atm:\n",
                   T, P);
  cout << node0 << "                   J/(mol*K)    cal/(mol*K)"<< endl;
  cout << node0
       << scprintf("  S_trans(T,P) = %9.4lf    %9.4lf\n",
                   S_trans, S_trans/4.184);
  cout << node0
       << scprintf("  S_rot(T)     = %9.4lf    %9.4lf\n", S_rot,S_rot/4.184);
  cout << node0
       << scprintf("  S_vib(T)     = %9.4lf    %9.4lf\n", S_vib,S_vib/4.184);
  cout << node0
       << scprintf("  S_el         = %9.4lf    %9.4lf\n", S_el,S_el/4.184);
  cout << node0 << scprintf("  Total entropy:\n");
  cout << node0 << scprintf("  S_total(T,P) = %9.4lf    %9.4lf\n", S, S/4.184);
  cout << node0 << endl;

  cout << node0 << "Various data used for thermodynamic analysis:" << endl;
  cout << node0 << endl;
  if (linear) cout << node0 << "Linear molecule" << endl;
  else cout << node0 << "Nonlinear molecule" << endl;
  cout << node0 << scprintf("Principal moments of inertia (amu*angstrom^2):"
          " %.5lf, %.5lf, %.5lf\n", pmi[0], pmi[1], pmi[2]);
  cout << node0 << "Point group: " << displacement_point_group_->symbol()
       << endl;
  cout << node0 << "Order of point group: " << ct.order() << endl;
  cout << node0 << "Rotational symmetry number: " << sigma << endl;
  if (linear) {
      cout << node0
           << scprintf("Rotational temperature (K): %.4lf\n", theta[1]);
    }
  else {
      cout << node0
           << scprintf("Rotational temperatures (K): %.4lf, %.4lf, %.4lf\n",
                       theta[0], theta[1], theta[2]);
    }
  cout << node0 << "Electronic degeneracy: " << degeneracy << endl;
  cout << node0 << endl;
}

void
MolecularFrequencies::animate(const RefRender& render,
                              const RefMolFreqAnimate& anim)
{
  int i,j;
  for (i=0; i<nirrep_; i++) {
      for (j=0; j<disym_->blocks()->size(i); j++) {
          char name[128];
          sprintf(name,"%s.%02d",
                displacement_point_group_->char_table().gamma(i).symbol(), j);
          anim->set_name(name);
          anim->set_mode(i,j);
          render->animate(anim);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// MolFreqAnimate

#define CLASSNAME MolFreqAnimate
#define PARENTS public AnimatedObject
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
MolFreqAnimate::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = AnimatedObject::_castdown(cd);
  return do_castdowns(casts,cd);
}

MolFreqAnimate::MolFreqAnimate(const RefKeyVal &keyval):
  AnimatedObject(keyval)
{
  renmol_ = keyval->describedclassvalue("rendered");
  molfreq_ = keyval->describedclassvalue("freq");
  dependent_mole_ = keyval->describedclassvalue("dependent_mole");
  irrep_ = keyval->intvalue("irrep");
  mode_ = keyval->intvalue("mode");
  nframe_ = keyval->intvalue("nframe");
  if (keyval->error() != KeyVal::OK) {
      nframe_ = 10;
    }
}

MolFreqAnimate::~MolFreqAnimate()
{
}

int
MolFreqAnimate::nobject()
{
  return nframe_;
}

RefRenderedObject
MolFreqAnimate::object(int iobject)
{
  BlockedSCMatrix *normco
      = BlockedSCMatrix::castdown(molfreq_->normal_coordinates());
  RefMolecule mol = renmol_->molecule();

  double scale = 0.2 * cos(M_PI*(iobject+0.5)/(double)nframe_);

  molfreq_->original_geometry();

  RefSCMatrix irrepblock = normco->block(irrep_);
  int ixyz, iatom, icoor=0;
  for (iatom=0; iatom<mol->natom(); iatom++) {
      for (ixyz=0; ixyz<3; ixyz++, icoor++) {
          mol->r(iatom,ixyz) += scale
                                   * irrepblock->get_element(icoor,mode_);
        }
    }

  if (dependent_mole_.nonnull()) dependent_mole_->obsolete();
  renmol_->init();

  char name[64];
  sprintf(name,"%02d",iobject);
  renmol_->set_name(name);

  molfreq_->original_geometry();

  return renmol_;
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
