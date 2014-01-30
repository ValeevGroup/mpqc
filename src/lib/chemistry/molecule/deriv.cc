//
// deriv.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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

#include <stdlib.h>
#include <fstream>

#include <util/misc/scexception.h>
#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <util/state/stateio.h>
#include <math/scmat/blocked.h>
#include <chemistry/molecule/deriv.h>
#include <chemistry/molecule/molfreq.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////
// MolecularHessian

static ClassDesc MolecularHessian_cd(
  typeid(MolecularHessian),"MolecularHessian",1,"public SavableState",
  0, 0, 0);

MolecularHessian::MolecularHessian()
{
  matrixkit_ = SCMatrixKit::default_matrixkit();
}

MolecularHessian::MolecularHessian(const Ref<KeyVal>&keyval)
{
  mol_ << keyval->describedclassvalue("molecule");
  desired_accuracy_ = keyval->doublevalue("accuracy", KeyValValuedouble(1e-4));
  matrixkit_ = SCMatrixKit::default_matrixkit();
}

MolecularHessian::MolecularHessian(StateIn&s):
  SavableState(s)
{
  mol_ << SavableState::restore_state(s);
  d3natom_ << SavableState::restore_state(s);
  matrixkit_ = SCMatrixKit::default_matrixkit();
}

MolecularHessian::~MolecularHessian()
{
}

void
MolecularHessian::save_data_state(StateOut&s)
{
  SavableState::save_state(mol_.pointer(),s);
  SavableState::save_state(d3natom_.pointer(),s);
}

RefSCDimension
MolecularHessian::d3natom()
{
  if (d3natom_.null()) d3natom_ = new SCDimension(mol_->natom()*3);
  return d3natom_;
}

RefSCMatrix
MolecularHessian::cartesian_to_symmetry(const Ref<Molecule> &mol,
                                        Ref<PointGroup> pg,
                                        Ref<SCMatrixKit> kit)
{
  int i;

  if (pg.null()) pg = mol->point_group();
  if (kit.null()) kit = SCMatrixKit::default_matrixkit();

  // create the character table for the point group
  CharacterTable ct = pg->char_table();

  int ng = ct.order();
  int nirrep = ct.nirrep();
  int natom = mol->natom();
  RefSCDimension d3natom = new SCDimension(3*natom);

  // Form the matrix of basis vectors in cartesian coordinates
  RefSCMatrix cartbasis(d3natom,d3natom,kit);
  cartbasis.assign(0.0);
  for (i=0; i<3*natom; i++) {
      cartbasis(i,i) = 1.0;
    }

  // Project out translations and rotations
  RefSCDimension dext(new SCDimension(6));
  // form a basis for the translation and rotation coordinates
  RefSCMatrix externalbasis(d3natom,dext,kit);
  externalbasis.assign(0.0);
  for (i=0; i<natom; i++) {
      SCVector3 atom(mol->r(i));
      for (int j=0; j<3; j++) {
          externalbasis(i*3 + j,j) = 1.0;
        }
      externalbasis(i*3 + 1, 3 + 0) =  atom[2];
      externalbasis(i*3 + 2, 3 + 0) = -atom[1];
      externalbasis(i*3 + 0, 3 + 1) =  atom[2];
      externalbasis(i*3 + 2, 3 + 1) = -atom[0];
      externalbasis(i*3 + 0, 3 + 2) =  atom[1];
      externalbasis(i*3 + 1, 3 + 2) = -atom[0];
    }
  // do an SVD on the external basis
  RefSCMatrix Uext(d3natom,d3natom,kit);
  RefSCMatrix Vext(dext,dext,kit);
  RefSCDimension min;
  if (d3natom.n()<dext.n()) min = d3natom;
  else min = dext;
  int nmin = min.n();
  RefDiagSCMatrix sigmaext(min,kit);
  externalbasis.svd(Uext,sigmaext,Vext);
  // find the epsilon rank
  const double epsilonext = 1.0e-4;
  int rankext = 0;
  for (i=0; i<nmin; i++) {
      if (sigmaext(i) > epsilonext) rankext++;
    }
  //ExEnv::out0() << indent << "The external rank is " << rankext << endl;
  // find the projection onto the externalbasis perp space
  if (rankext) {
      RefSCDimension drankext_tilde = new SCDimension(d3natom.n() - rankext);
      RefSCMatrix Uextr_tilde(d3natom,drankext_tilde,kit);
      Uextr_tilde.assign_subblock(Uext,
                                  0, d3natom.n()-1,
                                  0, drankext_tilde.n()-1,
                                  0, rankext);
      RefSymmSCMatrix projext_perp(d3natom, kit);
      projext_perp.assign(0.0);
      projext_perp.accumulate_symmetric_product(Uextr_tilde);
      cartbasis = projext_perp * cartbasis;
    }

  // Form the mapping of atom numbers to transformed atom number
  int **atom_map = new int*[natom];
  for (i=0; i < natom; i++) atom_map[i] = new int[ng];
  // loop over all centers
  for (i=0; i < natom; i++) {
      SCVector3 ac(mol->r(i));
      // then for each symop in the pointgroup, transform the coordinates of
      // center "i" and see which atom it maps into
      for (int g=0; g < ng; g++) {
          double np[3];
          SymmetryOperation so = ct.symm_operation(g);
          for (int ii=0; ii < 3; ii++) {
              np[ii] = 0;
              for (int jj=0; jj < 3; jj++) np[ii] += so(ii,jj) * ac[jj];
            }
          atom_map[i][g] = mol->atom_at_position(np, 0.05);
          if (atom_map[i][g] < 0) {
              throw ProgrammingError("atom mapping bad",
                                     __FILE__, __LINE__);
            }
        }
    }

  int *dims = new int[nirrep];

  RefSCMatrix *symmbasis = new RefSCMatrix[nirrep];

  // Project the cartesian basis into each irrep
  SymmetryOperation so;
  for (i=0; i<nirrep; i++) {
      IrreducibleRepresentation irrep = ct.gamma(i);
      RefSCMatrix *components = new RefSCMatrix[irrep.degeneracy()];
      // loop over the components of this irrep
      int j;
      for (j=0; j<irrep.degeneracy(); j++) {
          // form the projection matrix for this component of this irrep
          RefSCMatrix projmat(d3natom,d3natom,kit);
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
          RefSCMatrix cartbasis_ij = projmat * cartbasis;
          RefSCMatrix U(d3natom, d3natom, kit);
          RefSCMatrix V(d3natom, d3natom, kit);
          RefDiagSCMatrix sigma(d3natom, kit);
          cartbasis_ij.svd(U, sigma, V);
          // Compute the epsilon rank of cartbasis ij
          const double epsilon = 1.0e-3;
          int k, rank = 0;
          for (k=0; k<3*natom; k++) {
              if (sigma(k) > epsilon) rank++;
            }
          if (!rank) continue;
          // Find an orthogonal matrix that spans the range of cartbasis ij
          RefSCDimension drank = new SCDimension(rank);
          RefSCMatrix Ur(d3natom,drank,kit);
          Ur.assign_subblock(U,0, d3natom.n()-1, 0, drank.n()-1, 0, 0);
          // Reassign cartbasis_ij to the orthonormal basis
          cartbasis_ij = Ur;
          components[j] = cartbasis_ij;
        }
      int nbasisinirrep = 0;
      for (j=0; j<irrep.degeneracy(); j++) {
          nbasisinirrep += components[j].ncol();
        }

      dims[i] = nbasisinirrep;
      RefSCDimension dirrep = new SCDimension(nbasisinirrep);
      symmbasis[i] = kit->matrix(d3natom,dirrep);
      int offset = 0;
      for (j=0; j<irrep.degeneracy(); j++)
        if (components[j]) {
          symmbasis[i]->assign_subblock(
              components[j],
              0, d3natom.n()-1,
              offset, offset+components[j].ncol()-1,
              0, 0);
          offset += components[j].ncol();
        }
      delete[] components;
    }

  int total = 0;
  for (i=0; i<nirrep; i++) {
      total += dims[i];
    }
  Ref<SCBlockInfo> bi = new SCBlockInfo(total, nirrep, dims);
  for (i=0; i<nirrep; i++) {
      bi->set_subdim(i, symmbasis[i]->coldim());
    }
  RefSCDimension dsym = new SCDimension(bi);

  RefSCDimension bd3natom = new SCDimension(3*mol->natom());
  bd3natom->blocks()->set_subdim(0,d3natom);

  Ref<SCMatrixKit> symkit = new BlockedSCMatrixKit(kit);
  RefSCMatrix result(dsym, bd3natom, symkit);
  BlockedSCMatrix *bresult = dynamic_cast<BlockedSCMatrix*>(result.pointer());

  // put the symmetric basis in the result matrix
  for (i=0; i<nirrep; i++) {
      if (dims[i]>0) bresult->block(i).assign(symmbasis[i].t());
    }

  delete[] symmbasis;

  for (i=0; i<natom; i++) delete[] atom_map[i];
  delete[] atom_map;

  delete[] dims;
  return result;
}

void
MolecularHessian::set_energy(const Ref<MolecularEnergy> &)
{
}

MolecularEnergy*
MolecularHessian::energy() const
{
  return 0;
}

void
MolecularHessian::write_cartesian_hessian(const char *filename,
                                          const Ref<Molecule> &mol,
                                          const RefSymmSCMatrix &hess)
{
  int ntri = (3*mol->natom()*(3*mol->natom()+1))/2;
  double *hessv = new double[ntri];
  hess->convert(hessv);
  if (MessageGrp::get_default_messagegrp()->me() == 0) {
      int i,j;
      ofstream out(filename);
      // file format is version text 1
      out << "Hessian VT1" << endl;
      out << mol->natom() << " atoms" << endl;
      for (i=0; i<mol->natom(); i++) {
          out << scprintf("%2d % 15.12f % 15.12f % 15.12f",
                          mol->Z(i), mol->r(i,0), mol->r(i,1), mol->r(i,2))
              << endl;
          
        }
      const int nrow = 5;
      for (i=0; i<ntri; ) {
          for (j=0; j<nrow && i<ntri; j++,i++) {
              if (j>0) out << " ";
              out << scprintf("% 15.12f", hessv[i]);
            }
          out << endl;
        }
      out << "End Hessian" << endl;
    }
  delete[] hessv;
}

void
MolecularHessian::read_cartesian_hessian(const char *filename,
                                         const Ref<Molecule> &mol,
                                         const RefSymmSCMatrix &hess)
{
  int ntri = (3*mol->natom()*(3*mol->natom()+1))/2;
  vector<double> hessv(ntri);
  Ref<MessageGrp> grp = MessageGrp::get_default_messagegrp();
  if (grp->me() == 0) {
      int i;
      ifstream in(filename);
      const int nline = 100;
      char linebuf[nline];
      in.getline(linebuf, nline);
      if (strcmp(linebuf,"Hessian VT1")) {
          throw FileOperationFailed("not given a hessian file",
                                    __FILE__, __LINE__, filename,
                                    FileOperationFailed::Corrupt);
        }
      int natom;
      in >> natom;
      if (natom != mol->natom()) {
          throw FileOperationFailed("wrong number of atoms in hessian file",
                                    __FILE__, __LINE__, filename,
                                    FileOperationFailed::Corrupt);
        }
      in.getline(linebuf,nline);
      //ExEnv::outn() << "READ: should be atoms: " << linebuf << endl;
      for (i=0; i<mol->natom(); i++) {
          int Z;
          double x, y, z;
          in >> Z >> x >> y >> z;
          //ExEnv::outn() << "READ: " << Z << " " << x << " " << y << " " << z << endl;
        }
      for (i=0; i<ntri; i++) {
          in >> hessv[i];
          //ExEnv::outn() << "READ: hess[" << i << "] = " << hessv[i] << endl;
        }
      in.getline(linebuf, nline);
      //ExEnv::outn() << "READ: last line = " << linebuf << endl;
      if (strcmp(linebuf,"End Hessian")) {
          // try once more since there could be a left over new line
          in.getline(linebuf, nline);
          if (strcmp(linebuf,"End Hessian")) {
              //ExEnv::outn() << "READ: last line = " << linebuf << endl;
              throw FileOperationFailed("hessian file seems to be truncated",
                                        __FILE__, __LINE__, filename,
                                        FileOperationFailed::Corrupt);
            }
        }
    }
  grp->bcast(&(hessv[0]),ntri);
  hess->assign(&(hessv[0]));
}

void
MolecularHessian::set_desired_accuracy(double acc) {
  desired_accuracy_ = acc;
}

double
MolecularHessian::desired_accuracy() const {
  return desired_accuracy_;
}

/////////////////////////////////////////////////////////////////
// ReadMolecularHessian

static ClassDesc ReadMolecularHessian_cd(
  typeid(ReadMolecularHessian),"ReadMolecularHessian",1,"public MolecularHessian",
  0, create<ReadMolecularHessian>, create<ReadMolecularHessian>);

ReadMolecularHessian::ReadMolecularHessian(const Ref<KeyVal>&keyval):
  MolecularHessian(keyval)
{
  KeyValValuestring default_filename(SCFormIO::fileext_to_filename(".hess"));
  filename_ = keyval->stringvalue("filename", default_filename);
}

ReadMolecularHessian::ReadMolecularHessian(StateIn&s):
  SavableState(s),
  MolecularHessian(s)
{
  s.get(filename_);
}

ReadMolecularHessian::~ReadMolecularHessian()
{
}

void
ReadMolecularHessian::save_data_state(StateOut&s)
{
  MolecularHessian::save_data_state(s);
  s.put(filename_);
}

RefSymmSCMatrix
ReadMolecularHessian::cartesian_hessian()
{
  RefSymmSCMatrix hess = matrixkit()->symmmatrix(d3natom());
  read_cartesian_hessian(filename_.c_str(), mol_, hess);
  return hess;
}

/////////////////////////////////////////////////////////////////
// GuessMolecularHessian

static ClassDesc GuessMolecularHessian_cd(
  typeid(GuessMolecularHessian),"GuessMolecularHessian",1,"public MolecularHessian",
  0, create<GuessMolecularHessian>, create<GuessMolecularHessian>);

GuessMolecularHessian::GuessMolecularHessian(const Ref<KeyVal>&keyval):
  MolecularHessian(keyval)
{
  coor_ << keyval->describedclassvalue("coor");
  if (mol_.null()) mol_ = coor_->molecule();
}

GuessMolecularHessian::GuessMolecularHessian(StateIn&s):
  SavableState(s),
  MolecularHessian(s)
{
  coor_ << SavableState::restore_state(s);
}

GuessMolecularHessian::~GuessMolecularHessian()
{
}

void
GuessMolecularHessian::save_data_state(StateOut&s)
{
  MolecularHessian::save_data_state(s);
  SavableState::save_state(coor_.pointer(),s);
}

RefSymmSCMatrix
GuessMolecularHessian::cartesian_hessian()
{
  RefSymmSCMatrix hessian(coor_->dim(), coor_->matrixkit());
  coor_->guess_hessian(hessian);
  RefSymmSCMatrix xhessian(coor_->dim_natom3(), coor_->matrixkit());
  coor_->to_cartesian(xhessian,hessian);
  return xhessian;
}

/////////////////////////////////////////////////////////////////
// DiagMolecularHessian

static ClassDesc DiagMolecularHessian_cd(
  typeid(DiagMolecularHessian),"DiagMolecularHessian",1,"public MolecularHessian",
  0, create<DiagMolecularHessian>, create<DiagMolecularHessian>);

DiagMolecularHessian::DiagMolecularHessian(const Ref<KeyVal>&keyval):
  MolecularHessian(keyval)
{
  diag_ = keyval->doublevalue("diag",KeyValValuedouble(1.0));
}

DiagMolecularHessian::DiagMolecularHessian(StateIn&s):
  SavableState(s),
  MolecularHessian(s)
{
  s.get(diag_);
}

DiagMolecularHessian::~DiagMolecularHessian()
{
}

void
DiagMolecularHessian::save_data_state(StateOut&s)
{
  MolecularHessian::save_data_state(s);
  s.put(diag_);
}

RefSymmSCMatrix
DiagMolecularHessian::cartesian_hessian()
{
  RefSymmSCMatrix xhessian(d3natom(), matrixkit());
  xhessian->assign(0.0);
  xhessian->shift_diagonal(diag_);
  return xhessian;
}

/////////////////////////////////////////////////////////////////////////////
// MolecularGradient

static ClassDesc MolecularGradient_cd(
  typeid(MolecularGradient),"MolecularGradient",1,"public SavableState",
  0, 0, 0);

MolecularGradient::MolecularGradient()
{
  matrixkit_ = SCMatrixKit::default_matrixkit();
}

MolecularGradient::MolecularGradient(const Ref<KeyVal>&keyval)
{
  mol_ << keyval->describedclassvalue("molecule");
  desired_accuracy_ = keyval->doublevalue("accuracy", KeyValValuedouble(1e-4));
  matrixkit_ = SCMatrixKit::default_matrixkit();
}

MolecularGradient::MolecularGradient(StateIn&s):
  SavableState(s)
{
  mol_ << SavableState::restore_state(s);
  d3natom_ << SavableState::restore_state(s);
  matrixkit_ = SCMatrixKit::default_matrixkit();
  s.get(desired_accuracy_);
}

MolecularGradient::~MolecularGradient()
{
}

void
MolecularGradient::save_data_state(StateOut&s)
{
  SavableState::save_state(mol_.pointer(),s);
  SavableState::save_state(d3natom_.pointer(),s);
  s.put(desired_accuracy_);
}

RefSCDimension
MolecularGradient::d3natom()
{
  if (d3natom_.null()) d3natom_ = new SCDimension(mol_->natom()*3);
  return d3natom_;
}

void
MolecularGradient::set_energy(const Ref<MolecularEnergy> &)
{
}

MolecularEnergy*
MolecularGradient::energy() const
{
  return 0;
}

void
MolecularGradient::write_cartesian_gradient(const char *filename,
                                            const Ref<Molecule> &mol,
                                            const RefSCVector &grad)
{
  const int ncoord = 3 * mol->natom();
  double *gradv = new double[ncoord];
  grad->convert(gradv);
  if (MessageGrp::get_default_messagegrp()->me() == 0) {
      int i,j;
      ofstream out(filename);
      // file format is version text 1
      out << "Gradient VT1" << endl;
      out << mol->natom() << " atoms" << endl;
      for (i=0; i<mol->natom(); i++) {
          out << scprintf("%2d % 15.12f % 15.12f % 15.12f",
                          mol->Z(i), mol->r(i,0), mol->r(i,1), mol->r(i,2))
              << endl;

        }
      const int nrow = 5;
      for (i=0; i<ncoord; ) {
          for (j=0; j<nrow && i<ncoord; j++,i++) {
              if (j>0) out << " ";
              out << scprintf("% 20.15f", gradv[i]);
            }
          out << endl;
        }
      out << "End Gradient" << endl;
    }
  delete[] gradv;
}

void
MolecularGradient::read_cartesian_gradient(const char *filename,
                                           const Ref<Molecule> &mol,
                                           const RefSCVector &grad)
{
  const int ncoord = 3 * mol->natom();
  vector<double> gradv(ncoord);
  Ref<MessageGrp> grp = MessageGrp::get_default_messagegrp();
  if (grp->me() == 0) {
      int i;
      ifstream in(filename);
      const int nline = 100;
      char linebuf[nline];
      in.getline(linebuf, nline);
      if (strcmp(linebuf,"Gradient VT1")) {
          throw FileOperationFailed("not given a gradient file",
                                    __FILE__, __LINE__, filename,
                                    FileOperationFailed::Corrupt);
        }
      int natom;
      in >> natom;
      if (natom != mol->natom()) {
          throw FileOperationFailed("wrong number of atoms in gradient file",
                                    __FILE__, __LINE__, filename,
                                    FileOperationFailed::Corrupt);
        }
      in.getline(linebuf,nline);
      //ExEnv::outn() << "READ: should be atoms: " << linebuf << endl;
      for (i=0; i<mol->natom(); i++) {
          int Z;
          double x, y, z;
          in >> Z >> x >> y >> z;
          //ExEnv::outn() << "READ: " << Z << " " << x << " " << y << " " << z << endl;
        }
      for (i=0; i<ncoord; i++) {
          in >> gradv[i];
          //ExEnv::outn() << "READ: grad[" << i << "] = " << gradv[i] << endl;
        }
      in.getline(linebuf, nline);
      //ExEnv::outn() << "READ: last line = " << linebuf << endl;
      if (strcmp(linebuf,"End Gradient")) {
          // try once more since there could be a left over new line
          in.getline(linebuf, nline);
          if (strcmp(linebuf,"End Gradient")) {
              //ExEnv::outn() << "READ: last line = " << linebuf << endl;
              throw FileOperationFailed("gradient file seems to be truncated",
                                        __FILE__, __LINE__, filename,
                                        FileOperationFailed::Corrupt);
            }
        }
    }
  grp->bcast(&(gradv[0]),ncoord);
  grad->assign(&(gradv[0]));
}

void
MolecularGradient::set_desired_accuracy(double acc) {
  desired_accuracy_ = acc;
}

double
MolecularGradient::desired_accuracy() const {
  return desired_accuracy_;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
