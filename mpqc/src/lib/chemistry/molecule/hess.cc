//
// hess.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdlib.h>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <math/scmat/blocked.h>
#include <chemistry/molecule/hess.h>
#include <chemistry/molecule/molfreq.h>

/////////////////////////////////////////////////////////////////
// MolecularHessian

SavableState_REF_def(MolecularHessian);

#define CLASSNAME MolecularHessian
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
MolecularHessian::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

MolecularHessian::MolecularHessian()
{
  matrixkit_ = SCMatrixKit::default_matrixkit();
}

MolecularHessian::MolecularHessian(const RefKeyVal&keyval)
{
  mol_ = keyval->describedclassvalue("molecule");
  matrixkit_ = SCMatrixKit::default_matrixkit();
}

MolecularHessian::MolecularHessian(StateIn&s):
  SavableState(s)
{
  mol_.restore_state(s);
  d3natom_.restore_state(s);
  matrixkit_ = SCMatrixKit::default_matrixkit();
}

MolecularHessian::~MolecularHessian()
{
}

void
MolecularHessian::save_data_state(StateOut&s)
{
  mol_.save_state(s);
  d3natom_.save_state(s);
}

RefSCDimension
MolecularHessian::d3natom()
{
  if (d3natom_.null()) d3natom_ = new SCDimension(mol_->natom()*3);
  return d3natom_;
}

RefSCMatrix
MolecularHessian::cartesian_to_symmetry(const RefMolecule &mol,
                                        RefPointGroup pg,
                                        RefSCMatrixKit kit)
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
  cout << node0 << indent << "The external rank is " << rankext << endl;
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
              cerr << node0 << indent
                   << "FinDispMolecularHessian: atom mapping bad" << endl;
              abort();
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
      for (j=0; j<irrep.degeneracy(); j++) {
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
  RefSCBlockInfo bi = new SCBlockInfo(total, nirrep, dims);
  for (i=0; i<nirrep; i++) {
      bi->set_subdim(i, symmbasis[i]->coldim());
    }
  RefSCDimension dsym = new SCDimension(bi);

  RefSCDimension bd3natom = new SCDimension(3*mol->natom());
  bd3natom->blocks()->set_subdim(0,d3natom);

  RefSCMatrixKit symkit = new BlockedSCMatrixKit(kit);
  RefSCMatrix result(dsym, bd3natom, symkit);
  BlockedSCMatrix *bresult = BlockedSCMatrix::castdown(result.pointer());

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
MolecularHessian::set_energy(const RefMolecularEnergy &)
{
}

MolecularEnergy*
MolecularHessian::energy() const
{
  return 0;
}

/////////////////////////////////////////////////////////////////
// GuessMolecularHessian

#define CLASSNAME GuessMolecularHessian
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public MolecularHessian
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
GuessMolecularHessian::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = MolecularHessian::_castdown(cd);
  return do_castdowns(casts,cd);
}

GuessMolecularHessian::GuessMolecularHessian(const RefKeyVal&keyval):
  MolecularHessian(keyval)
{
  coor_ = keyval->describedclassvalue("coor");
  if (mol_.null()) mol_ = coor_->molecule();
}

GuessMolecularHessian::GuessMolecularHessian(StateIn&s):
  MolecularHessian(s)
  maybe_SavableState(s)
{
  coor_.restore_state(s);
}

GuessMolecularHessian::~GuessMolecularHessian()
{
}

void
GuessMolecularHessian::save_data_state(StateOut&s)
{
  MolecularHessian::save_data_state(s);
  coor_.save_state(s);
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

#define CLASSNAME DiagMolecularHessian
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public MolecularHessian
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
DiagMolecularHessian::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = MolecularHessian::_castdown(cd);
  return do_castdowns(casts,cd);
}

DiagMolecularHessian::DiagMolecularHessian(const RefKeyVal&keyval):
  MolecularHessian(keyval)
{
  diag_ = keyval->doublevalue("coor",KeyValValuedouble(1.0));
}

DiagMolecularHessian::DiagMolecularHessian(StateIn&s):
  MolecularHessian(s)
  maybe_SavableState(s)
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

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
