//
// petite.cc --- implementation of the PetiteList class
//
// Copyright (C) 1996 Limit Point Systems, Inc.
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

#include <stdio.h>

#include <util/misc/formio.h>
#include <chemistry/molecule/localdef.h>

#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/shellrot.h>

using namespace std;
using namespace sc;

////////////////////////////////////////////////////////////////////////////

int **
sc::compute_atom_map(const Ref<GaussianBasisSet> &basis)
{
  // grab references to the Molecule and BasisSet for convenience
  GaussianBasisSet& gbs = *basis.pointer();
  Molecule& mol = *gbs.molecule().pointer();

  // create the character table for the point group
  CharacterTable ct = mol.point_group()->char_table();

  int natom = gbs.ncenter();
  int ng = ct.order();
  int **atom_map;
  atom_map = new int*[natom];
  for (int i=0; i < natom; i++) atom_map[i] = new int[ng];
  
  double np[3];
  SymmetryOperation so;

  // loop over all centers
  for (int i=0; i < natom; i++) {
    SCVector3 ac(mol.r(i));

    // then for each symop in the pointgroup, transform the coordinates of
    // center "i" and see which atom it maps into
    for (int g=0; g < ng; g++) {
      so = ct.symm_operation(g);

      for (int ii=0; ii < 3; ii++) {
        np[ii] = 0;
        for (int jj=0; jj < 3; jj++)
          np[ii] += so(ii,jj) * ac[jj];
      }

      atom_map[i][g] = mol.atom_at_position(np, 0.05);
      if (atom_map[i][g] < 0) {
        ExEnv::out0() << "ERROR: Symmetry operation " << g << " did not map atom "
             << i+1 << " to another atom:" << endl;
        ExEnv::out0() << indent << "Molecule:" << endl;
        ExEnv::out0() << incindent;
        mol.print();
        ExEnv::out0() << decindent;
        ExEnv::out0() << indent << "attempted to find atom at" << endl;
        ExEnv::out0() << incindent;
        ExEnv::out0() << indent << np[0] << " " << np[1] << " " << np[2] << endl;
        abort();
      }
    }
  }

  return atom_map;
}
  

void
sc::delete_atom_map(int **atom_map, const Ref<GaussianBasisSet> &basis)
{
  if (atom_map) {
    int natom = basis->ncenter();
    for (int i=0; i < natom; i++)
      delete[] atom_map[i];
    delete[] atom_map;
  }
}

int **
sc::compute_shell_map(int **atom_map, const Ref<GaussianBasisSet> &basis)
{
  int **shell_map;

  GaussianBasisSet& gbs = *basis.pointer();
  Molecule& mol = *gbs.molecule().pointer();

  // create the character table for the point group
  CharacterTable ct = mol.point_group()->char_table();

  int natom = gbs.ncenter();
  int ng = ct.order();

  int nshell = basis->nshell();
  shell_map = new int*[nshell];
  for (int i=0; i < nshell; i++)
    shell_map[i] = new int[ng];

  for (int i=0; i<natom; i++) {
    // hopefully, shells on equivalent centers will be numbered in the same
    // order
    for (int s=0; s < gbs.nshell_on_center(i); s++) {
      int shellnum = gbs.shell_on_center(i,s);
      for (int g=0; g < ng; g++) {
        shell_map[shellnum][g] = gbs.shell_on_center(atom_map[i][g],s);
      }
    }
  }

  return shell_map;
}

void
sc::delete_shell_map(int **shell_map, const Ref<GaussianBasisSet> &basis)
{
  int nshell = basis->nshell();
  if (shell_map) {
    for (int i=0; i < nshell; i++)
      delete[] shell_map[i];
    delete[] shell_map;
  }
}

////////////////////////////////////////////////////////////////////////////

PetiteList::PetiteList(const Ref<GaussianBasisSet> &gbs,
                       const Ref<Integral>& ints) :
  gbs_(gbs),
  ints_(ints->clone())
{
  ints_->set_basis(gbs_);
  init();
}

PetiteList::~PetiteList()
{
  if (p1_)
    delete[] p1_;

  if (lamij_)
    delete[] lamij_;

  if (nbf_in_ir_)
    delete[] nbf_in_ir_;
  
  if (atom_map_) {
    for (int i=0; i < natom_; i++)
      delete[] atom_map_[i];
    delete[] atom_map_;
  }

  if (shell_map_) {
    for (int i=0; i < nshell_; i++)
      delete[] shell_map_[i];
    delete[] shell_map_;
  }

  natom_=0;
  nshell_=0;
  ng_=0;
  nblocks_=0;
  nirrep_=0;
  p1_=0;
  atom_map_=0;
  shell_map_=0;
  lamij_=0;
  nbf_in_ir_=0;
}

void
PetiteList::init()
{
  int i;

  // grab references to the Molecule and BasisSet for convenience
  GaussianBasisSet& gbs = *gbs_.pointer();
  Molecule& mol = *gbs.molecule().pointer();

  // create the character table for the point group
  CharacterTable ct = mol.point_group()->char_table();
  
  // initialize private members
  c1_=0;
  ng_ = ct.order();
  natom_ = mol.natom();
  nshell_ = gbs.nshell();
  nirrep_ = ct.nirrep();

  // if point group is C1, then zero everything
  if (ng_==1) {
    c1_=1;
    nblocks_=1;

    p1_=0;
    atom_map_=0;
    shell_map_=0;
    lamij_=0;
    nbf_in_ir_=0;
    return;
  }
  
  // allocate storage for arrays
  p1_ = new char[nshell_];
  lamij_ = new char[i_offset(nshell_)];

  atom_map_ = new int*[natom_];
  for (i=0; i < natom_; i++)
    atom_map_[i] = new int[ng_];
  
  shell_map_ = new int*[nshell_];
  for (i=0; i < nshell_; i++)
    shell_map_[i] = new int[ng_];
  
  // set up atom and shell mappings
  double np[3];
  SymmetryOperation so;
  
  // loop over all centers
  for (i=0; i < natom_; i++) {
    SCVector3 ac(mol.r(i));

    // then for each symop in the pointgroup, transform the coordinates of
    // center "i" and see which atom it maps into
    for (int g=0; g < ng_; g++) {
      so = ct.symm_operation(g);

      for (int ii=0; ii < 3; ii++) {
        np[ii] = 0;
        for (int jj=0; jj < 3; jj++)
          np[ii] += so(ii,jj) * ac[jj];
      }

      atom_map_[i][g] = mol.atom_at_position(np, 0.05);
      if (atom_map_[i][g] < 0) {
        ExEnv::out0() << "ERROR: Symmetry operation " << g << " did not map atom "
             << i+1 << " to another atom:" << endl;
        ExEnv::out0() << indent << "Molecule:" << endl;
        ExEnv::out0() << incindent;
        mol.print();
        ExEnv::out0() << decindent;
        ExEnv::out0() << indent << "attempted to find atom at" << endl;
        ExEnv::out0() << incindent;
        ExEnv::out0() << indent << np[0] << " " << np[1] << " " << np[2] << endl;
        abort();
      }
    }

    // hopefully, shells on equivalent centers will be numbered in the same
    // order
    for (int s=0; s < gbs.nshell_on_center(i); s++) {
      int shellnum = gbs.shell_on_center(i,s);
      for (int g=0; g < ng_; g++) {
        shell_map_[shellnum][g] = gbs.shell_on_center(atom_map_[i][g],s);
      }
    }
  }

  memset(p1_,0,nshell_);
  memset(lamij_,0,i_offset(nshell_));
  
  // now we do p1_ and lamij_
  for (i=0; i < nshell_; i++) {
    int g;

    // we want the highest numbered shell in a group of equivalent shells
    for (g=0; g < ng_; g++)
      if (shell_map_[i][g] > i)
        break;
    
    if (g < ng_)
      continue;
    
    // i is in the group P1
    p1_[i] = 1;

    for (int j=0; j <= i; j++) {
      int ij = i_offset(i)+j;
      int nij = 0;

      // test to see if IJ is in the group P2, if it is, then set lambda(ij)
      // equal to the number of equivalent shell pairs.  This number is
      // just the order of the group divided by the number of times ij is
      // mapped into itself
      int gg;
      for (gg=0; gg < ng_; gg++) {
        int gi = shell_map_[i][gg];
        int gj = shell_map_[j][gg];
        int gij = ij_offset(gi,gj);
        if (gij > ij)
          break;
        else if (gij == ij)
          nij++;
      }

      if (gg < ng_)
        continue;

      lamij_[ij] = (char) (ng_/nij);
    }
  }

  // form reducible representation of the basis functions
  double *red_rep = new double[ng_];
  memset(red_rep,0,sizeof(double)*ng_);
  
  for (i=0; i < natom_; i++) {
    for (int g=0; g < ng_; g++) {
      so = ct.symm_operation(g);
      int j= atom_map_[i][g];

      if (i!=j)
        continue;
      
      for (int s=0; s < gbs.nshell_on_center(i); s++) {
        for (int c=0; c < gbs(i,s).ncontraction(); c++) {
          int am=gbs(i,s).am(c);

          if (am==0)
            red_rep[g] += 1.0;
          else {
            ShellRotation r(am,so,ints_,gbs(i,s).is_pure(c));
            red_rep[g] += r.trace();
          }
        }
      }
    }
  }

  // and then use projection operators to figure out how many SO's of each
  // symmetry type there will be
  nblocks_ = 0;
  nbf_in_ir_ = new int[nirrep_];
  for (i=0; i < nirrep_; i++) {
    double t=0;
    for (int g=0; g < ng_; g++)
      t += ct.gamma(i).character(g)*red_rep[g];

    nbf_in_ir_[i] = ((int) (t+0.5))/ng_;
    if (ct.gamma(i).complex()) {
      nblocks_++;
      nbf_in_ir_[i] *= 2;
    } else {
      nblocks_ += ct.gamma(i).degeneracy();
    }
  }

  delete[] red_rep;
}

RefSCDimension
PetiteList::AO_basisdim()
{
  if (c1_)
    return SO_basisdim();
  
  RefSCDimension dim = new SCDimension(gbs_->nbasis(),1);
  dim->blocks()->set_subdim(0, gbs_->basisdim());
  return dim;
}

RefSCDimension
PetiteList::SO_basisdim()
{
  int i, j, ii;
  
  // grab a reference to the basis set
  GaussianBasisSet& gbs = *gbs_.pointer();
  
  // create the character table for the point group
  CharacterTable ct = gbs.molecule()->point_group()->char_table();

  // ncomp is the number of symmetry blocks we have
  int ncomp=nblocks();
  
  // saoelem is the current SO in a block
  int *nao = new int [ncomp];
  memset(nao,0,sizeof(int)*ncomp);

  if (c1_)
    nao[0] = gbs.nbasis();
  else {
    for (i=ii=0; i < nirrep_; i++) {
      int je = ct.gamma(i).complex() ? 1 : ct.gamma(i).degeneracy();
      for (j=0; j < je; j++,ii++)
        nao[ii] = nbf_in_ir_[i];
    }
  }

  RefSCDimension ret = new SCDimension(gbs.nbasis(),ncomp,nao);
  delete[] nao;

  for (i=ii=0; i < nirrep_; i++) {
    int nbas=(c1_) ? gbs.nbasis() : nbf_in_ir_[i];
    int je = ct.gamma(i).complex() ? 1 : ct.gamma(i).degeneracy();
    for (j=0; j < je; j++,ii++) {
      ret->blocks()->set_subdim(ii, new SCDimension(nbas));
    }

  }

  return ret;
}

void
PetiteList::print(ostream& o, int verbose)
{
  int i;

  o << indent << "PetiteList:" << endl << incindent;

  if (c1_) {
    o << indent << "is c1\n" << decindent;
    return;
  }
  
  if (verbose) {
    o
      << indent << "natom_ = " << natom_ << endl
      << indent << "nshell_ = " << nshell_ << endl
      << indent << "ng_ = " << ng_ << endl
      << indent << "nirrep_ = " << nirrep_ << endl << endl
      << indent << "atom_map_ =" << endl << incindent;

    for (i=0; i < natom_; i++) {
      o << indent;
      for (int g=0; g < ng_; g++)
        o << scprintf("%5d ",atom_map_[i][g]);
      o << endl;
    }

    o << endl << decindent
      << indent << "shell_map_ =" << endl << incindent;
    for (i=0; i < nshell_; i++) {
      o << indent;
      for (int g=0; g < ng_; g++)
        o << scprintf("%5d ",shell_map_[i][g]);
      o << endl;
    }

    o << endl << decindent
      << indent << "p1_ =" << endl << incindent;
    for (i=0; i < nshell_; i++)
      o << indent << scprintf("%5d\n",p1_[i]);

    o << decindent
      << indent << "lamij_ =" << endl << incindent;
    for (i=0; i < nshell_; i++) {
      o << indent;
      for (int j=0; j <= i; j++)
        o << scprintf("%5d ",lamij_[i_offset(i)+j]);
      o << endl;
    }
    o << endl << decindent;
  }

  CharacterTable ct = gbs_->molecule()->point_group()->char_table();
  for (i=0; i < nirrep_; i++)
    o << indent 
      << scprintf("%5d functions of %s symmetry\n",
                  nbf_in_ir_[i], ct.gamma(i).symbol());
}

// forms the basis function rotation matrix for the g'th symmetry operation
// in the point group
RefSCMatrix
PetiteList::r(int g)
{
  // grab a reference to the basis set
  GaussianBasisSet& gbs = *gbs_.pointer();
  
  SymmetryOperation so =
    gbs.molecule()->point_group()->char_table().symm_operation(g);

  RefSCMatrix ret = gbs.matrixkit()->matrix(gbs.basisdim(), gbs.basisdim());
  ret.assign(0.0);
  
  // this should be replaced with an element op at some point
  if (c1_) {
    for (int i=0; i < gbs.nbasis(); i++)
      ret.set_element(i,i,1.0);
    return ret;

  } else {
    for (int i=0; i < natom_; i++) {
      int j = atom_map_[i][g];

      for (int s=0; s < gbs.nshell_on_center(i); s++) {
        int func_i = gbs.shell_to_function(gbs.shell_on_center(i,s));
        int func_j = gbs.shell_to_function(gbs.shell_on_center(j,s));
      
        for (int c=0; c < gbs(i,s).ncontraction(); c++) {
          int am=gbs(i,s).am(c);

          if (am==0) {
            ret.set_element(func_j,func_i,1.0);
          } else {
            ShellRotation rr(am,so,ints_,gbs(i,s).is_pure(c));
            for (int ii=0; ii < rr.dim(); ii++)
              for (int jj=0; jj < rr.dim(); jj++)
                ret.set_element(func_j+jj,func_i+ii,rr(ii,jj));
          }

          func_i += gbs(i,s).nfunction(c);
          func_j += gbs(i,s).nfunction(c);
        }
      }
    }
  }
  return ret;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
