//
// molecule.cc
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
#include <string.h>
#include <stdio.h>

#include <util/misc/formio.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/formula.h>
#include <chemistry/molecule/localdef.h>
#include <math/scmat/cmatrix.h>

SavableState_REF_def(Molecule);

#define CLASSNAME Molecule
#define VERSION 5
#define PARENTS public SavableState
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
Molecule::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

Molecule::Molecule():
  r_(0), natoms_(0), Z_(0), mass_(0), labels_(0), charges_(0)
{
  pg_ = new PointGroup;
  atominfo_ = new AtomInfo();
  geometry_units_ = new Units("bohr");
  nuniq_ = 0;
  equiv_ = 0;
  nequiv_ = 0;
  atom_to_uniq_ = 0;
  init_symmetry_info();
}

Molecule::Molecule(const Molecule& mol):
 r_(0), natoms_(0), Z_(0), mass_(0), labels_(0), charges_(0)
{
  nuniq_ = 0;
  equiv_ = 0;
  nequiv_ = 0;
  atom_to_uniq_ = 0;
  *this=mol;
}

Molecule::~Molecule()
{
  clear();
}

void
Molecule::clear()
{
  if (r_) {
      delete[] r_[0];
      delete[] r_;
      r_ = 0;
    }
  if (labels_) {
      for (int i=0; i<natoms_; i++) {
          delete[] labels_[i];
        }
      delete[] labels_;
      labels_ = 0;
    }
  delete[] charges_;
  charges_ = 0;
  delete[] mass_;
  mass_ = 0;
  delete[] Z_;
  Z_ = 0;

  clear_symmetry_info();
}

Molecule::Molecule(const RefKeyVal&input):
 r_(0), natoms_(0), Z_(0), mass_(0), labels_(0), charges_(0)
{
  atominfo_ = input->describedclassvalue("atominfo");
  if (atominfo_.null()) atominfo_ = new AtomInfo;
  if (input->exists("pdb_file")) {
      geometry_units_ = new Units("angstrom");
      double ang_to_bohr = geometry_units_->to_atomic_units();
      const int LineLength = 85;
      char line[LineLength];
      char* filename = input->pcharvalue("pdb_file");
      FILE*fp = fopen(filename,"r");
      if (!fp) {
          cerr << node0 << indent
               << "Molecule::Molecule(const RefKeyVal&input): "
               << scprintf("pdb file not found: \"%s\"\n", filename);
          abort();
        }
      while(fgets(line,LineLength,fp)) {
          if (!strncmp(line,"HETA",4) || !strncmp(line,"ATOM",4)) {
              char atomsym[3];
              // find the atomic symbol
              int symletter=0, offset;
              for (offset=12; offset<16 && symletter<2; offset++) {
                  if (line[offset] != ' '
                      && (line[offset] < '0' || line[offset] > '9')) {
                      atomsym[symletter] = line[offset];
                      symletter++;
                    }
                }
              atomsym[symletter] = '\0';
              // skip dummy atoms
              if (!strcmp(atomsym,"Q")) continue;
              char position[9];
              position[8] = '\0';
              // x
              strncpy(position,&line[30],8);
              double x = atof(position);
              // y
              strncpy(position,&line[38],8);
              double y = atof(position);
              // z
              strncpy(position,&line[46],8);
              double z = atof(position);
              add_atom(AtomInfo::string_to_Z(atomsym),
                       x*ang_to_bohr, y*ang_to_bohr, z*ang_to_bohr);
            }
        }
      fclose(fp);
    }
  else {
      // check for old style units input first
      if (input->booleanvalue("angstrom")
          ||input->booleanvalue("angstroms")
          ||input->booleanvalue("aangstrom")
          ||input->booleanvalue("aangstroms")) {
          geometry_units_ = new Units("angstrom");
        }
      // check for new style units input
      else {
          geometry_units_ = new Units(input->pcharvalue("unit"),
                                      Units::Steal);
        }
        
      double conv = geometry_units_->to_atomic_units();

      // get the number of atoms and make sure that the geometry and the
      // atoms array have the same number of atoms.
      // right now we read in the unique atoms...then we will symmetrize.
      // the length of atoms must still equal the length of geometry, but
      // we'll try to set up atom_labels such that different lengths are
      // possible
      int natom = input->count("geometry");
      if (natom != input->count("atoms")) {
          cerr << node0 << indent
               << "Molecule: size of \"geometry\" != size of \"atoms\"\n";
          abort();
        }

      int i;
      for (i=0; i<natom; i++) {
          char *name, *label;
          int ghost = input->booleanvalue("ghost",i);
          double charge = input->doublevalue("charge",i);
          int have_charge = input->error() == KeyVal::OK;
          if (ghost) {
              have_charge = 1;
              charge = 0.0;
            }
          add_atom(AtomInfo::string_to_Z(name = input->pcharvalue("atoms",i)),
                   input->doublevalue("geometry",i,0)*conv,
                   input->doublevalue("geometry",i,1)*conv,
                   input->doublevalue("geometry",i,2)*conv,
                   label = input->pcharvalue("atom_labels",i),
                   input->doublevalue("mass",i),
                   have_charge, charge
              );
          delete[] name;
          delete[] label;
        }
    }

  char *symmetry = input->pcharvalue("symmetry");
  double symtol = input->doublevalue("symmetry_tolerance",
                                     KeyValValuedouble(1.0e-4));
  nuniq_ = 0;
  equiv_ = 0;
  nequiv_ = 0;
  atom_to_uniq_ = 0;
  if (symmetry && !strcmp(symmetry,"auto")) {
      set_point_group(highest_point_group(symtol), symtol*10.0);
    }
  else {
      pg_ = new PointGroup(input);

      // translate to the origin of the symmetry frame
      double r[3];
      for (int i=0; i<3; i++) {
          r[i] = -pg_->origin()[i];
          pg_->origin()[i] = 0;
        }
      translate(r);

      if (input->booleanvalue("redundant_atoms")) {
          init_symmetry_info();
          cleanup_molecule(symtol);
        }
      else {
          symmetrize();
        }
    }
  delete[] symmetry;
}

Molecule&
Molecule::operator=(const Molecule& mol)
{
  clear();

  pg_ = new PointGroup(*(mol.pg_.pointer()));
  atominfo_ = mol.atominfo_;
  geometry_units_ = new Units(mol.geometry_units_->string_rep());

  natoms_ = mol.natoms_;

  if (natoms_) {
      if (mol.mass_) {
          mass_ = new double[natoms_];
          memcpy(mass_,mol.mass_,natoms_*sizeof(double));
        }
      if (mol.charges_) {
          charges_ = new double[natoms_];
          memcpy(charges_,mol.charges_,natoms_*sizeof(double));
        }
      if (mol.labels_) {
          labels_ = new char *[natoms_];
          for (int i=0; i<natoms_; i++) {
              if (mol.labels_[i]) {
                  labels_[i] = strcpy(new char[strlen(mol.labels_[i])+1],
                                      mol.labels_[i]);
                }
              else labels_[i] = 0;
            }
        }
      r_ = new double*[natoms_];
      r_[0] = new double[natoms_*3];
      for (int i=0; i<natoms_; i++) {
          r_[i] = &(r_[0][i*3]);
        }
      memcpy(r_[0], mol.r_[0], natoms_*3*sizeof(double));
      Z_ = new int[natoms_];
      memcpy(Z_, mol.Z_, natoms_*sizeof(int));
    }

  init_symmetry_info();

  return *this;
}

void
Molecule::add_atom(int Z,double x,double y,double z,
                   const char *label,double mass,
                   int have_charge, double charge)
{
  int i;

  // allocate new arrays
  int *newZ = new int[natoms_+1];
  double **newr = new double*[natoms_+1];
  double *newr0 = new double[(natoms_+1)*3];
  char **newlabels = 0;
  if (label || labels_) {
      newlabels = new char*[natoms_+1];
    }
  double *newcharges = 0;
  if (have_charge || charges_) {
      newcharges = new double[natoms_+1];
    }
  double *newmass = 0;
  if (mass_ || mass != 0.0) {
      newmass = new double[natoms_+1];
    }

  // setup the r_ pointers
  for (i=0; i<=natoms_; i++) {
      newr[i] = &(newr0[i*3]);
    }

  // copy old data to new arrays
  if (natoms_) {
      memcpy(newZ,Z_,sizeof(int)*natoms_);
      memcpy(newr0,r_[0],sizeof(double)*natoms_*3);
      if (labels_) {
          memcpy(newlabels,labels_,sizeof(char*)*natoms_);
        }
      else if (newlabels) {
          memset(newlabels,0,sizeof(char*)*natoms_);
        }
      if (charges_) {
          memcpy(newcharges,charges_,sizeof(double)*natoms_);
        }
      else if (newcharges) {
          for (i=0; i<natoms_; i++) newcharges[i] = Z_[i];
        }
      if (mass_) {
          memcpy(newmass,mass_,sizeof(double)*natoms_);
        }
      else if (newmass) {
          memset(newmass,0,sizeof(double)*natoms_);
        }
    }

  // delete old data
  delete[] Z_;
  if (r_) {
      delete[] r_[0];
      delete[] r_;
    }
  delete[] labels_;
  delete[] charges_;
  delete[] mass_;

  // setup new pointers
  Z_ = newZ;
  r_ = newr;
  labels_ = newlabels;
  charges_ = newcharges;
  mass_ = newmass;

  // copy info for this atom into arrays
  Z_[natoms_] = Z;
  r_[natoms_][0] = x;
  r_[natoms_][1] = y;
  r_[natoms_][2] = z;
  if (mass_) mass_[natoms_] = mass;
  if (label) {
      labels_[natoms_] = strcpy(new char[strlen(label)+1],label);
    }
  else if (labels_) {
      labels_[natoms_] = 0;
    }
  if (have_charge) {
      charges_[natoms_] = charge;
    }
  else if (charges_) {
      charges_[natoms_] = Z;
    }

  natoms_++;
}

void
Molecule::print(ostream& os)
{
  int i;

  double conv = geometry_units_->from_atomic_units();

  // Be careful here, since MolecularFormula requires a smart
  // pointer that might try to delete this if i'm not referenced.
  reference();
  MolecularFormula *mf = new MolecularFormula(this);
  os << node0 << indent
     << "Molecular formula: " << mf->formula() << endl;
  delete mf;
  dereference();

  os << node0 << indent << "molecule<Molecule>: (" << endl;
  os << incindent;
  pg_->print(os);
  if (geometry_units_->string_rep()) {
      os << node0 << indent
         << "unit = \"" << geometry_units_->string_rep() << "\""
         << endl;
    }
  os << node0 << indent
     << scprintf("{%3s", "n")
     << scprintf(" %5s", "atoms");
  if (labels_) os << node0 << scprintf(" %11s", "atom_labels");
  int int_charges = 1;
  if (charges_) {
      for (i=0;i<natom();i++) if (charges_[i]!=(int)charges_[i]) int_charges=0;
      if (int_charges) {
          os << node0 << scprintf(" %7s", "charges");
        }
      else {
          os << node0 << scprintf(" %17s", "charges");
        }
    }
  os << node0
     << scprintf("  %16s", "")
     << scprintf(" %16s", "geometry   ")
     << scprintf(" %16s ", "");
  os << node0 << "}={" << endl;
  for (i=0; i<natom(); i++) {
      os << node0 << indent
         << scprintf(" %3d", i+1)
         << scprintf(" %5s", AtomInfo::symbol(Z_[i]));
      if (labels_) {
          const char *lab = labels_[i];
          if (lab == 0) lab = "";
          char  *qlab = new char[strlen(lab)+3];
          strcpy(qlab,"\"");
          strcat(qlab,lab);
          strcat(qlab,"\"");
          os << node0
             << scprintf(" %11s",qlab);
          delete[] qlab;
        }
      if (charges_) {
          if (int_charges) os << node0 << scprintf(" %7.4f", charges_[i]);
          else os << node0 << scprintf(" %17.15f", charges_[i]);
        }
      os << node0
         << scprintf(" [% 16.10f", conv * r(i,0))
         << scprintf(" % 16.10f", conv * r(i,1))
         << scprintf(" % 16.10f]", conv * r(i,2))
         << endl;
    }
  os << node0 << indent << "}" << endl;
  os << decindent;
  os << node0 << indent << ")" << endl;

  os << node0 << indent << "Atomic Masses:" << endl;
  for (i=0; i<natom(); i+=5) {
      os << node0 << indent;
      for (int j=i; j<i+5 && j<natom(); j++) {
          os << node0 << scprintf(" %10.5f", mass(j));
        }
      os << node0 << endl;
    }
}

int
Molecule::atom_label_to_index(const char *l) const
{
  int i;
  for (i=0; i<natom(); i++) {
      if (label(i) && !strcmp(l,label(i))) return i;
    }
  return -1;
}

double*
Molecule::charges() const
{
  double*result = new double[natoms_];
  if (charges_) {
      memcpy(result, charges_, sizeof(double)*natom());
    }
  else {
      for (int i=0; i<natom(); i++) result[i] = Z_[i];
    }
  return result;
}

double
Molecule::charge(int iatom) const
{
  if (charges_) return charges_[iatom];
  return Z_[iatom];
}

double
Molecule::nuclear_charge() const
{
  int i;
  double c = 0.0;
  for (i=0; i<natom(); i++) {
      c += charge(i);
    }
  return c;
}

void Molecule::save_data_state(StateOut& so)
{
  so.put(natoms_);
  pg_.save_state(so);
  geometry_units_.save_state(so);
  atominfo_.save_state(so);
  if (natoms_) {
      so.put(Z_, natoms_);
      so.put_array_double(r_[0], natoms_*3);
      so.put(charges_,natoms_);
    }
  if (mass_) {
      so.put(1);
      so.put_array_double(mass_, natoms_);
    }
  else {
      so.put(0);
    }
  if (labels_){
      so.put(1);
      for (int i=0; i<natoms_; i++) {
          so.putstring(labels_[i]);
        }
    }
  else {
      so.put(0);
    }
}

Molecule::Molecule(StateIn& si):
  r_(0), natoms_(0), Z_(0), mass_(0), labels_(0),
  SavableState(si)
{
  if (si.version(static_class_desc()) < 4) {
      cerr << "Molecule: cannot restore from old molecules" << endl;
      abort();
    }
  si.get(natoms_);
  pg_.restore_state(si);
  geometry_units_.restore_state(si);
  atominfo_.restore_state(si);
  if (natoms_) {
      si.get(Z_);
      r_ = new double*[natoms_];
      r_[0] = new double[natoms_*3];
      si.get_array_double(r_[0],natoms_*3);
      for (int i=1; i<natoms_; i++) {
          r_[i] = &(r_[0][i*3]);
        }
      if (si.version(static_class_desc()) > 4) {
          si.get(charges_);
        }
      else {
          charges_ = 0;
        }
    }
  int test;
  si.get(test);
  if (test) {
      mass_ = new double[natoms_];
      si.get_array_double(mass_, natoms_);
    }
  si.get(test);
  if (test){
      labels_ = new char*[natoms_];
      for (int i=0; i<natoms_; i++) {
          si.getstring(labels_[i]);
        }
    }

  nuniq_ = 0;
  equiv_ = 0;
  nequiv_ = 0;
  atom_to_uniq_ = 0;
  init_symmetry_info();
}

void
Molecule::set_point_group(const RefPointGroup&ppg, double tol)
{
  cout << node0 << indent
       << "Molecule: setting point group to " << ppg->symbol()
       << endl;
  pg_ = new PointGroup(*ppg.pointer());

  double r[3];
  for (int i=0; i<3; i++) {
      r[i] = -pg_->origin()[i];
      pg_->origin()[i] = 0;
    }
  translate(r);

  clear_symmetry_info();
  init_symmetry_info();

  cleanup_molecule(tol);
}

RefPointGroup
Molecule::point_group() const
{
  return pg_;
}

SCVector3
Molecule::center_of_mass() const
{
  SCVector3 ret;
  double M;

  ret = 0.0;
  M = 0.0;

  for (int i=0; i < natom(); i++) {
    double m = mass(i);
    ret += m * SCVector3(r(i));
    M += m;
  }

  ret *= 1.0/M;

  return ret;
}

double
Molecule::nuclear_repulsion_energy()
{
  int i, j;
  double e=0.0;

  for (i=1; i < natoms_; i++) {
    SCVector3 ai(r(i));
    double Zi = charge(i);
    
    for (j=0; j < i; j++) {
        SCVector3 aj(r(j));
        e += Zi * charge(j) / ai.dist(aj);
      }
    }

  return e;
}

void
Molecule::nuclear_repulsion_1der(int center, double xyz[3])
{
  int i,j,k;
  double rd[3],r2;
  double factor;

  xyz[0] = 0.0;
  xyz[1] = 0.0;
  xyz[2] = 0.0;
  for (i=0; i < natoms_; i++) {
      SCVector3 ai(r(i));
      double Zi = charge(i);

      for (j=0; j < i; j++) {
          if (center==i || center==j) {
              SCVector3 aj(r(j));

              r2 = 0.0;
              for (k=0; k < 3; k++) {
                  rd[k] = ai[k] - aj[k];
                  r2 += rd[k]*rd[k];
                }
        
              factor = - Zi * charge(j) * pow(r2,-1.5);
              if (center==j) factor = -factor;
              for (k=0; k<3; k++) {
                  xyz[k] += factor * rd[k];
                }
            }
        }
    }
}

void
Molecule::nuclear_efield(const double *position, double *efield)
{
  int i,j;
  double tmp;
  double rd[3];

  for (i=0; i<3; i++) efield[i] = 0.0;

  for (i=0; i<natoms_; i++) {
      SCVector3 a(r(i));
      tmp = 0.0;
      for (j=0; j<3; j++) {
          rd[j] = position[j] - a[j];
          tmp += rd[j]*rd[j];
        }
      tmp = charge(i)/(tmp*sqrt(tmp));
      for (j=0; j<3; j++) {
          efield[j] +=  rd[j] * tmp;
        }
    }
}

int
Molecule::atom_at_position(double *v, double tol) const
{
  SCVector3 p(v);
  for (int i=0; i < natom(); i++) {
      SCVector3 ai(r(i));
      if (p.dist(ai) < tol) return i;
    }
  return -1;
}

// We are given a molecule which may or may not have just the symmetry
// distinct atoms in it.  We have to go through the existing set of atoms,
// perform each symmetry operation in the point group on each of them, and
// then add the new atom if it isn't in the list already

void
Molecule::symmetrize(double tol)
{
  // if molecule is c1, don't do anything
  if (!strcmp(this->point_group()->symbol(),"c1")) {
    init_symmetry_info();
    return;
    }

  clear_symmetry_info();

  Molecule *newmol = new Molecule(*this);
  
  CharacterTable ct = this->point_group()->char_table();

  SCVector3 np;
  SymmetryOperation so;
  
  for (int i=0; i < natom(); i++) {
    SCVector3 ac(r(i));

    for (int g=0; g < ct.order(); g++) {
      so = ct.symm_operation(g);
      for (int ii=0; ii < 3; ii++) {
        np[ii]=0;
        for (int jj=0; jj < 3; jj++) np[ii] += so(ii,jj) * ac[jj];
      }

      int atom = newmol->atom_at_position(np.data(), tol);
      if (atom < 0) {
        newmol->add_atom(Z_[i],np[0],np[1],np[2],label(i));
      }
      else {
        if (Z(i) != newmol->Z(atom)
            || fabs(mass(i)-newmol->mass(atom))>1.0e-10) {
            cerr << node0 << "Molecule: symmetrize: atom mismatch" << endl;
            abort();
        }
      }
    }
  }
  
  RefUnits saved_units = geometry_units_;
  *this = *newmol;
  geometry_units_ = saved_units;
  delete newmol;

  init_symmetry_info();
}

void
Molecule::translate(const double *r)
{
  for (int i=0; i < natom(); i++) {
    r_[i][0] += r[0];
    r_[i][1] += r[1];
    r_[i][2] += r[2];
  }
}

// move the molecule to the center of mass
void
Molecule::move_to_com()
{
  SCVector3 com = -center_of_mass();
  translate(com.data());
}

// find the 3 principal coordinate axes, and rotate the molecule to be 
// aligned along them.  also rotate the symmetry frame contained in point_group
void
Molecule::transform_to_principal_axes(int trans_frame)
{
  // mol_move_to_com(mol);

  double *inert[3], inert_dat[9], *evecs[3], evecs_dat[9];
  double evals[3];

  int i;
  for (i=0; i < 3; i++) {
    inert[i] = &inert_dat[i*3];
    evecs[i] = &evecs_dat[i*3];
  }
  memset(inert_dat,'\0',sizeof(double)*9);
  memset(evecs_dat,'\0',sizeof(double)*9);

  for (i=0; i < natom(); i++) {
    SCVector3 ac(r(i));
    double m=mass(i);
    inert[0][0] += m * (ac[1]*ac[1] + ac[2]*ac[2]);
    inert[1][0] -= m * ac[0]*ac[1];
    inert[1][1] += m * (ac[0]*ac[0] + ac[2]*ac[2]);
    inert[2][0] -= m * ac[0]*ac[2];
    inert[2][1] -= m * ac[1]*ac[2];
    inert[2][2] += m * (ac[0]*ac[0] + ac[1]*ac[1]);
  }

  inert[0][1] = inert[1][0];
  inert[0][2] = inert[2][0];
  inert[1][2] = inert[2][1];

 // cleanup inert
  for (i=0; i < 3; i++) {
    for (int j=0; j <= i; j++) {
      if (fabs(inert[i][j]) < 1.0e-5) {
        inert[i][j]=inert[j][i]=0.0;
      }
    }
  }

  cmat_diag(inert, evals, evecs, 3, 1, 1e-14);

 // cleanup evecs
  for (i=0; i < 3; i++) {
    for (int j=0; j < 3; j++) {
      if (fabs(evecs[i][j]) < 1.0e-5) {
        evecs[i][j]=0.0;
      }
    }
  }

  double x,y,z;
  for (i=0; i < natom(); i++) {
    x = r(i,0); y = r(i,1); z = r(i,2);

    r_[i][0] = evecs[0][0]*x + evecs[1][0]*y + evecs[2][0]*z;
    r_[i][1] = evecs[0][1]*x + evecs[1][1]*y + evecs[2][1]*z;
    r_[i][2] = evecs[0][2]*x + evecs[1][2]*y + evecs[2][2]*z;
  }

  if (!trans_frame) return;
  
  SymmetryOperation tso=point_group()->symm_frame();

  for (i=0; i < 3; i++) {
    for (int j=0; j < 3; j++) {
      double t=0;
      for (int k=0; k < 3; k++) t += tso[i][k]*evecs[k][j];
      if (fabs(t) < 0.5)
        t = 0;
      else if (fabs(t) >= .5)
        t = 1;
      
      pg_->symm_frame()[i][j] = t;
    }
  }
}

// given a molecule, make sure that equivalent centers have coordinates
// that really map into each other

void
Molecule::cleanup_molecule(double tol)
{
  // if symmetry is c1, do nothing else
  if (!strcmp(point_group()->symbol(),"c1")) return;

  int i;
  SCVector3 up,np,ap;
  SymmetryOperation so;
  CharacterTable ct = point_group()->char_table();

  // first clean up the unique atoms by replacing each coordinate with the
  // average of coordinates obtained by applying all symmetry operations to
  // the original atom, iff the new atom ends up near the original atom
  for (i=0; i < nunique(); i++) {
      // up will store the original coordinates of unique atom i
      up = r(unique(i));
      // ap will hold the average coordinate (times the number of coordinates)
      // initialize it to the E result
      ap = up;
      int ncoor = 1;
      // loop through all sym ops except E
      for (int g=1; g < ct.order(); g++) {
          so = ct.symm_operation(g);
          for (int ii=0; ii < 3; ii++) {
              np[ii]=0;
              for (int jj=0; jj < 3; jj++) np[ii] += so(ii,jj) * up[jj];
            }
          if (np.dist(up) < 0.1) {
              for (int jj=0; jj < 3; jj++) ap[jj] += np[jj];
              ncoor++;
            }
        }
      // replace the unique coordinate with the average coordinate
      r_[unique(i)][0] = ap[0] / ncoor;
      r_[unique(i)][1] = ap[1] / ncoor;
      r_[unique(i)][2] = ap[2] / ncoor;
    }

  // find the atoms equivalent to each unique atom and eliminate
  // numerical errors that may be in the equivalent atom's coordinates

  // loop through unique atoms
  for (i=0; i < nunique(); i++) {
      // up will store the coordinates of unique atom i
      up = r(unique(i));

      // loop through all sym ops except E
      for (int g=1; g < ct.order(); g++) {
          so = ct.symm_operation(g);
          for (int ii=0; ii < 3; ii++) {
              np[ii]=0;
              for (int jj=0; jj < 3; jj++) np[ii] += so(ii,jj) * up[jj];
            }

          // loop through equivalent atoms
          int found = 0;
          for (int j=0; j < natom(); j++) {
              // see if j is generated from i
              if (np.dist(SCVector3(r(j))) < tol) {
                  r_[j][0] = np[0];
                  r_[j][1] = np[1];
                  r_[j][2] = np[2];
                  found = 1;
                }
            }
          if (!found) {
              cerr << node0
                   << "Molecule: cleanup: couldn't find atom at " << np << endl
                   << "  transforming uniq atom " << i << " at " << up << endl
                   << "  with symmetry op " << g << ":" << endl;
              so.print(cerr << node0);
              abort();
            }
        }
    }

}

///////////////////////////////////////////////////////////////////
// Compute the principal axes and the principal moments of inertia
///////////////////////////////////////////////////////////////////

void
Molecule::principal_moments_of_inertia(double *evals, double **evecs) const
{

  // The principal moments of inertia are computed in amu*angstrom^2
  // evals: principal moments of inertia
  // evecs: principal axes (optional argument)

  const double au_to_angs = 0.2800283608302436; // for moments of inertia

  double *inert[3];  // inertia tensor

  int i, j;
  int delete_evecs = 0;

  // (allocate and) initialize evecs, evals, and inert
  if (!evecs) {
    evecs = new double*[3];
    for (i=0; i<3; i++) evecs[i] = new double[3];
    delete_evecs = 1;
    }
  for (i=0; i<3; i++) {
    inert[i] = new double[3];
    memset(inert[i],'\0',sizeof(double)*3);
    memset(evecs[i],'\0',sizeof(double)*3);
    }
  memset(evals,'\0',sizeof(double)*3);

  SCVector3 com = center_of_mass();

  // compute inertia tensor
  SCVector3 ac;
  for (i=0; i<natom(); i++) {
    ac = r(i);
    // compute moments of inertia wrt center of mass
    for (j=0; j<3; j++) ac(j) -= com(j);
    double m=au_to_angs*mass(i);
    inert[0][0] += m * (ac[1]*ac[1] + ac[2]*ac[2]);
    inert[1][0] -= m * ac[0]*ac[1];
    inert[1][1] += m * (ac[0]*ac[0] + ac[2]*ac[2]);
    inert[2][0] -= m * ac[0]*ac[2];
    inert[2][1] -= m * ac[1]*ac[2];
    inert[2][2] += m * (ac[0]*ac[0] + ac[1]*ac[1]);
    }
  inert[0][1] = inert[1][0];
  inert[0][2] = inert[2][0];
  inert[1][2] = inert[2][1];

  cmat_diag(inert, evals, evecs, 3, 1, 1e-14);

  if (delete_evecs) {
    for (i=0; i<3; i++) delete[] evecs[i];
    delete[] evecs;
    }
  for (i=0; i<3; i++) {
    delete[] inert[i];
    }
}

int
Molecule::n_core_electrons()
{
  int i,n=0;
  for (i=0; i<natom(); i++) {
      int z = Z_[i];
      if (z > 2) n += 2;
      if (z > 10) n += 8;
      if (z > 18) n += 8;
      if (z > 30) n += 10;
      if (z > 36) n += 8;
      if (z > 48) n += 10;
      if (z > 54) n += 8;
      if (z > 72) {
          cerr << "Molecule::n_core_electrons: atomic number too large"
               << endl;
          abort();
        }
    }
  return n;
}

int
Molecule::max_z()
{
  int i, maxz=0;
  for (i=0; i<natom(); i++) {
      int z = Z_[i];
      if (z>maxz) maxz = z;
    }
  return maxz;
}

void
Molecule::print_pdb(ostream& os, char *title)
{
  RefUnits u = new Units("angstrom");
  double bohr = u->from_atomic_units();

  if (title)
      os << node0 << scprintf("%-10s%-60s\n","COMPND",title);
  else
      os << node0 << scprintf("%-10s%-60s\n","COMPND","Title");

  if (title)
      os << node0 << scprintf("REMARK   %s\n", title);

  int i;
  for (i=0; i < natom(); i++) {
    char symb[4];
    sprintf(symb,"%s1",AtomInfo::symbol(Z_[i]));

    os << node0 << scprintf(
        "HETATM%5d  %-3s UNK %5d    %8.3f%8.3f%8.3f  0.00  0.00   0\n",
        i+1, symb, 0, r(i,0)*bohr, r(i,1)*bohr, r(i,2)*bohr);
  }

  for (i=0; i < natom(); i++) {
    double at_rad_i = atominfo_->atomic_radius(Z_[i]);
    SCVector3 ai(r(i));

    os << node0 << scprintf("CONECT%5d",i+1);

    for (int j=0; j < natom(); j++) {

      if (j==i) continue;

      double at_rad_j = atominfo_->atomic_radius(Z_[j]);
      SCVector3 aj(r(j));

      if (ai.dist(aj) < 1.1*(at_rad_i+at_rad_j))
          os << node0 << scprintf("%5d",j+1);
    }

    os << node0 << endl;
  }

  os << node0 << "END" << endl;
  os.flush();
}

double
Molecule::mass(int atom) const
{
  if (!mass_ || mass_[atom] == 0) {
      return atominfo_->mass(Z_[atom]);
    }
  return mass_[atom];
}

const char *
Molecule::label(int atom) const
{
  if (!labels_) return 0;
  return labels_[atom];
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
