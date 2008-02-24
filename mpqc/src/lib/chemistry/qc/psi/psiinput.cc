//
// psiinput.cc
//
// Copyright (C) 2002 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <iostream>
#include <sstream>
#include <cmath>

#include <util/misc/string.h>
#include <util/misc/formio.h>
#include <util/class/scexception.h>
#include <math/symmetry/corrtab.h>
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/atominfo.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/psi/psiexenv.h>
#include <chemistry/qc/psi/psiinput.h>
#include <psifiles.h>

using namespace std;

namespace sc {
  
  PsiInput::PsiInput(const string& name) :
    file_() {
    filename_ = string(name);
    indentation_ = 0;
  }
  
  PsiInput::~PsiInput() {
  }
  
  void PsiInput::open() {
    file_.open(filename_.c_str(), ios::out);
    indentation_ = 0;
  }
  
  void PsiInput::close() {
    file_.close();
    indentation_ = 0;
  }
  
  void PsiInput::write_indent() {
    for (int i=0; i<indentation_; i++)
      file_ << " ";
  }
  
  void PsiInput::incindent(int i) {
    if (i > 0)
      indentation_ += i;
  }
  
  void PsiInput::decindent(int i) {
    if (i > 0)
      indentation_ -= i;
  }
  
  void PsiInput::begin_section(const char * s) {
    write_indent();
    file_ << s << ":("<< endl;
    incindent(2);
  }
  
  void PsiInput::end_section(void) {
    decindent(2);
    write_indent();
    file_ << ")"<< endl;
    write_string("\n");
  }
  
  void PsiInput::write_keyword(const char *keyword, const char *value) {
    write_indent();
    file_ << scprintf("%s = %s", keyword, value) << endl;
  }
  
  void PsiInput::write_keyword(const char *keyword, int value) {
    write_indent();
    file_ << scprintf("%s = %d", keyword, value) << endl;
  }
  
  void PsiInput::write_keyword(const char *keyword, bool value) {
    write_indent();
    file_ << scprintf("%s = %s", keyword, value ? "true" : "false") << endl;
  }
  
  void PsiInput::write_keyword(const char *keyword, double value) {
    write_indent();
    file_ << scprintf("%s = %20.15lf", keyword, value) << endl;
  }
  
  void PsiInput::write_keyword_array(const char *keyword, int num, int *values) {
    write_indent();
    file_ << scprintf("%s = (", keyword);
    for (int i=0; i<num; i++) {
      file_ << scprintf(" %d", values[i]);
    }
    file_ << ")"<< endl;
  }
  
  void PsiInput::write_keyword_array(const char *keyword, int num,
                                     double *values) {
    write_indent();
    file_ << scprintf("%s = (", keyword);
    for (int i=0; i<num; i++) {
      file_ << scprintf(" %20.15lf", values[i]);
    }
    file_ << ")"<< endl;
  }
  
  void PsiInput::write_string(const char *s) {
    write_indent();
    file_ << s;
  }
  
  void PsiInput::write_key_wq(const char *keyword, const char *value) {
    write_indent();
    file_ << scprintf("%s = \"%s\"", keyword, value) << endl;
  }
  
  void PsiInput::write_geom(const Ref<Molecule>& m) {
    // copy the molecule and rotate to the symmetry frame
    Ref<Molecule> mol = new Molecule(*m);
    mol->transform_to_symmetry_frame();
    
    // If the highest symmetry group is not the actual group - use subgroup keyword
    if (!mol->point_group()->equiv(mol->highest_point_group())) {
      write_keyword("subgroup", mol->point_group()->symbol());
    }
    
    write_keyword("units", "bohr");
    write_string("geometry = (\n");
    for (int i=0; i < mol->natom(); i++) {
      write_string("  (");
      char *s;
      file_ << mol->atom_symbol(i) <<scprintf(" %14.12lf %14.12lf %14.12lf",
                                              mol->r(i, 0), mol->r(i, 1),
                                              mol->r(i, 2))<< ")"<< endl;
    }
    write_string(")\n");
  }
  
  namespace {
    /// return false if bs has a mixture of pure and cart functions of am >= 1
    bool consistent_puream(const Ref<GaussianBasisSet>& bs) {
      bs->print();
      const int has_pure = bs->has_pure();
      const unsigned int nshell = bs->nshell();
      for (unsigned int s=0; s<nshell; ++s) {
        const GaussianShell& shell = bs->shell(s);
        const unsigned int ncon = shell.ncontraction();
        for (unsigned int c=0; c<ncon; ++c) {
          if (shell.am(c) >= 1&& has_pure != shell.is_pure(c))
            return false;
        }
      }
      return true;
    }
  }
  
  void PsiInput::write_basis(const Ref<GaussianBasisSet>& basis) {
    if (!consistent_puream(basis))
      throw FeatureNotImplemented("PsiInput::write_basis() -- Psi cannot handle yet basis sets that mix cartesian and harmonics functions",__FILE__,__LINE__);
    
    Ref<Molecule> molecule = basis->molecule();
    int natom = molecule->natom();
    
    write_string("basis = (\n");
    incindent(2);
    for (int atom=0; atom<natom; atom++) {
      int uatom = molecule->atom_to_unique(atom);
      
      // Replace all spaces with underscores in order for Psi libipv1 to parse properly
      char *name = strdup(basis->name());
      int len = strlen(name);
      for (int i=0; i<len; i++)
        if (name[i] == ' ')
          name[i] = '_';
      
      char *basisname = new char[strlen(basis->name()) + ((int)ceil(log10((long double)uatom+2))) + 5];
      sprintf(basisname, "\"%s%d\" \n", name, uatom);
      write_string(basisname);
      free(name);
    }
    decindent(2);
    write_string(")\n");
  }
  
  void PsiInput::write_basis_sets(const Ref<GaussianBasisSet>& basis) {
    begin_section("basis");
    Ref<Molecule> molecule = basis->molecule();
    Ref<AtomInfo> atominfo = basis->molecule()->atominfo();
    int nunique = molecule->nunique();
    
    for (int uatom=0; uatom<nunique; uatom++) {
      int atom = molecule->unique(uatom);
      std::string atomname = atominfo->name(molecule->Z(atom));
      
      // Replace all spaces with underscores in order for Psi libipv1 to parse properly
      char *name = strdup(basis->name());
      int len = strlen(name);
      for (int i=0; i<len; i++)
        if (name[i] == ' ')
          name[i] = '_';
      
      char *psibasisname = new char[atomname.size() + strlen(basis->name()) + ((int)ceil(log10((long double)uatom+2))) + 9];
      sprintf(psibasisname, "%s:\"%s%d\" = (\n", atomname.c_str(), name, uatom);
      write_string(psibasisname);
      free(name);
      incindent(2);
      int nshell = basis->nshell_on_center(atom);
      for (int sh=0; sh<nshell; sh++) {
        int shell = basis->shell_on_center(atom, sh);
        GaussianShell& Shell = basis->shell(shell);
        int ncon = Shell.ncontraction();
        int nprim = Shell.nprimitive();
        for (int con=0; con<ncon; con++) {
          char amstring[4];
          sprintf(amstring, "(%c\n", Shell.amchar(con));
          write_string(amstring);
          incindent(2);
          for (int prim=0; prim<nprim; prim++) {
            char primstring[50];
            sprintf(primstring, "(%20.10lf    %20.10lf)\n",
                    Shell.exponent(prim), Shell.coefficient_norm(con, prim));
            write_string(primstring);
          }
          decindent(2);
          write_string(")\n");
        }
      }
      decindent(2);
      write_string(")\n");
      delete[] psibasisname;
    }
    end_section();
  }
  
  void PsiInput::write_defaults(const Ref<PsiExEnv>& exenv, const char *dertype) {
    begin_section("psi");
    
    write_key_wq("label", " ");
    write_keyword("dertype", dertype);
    begin_section("files");
    begin_section("default");
    write_key_wq("name", (exenv->get_fileprefix()).c_str());
    int nscratch = exenv->get_nscratch();
    write_keyword("nvolume", nscratch);
    char *scrname;
    scrname = new char[10];
    for (int i=0; i<nscratch; i++) {
      sprintf(scrname, "volume%d", i+1);
      write_key_wq(scrname, (exenv->get_scratch(i)).c_str());
    }
    delete[] scrname;
    end_section();
    {
      ostringstream oss;
      oss << "file"<< PSIF_CHKPT << ": ( nvolume = 1 volume1 = \"./\" )"<< endl;
      write_string(oss.str().c_str());
    }
    end_section();
    
    end_section();
  }
  
  void PsiInput::print(std::ostream& o) {
    o << indent << "PsiInput:" << std::endl;
    o << indent; for(int i=0; i<7; ++i) { o << "----------"; }  o << endl;
    o << sc::incindent;
    std::ifstream f; f.open(filename_.c_str(), ios::in);
    while(!f.eof()) {
      char buf[256];
      f.getline(buf,256);
      o << indent << buf << endl;
    }
    f.close();
    o << sc::decindent << endl;
    o << indent; for(int i=0; i<7; ++i) { o << "----------"; } o << endl;
  }

}
