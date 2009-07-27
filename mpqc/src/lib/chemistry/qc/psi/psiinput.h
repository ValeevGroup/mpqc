//
// psiinput.h
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
#pragma interface
#endif

#ifndef _chemistry_qc_psi_input_h
#define _chemistry_qc_psi_input_h

using namespace std;

#include <fstream>
#include <string>
#include<util/ref/ref.h>
#include<chemistry/molecule/molecule.h>
#include<chemistry/qc/basis/basis.h>

namespace sc {

class PsiExEnv;
class CorrelationTable;

///////////////////////////////////////////////////
/// PsiInput is a Psi input file

class PsiInput: public RefCount {

  string filename_;
  std::ofstream file_;
  int me_;                // task id

  int indentation_;
  
  // No default constructor
  PsiInput() {};

  // can run on me_?
  bool can_run_on_me() { return me_ == 0; }
  
  public:
    PsiInput(const string& name);
    ~PsiInput();
    void open();
    void close();
    void print(std::ostream&o=ExEnv::out0());

    void begin_section(const char * s);
    void end_section();
    void write_indent();
    void incindent(int);
    void decindent(int);
    void write_comment(const char *);
    void write_keyword(const char *, const char *);
    void write_keyword(const char *, bool);
    void write_keyword(const char *, int);
    void write_keyword(const char *, double);
    template <typename T> void write_keyword_array(const char *, const std::vector<T>&);
    void write_keyword_array(const char *, int, int *);
    void write_keyword_array(const char *, int, double *);
    void write_string(const char *);
    void write_key_wq(const char *, const char *);

    /// Construct the "basis" keyword for input. All functions with angular momentum >= 1 must be Cartesian or all must be sph. harm.
    void write_basis(const Ref<GaussianBasisSet>&);
    /// Write basis sets explicitly
    void write_basis_sets(const Ref<GaussianBasisSet>&);
    void write_geom(const Ref<Molecule>&);
    
    void write_defaults(const Ref<PsiExEnv>&, const char *dertype);
};

}

#endif
