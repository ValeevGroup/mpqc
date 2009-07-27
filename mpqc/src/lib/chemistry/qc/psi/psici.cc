//
// psici.cc
//
// Copyright (C) 2008 Martin Torheyden
//
// Author: Martin Torheyden <mtorhey@vt.edu>
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

#include <assert.h>
#include <psifiles.h>
#include <ccfiles.h>

#include <chemistry/qc/mbptr12/pairiter.impl.h>
#include <chemistry/qc/mbptr12/print.h>
#include <chemistry/qc/psi/psiwfn.h>

#include <chemistry/qc/psi/psici.h>

using namespace std;
//using namespace sc;


namespace sc {

  static ClassDesc PsiCI_cd(typeid(PsiCI), "PsiCI", 1,
                            "public PsiCorrWavefunction", 0, create<PsiCI>,
                            create<PsiCI>);

  PsiCI::PsiCI(const Ref<KeyVal> &keyval) :
    PsiCorrWavefunction(keyval) {

    if(keyval->exists("opdm")) {
      eval_opdm_ = keyval->booleanvalue("opdm");
    }
    else {
      eval_opdm_ = false;
    }

    if(keyval->exists("opdm_print")) {
      opdm_print_ = keyval->booleanvalue("opdm_print");
    }
    else {
      opdm_print_ = false;
    }

    if(keyval->exists("tpdm")) {
      eval_tpdm_ = keyval->booleanvalue("tpdm");
    }
    else {
      eval_tpdm_ = false;
    }

    if(keyval->exists("tpdm_print")) {
      tpdm_print_ = keyval->booleanvalue("tpdm_print");
    }
    else {
      tpdm_print_ = false;
    }

    if(keyval->exists("root")) {
      root_ = keyval->intvalue("root");
    }
    else {
      root_ = 1;   /// corresponds to ground state.
    }

    if(keyval->exists("num_roots")) {
      num_roots_ = keyval->intvalue("num_roots");
      if(num_roots_<root_) {
        throw InputError("PsiCI_PT2R12::PsiCI_PT2R12(const Ref<KeyVal> &) -- num_roots should not be smaller than root.",__FILE__,__LINE__);
      }
    }
    else {
      num_roots_ = root_;
    }

    if(keyval->exists("ex_lvl")) {
      ex_lvl_ = keyval->intvalue("ex_lvl");
    }
    else {
      ex_lvl_ = 2;
    }

    if(keyval->exists("repl_otf")) {
      repl_otf_ = keyval->booleanvalue("repl_otf");
    }
    else {
      repl_otf_ = false;
    }

    if(keyval->exists("maxiter")) {
      maxiter_ = keyval->intvalue("maxiter");
    }
    else {
      maxiter_ = 12;
    }
  }

  PsiCI::PsiCI(StateIn& s) :
    PsiCorrWavefunction(s) {
    int eval_opdm; s.get(eval_opdm); eval_opdm_ = (bool)eval_opdm;
    int opdm_print; s.get(opdm_print); opdm_print_ = (bool)opdm_print;
    int eval_tpdm; s.get(eval_tpdm); eval_tpdm_ = (bool)eval_tpdm;
    int tpdm_print; s.get(tpdm_print); tpdm_print_ = (bool)tpdm_print;
    s.get(root_);
    s.get(num_roots_);
    s.get(ex_lvl_);
    int repl_otf; s.get(repl_otf); repl_otf_ = (bool)repl_otf;
    s.get(maxiter_);
  }

  PsiCI::~PsiCI(){}

  void PsiCI::save_data_state(StateOut &s){
    PsiCorrWavefunction::save_data_state(s);
    s.put((int)eval_opdm_);
    s.put((int)opdm_print_);
    s.put((int)eval_tpdm_);
    s.put((int)tpdm_print_);
    s.put(root_);
    s.put(num_roots_);
    s.put(ex_lvl_);
    s.put((int)repl_otf_);
    s.put(maxiter_);
  }

  void PsiCI::write_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->open();
    PsiCorrWavefunction::write_input(convergence);
    input->write_keyword("psi:wfn","detci");
    if(eval_opdm_==true){
      input->write_keyword("detci:opdm","true");
    }
    if(opdm_print_==true) {
      input->write_keyword("detci:opdm_print","true");
    }

    if(eval_tpdm_==true){
      input->write_keyword("detci:tpdm","true");
    }
    if(tpdm_print_==true) {
      input->write_keyword("detci:tpdm_print","true");
    }

    input->write_keyword("detci:root",root_);

    input->write_keyword("detci:num_roots",num_roots_);

    input->write_keyword("detci:ex_lvl",ex_lvl_);

    if(repl_otf_==true) {
      input->write_keyword("detci:repl_otf","true");
    }

    input->write_keyword("detci:maxiter",maxiter_);

    input->close();
  }

  void PsiCI::compute(){
    RefSymmSCMatrix opdm;
    RefSymmSCMatrix tpdm;
    PsiWavefunction::compute();
    if(eval_opdm_==true){
      opdm=onepdm();

      FILE *testf = fopen("opdm_test","w");
      print_onepdm_mat(testf,opdm,1.0e-6);
      fclose(testf);
    }
    if(eval_tpdm_==true){
      tpdm=twopdm();

      FILE *testf = fopen("tpdm_test","w");
      print_twopdm_mat(testf,tpdm,1.0e-6);
      fclose(testf);
    }
  }

}
