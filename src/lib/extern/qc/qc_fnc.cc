//
// qc_fnc.cc
//
// Copyright (C) 2012 Edward Valeev
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

#include <cassert>
#include <iostream>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/split.h>
#include <extern/qc/qc_fnc.h>

using namespace sc;

namespace {
  Ref<Molecule> qc_mol;

  Ref<GaussianBasisSet> qc_obs;
  struct gbs_data {
      gbs_data() {}
      gbs_data(const Ref<GaussianBasisSet>& gbs) {
        const int nshell = gbs->nshell();
        shell_to_center.resize(nshell);
        for(int s=0; s<nshell; ++s) {
          shell_to_center[s] = gbs->shell_to_center(s);
        }
        shell_to_nfunction.resize(nshell);
        for(int s=0; s<nshell; ++s) {
          shell_to_nfunction[s] = gbs->shell(s).nfunction();
        }
        shell_to_function.resize(nshell);
        for(int s=0; s<nshell; ++s) {
          shell_to_function[s] = gbs->shell_to_function(s);
        }
      }
      std::vector<int> shell_to_center;
      std::vector<int> shell_to_nfunction;
      std::vector<int> shell_to_function;
  };
  gbs_data qc_obs_data;

  Ref<Integral> qc_integral;
  Ref<OneBodyInt> qc_overlap;
  const double* qc_overlap_buf = 0;
  double* qc_overlap_extbuf = 0;
  Ref<OneBodyInt> qc_hcore;
  const double* qc_hcore_buf = 0;
  double* qc_hcore_extbuf = 0;
  Ref<TwoBodyInt> qc_twoecoulomb;
  const double* qc_twoecoulomb_buf = 0;
  double* qc_twoecoulomb_extbuf = 0;

  void die(const char* msg) {
    std::cerr << "error: " << msg << std::endl;
    abort();
  }
}

void init_molecule_(int natoms, const double* Z, const double* xyz, int use_symmetry) {
  qc_mol = new Molecule;
  for(int a=0; a<natoms; ++a) {
    qc_mol->add_atom(static_cast<int>(Z[a]), xyz[a*3], xyz[a*3+1], xyz[a*3+2]);
  }
  if (use_symmetry == 0) {
    Ref<PointGroup> c1_ptgrp = new PointGroup("C1");
    qc_mol->set_point_group(c1_ptgrp);
  }
  else {
    qc_mol->symmetrize(qc_mol->highest_point_group(1e-4));
  }

  ExEnv::out0() << std::endl;
  ExEnv::out0() << indent << "constructed Molecule object:" << std::endl;
  qc_mol->print(ExEnv::out0());
  ExEnv::out0() << std::endl;
}

void init_molecule_xyz_(const char* fname, int use_symmetry, int fname_nchar) {
  MPQC_ASSERT(false);
}

void init_basis_set_(const char* basis_name, int basis_name_nchar) {
  if (qc_mol == 0) {
    die("called init_basis_set before calling init_molecule!");
  }

  char* basis_name_cstr = new char[basis_name_nchar+1];
  strncpy(basis_name_cstr, basis_name, basis_name_nchar);
  basis_name_cstr[basis_name_nchar] = '\0';

  Ref<AssignedKeyVal> tmpkv = new AssignedKeyVal;
  tmpkv->assign("molecule", qc_mol.pointer());
  tmpkv->assign("name", basis_name_cstr);
  Ref<KeyVal> kv = tmpkv;
  qc_obs = new GaussianBasisSet(kv);

  // get rid of general constractions for simplicity
  if (qc_obs->max_ncontraction() > 1) {
    Ref<GaussianBasisSet> split_basis = new SplitBasisSet(qc_obs, qc_obs->name());
    qc_obs = split_basis;
  }

  qc_obs_data = gbs_data(qc_obs);

  ExEnv::out0() << std::endl << indent << "constructed GaussianBasisSet object:" << std::endl;
  qc_obs->print(ExEnv::out0());
  ExEnv::out0() << std::endl;
}

void init_basis_set_g94_(const char* fname, int fname_nchar) {
  MPQC_ASSERT(false);
  qc_obs_data = gbs_data(qc_obs);
}

int basis_set_nbasis_() {
  return qc_obs->nbasis();
}

int basis_set_nshell_() {
  return qc_obs->nshell();
}

int basis_set_max_nfunction_in_shell_() {
  return qc_obs->max_nfunction_in_shell();
}

int* basis_set_shell_to_nfunction_() {
  return &(qc_obs_data.shell_to_nfunction[0]);
}

void basis_set_shell_to_nfunction_subrt_(int* s2nf) {
  std::copy(qc_obs_data.shell_to_nfunction.begin(), qc_obs_data.shell_to_nfunction.end(), s2nf);
}

int* basis_set_shell_to_function_() {
  return &qc_obs_data.shell_to_function[0];
}

void basis_set_shell_to_function_subrt_(int* s2f) {
  std::copy(qc_obs_data.shell_to_function.begin(), qc_obs_data.shell_to_function.end(), s2f);
}

int* basis_set_shell_to_center_() {
  return &qc_obs_data.shell_to_center[0];
}

void basis_set_shell_to_center_subrt_(int* s2c) {
  std::copy(qc_obs_data.shell_to_center.begin(), qc_obs_data.shell_to_center.end(), s2c);
}

void init_integrals_() {
  if (qc_obs == 0) {
    die("called init_integrals before calling init_basis_set!");
  }
  int izero = 0;
  char** dptrzero = 0;
  Ref<Integral> integral = Integral::initial_integral(izero,dptrzero);
  if (integral) Integral::set_default_integral(integral);
  qc_integral = Integral::get_default_integral()->clone();
  qc_integral->set_basis(qc_obs);
}

void done_integrals_() {
  qc_integral = 0;
  done_overlap_integrals_();
  done_hcore_integrals_();
  done_twoecoulomb_integrals_();
}

void init_overlap_integrals_() {
  if (qc_integral == 0) {
    die("called init_overlap_integrals before calling init_integrals!");
  }
  qc_integral->set_basis(qc_obs, qc_obs);
  qc_overlap = qc_integral->overlap();
  qc_overlap_buf = qc_overlap->buffer();
}

const double* overlap_integrals_buffer_() {
  if (qc_overlap == 0) {
    die("called overlap_integrals_buffer before calling init_overlap_integrals!");
  }
  return qc_overlap_buf;
}

void set_overlap_integrals_buffer_(double* buf) {
  qc_overlap_extbuf = buf;
}

void compute_overlap_shell_(int bra, int ket) {
  qc_overlap->compute_shell(bra, ket);
  if (qc_overlap_extbuf) {
    const int nints = qc_obs->shell(bra).nfunction() * qc_obs->shell(ket).nfunction();
    std::copy(qc_overlap_buf, qc_overlap_buf+nints, qc_overlap_extbuf);
  }
}

void done_overlap_integrals_() {
  qc_overlap = 0;
  qc_overlap_buf = 0;
  qc_overlap_extbuf = 0;
}

void init_hcore_integrals_() {
  if (qc_integral == 0) {
    die("called init_hcore_integrals before calling init_integrals!");
  }
  qc_integral->set_basis(qc_obs, qc_obs);
  qc_hcore = qc_integral->hcore();
  qc_hcore_buf = qc_hcore->buffer();
}

const double* hcore_integrals_buffer_() {
  if (qc_hcore == 0) {
    die("called hcore_integrals_buffer before calling init_hcore_integrals!");
  }
  return qc_hcore_buf;
}

void set_hcore_integrals_buffer_(double* buf) {
  qc_hcore_extbuf = buf;
}

void compute_hcore_shell_(int bra, int ket) {
  qc_hcore->compute_shell(bra, ket);
  if (qc_hcore_extbuf) {
    const int nints = qc_obs->shell(bra).nfunction() * qc_obs->shell(ket).nfunction();
    std::copy(qc_hcore_buf, qc_hcore_buf+nints, qc_hcore_extbuf);
  }
}

void done_hcore_integrals_() {
  qc_hcore = 0;
  qc_hcore_buf = 0;
  qc_hcore_extbuf = 0;
}

void init_twoecoulomb_integrals_() {
  if (qc_integral == 0) {
    die("called init_twoecoulomb_integrals before calling init_integrals!");
  }
  qc_integral->set_basis(qc_obs, qc_obs, qc_obs, qc_obs);
  qc_twoecoulomb = qc_integral->electron_repulsion();
  qc_twoecoulomb_buf = qc_twoecoulomb->buffer();
}

const double* twoecoulomb_integrals_buffer_() {
  if (qc_twoecoulomb == 0) {
    die("called twoecoulomb_integrals_buffer before calling init_twoecoulomb_integrals!");
  }
  return qc_twoecoulomb_buf;
}

void set_twoecoulomb_integrals_buffer_(double *buf) {
  qc_twoecoulomb_extbuf = buf;
}

void compute_twoecoulomb_shell_(int bra1, int ket1, int bra2, int ket2) {
  qc_twoecoulomb->compute_shell(bra1, ket1, bra2, ket2);
  if (qc_twoecoulomb_extbuf) {
    const int nints = qc_obs->shell(bra1).nfunction() * qc_obs->shell(ket1).nfunction() *
                      qc_obs->shell(bra2).nfunction() * qc_obs->shell(ket2).nfunction();
    std::copy(qc_twoecoulomb_buf, qc_twoecoulomb_buf+nints, qc_twoecoulomb_extbuf);
  }
}

void done_twoecoulomb_integrals_() {
  qc_twoecoulomb = 0;
  qc_twoecoulomb_buf = 0;
  qc_twoecoulomb_extbuf = 0;
}




/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
