//
// vxb_eval_info.cc
//
// Copyright (C) 2003 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#include <stdexcept>
#include <stdlib.h>
#include <util/misc/formio.h>
#include <util/ref/ref.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/moindexspace.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

/*---------------
  R12IntEvalInfo
 ---------------*/
static ClassDesc R12IntEvalInfo_cd(
  typeid(R12IntEvalInfo),"R12IntEvalInfo",4,"virtual public SavableState",
  0, 0, create<R12IntEvalInfo>);

R12IntEvalInfo::R12IntEvalInfo(MBPT2_R12* mbptr12)
{
  // Default values
  memory_ = DEFAULT_SC_MEMORY;
  debug_ = 0;
  dynamic_ = false;
  print_percent_ = 10.0;

  wfn_ = mbptr12;
  ref_ = mbptr12->ref();
  integral_ = mbptr12->integral();
  bs_ = mbptr12->basis();
  bs_aux_ = mbptr12->aux_basis();

  matrixkit_ = SCMatrixKit::default_matrixkit();
  mem_ = MemoryGrp::get_default_memorygrp();
  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();

  integral_->set_basis(bs_);
  Ref<PetiteList> plist = integral_->petite_list();
  RefSCDimension oso_dim = plist->SO_basisdim();
  nocc_ = 0;
  for (int i=0; i<oso_dim->n(); i++) {
    if (ref_->occupation(i) == 2.0) nocc_++;
  }
  nfzc_ = mbptr12->nfzcore();
  nfzv_ = mbptr12->nfzvirt();
  nocc_act_ = nocc_ - nfzc_;
  noso_ = oso_dim.n();

  ints_method_ = mbptr12->r12ints_method();
  ints_file_ = mbptr12->r12ints_file();

  orbsym_ = 0;
  //eigen_(evals_,scf_vec_,occs_,orbsym_);
  eigen2_(evals_,scf_vec_,orbsym_);

  abs_method_ = mbptr12->abs_method();
  construct_ri_basis_(false);

  tfactory_ = new MOIntsTransformFactory(integral_,obs_space_);
  tfactory_->set_memory(memory_);
}

R12IntEvalInfo::R12IntEvalInfo(StateIn& si) : SavableState(si)
{
  wfn_ = require_dynamic_cast<Wavefunction*>(SavableState::restore_state(si),
					     "R12IntEvalInfo::R12IntEvalInfo");
  ref_ << SavableState::restore_state(si);
  integral_ << SavableState::restore_state(si);
  bs_ << SavableState::restore_state(si);
  bs_aux_ << SavableState::restore_state(si);
  bs_ri_ << SavableState::restore_state(si);

  matrixkit_ = SCMatrixKit::default_matrixkit();
  mem_ = MemoryGrp::get_default_memorygrp();
  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();

  si.get(nocc_);
  si.get(nocc_act_);
  si.get(nfzc_);
  si.get(nfzv_);
  si.get(noso_);

  int ints_method; si.get(ints_method); ints_method_ = (StoreMethod) ints_method;
  si.getstring(ints_file_);

  double memory; si.get(memory); memory_ = (size_t) memory;
  si.get(debug_);
  int dynamic; si.get(dynamic); dynamic_ = (bool) dynamic;

  if (si.version(::class_desc<R12IntEvalInfo>()) >= 2) {
    si.get(print_percent_);
  }

  if (si.version(::class_desc<R12IntEvalInfo>()) >= 3) {
    int absmethod; si.get(absmethod); abs_method_ = (LinearR12::ABSMethod) absmethod;
  }

  if (si.version(::class_desc<R12IntEvalInfo>()) >= 4) {
    obs_space_ << SavableState::restore_state(si);
    abs_space_ << SavableState::restore_state(si);
    ribs_space_ << SavableState::restore_state(si);
    act_occ_space_ << SavableState::restore_state(si);
    occ_space_ << SavableState::restore_state(si);
    act_vir_space_ << SavableState::restore_state(si);
    vir_space_ << SavableState::restore_state(si);
    tfactory_ << SavableState::restore_state(si);
  }

  orbsym_ = 0;
  //eigen_(evals_,scf_vec_,occs_,orbsym_);
  eigen2_(evals_,scf_vec_,orbsym_);
}

R12IntEvalInfo::~R12IntEvalInfo()
{
  free(ints_file_);
  delete[] orbsym_;
}

void R12IntEvalInfo::save_data_state(StateOut& so)
{
  SavableState::save_state(wfn_,so);
  SavableState::save_state(ref_.pointer(),so);
  SavableState::save_state(integral_.pointer(),so);
  SavableState::save_state(bs_.pointer(),so);
  SavableState::save_state(bs_aux_.pointer(),so);
  SavableState::save_state(bs_ri_.pointer(),so);

  so.put(nocc_);
  so.put(nocc_act_);
  so.put(nfzc_);
  so.put(nfzv_);
  so.put(noso_);

  so.put((int)ints_method_);
  so.putstring(ints_file_);

  so.put((double)memory_);
  so.put(debug_);
  so.put((int)dynamic_);
  so.put(print_percent_);
  so.put((int)abs_method_);
  
  SavableState::save_state(obs_space_.pointer(),so);
  SavableState::save_state(abs_space_.pointer(),so);
  SavableState::save_state(ribs_space_.pointer(),so);
  SavableState::save_state(act_occ_space_.pointer(),so);
  SavableState::save_state(occ_space_.pointer(),so);
  SavableState::save_state(act_vir_space_.pointer(),so);
  SavableState::save_state(vir_space_.pointer(),so);
  SavableState::save_state(tfactory_.pointer(),so);
}

char* R12IntEvalInfo::ints_file() const
{
  return strdup(ints_file_);
}

void
R12IntEvalInfo::set_memory(const size_t memory)
{
  if (memory > 0)
    memory_ = memory;
  tfactory_->set_memory(memory_);
}


void
R12IntEvalInfo::set_absmethod(LinearR12::ABSMethod abs_method)
{
  if (abs_method != abs_method_) {
    abs_method = abs_method_;
    construct_ri_basis_(false);
  }
}

RefSCMatrix
R12IntEvalInfo::orthog_ri() {
  
  if (orthog_ri_.null())
    construct_orthog_ri_();
  return orthog_ri_;
 };


/////////////////////////////////////////////////////////////////
// Function dquicksort performs a quick sort (smaller -> larger) 
// of the double data in item by the integer indices in index;
// data in item remain unchanged
//
// Both functions borrowed from lib/chemistry/qc/mbpt/mbpt.cc
//
/////////////////////////////////////////////////////////////////
static void
dqs(double *item,int *index,int left,int right)
{
  int i,j;
  double x;
  int y;
  
  i=left; j=right;
  x=item[index[(left+right)/2]];
  
  do {
    while(item[index[i]]<x && i<right) i++;
    while(x<item[index[j]] && j>left) j--;
    
    if (i<=j) {
      if (item[index[i]] != item[index[j]]) {
        y=index[i];
        index[i]=index[j];
        index[j]=y;
      }
      i++; j--;
    }
  } while(i<=j);
  
  if (left<j) dqs(item,index,left,j);
  if (i<right) dqs(item,index,i,right);
}

static void
dquicksort(double *item,int *index,int n)
{
  int i;
  if (n<=0) return;
  for (i=0; i<n; i++) {
    index[i] = i;
  }
  dqs(item,index,0,n-1);
}


// This function is obsoleve and will disappear soon
void R12IntEvalInfo::eigen_(RefDiagSCMatrix &vals, RefSCMatrix &vecs, RefDiagSCMatrix &occs, int*& orbsym)
{
  Ref<Molecule> molecule = bs_->molecule();
  Ref<SCMatrixKit> so_matrixkit = bs_->so_matrixkit();
  Ref<PetiteList> plist = ref_->integral()->petite_list();
  RefSCDimension oso_dim = plist->SO_basisdim();
  // Special dimension for MOs sorted by energy
  int* nfunc_per_block = new int[1];
  nfunc_per_block[0] = oso_dim.n();
  RefSCDimension mo_sorted_dim = new SCDimension(oso_dim.n(), 1, nfunc_per_block, "MO (sorted by energy)");
  mo_sorted_dim->blocks()->set_subdim(0, new SCDimension(nfunc_per_block[0]));

  int nbasis = bs_->nbasis();

  if (debug_) ExEnv::out0() << indent << "R12IntEvalInfo: eigen_" << endl;
  if (debug_) ExEnv::out0() << indent << "getting fock matrix" << endl;
  // get the closed shell AO fock matrices -- this is where SCF is computed
  RefSymmSCMatrix fock_c_so = ref_->fock(0);
  ExEnv::out0() << endl;
  
  // transform the AO fock matrices to the MO basis
  RefSymmSCMatrix fock_c_mo1 = so_matrixkit->symmmatrix(oso_dim);
  RefSCMatrix vecs_so_mo1 = ref_->eigenvectors();
  
  fock_c_mo1.assign(0.0);
  fock_c_mo1.accumulate_transform(vecs_so_mo1.t(), fock_c_so);
  fock_c_so = 0;
  
  if (debug_) ExEnv::out0() << indent << "diagonalizing" << endl;
  // diagonalize the fock matrix
  vals = fock_c_mo1.eigvals();
  
  // compute the AO to new MO scf vector
  if (debug_) ExEnv::out0() << indent << "AO to MO" << endl;
  RefSCMatrix so_ao = plist->sotoao();
  vecs = vecs_so_mo1.t() * so_ao;

  // fill in the occupations
  occs = matrixkit()->diagmatrix(vals.dim());
  for (int i=0; i<oso_dim->n(); i++) {
    occs(i) = ref_->occupation(i);
  }
  // allocate storage for symmetry information
  if (!orbsym) orbsym = new int[nbasis];
  // Check for degenerate eigenvalues.  Use unsorted eigenvalues since it
  // only matters if the degeneracies occur within a given irrep.
  BlockedDiagSCMatrix *bvals = dynamic_cast<BlockedDiagSCMatrix*>(vals.pointer());
  for (int i=0; i<bvals->nblocks(); i++) {
    int done = 0;
    RefDiagSCMatrix valsi = bvals->block(i);
    for (int j=1; j<valsi.n(); j++) {
      if (fabs(valsi(j)-valsi(j-1)) < 1.0e-7) {
	ExEnv::out0() << indent
		      << "NOTE: There are degenerate orbitals within an irrep."
		      << "  This will make"
		      << endl
		      << indent
		      << "      some diagnostics, such as the largest amplitude,"
		      << " nonunique."
		      << endl;
	done = 1;
	break;
      }
      if (done) break;
    }
  }
  // sort the eigenvectors and values if symmetry is not c1
  if (molecule->point_group()->char_table().order() != 1) {
    if (debug_) ExEnv::out0() << indent << "sorting eigenvectors" << endl;
    double *evals = new double[noso_];
    vals->convert(evals);
    int *indices = new int[noso_];
    dquicksort(evals,indices,noso_);
    delete[] evals;
    // make sure all nodes see the same indices and evals
    msg_->bcast(indices,noso_);
    RefSCMatrix newvecs(mo_sorted_dim, vecs.coldim(), matrixkit());
    RefDiagSCMatrix newvals(mo_sorted_dim, matrixkit());
    RefDiagSCMatrix newoccs(mo_sorted_dim, matrixkit());
    for (int i=0; i<noso_; i++) {
      newvals(i) = vals(indices[i]);
      newoccs(i) = occs(indices[i]);
      for (int j=0; j<nbasis; j++) {
	newvecs(i,j) = vecs(indices[i],j);
      }
    }
    occs = newoccs;
    vecs = newvecs;
    vals = newvals;
    
    // compute orbital symmetry information
    CharacterTable ct = molecule->point_group()->char_table();
    int orbnum = 0;
    int *tmp_irrep = new int[noso_];
    int *tmp_num = new int[noso_];
    for (int i=0; i<oso_dim->blocks()->nblock(); i++) {
      for (int j=0; j<oso_dim->blocks()->size(i); j++, orbnum++) {
	tmp_irrep[orbnum] = i;
	tmp_num[orbnum] = j;
      }
    }
    for (int i=0; i<noso_; i++) {
      orbsym[i] = tmp_irrep[indices[i]];
    }
    delete[] tmp_irrep;
    delete[] tmp_num;
    
    delete[] indices;
  }
  else {
    // compute orbital symmetry information for c1
    for (int i=0; i<noso_; i++) {
      orbsym[i] = 0;
    }
  }
  // check the splitting between frozen and nonfrozen orbitals
  if (nfzc_ && nfzc_ < noso_) {
    double split = vals(nfzc_) - vals(nfzc_-1);
    if (split < 0.2) {
      ExEnv::out0() << endl
		    << indent << "WARNING: "
		    << "R12IntEvalInfo: gap between frozen and active occupied orbitals is "
		    << split << " au" << endl << endl;
    }
  }
  if (nfzv_ && noso_-nfzv_-1 >= 0) {
    double split = vals(nbasis-nfzv_) - vals(nbasis-nfzv_-1);
    if (split < 0.2) {
      ExEnv::out0() << endl
		    << indent << "WARNING: "
		    << "R12IntEvalInfo: gap between frozen and active virtual orbitals is "
		    << split << " au" << endl << endl;
    }
  }
  if (debug_) ExEnv::out0() << indent << "R12IntEvalInfo: eigen_ done" << endl;
}


void R12IntEvalInfo::eigen2_(RefDiagSCMatrix &vals, RefSCMatrix &vecs, int*& orbsym)
{
  Ref<Molecule> molecule = bs_->molecule();
  Ref<SCMatrixKit> so_matrixkit = bs_->so_matrixkit();
  Ref<PetiteList> plist = ref_->integral()->petite_list();
  RefSCDimension oso_dim = plist->SO_basisdim();

  if (debug_) ExEnv::out0() << indent << "R12IntEvalInfo: eigen_" << endl;
  if (debug_) ExEnv::out0() << indent << "getting fock matrix" << endl;
  // get the closed shell AO fock matrices -- this is where SCF is computed
  RefSymmSCMatrix fock_c_so = ref_->fock(0);
  ExEnv::out0() << endl;
  
  // transform the AO fock matrices to the MO basis
  RefSymmSCMatrix fock_c_mo1 = so_matrixkit->symmmatrix(oso_dim);
  RefSCMatrix vecs_so_mo1 = ref_->eigenvectors();
  
  fock_c_mo1.assign(0.0);
  fock_c_mo1.accumulate_transform(vecs_so_mo1.t(), fock_c_so);
  fock_c_so = 0;
  
  if (debug_) ExEnv::out0() << indent << "diagonalizing" << endl;
  // diagonalize the fock matrix
  vals = fock_c_mo1.eigvals();
  
  // compute the AO to new MO scf vector
  if (debug_) ExEnv::out0() << indent << "AO to MO" << endl;
  RefSCMatrix so_ao = plist->sotoao();
  vecs = so_ao.t() * vecs_so_mo1;

  obs_space_ = new MOIndexSpace("MOs sorted by energy", vecs, bs_, vals, 0, 0);
  occ_space_ = new MOIndexSpace("occupied MOs sorted by energy", vecs, bs_, vals, 0, noso_ - nocc_);
  if (nfzc_ == 0)
    act_occ_space_ = occ_space_;
  else
    act_occ_space_ = new MOIndexSpace("active occupied MOs sorted by energy", vecs, bs_, vals, nfzc_, noso_ - nocc_);
  vir_space_ = new MOIndexSpace("unoccupied MOs sorted by energy", vecs, bs_, vals, nocc_, 0);
  if (nfzv_ == 0)
    act_vir_space_ = vir_space_;
  else
    act_vir_space_ = new MOIndexSpace("active unoccupied MOs sorted by energy", vecs, bs_, vals, nocc_, nfzv_);

  vecs = obs_space_->coefs().t();
  vals = obs_space_->evals();
  vector<int> mosym = obs_space_->mosym();
  int nmo = mosym.size();
  if (nmo != noso_)
    throw std::runtime_error("R12IntEvalInfo::eigen2_() -- logic error in MOIndexSpace");
  if (orbsym == 0)
    orbsym = new int[nmo];
  for(int i=0; i<nmo; i++)
    orbsym[i] = mosym[i];

  if (debug_) ExEnv::out0() << indent << "R12IntEvalInfo: eigen2_ done" << endl;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
