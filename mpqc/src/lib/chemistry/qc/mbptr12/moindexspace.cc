//
// moindexspace.cc
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
#include <algorithm>
#include <stdlib.h>
#include <util/misc/formio.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/mbptr12/moindexspace.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

/*---------------
  MOIndexSpace
 ---------------*/
static ClassDesc MOIndexSpace_cd(
  typeid(MOIndexSpace),"MOIndexSpace",1,"virtual public SavableState",
  0, 0, create<MOIndexSpace>);

MOIndexSpace::MOIndexSpace(std::string name, const RefSCMatrix& full_coefs, const Ref<GaussianBasisSet> basis,
                           const vector<int>& offsets, const vector<int>& nmopi, IndexOrder moorder,
                           const RefDiagSCMatrix& evals) :
  name_(name), basis_(basis), offsets_(offsets), nmo_(nmopi), full_rank_(full_coefs.coldim().n()),
  moorder_(moorder)
{
  full_coefs_to_coefs(full_coefs, evals);
  
  init();
}

MOIndexSpace::MOIndexSpace(std::string name, const RefSCMatrix& full_coefs, const Ref<GaussianBasisSet> basis,
                           const RefDiagSCMatrix& evals, int nfzc, int nfzv, IndexOrder moorder) :
  name_(name), basis_(basis), full_rank_(full_coefs.coldim().n()), moorder_(moorder)
{
  if (evals.null())
    throw std::runtime_error("MOIndexSpace::MOIndexSpace() -- null eigenvalues matrix");
  if (nfzc < 0 || nfzc > full_coefs.coldim().n())
    throw std::runtime_error("MOIndexSpace::MOIndexSpace() -- invalid nfzc");
  if (nfzc + nfzv > full_coefs.coldim().n())
    throw std::runtime_error("MOIndexSpace::MOIndexSpace() -- invalid nfzc+nfzv");
  frozen_to_blockinfo(nfzc,nfzv,evals);
  full_coefs_to_coefs(full_coefs, evals);
  
  init();
}

MOIndexSpace::MOIndexSpace(std::string name, const RefSCMatrix& full_coefs, const Ref<GaussianBasisSet> basis) :
  name_(name), basis_(basis), full_rank_(full_coefs.coldim().n()), moorder_(symmetry)
{
  Ref<SCBlockInfo> modim_blocks = full_coefs.coldim()->blocks();
  int nb = modim_blocks->nblock();
  offsets_.resize(nb);
  nmo_.resize(nb);
  for(int i=0; i<nb; i++) {
    offsets_[i] = 0;
    nmo_[i] = modim_blocks->size(i);
  }
  
  full_coefs_to_coefs(full_coefs, 0);
  
  init();
}

/*MOIndexSpace::MOIndexSpace(std::string name, const RefSCMatrix& coefs, const Ref<GaussianBasisSet> basis,
                           IndexOrder moorder, const vector<int>& mosym, const RefDiagSCMatrix& evals) :
  name_(name), coefs_(coefs), evals_(evals), basis_(basis), mosym_(mosym), full_rank_(coefs_.coldim().n()),
  rank_(full_rank_), moorder_(moorder)
{
  check_mosym();
  
  init();
}*/

MOIndexSpace::MOIndexSpace(std::string name, const Ref<MOIndexSpace>& orig_space, const RefSCMatrix& new_coefs,
                           const Ref<GaussianBasisSet>& new_basis) :
  name_(name), mosym_(orig_space->mosym_), evals_(orig_space->evals_),
  rank_(orig_space->rank_), full_rank_(orig_space->full_rank_), nblocks_(orig_space->nblocks_),
  offsets_(orig_space->offsets_), nmo_(orig_space->nmo_), map_to_full_space_(orig_space->map_to_full_space_),
  moorder_(orig_space->moorder_)
{
  if (rank_ != new_coefs.coldim()->n())
    throw std::runtime_error("MOIndexSpace::MOIndexSpace() -- new_coefs have different number of orbitals");
  coefs_ = new_coefs;
  basis_ = new_basis;
  init();
}

MOIndexSpace::MOIndexSpace(StateIn& si) : SavableState(si)
{
  si.get(name_);

  coefs_.restore(si);
  evals_.restore(si);
  basis_ << SavableState::restore_state(si);
  si.get(mosym_);

  si.get(rank_);
  si.get(full_rank_);
  si.get(nblocks_);
  si.get(offsets_);
  si.get(nmo_);
  si.get(map_to_full_space_);

  int moorder; si.get(moorder); moorder = (int) moorder_;
  
  init();
}

MOIndexSpace::~MOIndexSpace()
{
}

void
MOIndexSpace::save_data_state(StateOut& so)
{
  so.put(name_);

  coefs_.save(so);
  evals_.save(so);
  modim_ = evals_.dim();
  SavableState::save_state(basis_.pointer(),so);
  so.put(mosym_);

  so.put(rank_);
  so.put(full_rank_);
  so.put(nblocks_);
  so.put(offsets_);
  so.put(nmo_);
  so.put(map_to_full_space_);

  so.put((int)moorder_);
}

const std::string&
MOIndexSpace::name() const { return name_; }

const Ref<GaussianBasisSet>&
MOIndexSpace::basis() const { return basis_; }

const RefSCMatrix&
MOIndexSpace::coefs() const { return coefs_; }

const RefDiagSCMatrix&
MOIndexSpace::evals() const { return evals_; }

const vector<int>&
MOIndexSpace::mosym() const { return mosym_; }

MOIndexSpace::IndexOrder
MOIndexSpace::moorder() const { return moorder_; }

int
MOIndexSpace::rank() const { return rank_; }

int
MOIndexSpace::full_rank() const { return full_rank_; }

int
MOIndexSpace::nblocks() const { return nblocks_; }

const vector<int>&
MOIndexSpace::nmo() const { return nmo_; }

const vector<int>&
MOIndexSpace::offsets() const { return offsets_; }

int
MOIndexSpace::to_full_space(const int i) const { return map_to_full_space_.at(i); }

void
MOIndexSpace::check_mosym() const
{
  int ng = basis_->molecule()->point_group()->char_table().order();
  
  for(vector<int>::const_iterator p=mosym_.begin(); p != mosym_.end(); ++p) {
    if (*p < 0 || *p >= ng)
      throw std::runtime_error("MOIndexSpace::check_mosym() -- invalid value in the list of orbital irreps");
  }
}


void
MOIndexSpace::frozen_to_blockinfo(const int nfzc, const int nfzv, const RefDiagSCMatrix& evals)
{
  int rank = evals.dim().n();

  int nb = evals.dim()->blocks()->nblock();
  nmo_.resize(nb);
  offsets_.resize(nb);
  for(int b=0; b<nb; b++) {
    nmo_[b] = evals.dim()->blocks()->size(b);
    offsets_[b] = 0;
  }
  
  // Get the energies of the orbitals in this space
  double* energy = new double[rank];
  int* index_map = new int[rank];
  vector<int> blocked_index_to_irrep(rank);
  int ii = 0;     // blocked index to this space
  int offset = 0;
  for(int b=0; b<nb; b++) {
    for(int i=0; i<nmo_[b]; i++, ii++) {
      energy[ii] = evals.get_element(i+offset);
      blocked_index_to_irrep[ii] = b;
    }
    offset += nmo_[b];
  }
    
  // Do the sort
  dquicksort(energy,index_map,rank);
  delete[] energy;
  
  // Get rid of nfzc lowest orbitals
  for(int i=0; i<nfzc; i++) {
    int b = blocked_index_to_irrep[index_map[i]];
    ++offsets_[b];
    --nmo_[b];
  }

  // Get rid of nfzv highest orbitals
  for(int i=rank-nfzv; i<rank; i++) {
    int b = blocked_index_to_irrep[index_map[i]];
    --nmo_[b];
  }
  
  delete[] index_map;
}

void
MOIndexSpace::full_coefs_to_coefs(const RefSCMatrix& full_coefs, const RefDiagSCMatrix& evals)
{
  // compute the rank of this
  rank_ = 0;
  for(vector<int>::const_iterator p=nmo_.begin(); p != nmo_.end(); ++p) {
    rank_ += *p;
  }

  mosym_.resize(rank_);
  RefSCDimension modim = full_coefs.coldim();  // the dimension of the full space
  
  // In general vectors are ordered differently from the original
  int* index_map = new int[rank_];                      // maps index in this (sorted) space to this (blocked) space
  vector<int> blocked_subindex_to_full_index(rank_);    // maps index from this space(in blocked form) into the full space
  vector<int> blocked_subindex_to_irrep(rank_);         // maps index from this space(in blocked form) to the irrep
  if (moorder_ == symmetry) {
    // coefs_ has the same number of blocks as full_coefs_
    int nb = modim->blocks()->nblock();
    int* nfunc_per_block = new int[nb];
    for(int i=0; i<nb; i++)
      nfunc_per_block[i] = nmo_[i];
    modim_ = new SCDimension(rank_, nb, nfunc_per_block, ("MO(" + name_ + ")").c_str());
    for(int i=0; i<nb; i++)
      modim_->blocks()->set_subdim(i, new SCDimension(nfunc_per_block[i]));
    delete[] nfunc_per_block;
  
    // The sorted->blocked reordering array is trivial when no resorting is done
    for(int i=0; i<rank_; i++) {
      index_map[i] = i;
    }
    
    int ii = 0;     // blocked index to this space
    int offset = 0;
    for(int b=0; b<nb; b++) {
      for(int i=0; i<nmo_[b]; i++, ii++) {
        blocked_subindex_to_full_index[ii] = i+offsets_[b]+offset;
        blocked_subindex_to_irrep[ii] = b;
      }
      offset += modim->blocks()->size(b);
    }    
  }
  else if (moorder_ == energy) {
    //
    // Sort vectors by their energy
    //
    
    // Get the energies of the orbitals in this space
    double* energy = new double[rank_];
    int nb = nmo_.size();
    int ii = 0;     // blocked index to this space
    int offset = 0;
    for(int b=0; b<nb; b++) {
      for(int i=0; i<nmo_[b]; i++, ii++) {
        energy[ii] = evals.get_element(i+offsets_[b]+offset);
        blocked_subindex_to_full_index[ii] = i+offsets_[b]+offset;
        blocked_subindex_to_irrep[ii] = b;
      }
      offset += modim->blocks()->size(b);
    }
    
    // Do the sort
    dquicksort(energy,index_map,rank_);
    
    // coefs_ has 1 block
    int* nfunc_per_block = new int[1];
    nfunc_per_block[0] = rank_;
    modim_ = new SCDimension(rank_, 1, nfunc_per_block, ("MO(" + name_ + ")").c_str());
    if (rank_)
      modim_->blocks()->set_subdim(0, new SCDimension(nfunc_per_block[0]));
    
    // Recompute offsets_ and nmo_ to conform the energy ordering
    offset = 0;
    for(vector<int>::const_iterator p=offsets_.begin(); p != offsets_.end(); ++p) {
      offset += *p;
    }
    offsets_.resize(1);
    offsets_[0] = offset;
    nmo_.resize(1);
    nmo_[0] = rank_;
    
    delete[] energy;
    delete[] nfunc_per_block;
  }
  else
    throw std::runtime_error("MOIndexSpace::full_coefs_to_coefs() -- moorder should be either energy or symmetry");
  
  // Copy required columns of full_coefs_ into coefs_
  RefSCDimension aodim = full_coefs.rowdim();
  Ref<SCMatrixKit> so_matrixkit = basis_->so_matrixkit();
  coefs_ = so_matrixkit->matrix(aodim, modim_);
  evals_ = so_matrixkit->diagmatrix(modim_);
  for (int i=0; i<rank_; i++) {
    int ii = blocked_subindex_to_full_index[index_map[i]];
    mosym_[i] = blocked_subindex_to_irrep[index_map[i]];
    for (int j=0; j<aodim.n(); j++) {
      coefs_(j,i) = full_coefs(j,ii);
    }
  }
  if (evals.nonnull())
    for (int i=0; i<rank_; i++) {
      int ii = blocked_subindex_to_full_index[index_map[i]];
      evals_(i) = evals(ii);
    }
  else
    evals_.assign(0.0);

  nblocks_ = modim_->blocks()->nblock();

  // Compute the map to the full space
  map_to_full_space_.resize(rank_);
  for (int i=0; i<rank_; i++) {
    map_to_full_space_[i] = blocked_subindex_to_full_index[index_map[i]];
  }

  delete[] index_map;  
}

void
MOIndexSpace::init()
{
}


size_t
MOIndexSpace::memory_in_use() const
{
  size_t memory = (size_t)basis_->nbasis() * rank_ * sizeof(double);
  return memory;
}

void
MOIndexSpace::print(ostream&o) const
{
  o << indent << "MOIndexSpace \"" << name_ << "\":" << endl;
  o << incindent;
  o << indent << "Basis Set:" << endl;
  o << incindent; basis_->print(o); o << decindent << endl;
  o << decindent;
}

void
MOIndexSpace::print_summary(ostream& o) const
{
  o << indent << "MOIndexSpace \"" << name_ << "\":" << endl;
  o << incindent;
  o << indent << "GaussianBasisSet \"" << basis_->name() << "\""<< endl;
  o << indent << "  rank  nbasis  nshell  nfuncmax" << endl;
  o << indent << scprintf("  %-6i %-6i  %-6i   %-6i",
                          rank_,
                          basis_->nbasis(),
                          basis_->nshell(),
                          basis_->max_nfunction_in_shell()) << endl;
  o << decindent;

}

/////////////////////////////////////////////////////////////////
// Function dquicksort performs a quick sort (smaller -> larger) 
// of the double data in item by the integer indices in index;
// data in item remain unchanged
/////////////////////////////////////////////////////////////////
namespace {
  // use this to compute permutation corresponding to a sort
  class IndexedValue {
    int index_;
    double value_;
  public:
    IndexedValue(int index, double value) : index_(index), value_(value) {}
    int index() const { return index_; }
    double value() const { return value_; }

    bool operator<(const IndexedValue& a) const {
      const double small_diff = 1.0e-12;
      if (fabs(value_-a.value_) < small_diff)
        return false;
      else
        return value_ < a.value_;
    }
  };

};


void
MOIndexSpace::dquicksort(double *item,int *index,int n)
{
  if (n<=0) return;
  typedef vector<IndexedValue> vectype;
  typedef vector<IndexedValue>::iterator iter;
  vector<IndexedValue> vals;
  for (int i=0; i<n; i++) {
    IndexedValue val(i,item[i]);
    vals.push_back(val);
  }
  stable_sort(vals.begin(),vals.end());
  for (int i=0; i<n; i++) {
    index[i] = vals.at(i).index();
  }
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
