//
// orbitalspace.cc
//
// Copyright (C) 2003 Edward Valeev
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

#include <stdexcept>
#include <algorithm>
#include <stdlib.h>
#include <cfloat>
#include <cassert>
#include <util/misc/formio.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/wfn/orbitalspace.h>

#define DEBUG_DENSE_TRANSFORM 0

using namespace std;
using namespace sc;

inline int max(int a, int b) {
  return (a > b) ? a : b;
}

ParticleHoleOrbitalAttributes ParticleHoleOrbitalAttributes::None(0x0000);
ParticleHoleOrbitalAttributes ParticleHoleOrbitalAttributes::Hole(0x0001);
ParticleHoleOrbitalAttributes ParticleHoleOrbitalAttributes::Particle(0x0002);
ParticleHoleOrbitalAttributes ParticleHoleOrbitalAttributes::Any(0x0001 | 0x0002);

/*---------------
 OrbitalSpace
 ---------------*/
static ClassDesc OrbitalSpace_cd(typeid(OrbitalSpace), "OrbitalSpace", 1,
                                 "virtual public SavableState", 0, 0, create<OrbitalSpace> );

OrbitalSpace::OrbitalSpace() {}

OrbitalSpace::OrbitalSpace(const OrbitalSpace& O) :
  id_(O.id_), name_(O.name_), basis_(O.basis_), integral_(O.integral_),
      block_sizes_(O.block_sizes_) {
  init();
}

OrbitalSpace::OrbitalSpace(const std::string& id, const std::string& name,
                           const RefSCMatrix& full_coefs,
                           const RefDiagSCMatrix& evals,
                           const Ref<GaussianBasisSet>& basis,
                           const Ref<Integral>& integral,
                           const std::vector<unsigned int>& block_offsets,
                           const std::vector<unsigned int>& block_sizes) :
  id_(id), name_(name), basis_(basis), integral_(integral),
  block_sizes_(block_sizes) {
  full_coefs_to_coefs(full_coefs, evals, block_offsets, general);

  init();
}


OrbitalSpace::OrbitalSpace(const std::string& id, const std::string& name,
                           const RefSCMatrix& full_coefs,
                           const Ref<GaussianBasisSet>& basis,
                           const Ref<Integral>& integral,
                           const std::vector<unsigned int>& block_offsets,
                           const std::vector<unsigned int>& block_sizes,
                           const IndexOrder& moorder,
                           const RefDiagSCMatrix& evals) :
  id_(id), name_(name), basis_(basis), integral_(integral),
      block_sizes_(block_sizes) {
  full_coefs_to_coefs(full_coefs, evals, block_offsets, moorder);

  init();
}

OrbitalSpace::OrbitalSpace(const std::string& id, const std::string& name,
                           const RefSCMatrix& full_coefs, const Ref<
                               GaussianBasisSet>& basis,
                           const Ref<Integral>& integral,
                           const RefDiagSCMatrix& evals, unsigned int nfzc,
                           unsigned int nfzv, const IndexOrder& moorder) :
  id_(id), name_(name), basis_(basis), integral_(integral) {
  if (evals == 0)
    throw std::runtime_error(
                             "OrbitalSpace::OrbitalSpace() -- null eigenvalues matrix");
  if (nfzc > full_coefs.coldim().n())
    throw std::runtime_error("OrbitalSpace::OrbitalSpace() -- invalid nfzc");
  if (nfzc + nfzv > full_coefs.coldim().n())
    throw std::runtime_error(
                             "OrbitalSpace::OrbitalSpace() -- invalid nfzc+nfzv");
  std::vector<unsigned int> block_offsets = frozen_to_blockinfo(nfzc, nfzv,
                                                                evals);
  full_coefs_to_coefs(full_coefs, evals, block_offsets, moorder);

  init();
}

OrbitalSpace::OrbitalSpace(const std::string& id, const std::string& name,
                           const RefSCMatrix& full_coefs,
                           const Ref<GaussianBasisSet>& basis,
                           const Ref<Integral>& integral) :
  id_(id), name_(name), basis_(basis), integral_(integral) {
  Ref<SCBlockInfo> modim_blocks = full_coefs.coldim()->blocks();
  unsigned int nb = modim_blocks->nblock();
  std::vector<unsigned int> block_offsets(nb);
  block_sizes_.resize(nb);
  for (int i = 0; i < nb; i++) {
    block_offsets[i] = 0;
    block_sizes_[i] = modim_blocks->size(i);
  }

  full_coefs_to_coefs(full_coefs, 0, block_offsets, symmetry);

  init();
}

OrbitalSpace::OrbitalSpace(const std::string& id, const std::string& name,
                           const Ref<OrbitalSpace>& orig_space,
                           const RefSCMatrix& new_coefs,
                           const Ref<GaussianBasisSet>& new_basis) :
  id_(id), name_(name), integral_(orig_space->integral()),
      orbsym_(orig_space->orbsym_), evals_(orig_space->evals_),
      block_sizes_(orig_space->block_sizes_) {
  if (orig_space->rank() != new_coefs.coldim()->n())
    throw std::runtime_error(
                             "OrbitalSpace::OrbitalSpace() -- new_coefs have different number of orbitals");
  coefs_ = new_coefs;
  basis_ = new_basis;
  init();
}

OrbitalSpace::OrbitalSpace(StateIn& si) :
  SavableState(si) {
  si.get(id_);
  si.get(name_);

  coefs_.restore(si);
  evals_.restore(si);
  basis_ << SavableState::restore_state(si);
  integral_ << SavableState::restore_state(si);
  si.get(orbsym_);
  si.get(block_sizes_);

  init();
}

OrbitalSpace::~OrbitalSpace() {
}

void OrbitalSpace::save_data_state(StateOut& so) {
  so.put(id_);
  so.put(name_);

  coefs_.save(so);
  evals_.save(so);
  SavableState::save_state(basis_.pointer(), so);
  SavableState::save_state(integral_.pointer(), so);
  so.put(orbsym_);
  so.put(block_sizes_);
}

const std::string&
OrbitalSpace::id() const {
  return id_;
}
const std::string&
OrbitalSpace::name() const {
  return name_;
}

const RefSCDimension&
OrbitalSpace::dim() const {
  return dim_;
}

const Ref<GaussianBasisSet>&
OrbitalSpace::basis() const {
  return basis_;
}

const Ref<Integral>&
OrbitalSpace::integral() const {
  return integral_;
}

const RefSCMatrix&
OrbitalSpace::coefs() const {
  return coefs_;
}

RefSCMatrix
OrbitalSpace::coefs_nb() const {
  RefSCMatrix result = coefs_;

  SCMatrixKit* kit = coefs_.kit().pointer();
  BlockedSCMatrixKit* bkit = dynamic_cast<BlockedSCMatrixKit*>(kit);
  if (bkit) { // if using blocked kit, rebuild using subkit (which should be nonblocked
    SCMatrixKit* subkit = bkit->subkit().pointer();
    BlockedSCMatrixKit* bsubkit = dynamic_cast<BlockedSCMatrixKit*>(subkit);
    while (bsubkit) {
      subkit = bsubkit->subkit().pointer();
      bsubkit = dynamic_cast<BlockedSCMatrixKit*>(subkit);
    }

    result = subkit->matrix(coefs_.rowdim(), coefs_.coldim());
    result->convert(coefs_);
  }
  // else do nothing

  return result;
}

const RefDiagSCMatrix&
OrbitalSpace::evals() const {
  return evals_;
}

const std::vector<unsigned int>&
OrbitalSpace::orbsym() const {
  return orbsym_;
}

unsigned int OrbitalSpace::rank() const {
  return dim_.n();
}

unsigned int OrbitalSpace::nblocks() const {
  return block_sizes_.size();
}

const std::vector<unsigned int>&
OrbitalSpace::block_sizes() const {
  return block_sizes_;
}

void OrbitalSpace::check_orbsym() const {
  const int ng = basis_->molecule()->point_group()->char_table().order();

  for (std::vector<unsigned int>::const_iterator p = orbsym_.begin();
       p != orbsym_.end();
       ++p) {
    if (*p >= ng)
      throw ProgrammingError("OrbitalSpace::check_orbsym() -- invalid value in the list of orbital irreps",
                             __FILE__,__LINE__,class_desc());
  }
}

std::vector<unsigned int> OrbitalSpace::frozen_to_blockinfo(unsigned int nfzc,
                                                            unsigned int nfzv,
                                                            const RefDiagSCMatrix& evals) {
  unsigned int rank = evals.dim().n();

  unsigned int nb = evals.dim()->blocks()->nblock();
  std::vector<unsigned int> offsets(nb);
  block_sizes_.resize(nb);
  offsets.resize(nb);
  for (unsigned int b = 0; b < nb; b++) {
    block_sizes_[b] = evals.dim()->blocks()->size(b);
    offsets[b] = 0;//initialize offsets vector
  }

  // Get the energies of the orbitals in this space
  double* energy = new double[rank];
  unsigned int* index_map = new unsigned int[rank];
  std::vector<unsigned int> blocked_index_to_irrep(rank);
  unsigned int ii = 0; // blocked index to this space
  unsigned int offset = 0;
  for (unsigned int b = 0; b < nb; b++) {
    for (unsigned int i = 0; i < block_sizes_[b]; i++, ii++) {
      energy[ii] = evals.get_element(i + offset);
      blocked_index_to_irrep[ii] = b; //therefore evals are ordered according to symmetry
    }
    offset += block_sizes_[b];
  }

  // Do the sort
  dquicksort(energy, index_map, rank);//index_map contains the orbital indices in energy ordering
  delete[] energy;

  // Get rid of nfzc lowest orbitals
  for (unsigned int i = 0; i < nfzc; i++) {
    unsigned int b = blocked_index_to_irrep[index_map[i]];
    ++offsets[b];    //tells the starting orbital in each irrep, excluding fzc; makes sense if the orbitals in each sym block are already in energy ordering
    --block_sizes_[b];
  }

  // Get rid of nfzv highest orbitals
  for (unsigned int i = rank - nfzv; i < rank; i++) {
    unsigned int b = blocked_index_to_irrep[index_map[i]];
    --block_sizes_[b];
  }

  delete[] index_map;
  return offsets;
}

void OrbitalSpace::full_coefs_to_coefs(const RefSCMatrix& full_coefs,
                                       const RefDiagSCMatrix& evals,
                                       const std::vector<unsigned int>& offsets,
                                       IndexOrder moorder) {
  // compute the rank of this
  unsigned int rank = 0;
  for (vector<unsigned int>::const_iterator p = block_sizes_.begin();
       p != block_sizes_.end();
       ++p) {
    rank += *p;
  }//rank equals # of active orbitals

  orbsym_.resize(rank);
  RefSCDimension modim = full_coefs.coldim(); // the dimension of the full space

  // In general vectors are ordered differently from the original
  unsigned int* index_map = new unsigned int[rank]; // maps index in this (sorted) space to this (blocked) space
  std::vector<unsigned int> blocked_subindex_to_full_index(rank); // maps index from this space(in blocked form) into the full space
  std::vector<unsigned int> blocked_subindex_to_irrep(rank); // maps index from this space(in blocked form) to the irrep
  if (moorder == symmetry) {
    // coefs_ has the same number of blocks as full_coefs_
    const unsigned int nb = modim->blocks()->nblock();
    int* nfunc_per_block = new int[nb];
    for (unsigned int i = 0; i < nb; i++)
      nfunc_per_block[i] = block_sizes_[i];
    dim_ = new SCDimension(rank, nb, nfunc_per_block,
                           ("MO(" + name_ + ")").c_str());
    if (rank) {
      for (unsigned int i = 0; i < nb; i++)
        dim_->blocks()->set_subdim(i, new SCDimension(nfunc_per_block[i]));
    }
    delete[] nfunc_per_block;

    // The sorted->blocked reordering array is trivial when no resorting is done
    for (unsigned int i = 0; i < rank; i++) {
      index_map[i] = i;
    }//it will change for energy ordering

    unsigned int ii = 0; // blocked index to this space
    unsigned int offset = 0;
    for (unsigned int b = 0; b < nb; b++) {
      for (unsigned int i = 0; i < block_sizes_[b]; i++, ii++) { //only active (non-virtual) orbs are considered (due to block_sizes_, see frozen_to_blockinfo())
        blocked_subindex_to_full_index[ii] = i + offsets[b] + offset;
        blocked_subindex_to_irrep[ii] = b;
      }
      offset += modim->blocks()->size(b);
    }
  } else if (moorder == energy) {
    //
    // Sort vectors by their energy
    //

    // Get the energies of the orbitals in this space
    double* energy = new double[rank];
    const unsigned int nb = block_sizes_.size();
    unsigned int ii = 0; // blocked index to this space
    unsigned int offset = 0;
    for (unsigned int b = 0; b < nb; b++) {
      for (unsigned int i = 0; i < block_sizes_[b]; i++, ii++) {//only active (non-virtual) orbs are considered (due to block_sizes_, see frozen_to_blockinfo())
        energy[ii] = evals.get_element(i + offsets[b] + offset);
        blocked_subindex_to_full_index[ii] = i + offsets[b] + offset;
        blocked_subindex_to_irrep[ii] = b;
      }
      offset += modim->blocks()->size(b);
    }

    // Do the sort
    dquicksort(energy, index_map, rank);//index_map is for active orbitals

    // coefs_ has 1 block
    int* nfunc_per_block = new int[1];
    nfunc_per_block[0] = rank;
    dim_ = new SCDimension(rank, 1, nfunc_per_block,
                           ("MO(" + name_ + ")").c_str());//one block
    if (rank)
      dim_->blocks()->set_subdim(0, new SCDimension(nfunc_per_block[0]));

    // Recompute block_sizes_ to conform the energy ordering
    block_sizes_.resize(1);
    block_sizes_[0] = rank;

    delete[] energy;
    delete[] nfunc_per_block;
  } else
    throw std::runtime_error("OrbitalSpace::full_coefs_to_coefs() -- moorder should be either energy or symmetry");

  // Copy required columns of full_coefs_ into coefs_
  RefSCDimension aodim = full_coefs.rowdim();
  Ref<SCMatrixKit> so_matrixkit = basis_->so_matrixkit();
  coefs_ = so_matrixkit->matrix(aodim, dim_);
  evals_ = so_matrixkit->diagmatrix(dim_);
  for (unsigned int i = 0; i < rank; i++) {
    const unsigned int ii = blocked_subindex_to_full_index[index_map[i]];
    orbsym_[i] = blocked_subindex_to_irrep[index_map[i]];
    for (unsigned int j = 0; j < aodim.n(); j++) {
      coefs_(j, i) = full_coefs(j, ii);
    }
  }
  if (evals)
    for (unsigned int i = 0; i < rank; i++) {
      const unsigned int ii = blocked_subindex_to_full_index[index_map[i]];
      evals_( i) = evals(ii);
    }
  else
    evals_.assign(0.0);

  delete[] index_map;
}

void OrbitalSpace::init() {
  if (id_.size() > max_id_length)
    throw ProgrammingError(
                           "OrbitalSpace constructed with id longer than allowed",
                           __FILE__,__LINE__);

  dim_ = evals_.dim();
}

void OrbitalSpace::init(const std::string& id, const std::string& name,
                        const Ref<GaussianBasisSet>& basis, const Ref<Integral>& integral,
                        const RefSCMatrix& coefs,
                        const RefDiagSCMatrix& evals,
                        const std::vector<unsigned int>& orbsyms,
                        unsigned int nblocks,
                        const std::vector<BlockedOrbital>& indexmap) {
  id_ = id;
  name_ = name;
  basis_ = basis;
  integral_ = integral;

  const unsigned int norbs = indexmap.size();
  const bool debug = false;
  // compute block sizes
  block_sizes_.resize(nblocks);  std::fill(block_sizes_.begin(), block_sizes_.end(), 0);
  for(int b=0; b<nblocks; ++b) {
    for(size_t o=0; o<norbs; ++o) {
      if (indexmap[o].attr() == b)
        ++block_sizes_[b];
    }
  }
  // compute block offsets
  block_offsets_.resize(nblocks);  std::fill(block_offsets_.begin(), block_offsets_.end(), 0);
  for(int b=1; b<nblocks; ++b)
    block_offsets_[b] = block_offsets_[b-1] + block_sizes_[b-1];

  // build new blocked dimension
  int* nfunc_per_block = new int[nblocks];
  for (unsigned int i = 0; i < nblocks; i++)
    nfunc_per_block[i] = block_sizes_[i];
  dim_ = new SCDimension(norbs, nblocks, nfunc_per_block, id.c_str());
  for (unsigned int i = 0; i < nblocks; i++)
    dim_->blocks()->set_subdim(i, new SCDimension(nfunc_per_block[i]));
  delete[] nfunc_per_block;

  // Map coefficients, eigenvalues, etc.
  RefSCDimension aodim = coefs.rowdim();
  Ref<SCMatrixKit> so_matrixkit = basis_->so_matrixkit();
  coefs_ = so_matrixkit->matrix(aodim, dim_);
  evals_ = so_matrixkit->diagmatrix(dim_);
  orbsym_.resize(norbs);
  if(debug)
  {
    coefs->print("debug orbspace: print coefs");
    evals->print("debug orbspace: print evals");
  }
  for (unsigned int i = 0; i < norbs; ++i) {
    const unsigned int ii = indexmap[i].index();
    orbsym_[i] = orbsyms[ii];
    for (unsigned int j = 0; j < aodim.n(); j++)
    {
      coefs_(j, i) = coefs(j, ii);
    }
    evals_(i) = evals(ii);
  }
  if(debug)
  {
    coefs_->print("debug orbspace: print coefs");
    evals_->print("debug orbspace: print evals");
  }

  init();
}

size_t OrbitalSpace::memory_in_use() const {
  size_t memory = (size_t) basis_->nbasis() * rank() * sizeof(double);
  return memory;
}

void OrbitalSpace::print(ostream&o) const {
  o << indent << "OrbitalSpace \"" << name_ << "\":" << endl;
  o << incindent;
  o << indent << "Basis Set:" << endl;
  o << incindent;
  basis_->print(o);
  o << decindent << endl;
  o << decindent;
}

void OrbitalSpace::print_detail(ostream&o) const {
  o << indent << "OrbitalSpace \"" << name_ << "\":" << endl;
  o << incindent;
  o << indent << "id = " << id_ << endl;
  o << indent << "Basis Set:" << endl;
  o << incindent;
  basis_->print(o);
  o << decindent << endl;
  evals_.print("Eigenvalues");
  coefs_.print("Coefficients");
  o << decindent;
}

void OrbitalSpace::print_summary(ostream& o) const {
  o << indent << "OrbitalSpace \"" << name_ << "\" (id = " << id_ << "):" << endl;
  o << incindent;
  o << indent << "GaussianBasisSet \"" << basis_->name() << "\"" << endl;
  o << indent << "  rank  nbasis  nshell  nfuncmax" << endl;
  o << indent << scprintf("  %-6i %-6i  %-6i   %-6i", rank(), basis_->nbasis(),
                          basis_->nshell(), basis_->max_nfunction_in_shell())
      << endl;
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
      IndexedValue(int index, double value) :
        index_(index), value_(value) {
      }
      int index() const {
        return index_;
      }
      // not used at the moment
#if 0
      double value() const {
        return value_;
      }
#endif

      bool operator<(const IndexedValue& a) const {
        const double small_diff = 1.0e-12;
        if (fabs(value_ - a.value_) < small_diff)
          return false;
        else
          return value_ < a.value_;
      }
  };

}
;

void OrbitalSpace::dquicksort(double *item, unsigned int *index, unsigned int n) {
  typedef std::vector<IndexedValue> vectype;
  typedef std::vector<IndexedValue>::iterator iter;
  std::vector<IndexedValue> vals;
  for (unsigned int i = 0; i < n; i++) {
    IndexedValue val(i, item[i]);
    vals.push_back(val);
  }
  stable_sort(vals.begin(), vals.end());
  for (unsigned int i = 0; i < n; i++) {
    index[i] = vals.at(i).index();
  }
}

/////////////////////////////////////////////////////////////////////////////

static ClassDesc MaskedOrbitalSpace_cd(typeid(MaskedOrbitalSpace), "MaskedOrbitalSpace", 1,
                                       "public OrbitalSpace", 0, 0, create<MaskedOrbitalSpace> );

MaskedOrbitalSpace::MaskedOrbitalSpace(StateIn& si) :
  OrbitalSpace(si) {}

void
MaskedOrbitalSpace::save_data_state(StateOut& so) {
  OrbitalSpace::save_data_state(so);
}

MaskedOrbitalSpace::MaskedOrbitalSpace(const std::string& id,
                                       const std::string& name,
                                       const Ref<OrbitalSpace>& orig_space,
                                       const std::vector<bool>& mask) :
  OrbitalSpace() {

  // validate input
  const size_t num_orig_orbs = mask.size();
  if (mask.size() != orig_space->rank())
    throw ProgrammingError("MaskedOrbitalSpace::MaskedOrbitalSpace -- mask must have same size as the original OrbitalSpace",
                           __FILE__,__LINE__);

  // create vector of BlockedOrbitals
  std::vector<BlockedOrbital> blocked_orbs;
  const unsigned int nblocks = orig_space->nblocks();
  size_t block_offset = 0;
  const std::vector<unsigned int>& block_sizes = orig_space->block_sizes();
  for(unsigned int b=0; b<nblocks; ++b) {
    const size_t block_size = block_sizes[b];
    for(size_t o=0; o<block_size; ++o) {
      const size_t oo = o + block_offset;
      if (mask[oo]) {
        blocked_orbs.push_back( BlockedOrbital(oo, b) );
      }
    }
    block_offset += block_size;
  }

  init(id, name,
       orig_space->basis(),
       orig_space->integral(),
       orig_space->coefs(),
       orig_space->evals(),
       orig_space->orbsym(),
       nblocks, blocked_orbs);

}

/////////////////////////////////////////////////////////////////////////////
namespace {
  RefSCMatrix make_unit_matrix(const Ref<GaussianBasisSet>& basis,
                               const Ref<Integral>& integral) {
    const int nao = basis->nbasis();
    RefSCDimension obs_ao_dim = new SCDimension(nao,1);
    obs_ao_dim->blocks()->set_subdim(0,new SCDimension(nao));
    Ref<PetiteList> pl = integral->petite_list(basis);
    RefSCMatrix obs_ao_coefs = basis->so_matrixkit()->matrix(pl->AO_basisdim(),obs_ao_dim);
    obs_ao_coefs.assign(0.0);
    for(int ao=0; ao<nao; ++ao)
      obs_ao_coefs.set_element(ao,ao,1.0);
    return obs_ao_coefs;
  }
}

static ClassDesc AtomicOrbitalSpace_cd(typeid(AtomicOrbitalSpace), "AtomicOrbitalSpace", 1,
                                       "public OrbitalSpace", 0, 0, create<AtomicOrbitalSpace> );

AtomicOrbitalSpace::AtomicOrbitalSpace(StateIn& si) :
  OrbitalSpace(si) {}

void
AtomicOrbitalSpace::save_data_state(StateOut& so) {
  OrbitalSpace::save_data_state(so);
}

AtomicOrbitalSpace::AtomicOrbitalSpace(const std::string& id, const std::string& name,
                                       const Ref<GaussianBasisSet>& basis,
                                       const Ref<Integral>& integral) :
  OrbitalSpace(id, name, make_unit_matrix(basis,integral), basis, integral) {}

/////////////////////////////////////////////////////////////////////////////

static ClassDesc NonblockedOrbitalSpace_cd(typeid(NonblockedOrbitalSpace), "NonblockedOrbitalSpace", 1,
                                           "public OrbitalSpace", 0, 0, create<NonblockedOrbitalSpace> );

NonblockedOrbitalSpace::NonblockedOrbitalSpace(StateIn& si) :
  OrbitalSpace(si) {}

void
NonblockedOrbitalSpace::save_data_state(StateOut& so) {
  OrbitalSpace::save_data_state(so);
}

NonblockedOrbitalSpace::NonblockedOrbitalSpace(const std::string& id,
                                               const std::string& name,
                                               const Ref<OrbitalSpace>& orig_space) :
  OrbitalSpace() {

  // create vector of BlockedOrbitals
  std::vector<BlockedOrbital> blocked_orbs;
  const unsigned int nblocks = orig_space->nblocks();
  const std::vector<unsigned int>& block_sizes = orig_space->block_sizes();
  size_t block_offset = 0;
  for(unsigned int b=0; b<nblocks; ++b) {
    const size_t block_size = block_sizes[b];
    for(size_t o=0; o<block_size; ++o) {
      const size_t oo = o + block_offset;
      blocked_orbs.push_back( BlockedOrbital(oo, 0) );
    }
    block_offset += block_size;
  }

  init(id, name,
       orig_space->basis(),
       orig_space->integral(),
       orig_space->coefs(),
       orig_space->evals(),
       orig_space->orbsym(),
       1, blocked_orbs);

}

/////////////////////////////////////////////////////////////////////////////

static ClassDesc OrbitalSpaceUnion_cd(typeid(OrbitalSpaceUnion), "OrbitalSpace  ", 1,
                                           "public OrbitalSpace", 0, 0, create<OrbitalSpaceUnion> );

OrbitalSpaceUnion::OrbitalSpaceUnion(StateIn& si) :
  OrbitalSpace(si) {}

void
OrbitalSpaceUnion::save_data_state(StateOut& so) {
  OrbitalSpace::save_data_state(so);
}

OrbitalSpaceUnion::OrbitalSpaceUnion(const std::string& id, const std::string& name,
                                     const OrbitalSpace& s1, const OrbitalSpace& s2,
                                     bool merge_blocks) :
  OrbitalSpace() {

  // validate input
  if (s1.integral()->equiv(s2.integral()) == false)
    throw std::logic_error("OrbitalSpaceUnion::OrbitalSpaceUnion -- incompatible integral matrices are not yet supported");
  if (merge_blocks && s1.nblocks() != s2.nblocks())
    throw std::logic_error("OrbitalSpaceUnion::OrbitalSpaceUnion : merge_blocks is true but different number of blocks");
  const unsigned int nblocks = merge_blocks ? s1.nblocks() : s1.nblocks() + s2.nblocks();

  // if trying to add a space to itself and blocks to be merged, throw
  if (merge_blocks && &s1 == &s2)
    MPQC_ASSERT(false);

  // compute the basis set
  // try 3 cases:
  // s1.basis() is contained in s2
  // s2.basis() is contained in s1
  // neither, get s1.basis() + s2.basis()
  Ref<GaussianBasisSet> bs = 0;
  try { std::vector<unsigned int> tmp = *(s1.basis()) << *(s2.basis()); bs = s1.basis(); }
  catch (...) { }
  if (bs == 0) {
    try { std::vector<unsigned int> tmp = *(s2.basis()) << *(s1.basis()); bs = s2.basis(); }
    catch (...) { }
  }
  if (bs == 0) {
    bs = s1.basis() + s2.basis();
  }
  const std::vector<unsigned int> map1 = (*bs) << (*s1.basis());
  const std::vector<unsigned int> map2 = (*bs) << (*s2.basis());

  const unsigned int nao1 = s1.basis()->nbasis();
  const unsigned int norbs1 = s1.rank();
  const vector<unsigned int>& orbsym1 = s1.orbsym();
  const unsigned int nao2 = s2.basis()->nbasis();
  const unsigned int norbs2 = s2.rank();
  const vector<unsigned int>& orbsym2 = s2.orbsym();

  // merge coefficient, energy and symmetry data, don't pay attention to blocking yet
  Ref<SCMatrixKit> so_matrixkit = bs->so_matrixkit();
  const int norbs = norbs1 + norbs2;
  const int nbs = bs->nbasis();
  RefSCDimension dim = new SCDimension(norbs, 1, &norbs);
  dim->blocks()->set_subdim(0, new SCDimension(norbs) );
  Ref<Integral> intf = s1.integral()->clone();
  Ref<PetiteList> plist = intf->petite_list(bs);
  RefSCDimension aodim = plist->AO_basisdim();
  RefSCMatrix coefs = so_matrixkit->matrix(aodim, dim);  coefs->assign(0.0);
  RefDiagSCMatrix evals = so_matrixkit->diagmatrix(dim);
  std::vector<unsigned int> orbsym(dim.n());
  int ii = 0;
  for (unsigned int i=0; i<norbs1; ++i, ++ii) {
    orbsym[ii] = orbsym1[i];
    for (unsigned int j=0; j<nao1; j++) {
      const unsigned int jj = map1[j];
      coefs(jj, ii) = s1.coefs().get_element(j, i);
    }
    evals(ii) = s1.evals().get_element(i);
  }
  for (unsigned int i=0; i<norbs2; ++i, ++ii) {
    orbsym[ii] = orbsym2[i];
    for (unsigned int j=0; j<nao2; j++) {
      const unsigned int jj = map2[j];
      coefs(jj, ii) = s2.coefs().get_element(j, i);
    }
    evals(ii) = s2.evals().get_element(i);
  }

  // create vector of BlockedOrbitals
  std::vector<BlockedOrbital> blocked_orbs;
  const std::vector<unsigned int>& block_sizes1 = s1.block_sizes();
  const std::vector<unsigned int>& block_sizes2 = s2.block_sizes();
  // blocking info depends on merge_blocks
  if (!merge_blocks) {
    size_t block_offset = 0;
    unsigned int bb = 0;
    for (unsigned int b = 0; b < block_sizes1.size(); ++b, ++bb) {
      const size_t block_size = block_sizes1[b];
      for (size_t o = 0; o < block_size; ++o) {
        const size_t oo = o + block_offset;
        blocked_orbs.push_back(BlockedOrbital(oo, bb));
      }
      block_offset += block_size;
    }
    for (unsigned int b = 0; b < block_sizes2.size(); ++b, ++bb) {
      const size_t block_size = block_sizes2[b];
      for (size_t o = 0; o < block_size; ++o) {
        const size_t oo = o + block_offset;
        blocked_orbs.push_back(BlockedOrbital(oo, bb));
      }
      block_offset += block_size;
    }
  }
  else {
    // if merging blocks, get rid of duplicate orbitals?
    // Since merging may produce linear dependencies this doesn't solve all problems
    // hence the user should check the result for linear dependencies, etc.

    blocked_orbs.resize(norbs);
    // include orbitals from space1, block-by-block
    size_t block_offset = 0;
    int current = 0;
    for (unsigned int b = 0; b < block_sizes1.size(); ++b) {
      const size_t block_size1 = block_sizes1[b];
      const size_t block_size2 = block_sizes2[b];
      for (size_t o = 0; o < block_size1; ++o) {
        const size_t oo = o + block_offset;
        blocked_orbs[oo] = BlockedOrbital(current, b);
        ++current;
      }
      block_offset += block_size1 + block_size2;
    }
    // include orbitals from space2, block-by-block
    block_offset = 0;
    for (unsigned int b = 0; b < block_sizes2.size(); ++b) {
      const size_t block_size1 = block_sizes1[b];
      const size_t block_size2 = block_sizes2[b];
      block_offset += block_size1;
      for (size_t o = 0; o < block_size2; ++o) {
        const size_t oo = o + block_offset;
        blocked_orbs[oo] = BlockedOrbital(current, b);
        ++current;
      }
      block_offset += block_size2;
    }
  }

  init(id, name,
       bs,
       s1.integral(),
       coefs,
       evals,
       orbsym,
       nblocks, blocked_orbs);

}

/////////////////////////////////////////////////////////////////////////////

ClassDesc EmptyOrbitalSpace::class_desc_(
          typeid(EmptyOrbitalSpace), "EmptyOrbitalSpace", 1,
          "public OrbitalSpace", 0, 0, create<EmptyOrbitalSpace> );

EmptyOrbitalSpace::EmptyOrbitalSpace(StateIn& si) :
  OrbitalSpace(si) {}

void
EmptyOrbitalSpace::save_data_state(StateOut& so) {
  OrbitalSpace::save_data_state(so);
}

EmptyOrbitalSpace::EmptyOrbitalSpace(const std::string& id,
                                     const std::string& name,
                                     const Ref<GaussianBasisSet>& bs,
                                     const Ref<Integral>& ints,
                                     const IndexOrder& moorder) :
  OrbitalSpace() {

  const unsigned int nblocks = (moorder == symmetry) ? bs->molecule()->point_group()->char_table().nirrep()
                                                     : 1;
  std::vector<BlockedOrbital> blocked_orbs;
  std::vector<int> func_per_blk(nblocks, 0);
  RefSCDimension emptydim = new SCDimension(0, nblocks, &(func_per_blk[0]), "");
  for(int b=0; b<nblocks; ++b)  emptydim->blocks()->set_subdim(b, new SCDimension(0));
  const int nao = bs->nbasis();
  RefSCDimension aodim = new SCDimension(nao,1);
  aodim->blocks()->set_subdim(0,new SCDimension(nao));
  RefSCMatrix coefs = bs->so_matrixkit()->matrix(aodim, emptydim);
  RefDiagSCMatrix evals = bs->so_matrixkit()->diagmatrix(emptydim);
  std::vector<unsigned int> orbsym;

  init(id, name,
       bs,
       ints,
       coefs,
       evals,
       orbsym,
       nblocks, blocked_orbs);
}

/////////////////////////////////////////////////////////////////////////////

MOIndexMap sc::operator<<(const OrbitalSpace& s2, const OrbitalSpace& s1) {
  typedef std::vector<int> OrbitalMap;
  const bool must_use_same_basis = true;
  OrbitalMap int_map = map(s2, s1, must_use_same_basis);
  OrbitalMap::iterator i = std::find(int_map.begin(), int_map.end(), -1);
  if (i != int_map.end())
    throw CannotConstructMap();

  MOIndexMap result(int_map.size());
  std::copy(int_map.begin(), int_map.end(), result.begin());
  return result;
}

std::vector<int>
sc::map(const OrbitalSpace& s2, const OrbitalSpace& s1, bool expect_same_basis) {
  const unsigned int rank1 = s1.rank();
  const unsigned int rank2 = s2.rank();
  if (!s1.integral()->equiv(s2.integral()))
    throw CannotConstructMap();
  const bool same_basis = (s1.basis() == s2.basis());
  if (expect_same_basis && !same_basis)
    throw CannotConstructMap();
  std::vector<int> bfmap;
  if (!same_basis) {
    try {
      bfmap = sc::map(*(s2.basis()),*(s1.basis()));
    }
    catch (...) {
      throw CannotConstructMap();
    }
  }

  const RefSCMatrix& c1 = s1.coefs().t();
  const RefSCMatrix& c2 = s2.coefs().t();
#if 0
  c1.print("map(OrbitalSpace,OrbitalSpace): c1");
  c2.print("map(OrbitalSpace,OrbitalSpace): c2");
#endif
  const unsigned int nao1 = c1.coldim().n();
  const unsigned int nao2 = c2.coldim().n();

  typedef std::vector<int> maptype;
  maptype result(rank1);

  // if objects are the same, map is trivial
  if (&s1 == &s2) {
    for (unsigned int mo1 = 0; mo1 < rank1; mo1++)
      result[mo1] = mo1;
    return result;
  }

  // for each MO in 1 find vector in 2
  for (unsigned int mo1 = 0; mo1 < rank1; mo1++) {
    bool found_match = false;
    for (unsigned int mo2 = 0; mo2 < rank2; mo2++) {
      bool vectors_do_not_match = false;
      // compare vectors
      for (int ao1 = 0; ao1 < nao1; ++ao1) {
        const int ao2 = same_basis ? ao1 : bfmap[ao1];
        if (ao2 == -1) { // ao1 does not have a match in s2.basis() -- this is OK as long as its coefficient is 0
          if (fabs(c1.get_element(mo1, ao1)) > 1.0e-12) {
            vectors_do_not_match = true;
            break;
          }
        }
        else if (fabs(c1.get_element(mo1, ao1) - c2.get_element(mo2, ao2)) > 1.0e-12) {
          vectors_do_not_match = true;
#if 0
          ExEnv::out0() << "map(OrbitalSpace,OrbitalSpace): (mo1,mo2,ao) = "
          << mo1 << "," << mo2 << "," << ao << "  delta = "
          << fabs(c1.get_element(mo1,ao)-c2.get_element(mo2,ao)) << std::endl;
#endif
          break;
        }
      }
      // go to next mo1 if found match
      if (!vectors_do_not_match) {
        found_match = true;
#if 0
        ExEnv::out0() << "map(OrbitalSpace,OrbitalSpace): found match (mo1,mo2) = "
        << mo1 << "," << mo2 << std::endl;
#endif
        result[mo1] = mo2;
        mo2 = rank2;
      }
    }
    // if this mo1 doesn't match any mo2 -- punt
    if (!found_match)
      result[mo1] = -1;
  }

#if 0
  ExEnv::out0() << "map(OrbitalSpace,OrbitalSpace):" << endl;
  for(int i=0; i<rank1; ++i)
    ExEnv::out0() << result[i] << endl;
#endif

  return result;
}

sc::SparseMOIndexMap sc::sparsemap(const OrbitalSpace& s2,
                                   const OrbitalSpace& s1, double hardzero) {
  const unsigned int rank1 = s1.rank();
  const unsigned int rank2 = s2.rank();
  if (rank1 > rank2 || s1.basis() != s2.basis() || s1.integral()->class_desc()
      != s2.integral()->class_desc())
    throw CannotConstructMap();

  const RefSCMatrix& c1 = s1.coefs().t();
  const RefSCMatrix& c2 = s2.coefs().t();
#if 0
  c1.print("sparsemap(OrbitalSpace,OrbitalSpace): c1");
  c2.print("sparsemap(OrbitalSpace,OrbitalSpace): c2");
#endif
  const unsigned int nao = c1.coldim().n();

  typedef SparseMOIndexMap maptype;
  maptype map(rank1);

  // if objects are the same, map is trivial
  if (&s1 == &s2) {
    for (unsigned int mo1 = 0; mo1 < rank1; mo1++)
      map[mo1] = make_pair(mo1, 1.0);
    return map;
  }

  std::vector<int> has_been_mapped(rank2, -1);
  // for each MO in 1 find vector in 2
  for (unsigned int mo1 = 0; mo1 < rank1; mo1++) {
    bool found_match = false;
    for (unsigned int mo2 = 0; mo2 < rank2; mo2++) {
      // if mo2 is not yet mapped by one of MOs in 1
      if (has_been_mapped[mo2] != -1)
        continue;

      // compare vectors assuming same phase
      bool vectors_do_not_match = false;
      double phase = 1.0;
      for (unsigned int ao = 0; ao < nao; ao++) {
        if (fabs(c1.get_element(mo1, ao) - phase * c2.get_element(mo2, ao))
            > hardzero) {
          vectors_do_not_match = true;
#if 0
          ExEnv::out0() << "sparsemap(OrbitalSpace,OrbitalSpace): (mo1,mo2,ao,phase) = "
          << mo1 << "," << mo2 << "," << ao << "," << phase << "  delta = "
          << fabs(c1.get_element(mo1,ao)-phase*c2.get_element(mo2,ao)) << std::endl;
#endif
          break;
        }
      }

      // compare vectors assuming opposite phase
      if (vectors_do_not_match) {
        vectors_do_not_match = false;
        phase = -1.0;
        for (unsigned int ao = 0; ao < nao; ao++) {
          if (fabs(c1.get_element(mo1, ao) - phase * c2.get_element(mo2, ao))
              > hardzero) {
            vectors_do_not_match = true;
#if 0
            ExEnv::out0() << "sparsemap(OrbitalSpace,OrbitalSpace): (mo1,mo2,ao,phase) = "
            << mo1 << "," << mo2 << "," << ao << "," << phase << "  delta = "
            << fabs(c1.get_element(mo1,ao)-phase*c2.get_element(mo2,ao)) << std::endl;
#endif
            break;
          }
        }
      }

      // go to next mo1 if found match
      if (!vectors_do_not_match) {
        found_match = true;
#if 0
        ExEnv::out0() << "sparsemap(OrbitalSpace,OrbitalSpace): found match (mo1,mo2,phase) = "
        << mo1 << "," << mo2 << "," << phase << std::endl;
#endif
        map[mo1] = make_pair(mo2, phase);
        has_been_mapped[mo2] = 1;
        mo2 = rank2;
      }
    }
    // if this mo1 doesn't match any mo2 -- punt
    if (!found_match)
      throw CannotConstructMap();
  }

  return map;
}

RefSCMatrix sc::transform(const OrbitalSpace& s2, const OrbitalSpace& s1,
                          const Ref<SCMatrixKit>& kit) {
  const unsigned int rank1 = s1.rank();
  const unsigned int rank2 = s2.rank();
  // simple tests for whether the transform can be constructed
  if (rank1 < rank2 || s1.integral()->class_desc()
      != s2.integral()->class_desc())
    throw CannotConstructMap();
  Ref<Integral> integral = s1.integral()->clone();
  const Ref<GaussianBasisSet>& bs1 = s1.basis();
  const Ref<GaussianBasisSet>& bs2 = s2.basis();

  const RefSCMatrix& c1 = s1.coefs().t();
  const RefSCMatrix& c2 = s2.coefs().t();
#if DEBUG_DENSE_TRANSFORM
  c1.print("transform(OrbitalSpace,OrbitalSpace): c1");
  c2.print("transform(OrbitalSpace,OrbitalSpace): c2");
#endif

  RefSCMatrix tform(c2.rowdim(), c1.rowdim(), kit);
  tform.assign(0.0);

  // if objects are the same, the transform is trivial
  if (&s1 == &s2) {
    RefDiagSCMatrix unit(c2.rowdim(), kit);
    unit.assign(1.0);
    tform.accumulate(unit);
    return tform;
  }

  // compute the overlap matrix in AO basis
  RefSCMatrix S21 = overlap(s2, s1, kit);

  // the transform matrix is the inverse of S21
  RefSCMatrix U21 = S21.gi();

#if DEBUG_DENSE_TRANSFORM
  // test
  (S21).print("S21");
  (U21).print("U21");
  (U21 * S21).print("U21 * S21");
  (S21 * U21).print("S21 * U21");
#endif

  return U21;
}

RefSCMatrix sc::overlap(const OrbitalSpace& space1, const OrbitalSpace& space2,
                        const Ref<SCMatrixKit>& kit) {
  if (!space1.integral()->equiv(space2.integral()))
    throw ProgrammingError(
                           "sc::overlap(s1,s2) -- s1 and s2 are supported by incompatible Integral factories");
  const Ref<GaussianBasisSet> bs1 = space1.basis();
  const Ref<GaussianBasisSet> bs2 = space2.basis();
  const bool bs1_eq_bs2 = (bs1 == bs2);
  int nshell1 = bs1->nshell();
  int nshell2 = bs2->nshell();

  RefSCMatrix vec1t = space1.coefs().t();
  RefSCMatrix vec2 = space2.coefs();

  Ref<Integral> localints = space1.integral()->clone();
  localints->set_basis(bs1, bs2);
  Ref<OneBodyInt> ov_ints = localints->overlap();

  // form AO moment matrices
  RefSCDimension aodim1 = vec1t.coldim();
  RefSCDimension aodim2 = vec2.rowdim();
  RefSCMatrix s(aodim1, aodim2, vec1t.kit());
  s.assign(0.0);

  for (int sh1 = 0; sh1 < nshell1; sh1++) {
    int bf1_offset = bs1->shell_to_function(sh1);
    int nbf1 = bs1->shell(sh1).nfunction();

    int sh2max;
    if (bs1_eq_bs2)
      sh2max = sh1;
    else
      sh2max = nshell2 - 1;

    for (int sh2 = 0; sh2 <= sh2max; sh2++) {
      int bf2_offset = bs2->shell_to_function(sh2);
      int nbf2 = bs2->shell(sh2).nfunction();

      ov_ints->compute_shell(sh1, sh2);
      const double *ovintsptr = ov_ints->buffer();

      int bf1_index = bf1_offset;
      for (int bf1 = 0; bf1 < nbf1; bf1++, bf1_index++, ovintsptr += nbf2) {
        int bf2_index = bf2_offset;
        const double *ptr = ovintsptr;
        int bf2max;
        if (bs1_eq_bs2 && sh1 == sh2)
          bf2max = bf1;
        else
          bf2max = nbf2 - 1;
        for (int bf2 = 0; bf2 <= bf2max; bf2++, bf2_index++) {

          s.set_element(bf1_index, bf2_index, *(ptr++));

        }
      }
    }
  }

  // and clean up a bit
  ov_ints = 0;

  // Symmetrize matrices, if necessary
  if (bs1_eq_bs2) {

    const int nbasis = bs1->nbasis();

    for (int bf1 = 0; bf1 < nbasis; bf1++)
      for (int bf2 = 0; bf2 <= bf1; bf2++) {
        s(bf2, bf1) = s(bf1, bf2);
      }

  }

  // finally, transform
  RefSCMatrix s12 = vec1t * s * vec2;
  // and convert to S (kits may differ!)
  double* ss = new double[s12.rowdim().n() * s12.coldim().n()];
  s12.convert(ss);
  s12 = 0;
  RefSCMatrix S(vec1t.rowdim(), vec2.coldim(), kit);
  S.assign(ss);

  // and clean up a bit
  s = 0;

  return S;
}

bool sc::in(const OrbitalSpace& s1, const OrbitalSpace& s2) {
  bool result = true;
  try {
    s2 << s1;
  } catch (CannotConstructMap& e) {
    result = false;
  }
  return result;
}

std::pair<std::string, Ref<OrbitalSpace> > sc::make_keyspace_pair(const Ref<
    OrbitalSpace>& space, SpinCase1 spin) {
  return std::make_pair(ParsedOrbitalSpaceKey::key(space->id(), spin), space);
}

std::string sc::new_unique_key(const Ref<OrbitalSpaceRegistry> & oreg)
{
  static unsigned int index = 0;
  const std::string id = "Z";
  bool found_in_registry;
  do {
    ++index;
    found_in_registry = false;
    for(int s=AnySpinCase1; s!=InvalidSpinCase1; ++s) {
      const SpinCase1 spin = static_cast<SpinCase1>(s);
      std::ostringstream oss;
      oss << id << index;
      const std::string trial_key = ParsedOrbitalSpaceKey::key(oss.str(), spin);
      if (oreg->key_exists(trial_key)) {
        found_in_registry = true;
        break;
      }
    }
  } while (found_in_registry);
  std::ostringstream oss;
  oss << id << index;
  return oss.str();
}

bool sc::operator==(const OrbitalSpace& space1, const OrbitalSpace& space2) {
  if (&space1 == &space2) return true;
  if (!space1.integral()->equiv(space2.integral()) || space1.rank()
      != space2.rank() || space1.nblocks() != space2.nblocks()
      || space1.block_sizes() != space2.block_sizes()
      || !space1.basis()->equiv(space2.basis()) || (space1.coefs()
      - space2.coefs())->maxabs() > DBL_EPSILON)
    return false;
  else
    return true;
}

/////////////////////////////////////////////////////////////////////////////

namespace {
  // pop off str from beginning up to token.
  std::string pop_till_token(std::string& str, char token) {
    const size_t next_token_pos = str.find_first_of(token);
    std::string result;
    if (next_token_pos != std::string::npos) {
      result = str.substr(0, next_token_pos);
      str.erase(0, next_token_pos + 1);
    } else {
      result = str;
      str.clear();
    }
    return result;
  }
}

ParsedOrbitalSpaceKey::ParsedOrbitalSpaceKey(const std::string& key) :
  key_(key) {
  std::string keycopy(key);
  label_ = pop_till_token(keycopy, '[');
  if (!keycopy.empty()) {
    const std::string spinkey = pop_till_token(keycopy, ']');
    if (spinkey == std::string("A"))
      spin_ = Alpha;
    else if (spinkey == std::string("B"))
      spin_ = Beta;
    else
      throw exception();
  } else {
    spin_ = AnySpinCase1;
  }
}

std::string ParsedOrbitalSpaceKey::key(const std::string& label, SpinCase1 spin) {
  if (spin != AnySpinCase1) {
    std::ostringstream oss;
    oss << label << "[" << (spin == Alpha ? "A" : "B") << "]";
    return oss.str();
  }
  return label;
}

/////////////////////////////////////////////////////////////////////////////

ParsedTransformedOrbitalSpaceKey::ParsedTransformedOrbitalSpaceKey(const std::string& key) :
  key_(key) {
  std::string keycopy(key);

  const std::string tformed_key = pop_till_token(keycopy, '_');
  if (keycopy.empty()) throw exception();
  ParsedOrbitalSpaceKey parsed_tformed_key(tformed_key);
  label_ = parsed_tformed_key.label();
  spin_ = parsed_tformed_key.spin();

  const std::string oper_key = pop_till_token(keycopy, '(');
  if (keycopy.empty()) throw exception();
  transform_operator_ = OneBodyOper::invalid;
  for (int o = OneBodyOper::gamma; o != OneBodyOper::invalid; ++o) {
    OneBodyOper::type otype = static_cast<OneBodyOper::type> (o);
    if (oper_key == OneBodyOper::to_string(otype)) {
      transform_operator_ = otype;
      break;
    }
  }
  if (transform_operator_ == OneBodyOper::invalid)
    throw exception();

  const std::string orig_key = pop_till_token(keycopy, ')');
  if (not keycopy.empty()) throw exception();
  ParsedOrbitalSpaceKey parsed_orig_key(orig_key);
  support_label_ = parsed_orig_key.label();
  support_spin_ = parsed_orig_key.spin();
}

std::string
ParsedTransformedOrbitalSpaceKey::key(const std::string& label, SpinCase1 spin,
                                      const std::string& original_label,
                                      SpinCase1 original_spin, OneBodyOper::type oper) {
  std::string result(ParsedOrbitalSpaceKey::key(label,spin) + "_" + OneBodyOper::to_string(oper) +
      "(" + ParsedOrbitalSpaceKey::key(original_label,original_spin) + ")");
  return result;
}

bool
ParsedTransformedOrbitalSpaceKey::valid_key(const std::string& key) {
  std::string keycopy(key);

  const std::string tformed_key = pop_till_token(keycopy, '_');
  if (keycopy.empty()) return false;
  ParsedOrbitalSpaceKey parsed_tformed_key(tformed_key);

  const std::string oper_key = pop_till_token(keycopy, '(');
  if (keycopy.empty()) return false;
  OneBodyOper::type transform_operator = OneBodyOper::invalid;
  for (int o = OneBodyOper::gamma; o != OneBodyOper::invalid; ++o) {
    OneBodyOper::type otype = static_cast<OneBodyOper::type> (o);
    if (oper_key == OneBodyOper::to_string(otype)) {
      transform_operator = otype;
      break;
    }
  }
  if (transform_operator == OneBodyOper::invalid)
    return false;

  const std::string orig_key = pop_till_token(keycopy, ')');
  if (not keycopy.empty()) return false;
  ParsedOrbitalSpaceKey parsed_orig_key(orig_key);

  return true;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
