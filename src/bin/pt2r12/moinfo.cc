//
// moinfo.cc
//
// Copyright (C) 2011 Edward Valeev
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

#include <moinfo.h>
#include <iostream>
#include <numeric>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/split.h>

using namespace sc;

namespace {
  /// pop off str from beginning up to first apperance of separator sep
  /// the string plus ALL following instances of sep are erased from str
  std::string pop_till_sep(std::string& str, char sep) {
    while (!str.empty() && *(str.begin()) == sep) {
      str.erase(str.begin());
    }
    const size_t next_sep_pos = str.find_first_of(sep);
    std::string result;
    if (next_sep_pos != std::string::npos) {
      result = str.substr(0, next_sep_pos);
      str.erase(0, next_sep_pos + 1);
    } else {
      result = str;
      str.clear();
    }
    while (!str.empty() && *(str.begin()) == sep) {
      str.erase(str.begin());
    }
    return result;
  }

  /// converts str's characters to lowercase
  void tolower(std::string& str) {
    for(std::string::iterator i=str.begin();
        i!=str.end(); ++i)
      *i = std::tolower(*i);
  }

  /// splits str into a set of tokens, separated by one or more consecutive separator characters sep
  std::vector<std::string> split(std::string str,
                                 char sep) {
    std::vector<std::string> result;
    while (str.empty() == false) {
      std::string token = pop_till_sep(str, ' ');
      result.push_back(token);
    }
    return result;
  }

  /// parses str as a set of T
  template <typename T>
  std::vector<T> parse(std::string str) {
    std::istringstream iss(str);
    std::vector<T> result;
    while(!iss.eof()) {
      T t; iss >> t; result.push_back(t);
    }
    return result;
  }

  /// maps indices
  template <typename T, typename Index>
  std::vector<T> remap(const std::vector<T>& orig,
                              const std::vector<Index>& indexmap) {
    std::vector<T> result(orig.size());
    for(size_t i=0; i<orig.size(); ++i)
      result[indexmap[i]] = orig[i];
    return result;
  }

  /// read a line into std::string
  std::string readline(std::istream& in) {
    const size_t nline = 256;
    char linebuf[nline];
    in.getline(linebuf, nline);
    return linebuf;
  }
  /// skip EOL
  void skipeol(std::istream& in) {
    readline(in);
  }
}

ExternMOInfo::ExternMOInfo(const std::string & filename,
                           const Ref<Integral>& integral)
{
  std::ifstream in(filename.c_str());
  std::string strjunk;

  //////
  // parse molecular geometry
  //////
  Ref<Molecule> molecule = new Molecule;
  bool have_atoms = true;
  while (have_atoms) {
    double charge, x, y, z;
    in >> charge >> x >> y >> z;
    if (charge != -1.0) {
      //x /= 0.529177249;
      //y /= 0.529177249;
      //z /= 0.529177249;
      molecule->add_atom(charge, x, y, z);
    }
    else {
      have_atoms = false;
      // to use getline after now skip EOL
      const size_t nline = 256;
      char linebuf[nline];
      in.getline(linebuf, nline);
    }
  }

  //////
  // parse the point group info
  //////
  // read the Schoenflies symbol
  std::string pointgroup_symbol; in >> pointgroup_symbol;
  Ref<PointGroup> pg = new PointGroup(pointgroup_symbol.c_str());
  molecule->set_point_group(pg);
  molecule->print();
  // first map MPQC irrep label to irrep index
  const unsigned int nirrep = pg->order();
  std::map<std::string, unsigned int> irreplabel_to_index;
  for(unsigned int g=0; g<nirrep; ++g) {
    std::string glabel = pg->char_table().gamma(g).symbol();
    tolower(glabel);
    irreplabel_to_index[glabel] = g;
    std::cout << "irrep " << g << " " << glabel << std::endl;
  }
  // now get irrep labels and determine mapping from extern irrep index to mpqc irrep index
  std::vector<unsigned int> extern_to_mpqc_irrep_map;
  skipeol(in);
  std::string irrep_labels = readline(in);
  std::cout << "irrep_labels = " << irrep_labels << std::endl;
  if (irrep_labels != "-1") { // unless -1 follows the point group symbol, the irrep labels are given next
                              // map them to the MPQC labels
    while (irrep_labels.empty() == false) {
      std::string irreplabel = pop_till_sep(irrep_labels, ' ');
      tolower(irreplabel);
      extern_to_mpqc_irrep_map.push_back(irreplabel_to_index[irreplabel]);
      std::cout << "irreplabel = " << irreplabel << " (first char = '" << irreplabel[0] << "') irrep_labels = " << irrep_labels << std::endl;
    }
    assert(extern_to_mpqc_irrep_map.size() == nirrep);
    in >> strjunk;
    std::cout << strjunk << std::endl;
  } else { // irrep labels are not given? assume everything same as in MPQC
    for(int g=0; g<pg->order(); ++g)
      extern_to_mpqc_irrep_map.push_back(g);
  }

  //////
  // parse basis set
  //////
  Ref<GaussianBasisSet> basis;
  {
    Ref<AssignedKeyVal> tmpkv = new AssignedKeyVal;
    // add molecule
    tmpkv->assign("molecule", molecule.pointer());
    const unsigned int num_atoms = molecule->natom();
    bool have_gencon = false;

    // make name for atomic basis sets
    for (unsigned int a = 0; a < num_atoms; ++a) {
      std::ostringstream oss1, oss2;
      oss1 << "basis:" << a;
      oss2 << "custom" << a;
      tmpkv->assign(oss1.str().c_str(), oss2.str().c_str());
    }

    // read in first line
    std::string buffer;
    const size_t nline = 256;
    char linebuf[nline];
    in.getline(linebuf, nline);
    buffer = linebuf;

    // read shells for each atom
    int current_shell_index = -1;
    bool have_more_shells;
    bool have_more_prims;
    std::map<int, int> current_shell_index_on_atom;
    std::string shell_prefix;
    std::string am;
    int puream = 0;

    do { // shell loop
      int current_primitive_index = -1;
      do { // primitive loop

        const size_t first_not_whitespace = buffer.find_first_not_of(" ");
        const char first_char = buffer[first_not_whitespace];
        if (first_char == '-') { // done with the basis set specification
          have_more_shells = false;
          have_more_prims = false;
        } else { // have new shell or primitive? parse further
          const size_t first_whitespace = buffer.find_first_of(
              " ", first_not_whitespace);
          const size_t first_not_whitespace = buffer.find_first_not_of(
              " ", first_whitespace);
          const char first_char = buffer[first_not_whitespace];

          if (isalpha(first_char) == true) { // start of new shell
            have_more_shells = true;
            have_more_prims = false;
          } else { // have a primitive
            have_more_shells = false;
            have_more_prims = true;
            ++current_primitive_index;
          }
        }

        std::istringstream iss(buffer);

        // have new shell?
        if (have_more_shells) {

          int current_atom_index;
          iss >> current_atom_index >> am >> puream;
          tolower(am);
          --current_atom_index;
          if (current_shell_index_on_atom.find(current_atom_index)
              == current_shell_index_on_atom.end())
            current_shell_index_on_atom[current_atom_index] = -1;
          current_shell_index =
              ++current_shell_index_on_atom[current_atom_index];
          std::ostringstream oss;
          oss << ":basis:" << molecule->atom_name(current_atom_index)
              << ":custom" << current_atom_index << ":" << current_shell_index;
          shell_prefix = oss.str();

//          if (am == "SP") {
//            have_gencon = true;
//          } else if (am[0] != 'S' && am[0] != 'P') {
//            iss >> puream;
//          }

          tmpkv->assign((shell_prefix + ":type:0:am").c_str(),
                        std::string(1, am[0]));
          if (puream)
            tmpkv->assign((shell_prefix + ":type:0:puream").c_str(),
                          std::string("true"));
          if (am == "SP") {
            tmpkv->assign((shell_prefix + ":type:1:am").c_str(),
                          std::string(1, am[1]));
          }

        }

        // have new primitive?
        if (have_more_prims) {
          double exponent, coef0, coef1;
          iss >> exponent >> coef0;
          if (am == "SP") {
            have_gencon = true;
            iss >> coef1;
          }
          { // exponent
            std::ostringstream oss;
            oss << shell_prefix << ":exp:" << current_primitive_index;
            tmpkv->assign(oss.str().c_str(), exponent);
          }
          { // coef 0
            std::ostringstream oss;
            oss << shell_prefix << ":coef:0:" << current_primitive_index;
            tmpkv->assign(oss.str().c_str(), coef0);
          }
          if (am == "SP") { // coef 1
            std::ostringstream oss;
            oss << shell_prefix << ":coef:1:" << current_primitive_index;
            tmpkv->assign(oss.str().c_str(), coef1);
          }

        } // done with the current primitive

        // read next line
        if (have_more_prims || have_more_shells) {
          const size_t nline = 256;
          char linebuf[nline];
          in.getline(linebuf, nline);
          buffer = linebuf;
        }

      } while (have_more_prims); // loop over prims in this shell

    } while (have_more_shells); // loop over shells on this atom
    Ref<KeyVal> kv = tmpkv;
    basis = new GaussianBasisSet(kv);

    if (have_gencon) { // if have general contractions like SP shells in Pople-style basis sets split them so that IntegralLibint2 can handle them
      Ref<GaussianBasisSet> split_basis = new SplitBasisSet(basis);
      basis = split_basis;
    }
  }
  basis->print();

#if 0
  unsigned int nmo;  in >> nmo;
  unsigned int nao;  in >> nao;
  unsigned int ncore; in >> ncore;
  nfzc_ = ncore;
  in >> nfzv_;
  unsigned int nact; in >> nact;
  nocc_ = nact + ncore;
#endif

  ////////////////////////////////////////////////
  // read the number of MOs and the coefficients
  ////////////////////////////////////////////////
  unsigned int nmo;
  in >> nmo;

#if 0
  orbsym_.resize(nmo);
  { // compute orbital symmetries from irrep labels
    for(int o=0; o<nmo; ++o) {
  //    unsigned int irrep; in >> irrep;
  //    --irrep;
      std::string irreplabel;
      in >> irreplabel;
      // convert irrep_label to irrep index
      orbsym_[o] = irreplabel_to_index[irreplabel];
      std::cout << "mo " << o << " " << irreplabel << " " << orbsym_[o] << std::endl;
    }
  }
#endif

  // use properly blocked AO dimension
  Ref<Integral> localints = integral->clone();
  localints->set_basis(basis);
  RefSCDimension aodim = localints->petite_list()->AO_basisdim();

  // make a dummy MO dimension for now -- it will be recreated by OrderedOrbitalSpace
  RefSCDimension modim = new SCDimension(nmo, 1);
  modim->blocks()->set_subdim(0, new SCDimension(nmo));
  RefSCMatrix coefs_extern = basis->so_matrixkit()->matrix(aodim, modim);
  coefs_extern.assign(0.0);
  bool have_coefs = true;
  while(have_coefs) {
    int row, col;
    double value;
    in >> row >> col >> value;
    if (row != -1) {
      --row;  --col;
      coefs_extern.set_element(row,col,value);
    }
    else
      have_coefs = false;
  }

  ////////////////////////////////////////////
  // read the orbital info
  ////////////////////////////////////////////
  skipeol(in);
  std::string token = readline(in);
  std::vector<unsigned int> mopi = parse<unsigned int>(token);  assert(mopi.size() == pg->order());
  token = readline(in);
  fzcpi_ = parse<unsigned int>(token);  assert(fzcpi_.size() == pg->order());
  token = readline(in);
  fzvpi_ = parse<unsigned int>(token);  assert(fzvpi_.size() == pg->order());
  token = readline(in);
  inactpi_ = parse<unsigned int>(token);  assert(inactpi_.size() == pg->order());
  token = readline(in);
  actpi_ = parse<unsigned int>(token);  assert(actpi_.size() == pg->order());
  assert(std::accumulate(mopi.begin(), mopi.end(), 0) == nmo);
  unsigned int junk; in >> junk;

  ////////////////////////////////////////////
  // remap the coefficients and orbital info
  ////////////////////////////////////////////
  // pseudo occupation numbers -- frozen core and inactive are 2.0, active are 1.0, the rest are 0.0
  RefDiagSCMatrix pseudo_occnums = coefs_extern.kit()->diagmatrix(coefs_extern.coldim());
  pseudo_occnums.assign(0.0);
  unsigned int mo = 0;
  for(unsigned int irrep=0; irrep<pg->order(); ++irrep) {
    unsigned int mo_in_irrep = 0;
    for(unsigned i=0; i<fzcpi_[irrep]; ++i, ++mo, ++mo_in_irrep)
      pseudo_occnums.set_element(i, 2.0);
    for(unsigned i=0; i<inactpi_[irrep]; ++i, ++mo, ++mo_in_irrep)
      pseudo_occnums.set_element(i, 2.0);
    for(unsigned i=0; i<actpi_[irrep]; ++i, ++mo, ++mo_in_irrep)
      pseudo_occnums.set_element(i, 1.0);
    mo += mopi[irrep] - mo_in_irrep; // skip the rest of the orbitals in this block
  }

  // pseudo eigenvalues are the occupation numbers with changed sign
  RefDiagSCMatrix pseudo_evals = pseudo_occnums.copy();
  pseudo_evals.scale(-1.0);

  // compute map from MO to the mpqc irrep
  std::vector<unsigned int> orbsym;
  {
    for (unsigned int g = 0; g < mopi.size(); ++g) {
      const unsigned int g_mpqc = extern_to_mpqc_irrep_map[g];
      for (unsigned int mo_g = 0; mo_g < mopi[g]; ++mo_g)
        orbsym.push_back(g_mpqc);
    }
  }

  orbs_ = new OrdOrbitalSpace(
      std::string("p"), std::string("MOInfo orbitals"),
      basis, integral, coefs_extern, pseudo_evals,
      pseudo_occnums, orbsym, SymmetryMOOrder(pg->order()) );

  // remap all orbitals to MPQC irreps
  fzcpi_ = remap(fzcpi_, extern_to_mpqc_irrep_map);
  inactpi_ = remap(inactpi_, extern_to_mpqc_irrep_map);
  actpi_ = remap(actpi_, extern_to_mpqc_irrep_map);
  fzvpi_ = remap(fzvpi_, extern_to_mpqc_irrep_map);

  //////////////////////////////////////////////////////////////////////////////////
  // lastly determine the index maps from extern MO indices to the MPQC MO indices
  //////////////////////////////////////////////////////////////////////////////////
  {
    assert(indexmap_.empty());
    // use blockinfo to get MO offets for orbitals in MPQC order
    Ref<SCBlockInfo> blkinfo = orbs_->coefs()->coldim()->blocks();
    for (unsigned int g = 0; g < pg->order(); ++g) {
      const unsigned int g_mpqc = extern_to_mpqc_irrep_map[g];
      const unsigned int mpqc_mo_offset = blkinfo->start(g_mpqc);
      const unsigned int nmo_g = mopi[g];
      for (unsigned int mo_g = 0; mo_g < nmo_g; ++mo_g)
        indexmap_.push_back( mpqc_mo_offset + mo_g );
    }
  }
  {
    assert(occindexmap_.empty());
    // use blockinfo to get MO offets for orbitals in MPQC order
    Ref<SCBlockInfo> blkinfo = orbs_->coefs()->coldim()->blocks();
    for (unsigned int g = 0; g < pg->order(); ++g) {
      const unsigned int g_mpqc = extern_to_mpqc_irrep_map[g];
      const unsigned int mpqc_mo_offset = blkinfo->start(g_mpqc);
      const unsigned int nmo_g = fzcpi_[g] + inactpi_[g] + actpi_[g];
      for (unsigned int mo_g = 0; mo_g < nmo_g; ++mo_g)
        occindexmap_.push_back( mpqc_mo_offset + mo_g );
    }
  }

  coefs_extern.print("ExternMOInfo:: extern MO coefficients");
  orbs_->coefs().print("ExternMOInfo:: reordered MO coefficients");

  in.close();
}

const std::vector<unsigned int>& ExternMOInfo::indexmap() const
{
  return indexmap_;
}

const std::vector<unsigned int>& ExternMOInfo::occindexmap() const
{
  return occindexmap_;
}

//Ref<GaussianBasisSet> ExternMOInfo::basis() const
//{
//    return basis_;
//}
//
//RefSCMatrix ExternMOInfo::coefs() const
//{
//    return coefs_;
//}
//
//unsigned int ExternMOInfo::nfzc() const
//{
//  // return nfzc_;
//  return 0;
//}
//
//unsigned int ExternMOInfo::nfzv() const
//{
//  //return nfzv_;
//  return 0;
//}
//
//unsigned int ExternMOInfo::nocc() const
//{
//  //return nocc_;
//  return 7;
//}

//const std::vector<unsigned int>& ExternMOInfo::mopi() const
//{
//  return mopi_;
//}
//
const std::vector<unsigned int>& ExternMOInfo::fzcpi() const
{
  return fzcpi_;
}

const std::vector<unsigned int>& ExternMOInfo::fzvpi() const
{
  return fzvpi_;
}

const std::vector<unsigned int>& ExternMOInfo::actpi() const
{
  return actpi_;
}

const std::vector<unsigned int>& ExternMOInfo::inactpi() const
{
  return inactpi_;
}

//std::vector<unsigned int> ExternMOInfo::orbsym() const
//{
//  return orbsym;
//}
//
/////////////////////////////////////////////////////////////////////////////

ClassDesc ExternSpinFreeRDMOne::class_desc_(typeid(ExternSpinFreeRDMOne),
                                        "ExternSpinFreeRDMOne",
                                        1,               // version
                                        "public SpinFreeRDM<One>", // must match parent
                                        0,               // not DefaultConstructible
                                        0,               // not KeyValConstructible
                                        0                // not StateInConstructible
                                        );

ExternSpinFreeRDMOne::ExternSpinFreeRDMOne(const std::string & filename,
                                           const std::vector<unsigned int>& indexmap,
                                           const Ref<OrbitalSpace>& orbs) :
  SpinFreeRDM<One>(Ref<Wavefunction>()), orbs_(orbs)
{
  std::ifstream in(filename.c_str());
  if (in.is_open() == false) {
    std::ostringstream oss;
    oss << "ExternSpinFreeRDMOne: could not open file " << filename;
    throw std::runtime_error(oss.str().c_str());
  }
  rdm_ = orbs->coefs().kit()->symmmatrix(orbs->coefs().coldim()); rdm_.assign(0.0);
  bool have_coefs = true;
  while(have_coefs) {
    int row, col;
    double value;
    in >> row >> col >> value;
    if (row != -1) {
      --row;  --col;
      const unsigned int mbra = indexmap[row];
      const unsigned int mket = indexmap[col];
      rdm_.set_element(mbra, mket, value);
    }
    else
      have_coefs = false;
  }
  in.close();

  rdm_.print("ExternSpinFreeRDMOne:: MO density");
}

ExternSpinFreeRDMOne::ExternSpinFreeRDMOne(const RefSymmSCMatrix& rdm, const Ref<OrbitalSpace>& orbs) :
  SpinFreeRDM<One>(Ref<Wavefunction>()), rdm_(rdm), orbs_(orbs)
{
  rdm_.print("ExternSpinFreeRDMOne:: MO density");
}

ExternSpinFreeRDMOne::~ExternSpinFreeRDMOne()
{
}

/////////////////////////////////////////////////////////////////////////////

ClassDesc ExternSpinFreeRDMTwo::class_desc_(typeid(ExternSpinFreeRDMTwo),
                                        "ExternSpinFreeRDMTwo",
                                        1,               // version
                                        "public SpinFreeRDM<Two>", // must match parent
                                        0,               // not DefaultConstructible
                                        0,               // not KeyValConstructible
                                        0                // not StateInConstructible
                                        );

sc::ExternSpinFreeRDMTwo::ExternSpinFreeRDMTwo(const std::string & filename,
                                               const std::vector<unsigned int>& indexmap,
                                               const Ref<OrbitalSpace> & orbs) :
    SpinFreeRDM<Two>(Ref<Wavefunction>()), filename_(filename), orbs_(orbs)
{
  std::ifstream in(filename_.c_str());
  if (in.is_open() == false) {
    std::ostringstream oss;
    oss << "ExternSpinFreeRDMTwo: could not open file " << filename_;
    throw std::runtime_error(oss.str().c_str());
  }
  const unsigned int norbs = orbs->coefs().coldim().n();
  RefSCDimension dim = new SCDimension(norbs * norbs);
  dim->blocks()->set_subdim(0, new SCDimension(dim->n()));
  rdm_ = orbs->coefs().kit()->symmmatrix(dim);
  rdm_.assign(0.0);
  bool have_coefs = true;
  while (have_coefs) {
    int bra1, bra2, ket1, ket2;
    double value;
    in >> bra1 >> bra2 >> ket1 >> ket2 >> value;
    // these index orderings are hard to distinguish! Both give the same energy
    // can only distinguish if ALL elements are given, then for the above ordering
    // element 1 2 3 4 will be equivalent to 3 4 1 2, 4 3 2 1, and 2 1 4 3
    // for the ordering below same element will be equivalent to
    // 3 4 1 2, 4 3 2 1, and 2 3 4 1 !!!
    // in >> bra1 >> ket2 >> ket1 >> bra2 >> value;
    if (bra1 != -1) {
      --bra1;
      --bra2;
      --ket1;
      --ket2;
      
      const unsigned int mbra1 = indexmap[bra1];
      const unsigned int mbra2 = indexmap[bra2];
      const unsigned int mket1 = indexmap[ket1];
      const unsigned int mket2 = indexmap[ket2];

      rdm_.set_element(mbra1 * norbs + mbra2, mket1 * norbs + mket2, value);
      rdm_.set_element(mbra2 * norbs + mbra1, mket2 * norbs + mket1, value);
    } else
      have_coefs = false;
  }
  in.close();

  //rdm_.print("ExternSpinFreeRDMTwo:: MO density");
}

ExternSpinFreeRDMTwo::~ExternSpinFreeRDMTwo()
{
}

Ref< SpinFreeRDM<One> >
ExternSpinFreeRDMTwo::rdm_m_1() const
{
  RefSymmSCMatrix rdm1 = orbs_->coefs().kit()->symmmatrix(orbs_->coefs().coldim());
  rdm1.assign(0.0);
  const unsigned int norbs = rdm1.n();
  for (unsigned int b1 = 0; b1 < norbs; ++b1) {
    const unsigned b12_offset = b1 * norbs;
    for (unsigned int k1 = 0; k1 <= b1; ++k1) {
      const unsigned k12_offset = k1 * norbs;
      double value = 0.0;
      for (unsigned int i2 = 0; i2 < norbs; ++i2) {
        value += rdm_.get_element(b12_offset + i2, k12_offset + i2);
      }
      rdm1.set_element(b1, k1, value);
    }
  }

  // compute the number of electrons and scale 1-rdm:
  // trace of the 2-rdm is n(n-1)
  // trace of the 1-rdm should be n
  const double trace_2rdm = rdm1.trace();
  const double nelectron = (1.0 + std::sqrt(1.0 + 4.0 * trace_2rdm)) / 2.0;
  rdm1.scale(1.0 / (nelectron - 1.0));

  Ref < SpinFreeRDM<One> > result = new ExternSpinFreeRDMOne(rdm1, orbs_);
  return result;
}





/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
