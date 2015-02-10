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

#include <extern/moinfo/moinfo.h>
#include <iostream>
#include <numeric>
#include <cassert>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/split.h>
#include <chemistry/qc/basis/cart.h>
#include <chemistry/qc/lcao/fockbuilder.h>

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
    const size_t nline = 2<<15; // 2^16 = 64 kiB
    char linebuf[nline];
    in.getline(linebuf, nline);
    return linebuf;
  }
  /// skip EOL
  void skipeol(std::istream& in) {
    readline(in);
  }
}

////////////////////////////////
ClassDesc
ExternMOInfo::class_desc_(typeid(ExternMOInfo),
                     "ExternMOInfo",
                     1,               // version
                     "public DescribedClass",
                     0, // this class is not DefaultConstructible
                     0, // this class is not KeyValConstructible
                     0  // this class is not StateInConstructible
                     );

ExternMOInfo::ExternMOInfo(std::string filename,
                           Ref<Integral> integral,
                           std::string basisname)
{
  std::ifstream in(filename.c_str());
  if(!in.is_open())
  {
    throw InputError((std::string("Failed to open data file: ") + filename).c_str(), __FILE__, __LINE__);
  }
  std::string strjunk;

  //////
  // parse molecular geometry
  //////
  Ref<Molecule> molecule = new Molecule;
  bool have_atoms = true;
  while (have_atoms)
  {
    double charge, x, y, z;
    in >> charge >> x >> y >> z;
    if (charge != -1.0) {
      // geometries must be in atomic units
      molecule->add_atom(static_cast<int>(charge), x, y, z);
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
  // coordinates may be given in non-canonical frame for the given point group
  // detect the highest Abelian point group (hopefully same as given) and use its frame
  Ref<PointGroup> pg;
  {
    Ref<PointGroup> highestpg = molecule->highest_point_group();
    pg = new PointGroup(pointgroup_symbol.c_str(), highestpg->symm_frame(), highestpg->origin());
    try {
      molecule->set_point_group(pg);
    }
    catch(AlgorithmException& ex) {
      try {
        ex.elaborate() << "in ExternMOInfo ctor: could not detect point group " << pointgroup_symbol;
      }
      catch(...) {}
      throw ex;
    }
  molecule->print();
  ExEnv::out0() << indent << "nuclear repulsion energy = "
                << scprintf("%25.15f",molecule->nuclear_repulsion_energy()) << std::endl;
  // first map MPQC irrep label to irrep index
  const unsigned int nirrep = pg->order();
  std::map<std::string, unsigned int> irreplabel_to_index;
  for(unsigned int g=0; g<nirrep; ++g) {
    std::string glabel = pg->char_table().gamma(g).symbol();
    tolower(glabel);
    irreplabel_to_index[glabel] = g;
    //std::cout << "irrep " << g << " " << glabel << std::endl;
  }
  // some peculiar programs ... cough MOLCAS cough ... may use non-standard (i.e. non-Cotton) irrep orderings
  // if so, pick up extern irrep labels and determine mapping from extern irrep index to mpqc irrep index
  std::vector<unsigned int> extern_to_mpqc_irrep_map;
  skipeol(in);
  std::string irrep_labels = readline(in);
  //std::cout << "irrep_labels = " << irrep_labels << std::endl;
  if (irrep_labels != "-1") { // unless -1 follows the point group symbol, the irrep labels are given next
                              // map them to the MPQC labels
    while (irrep_labels.empty() == false) {
      std::string irreplabel = pop_till_sep(irrep_labels, ' ');
      tolower(irreplabel);
      extern_to_mpqc_irrep_map.push_back(irreplabel_to_index[irreplabel]);
      //std::cout << "irreplabel = " << irreplabel << " (first char = '" << irreplabel[0] << "') irrep_labels = " << irrep_labels << std::endl;
    }
    MPQC_ASSERT(extern_to_mpqc_irrep_map.size() == nirrep);
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
    // optional descriptive label may be also given for atomic basis sets
    if (!basisname.empty()) {
      tmpkv->assign("name", basisname.c_str());
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
          const bool first_char_is_a_letter = (first_char >= 'a' && first_char <= 'z') ||
                                              (first_char >= 'A' && first_char <= 'Z');

          if (first_char_is_a_letter) { // start of new shell
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
          puream = 0; // if puream not given, make the basis function cartesian
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

          tmpkv->assign((shell_prefix + ":type:0:am").c_str(),
                        std::string(1, am[0]));
          if (puream)
            tmpkv->assign((shell_prefix + ":type:0:puream").c_str(),
                          std::string("true"));
          if (am == "sp") {
            tmpkv->assign((shell_prefix + ":type:1:am").c_str(),
                          std::string(1, am[1]));
          }

        }

        // have new primitive?
        if (have_more_prims) {
          double exponent, coef0, coef1;
          iss >> exponent >> coef0;
          if (am == "sp") {
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
          if (am == "sp") { // coef 1
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
      Ref<GaussianBasisSet> split_basis = new SplitBasisSet(basis, basisname);
      basis = split_basis;
    }
  }
  basis->print();

  ////////////////////////////////////////////////
  // read the number of MOs and the coefficients
  ////////////////////////////////////////////////
  unsigned int nmo;
  in >> nmo;

  // use properly blocked AO dimension
  Ref<Integral> localints = integral->clone();
  localints->set_basis(basis);
  RefSCDimension aodim = localints->petite_list()->AO_basisdim();

  // make a dummy MO dimension for now -- it will be recreated by OrderedOrbitalSpace
  RefSCDimension modim = new SCDimension(nmo, 1);
  modim->blocks()->set_subdim(0, new SCDimension(nmo));

  // read coefficients
  RefSCMatrix coefs_extern = basis->so_matrixkit()->matrix(aodim, modim);
  // some wicked programs ... cough GAMESS cough ... report MO coefficients in *cartesian* basis always
#ifdef PT2R12GAMESS
  Ref<GaussianBasisSet> cbasis;
  if (basis->has_pure()) {
    // construct CartesianBasisSet and its AO dimension
    integral->set_basis(basis);
    cbasis = new CartesianBasisSet(basis, integral);
    localints->set_basis(cbasis);
    RefSCDimension aodim_cart = localints->petite_list()->AO_basisdim();
    coefs_extern = basis->so_matrixkit()->matrix(aodim_cart, modim);
  }
#endif
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

  // some wicked programs ... cough GAMESS cough ... report MO coefficients in *cartesian* basis always
#ifdef PT2R12GAMESS
  // if yes, may need to transform to the "real" basis?
  if (basis->has_pure()) {
    // now compute overlap between spherical and cartesian basis
    Ref<Integral> localints = integral->clone();
    localints->set_basis(basis,cbasis);
    RefSCMatrix S_sph_cart = detail::onebodyint_ao<&Integral::overlap>(basis, cbasis, localints);

    // SO basis is always blocked, so first make sure S_sph_cart is blocked
    {
      Ref<Integral> braints = integral->clone();  braints->set_basis(basis);
      Ref<PetiteList> brapl = braints->petite_list();
      Ref<Integral> ketints = integral->clone();  ketints->set_basis(cbasis);
      Ref<PetiteList> ketpl = ketints->petite_list();

      RefSCMatrix S_sph_cart_blk = basis->so_matrixkit()->matrix(brapl->AO_basisdim(),ketpl->AO_basisdim());
      S_sph_cart_blk->convert(S_sph_cart);
      S_sph_cart = S_sph_cart_blk;
    }

    // compute projector from cart to spherical basis
    // P(cart->sph) = S^-1 (sph/sph) * S (sph/cart)
    RefSymmSCMatrix S_sph_sph;
    {
      localints->set_basis(basis,basis);
      S_sph_sph = detail::onebodyint_ao<&Integral::overlap>(basis, localints);
    }
    Ref<OverlapOrthog> orthog = new OverlapOrthog(OverlapOrthog::Symmetric, S_sph_sph,
                                                  S_sph_sph.kit(), 1e-8, 0);
    coefs_extern = orthog->overlap_inverse() * S_sph_cart * coefs_extern;
  }
#endif

  ////////////////////////////////////////////
  // read the orbital info
  ////////////////////////////////////////////
#if PT2R12GAMESS // GAMESS reports orbitals in energy order

  skipeol(in);
  std::string token = readline(in);
  std::vector<unsigned int> orbsym = parse<unsigned int>(token);  MPQC_ASSERT(orbsym.size() == nmo);
  // map irrep indices to MPQC order; N.B. for MOs start from 1, hence subtraction of 1
  for(std::vector<unsigned int>::iterator i=orbsym.begin();
      i!=orbsym.end();
      ++i)
    *i = extern_to_mpqc_irrep_map[*i - 1];

  unsigned int nfzc;    in >> nfzc;
  unsigned int ninact;  in >> ninact;
  unsigned int nact;    in >> nact;
  unsigned int ncorr;    in >> ncorr;
  unsigned int nfzv;    in >> nfzv;
  const unsigned int nuocc = nmo - nfzc - ninact - nact - nfzv;
  mopi_.resize(pg->order());     std::fill(mopi_.begin(), mopi_.end(), 0u);
  fzcpi_.resize(pg->order());    std::fill(fzcpi_.begin(), fzcpi_.end(), 0u);
  inactpi_.resize(pg->order());  std::fill(inactpi_.begin(), inactpi_.end(), 0u);
  actpi_.resize(pg->order());    std::fill(actpi_.begin(), actpi_.end(), 0u);
  corrpi_.resize(pg->order());   std::fill(corrpi_.begin(), corrpi_.end(), 0u);
  fzvpi_.resize(pg->order());    std::fill(fzvpi_.begin(), fzvpi_.end(), 0u);

  RefDiagSCMatrix pseudo_occnums = coefs_extern.kit()->diagmatrix(coefs_extern.coldim());
  pseudo_occnums.assign(0.0);

  unsigned int mo = 0;
  for(unsigned int i=0; i<nfzc;   ++i, ++mo)  { ++fzcpi_[orbsym[mo]];   pseudo_occnums.set_element(mo, 2.0); }
  // inactive orbitals (i.e. doubly occupied) are correlated
  for(unsigned int i=0; i<ninact; ++i, ++mo)  { ++inactpi_[orbsym[mo]]; ++corrpi_[orbsym[mo]]; pseudo_occnums.set_element(mo, 2.0); }
  for(unsigned int i=0; i<nact;   ++i, ++mo)
  {
    ++actpi_[orbsym[mo]];
    pseudo_occnums.set_element(mo, 1.0);
    if(i<ncorr) ++corrpi_[orbsym[mo]];
  }
  mo += nuocc;
  for(unsigned int i=0; i<nfzv;   ++i, ++mo)  { ++fzvpi_[orbsym[mo]];                                        }
  for(unsigned int i=0; i<nmo;   ++i, ++mo)   { ++mopi_[orbsym[i]];                                          }

#else // the default is to report orbitals in symmetry order

  skipeol(in);
  std::string token;
  token = readline(in);  mopi_    = parse<unsigned int>(token); MPQC_ASSERT(mopi_.size() == pg->order());
  token = readline(in);  fzcpi_   = parse<unsigned int>(token); MPQC_ASSERT(fzcpi_.size() == pg->order());
  token = readline(in);  inactpi_ = parse<unsigned int>(token); MPQC_ASSERT(inactpi_.size() == pg->order());
  token = readline(in);  actpi_   = parse<unsigned int>(token); MPQC_ASSERT(actpi_.size() == pg->order());
  token = readline(in);  fzvpi_   = parse<unsigned int>(token); MPQC_ASSERT(fzvpi_.size() == pg->order());
  MPQC_ASSERT(std::accumulate(mopi_.begin(), mopi_.end(), 0) == nmo);
  unsigned int junk; in >> junk; // fzcpi_ etc are in molcas symmetry order

  // by default correlate all inactive and active orbitals
  corrpi_ = actpi_;
  for(size_t i=0; i<corrpi_.size(); ++i)
    corrpi_[i] += inactpi_[i];

  const unsigned int nfzc   = std::accumulate(fzcpi_.begin(),   fzcpi_.end(),   0u);
  const unsigned int ninact = std::accumulate(inactpi_.begin(), inactpi_.end(), 0u);
  const unsigned int nact   = std::accumulate(actpi_.begin(),   actpi_.end(),   0u);
  const unsigned int ncorr   = std::accumulate(corrpi_.begin(),   corrpi_.end(),   0u);
  const unsigned int nfzv   = std::accumulate(fzvpi_.begin(),   fzvpi_.end(),   0u);
  const unsigned int nuocc = nmo - nfzc - ninact - nact - nfzv;

  // compute MPQC irrep indices of the extern MOs
  std::vector<unsigned int> orbsym;
  {
    for (unsigned int g = 0; g < mopi_.size(); ++g) {
      const unsigned int g_mpqc = extern_to_mpqc_irrep_map[g];
      for (unsigned int mo_g = 0; mo_g < mopi_[g]; ++mo_g)
        orbsym.push_back(g_mpqc);
    }
  }

  // pseudo occupation numbers -- frozen core and inactive are 2.0, active are 1.0, the rest are 0.0
  RefDiagSCMatrix pseudo_occnums = coefs_extern.kit()->diagmatrix(coefs_extern.coldim());
  pseudo_occnums.assign(0.0);
  unsigned int mo = 0;
  for(unsigned int irrep=0; irrep<pg->order(); ++irrep) {
    unsigned int mo_in_irrep = 0;
    for(unsigned i=0; i<fzcpi_[irrep]; ++i,++mo_in_irrep)
      pseudo_occnums.set_element(mo+i, 2.0);
    mo += fzcpi_[irrep];
    for(unsigned i=0; i<inactpi_[irrep]; ++i, ++mo_in_irrep)
      pseudo_occnums.set_element(mo + i, 2.0);
    mo += inactpi_[irrep];
    for(unsigned i=0; i<actpi_[irrep]; ++i, ++mo_in_irrep)
      pseudo_occnums.set_element(mo +i, 1.0);
    mo += actpi_[irrep];
    mo += mopi_[irrep] - mo_in_irrep; // skip the rest of the orbitals in this block
  }

#endif

  in.close();

  // pseudo eigenvalues are the occupation numbers with changed sign; this is used to sort orbs
  RefDiagSCMatrix pseudo_evals = pseudo_occnums.copy();
  pseudo_evals.scale(-1.0);

  ////////////////////////////////////////////
  // remap the coefficients and orbital info
  ////////////////////////////////////////////

  // sym->E->occ order
  orbs_sb_ = new SymmOrbitalSpace(
      std::string("p"), std::string("symmetry-ordered MOInfo orbitals"),
      basis, integral, coefs_extern, pseudo_evals,
      pseudo_occnums, orbsym, SymmetryMOOrder(pg->order()) );
  // occ->sym->E
  orbs_ = new CorrOrbitalSpace(
      std::string("p~"), std::string("energy-ordered MOInfo orbitals"),
      basis, integral, coefs_extern, pseudo_evals,
      pseudo_occnums, orbsym, EnergyMOOrder<std::less<double> >() );

  // right now irrep->#oforbs arrays still use extern irrep ordering
  // remap irrep->#oforbs arrays to MPQC irrep ordering
  fzcpi_ = remap(fzcpi_, extern_to_mpqc_irrep_map);
  inactpi_ = remap(inactpi_, extern_to_mpqc_irrep_map);
  actpi_ = remap(actpi_, extern_to_mpqc_irrep_map);
  corrpi_ = remap(corrpi_, extern_to_mpqc_irrep_map);
  fzvpi_ = remap(fzvpi_, extern_to_mpqc_irrep_map);
  mopi_ = remap(mopi_, extern_to_mpqc_irrep_map);

  //////////////////////////////////////////////////////////////////////////////////
  // lastly determine the index maps from extern MO indices to the MPQC MO indices
  //////////////////////////////////////////////////////////////////////////////////

  // compute the maps using standard OrbitalSpace operators
  // first, compute extern spaces: full, active, and occupied
  // the goal is to keep the ordering, hence we'll just assign all orbitals to irrep 0
  std::vector<unsigned int> dummy_orbsym(nmo, 0u);
  Ref<OrbitalSpace> orbs_extern = new CorrOrbitalSpace(
      std::string("p(ext)"), std::string("extern MOInfo orbitals"),
      basis, integral, coefs_extern, pseudo_evals,
      pseudo_occnums, dummy_orbsym, EnergyMOOrder<std::less<double> >() );
  Ref<OrbitalSpace> orbs_occ_extern = new OrbitalSpace(std::string("i(ext)"), std::string("extern occupied MOInfo orbitals"),
                                                       orbs_extern->coefs(),
                                                       orbs_extern->basis(),
                                                       orbs_extern->integral(), orbs_extern->evals(),
                                                       0, nuocc + nfzv);
  Ref<OrbitalSpace> orbs_act_extern = new OrbitalSpace(std::string("m(ext)"),
                                                       std::string("extern active MOInfo orbitals"),
                                                       orbs_extern->coefs(),
                                                       orbs_extern->basis(),
                                                       orbs_extern->integral(), orbs_extern->evals(),
                                                       nfzc + ninact, nuocc + nfzv);

  // second, compute MPQC spaces: active, and occupied
  Ref<OrbitalSpace> orbs_sb_occ = new OrbitalSpace(std::string("i"),
                                                   std::string("symmetry-ordered occupied MOInfo orbitals"),
                                                   orbs_sb_->coefs(),
                                                   orbs_sb_->basis(),
                                                   orbs_sb_->integral(), orbs_sb_->evals(),
                                                   0, nuocc + nfzv,
                                                   OrbitalSpace::symmetry);

  indexmap_ =    (*orbs_sb_) << (*orbs_extern);
  actindexmap_occ_ = (*orbs_sb_occ) << (*orbs_act_extern);
  occindexmap_occ_ = (*orbs_sb_occ) << (*orbs_occ_extern);

  if (0) {
    coefs_extern.print("ExternMOInfo:: extern MO coefficients");
    orbs_sb_->coefs().print("ExternMOInfo:: reordered MO coefficients");
  }
}

const std::vector<unsigned int>& ExternMOInfo::indexmap() const
{
  return indexmap_;
}

const std::vector<unsigned int>& ExternMOInfo::occindexmap_occ() const
{
  return occindexmap_occ_;
}

const std::vector<unsigned int>& ExternMOInfo::actindexmap_occ() const
{
  return actindexmap_occ_;
}

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
const std::vector<unsigned int>& ExternMOInfo::corrpi() const
{
  return corrpi_;
}
const std::vector<unsigned int>& ExternMOInfo::inactpi() const
{
  return inactpi_;
}
const std::vector<unsigned int>& ExternMOInfo::mopi() const
{
  return mopi_;
}

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
  scmat_ = orbs->coefs().kit()->symmmatrix(orbs->coefs().coldim()); scmat_.assign(0.0);
  bool have_coefs = true;
  while(have_coefs) {
    int row, col;
    double value;
    in >> row >> col >> value;
    if (row != -1) {
      --row;  --col;
      const unsigned int mbra = indexmap[row];
      const unsigned int mket = indexmap[col];
      scmat_.set_element(mbra, mket, value);
    }
    else
      have_coefs = false;
  }
  in.close();

  //scmat_.print("ExternSpinFreeRDMOne:: MO density");
}

ExternSpinFreeRDMOne::ExternSpinFreeRDMOne(const RefSymmSCMatrix& rdm, const Ref<OrbitalSpace>& orbs) :
  SpinFreeRDM<One>(Ref<Wavefunction>()), orbs_(orbs)
{
  scmat_ = rdm;
  //scmat_.print("ExternSpinFreeRDMOne:: MO density");
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

ExternSpinFreeRDMTwo::ExternSpinFreeRDMTwo(const std::string & filename,
                                           const std::vector<unsigned int>& indexmap,
                                           const Ref<OrbitalSpace> & occ_orbs):
    SpinFreeRDM<Two>(Ref<Wavefunction>()), filename_(filename), orbs_(occ_orbs)
{
  // if the file reports all occupied indices
  if (occ_orbs->rank() == indexmap.size())
    this->init_from_rdm2_occspace(indexmap,occ_orbs);
  else // more complicated -- only active density is reported, need special handling
    this->init_from_rdm2_actspace(indexmap,occ_orbs);
}

void
ExternSpinFreeRDMTwo::init_from_rdm2_actspace(const std::vector<unsigned int>& indexmap,
                                              const Ref<OrbitalSpace> & occ_orbs)
{
  // goal: build 2-rdm in occupied space
  // the logic: build 1-rdm prototype -> build active space 2-rdm from reading file
  // -> build active space 1-rdm (from active space 2-rdm)
  // -> finish building 1-rdm in occupied space
  // -> build occupied space 2-rdm prototype by using 1-rdm
  // -> rectify the active section (all 4 indices are active) of 2-rdm by active space 2-rdm
  std::ifstream in(filename_.c_str());
  if (in.is_open() == false) {
    std::ostringstream oss;
    oss << "ExternSpinFreeRDMTwo: could not open file " << filename_;
    throw std::runtime_error(oss.str().c_str());
  }
  const unsigned int nocc = occ_orbs->coefs().coldim().n();
  const unsigned int nact = indexmap.size(); // # of active orbitals

  RefSCDimension dim = new SCDimension(nocc * nocc);
  RefSCDimension dim1 = new SCDimension(nocc);
  dim->blocks()->set_subdim(0, new SCDimension(dim->n()));
  dim1->blocks()->set_subdim(0, new SCDimension(dim1->n()));
  scmat_ = occ_orbs->coefs().kit()->symmmatrix(dim);
  scmat_.assign(0.0);
  RefSymmSCMatrix  rdm1 = occ_orbs->coefs().kit()->symmmatrix(dim1); // to store 1-rdm in obs
  rdm1.assign(0.0);

  // first fill "spin-free 1-rdm" with a diagonal mat ('2': alpha + beta)
  // the active block will be overwritten later.
  for (int i = 0; i < nocc; ++i)
  {
    rdm1.set_element(i,i, 2.0);
  }

  // now we construct 2-rdm in 'active' space.
  RefSCDimension act_dim2 = new SCDimension(nact * nact);
  act_dim2->blocks()->set_subdim(0, new SCDimension(act_dim2->n()));
  RefSymmSCMatrix ext_act_rdm2 = occ_orbs->coefs().kit()->symmmatrix(act_dim2);
  RefSCDimension act_dim1 = new SCDimension(nact);
  act_dim1->blocks()->set_subdim(0, new SCDimension(act_dim1->n()));
  RefSymmSCMatrix ext_act_rdm1 = occ_orbs->coefs().kit()->symmmatrix(act_dim1);
  ext_act_rdm2.assign(0.0);
  ext_act_rdm1.assign(0.0);
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
      ext_act_rdm2.set_element(bra1 * nact + bra2, ket1 * nact + ket2, value);
      ext_act_rdm2.set_element(bra2 * nact + bra1, ket2 * nact + ket1, value);
    } else
      have_coefs = false;
  }
  in.close();

  const double trace_act2 = ext_act_rdm2.trace(); // = n_act (n_act -1)
  const double nact_particle = (1.0 + std::sqrt(1.0 + 4.0 * trace_act2)) / 2.0;
  if(nact_particle == 1) // refer to the line "ext_act_rdm1.scale(1.0 /(nact_particle-1.0))"; basically for nact_particle ==1
                         // 2-rdm elements are all zero, so we can't construct 1-rdm from 2-rdm. This exception is practially
                         // useless, so let's just throw exception for now.
    throw AlgorithmException("The case of one-active electron is not dealt with properly: 1-rdm can not be constructed");

  // build active space 1-rdm from partial trace of 2-rdm in active space.
  for (unsigned int b1 = 0; b1 < nact; ++b1) {
    const unsigned b12_offset = b1 * nact;
    for (unsigned int k1 = 0; k1 <= b1; ++k1) {
      const unsigned k12_offset = k1 * nact;
      double value = 0.0;
      for (unsigned int i2 = 0; i2 < nact; ++i2) {
        value += ext_act_rdm2.get_element(b12_offset + i2, k12_offset + i2);
      }
      ext_act_rdm1.set_element(b1, k1, value);
    }
  }
  ext_act_rdm1.scale(1.0 /(nact_particle-1.0));

  // overwrite the active block of orignial 1-rdm in occupied space
  for (int i = 0; i < nact; ++i)
  {
    for (int j = 0; j < nact; ++j)
    {
      rdm1.set_element(indexmap[i], indexmap[j], ext_act_rdm1.get_element(i,j));
    }
  } // finish constructing rdm1 mat. It will be used to build rdm2


  // To build 2-rdm in occupied space: the trick is that unless all 4 indices are active,
  // any 2-rdm element Gamma^pq_rs can be decomposed as Gamma^p_r * Gamma^q_s - 0.5*Gamma^p_s * Gamma^q_r.
  // We use this to construct the 2-rdm and then overwrite the 'all active' section with act_rdm2
  for (int p = 0; p < nocc; ++p)
  {
    for (int q = 0; q < nocc; ++q)
    {
      const unsigned int uppind = p*nocc + q;
      for (int r = 0; r < nocc; ++r)
      {
        for (int s = 0; s < nocc; ++s)
        {
          const unsigned int lowind = r*nocc + s;
          if(uppind>=lowind)
            scmat_.set_element(uppind,lowind, rdm1.get_element(p,r)*rdm1.get_element(q,s)
                                            -0.5*rdm1.get_element(p,s)*rdm1.get_element(q,r));
        }
      }
    }
  }
  for (int p = 0; p < nact; ++p)
  {
    const unsigned int indp = indexmap[p];
    for (int q = 0; q < nact; ++q)
    {
      const unsigned int indq = indexmap[q];
      for (int r = 0; r < nact; ++r)
      {
        const unsigned int indr = indexmap[r];
        for (int s = 0; s < nact; ++s)
        {
          const unsigned int inds = indexmap[s];
          scmat_.set_element(indp *nocc + indq, indr *nocc+inds, ext_act_rdm2.get_element(p*nact+q, r*nact+s));
        }
      }
    }
  }
  //scmat_.print("ExternSpinFreeRDMTwo:: MO density");
}

void
ExternSpinFreeRDMTwo::init_from_rdm2_occspace(const std::vector<unsigned int>& indexmap,
                                              const Ref<OrbitalSpace> & occ_orbs)
{
  std::ifstream in(filename_.c_str());
  if (in.is_open() == false) {
    std::ostringstream oss;
    oss << "ExternSpinFreeRDMTwo: could not open file " << filename_;
    throw std::runtime_error(oss.str().c_str());
  }

  const unsigned int norbs = occ_orbs->coefs().coldim().n();
  RefSCDimension dim = new SCDimension(norbs * norbs);
  dim->blocks()->set_subdim(0, new SCDimension(dim->n()));
  scmat_ = occ_orbs->coefs().kit()->symmmatrix(dim);
  scmat_.assign(0.0);
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
    //in >> bra1 >> ket2 >> ket1 >> bra2 >> value;
    if (bra1 != -1) {
      --bra1;
      --bra2;
      --ket1;
      --ket2;
      const unsigned int mbra1 = indexmap[bra1];
      const unsigned int mbra2 = indexmap[bra2];
      const unsigned int mket1 = indexmap[ket1];
      const unsigned int mket2 = indexmap[ket2];

      scmat_.set_element(mbra1 * norbs + mbra2, mket1 * norbs + mket2, value);
      scmat_.set_element(mbra2 * norbs + mbra1, mket2 * norbs + mket1, value);
    } else
      have_coefs = false;
  }
  in.close();

  //scmat_.print("ExternReadRDMTwo:: MO density");
}

ExternSpinFreeRDMTwo::~ExternSpinFreeRDMTwo()
{
}

Ref< SpinFreeRDM<One> >
ExternSpinFreeRDMTwo::rdm_m_1() const
{
  if (rdm1_.null()) {
    RefSymmSCMatrix rdm1 = orbs_->coefs().kit()->symmmatrix(
        orbs_->coefs().coldim());
    rdm1.assign(0.0);
    const unsigned int norbs = rdm1.n();
    for (unsigned int b1 = 0; b1 < norbs; ++b1) {
      const unsigned b12_offset = b1 * norbs;
      for (unsigned int k1 = 0; k1 <= b1; ++k1) {
        const unsigned k12_offset = k1 * norbs;
        double value = 0.0;
        for (unsigned int i2 = 0; i2 < norbs; ++i2) {
          value += scmat_.get_element(b12_offset + i2, k12_offset + i2);
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

    rdm1_ = new ExternSpinFreeRDMOne(rdm1, orbs_);
  }

  return rdm1_;
}

const Ref<DistArray4>&
ExternSpinFreeRDMTwo::da4() const {
  if (da4_.null()) {

    const int n = orbs()->dim().n();
    da4_ = make_distarray4(1, n, n, n, n);
    da4_->activate();

    BlockedSymmSCMatrix* blocked_scmat = dynamic_cast<BlockedSymmSCMatrix*>(scmat_.pointer());

    std::vector<double> k1k2_buf(n*n);
    int b12 = 0;
    for(int b1=0; b1<n; ++b1) {
      for(int b2=0; b2<n; ++b2, ++b12) {

        if (blocked_scmat) {
          const size_t n2 = n*n;
          for(size_t k12=0; k12<n2; ++k12)
            k1k2_buf[k12] = scmat_.get_element(b12, k12);
        }
        else {
          RefSCVector b12_row = scmat_.get_row(b12);
          b12_row.convert(&k1k2_buf[0]);
        }

        da4_->store_pair_block(b1, b2, 0, &(k1k2_buf[0]));
      }
    }

    if (da4_->data_persistent()) da4_->deactivate();
  }

  return da4_;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
