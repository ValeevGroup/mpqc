//
// approx_pairs.cc
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Jun 6, 2014
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

#include "approx_pairs.h"

using namespace sc;
using namespace std;

//#define USE_OLD_SOLIDHARM_ORDERING 0
//#if not USE_OLD_SOLIDHARM_ORDERING // CCA ordering
//static inline int kpure(int l, int m) { return l+m; }
//static inline int ipure_molden(int l, int m) { return m<0?2*-m:(m==0?0:2*m-1); }
//static inline int off_molden(int ioff, int nbf) { int m = ioff - nbf/2; return (m<0?2*-m:(m==0?0:2*m-1)); }
//#endif

ClassDesc ApproximatePairWriter::cd_(
    typeid(ApproximatePairWriter),
    "ApproximatePairWriter",
    1, // verion number
    "public DescribedClass",
    0, // default constructor pointer
    create<ApproximatePairWriter>, // KeyVal constructor pointer
    0 // StateIn constructor pointer
);



ApproximatePairWriter::ApproximatePairWriter(
    const Ref<KeyVal>& keyval
)
{
  // Get the filename
  if (keyval->exists("filename")) {
    filename_ = keyval->stringvalue("filename");
  }
  else {
    filename_ = "-";
  }

  // Get the basis for expressing the exact pairs
  exbs_ << keyval->describedclassvalue("exact_fit_basis", KeyValValueRefDescribedClass(0));
  if(exbs_.null()){
    throw InputError("Approximate pair writer requires a large basis to express the 'exact' pairs in",
        __FILE__, __LINE__, "exact_fit_basis");
  }

  // Read the pairs to print as [[center, function_on_center], [center, function_on_center]]
  if(not keyval->exists("pairs")) {
    throw InputError("Approximate pair writer needs list of pairs to work with",
        __FILE__, __LINE__, "pairs");
  }
  int npair = keyval->count("pairs");

  for(int i = 0; i < npair; ++i) {
    if(keyval->count("pairs", i) != 2) throw InputError("Invalid pair specification", __FILE__, __LINE__, "pairs");
    if(keyval->count("pairs", i, 0) != 2) throw InputError("Invalid pair specification", __FILE__, __LINE__, "pairs");
    if(keyval->count("pairs", i, 1) != 2) throw InputError("Invalid pair specification", __FILE__, __LINE__, "pairs");

    int shell_1_center = keyval->Va_intvalue("pairs", 3, i, 0, 0);
    int shell_1_offset = keyval->Va_intvalue("pairs", 3, i, 0, 1);
    int shell_2_center = keyval->Va_intvalue("pairs", 3, i, 1, 0);
    int shell_2_offset = keyval->Va_intvalue("pairs", 3, i, 1, 1);
    pairs_tmp_.push_back(std::array<int,4>{{shell_1_center, shell_1_offset, shell_2_center, shell_2_offset}});

  }

  normalize_pairs_ = keyval->booleanvalue("normalize_pairs", KeyValValueboolean(normalize_pairs_));

}


void
ApproximatePairWriter::initialize()
{
  if(initialized_) return;

  wfn_->integral()->set_basis(exbs_, exbs_);
  ints2c_ex_ = wfn_->integral()->coulomb<2>();
  wfn_->integral()->set_basis(wfn_->gbs_, wfn_->gbs_, exbs_);
  ints3c_ex_ = wfn_->integral()->coulomb<3>();

  for(auto&& pary : pairs_tmp_) {
    pairs_.emplace_back(
        wfn_->gbs_->shell_on_center(pary[0], pary[1]),
        wfn_->gbs_->shell_on_center(pary[2], pary[3])
    );
  }

  //----------------------------------------//

  initialized_ = true;
}


void
ApproximatePairWriter::write_pairs()
{

  initialize();

  //----------------------------------------//

  std::ostream *out;
  if (filename_ == "-") {
      out = &(ExEnv::out0());
  }
  else {
      ExEnv::out0() << incindent << indent << "ApproximatePairWriter"
                    << " is writing its output to \"" << filename_ << "\""
                    << std::endl;
      ExEnv::out0() << decindent;
      out = new std::ofstream(filename_.c_str());
  }

  //----------------------------------------//

  // Always write the [Molden Format] tag

  *out << "[Molden Format]" << endl;

  //----------------------------------------//

  // Write various sections

  write_atoms_section(*out);
  write_gto_section(*out);
  write_mo_section(*out);

  //----------------------------------------//

  if (filename_ == "-") {
      *out << decindent;
    }
  else {
      delete out;
    }

}

void
ApproximatePairWriter::write_atoms_section(std::ostream &out)
{
  Ref<Molecule> mol = wfn_->molecule();
  out << "[Atoms]";
  out << " AU";
  out << endl;

  out.fill(' ');
  out << std::fixed << std::setprecision(8);
  for(int iatom = 0; iatom < mol->natom(); ++iatom){
    out << std::setw(4) << mol->atom_symbol(iatom)
        << std::setw(6) << (iatom + 1)
        << std::setw(6) << mol->Z(iatom)
        << setw(14) << mol->r(iatom, 0)
        << setw(14) << mol->r(iatom, 1)
        << setw(14) << mol->r(iatom, 2) << endl;
  }

}

void
ApproximatePairWriter::write_gto_section(std::ostream &out)
{
  Ref<GaussianBasisSet> gbs = wfn_->gbs_;
  Ref<GaussianBasisSet> dfbs = wfn_->dfbs_;
  //assert(gbs->has_pure());
  //assert(dfbs->has_pure());
  //assert(exbs_->has_pure());
  assert(exbs_->ncenter() == gbs->ncenter());
  assert(exbs_->ncenter() == dfbs->ncenter());

  out << "[GTO]" << endl;
  int offset = 0;
  for(int center = 0; center < exbs_->ncenter(); ++center){
    out << setw(4) << (center + 1)
        << setw(4) << 0
        << endl;
    //----------------------------------------//
    int max_am = 1;
    for(int am = 0; am <= max_am; ++am) {
      for(int ish = 0; ish < gbs->nshell_on_center(center); ++ish){
        GaussianShell sh = gbs->shell(center, ish);
        if(sh.max_am() > max_am) max_am = sh.am(0);
        if(sh.max_am() != am) continue;
        int shidx = gbs->shell_on_center(center, ish);
        shell_map_gbs_[offset] = shidx;
        if(sh.ncontraction() == 1){
          out << setw(4) << sh.amchar(0)
              << setw(6) << sh.nprimitive()
              << setw(6) << "1.00"
              << endl;
          for(int iprim = 0; iprim < sh.nprimitive(); ++iprim){
            out << std::scientific << std::uppercase << std::setprecision(10);
            out << setw(20) << sh.exponent(iprim)
                << setw(20) << sh.coefficient_norm(0, iprim)
                << endl;
          } // end loop over primatives
        } // end if ncontraction == 1
        else{
          throw FeatureNotImplemented("generally contracted basis sets", __FILE__, __LINE__, class_desc());
        }
        ++offset;
      }
      //----------------------------------------//
      for(int ish = 0; ish < dfbs->nshell_on_center(center); ++ish){
        GaussianShell sh = dfbs->shell(center, ish);
        if(sh.max_am() > max_am) max_am = sh.am(0);
        if(sh.max_am() != am) continue;
        int shidx = dfbs->shell_on_center(center, ish);
        shell_map_dfbs_[offset] = shidx;
        if(sh.ncontraction() == 1){
          out << setw(4) << sh.amchar(0)
              << setw(6) << sh.nprimitive()
              << setw(6) << "1.00"
              << endl;
          for(int iprim = 0; iprim < sh.nprimitive(); ++iprim){
            out << std::scientific << std::uppercase << std::setprecision(10);
            out << setw(20) << sh.exponent(iprim)
                << setw(20) << sh.coefficient_norm(0, iprim)
                << endl;
          } // end loop over primatives
        } // end if ncontraction == 1
        else{
          throw FeatureNotImplemented("generally contracted basis sets", __FILE__, __LINE__, class_desc());
        }
        ++offset;
      }
      //----------------------------------------//
      for(int ish = 0; ish < exbs_->nshell_on_center(center); ++ish){
        GaussianShell sh = exbs_->shell(center, ish);
        if(sh.max_am() > max_am) max_am = sh.am(0);
        if(sh.max_am() != am) continue;
        int shidx = exbs_->shell_on_center(center, ish);
        shell_map_exbs_[offset] = shidx;
        if(sh.ncontraction() == 1){
          out << setw(4) << sh.amchar(0)
              << setw(6) << sh.nprimitive()
              << setw(6) << "1.00"
              << endl;
          for(int iprim = 0; iprim < sh.nprimitive(); ++iprim){
            out << std::scientific << std::uppercase << std::setprecision(10);
            out << setw(20) << sh.exponent(iprim)
                << setw(20) << sh.coefficient_norm(0, iprim)
                << endl;
          } // end loop over primatives
        } // end if ncontraction == 1
        else{
          throw FeatureNotImplemented("generally contracted basis sets", __FILE__, __LINE__, class_desc());
        }
        ++offset;
      }
    }
    //----------------------------------------//
    // format calls for a blank line here
    out << endl;
  } // end loop over centers
}

void
ApproximatePairWriter::write_mo(
    std::ostream& out,
    const Eigen::VectorXd& coefs,
    int orb_type,
    const std::string& name
)
{
  decltype(shell_map_exbs_)& shell_map =
      (orb_type == GBS_orbs)  ? shell_map_gbs_ : (
       orb_type == DFBS_pair ?  shell_map_dfbs_ :
                                shell_map_exbs_
  );

  const int nshell_tot = wfn_->gbs_->nshell() + wfn_->dfbs_->nshell() + exbs_->nshell();
  Ref<GaussianBasisSet>& bs =
      (orb_type == GBS_orbs)  ? wfn_->gbs_ : (
       orb_type == DFBS_pair ?  wfn_->dfbs_ :
                                exbs_
  );

  out << setw(8) << "Sym="
      << name
      << endl
      << setw(8) << "Ene="
      << setw(16) << 0.0
      << endl
      << setw(8) << "Spin="
      << "Alpha"
      << endl
      << setw(8) << "Occup="
      << 0.0
      << endl;
  int fxn_offset = 0;
  for(int ish = 0; ish < nshell_tot; ++ish){
    Ref<GaussianBasisSet> idxbs;
    bool use_coefs = false;
    decltype(shell_map_exbs_)* idxshmap;
    if(shell_map_gbs_.find(ish) != shell_map_gbs_.end()) {
      idxbs = wfn_->gbs_;
      idxshmap = &shell_map_gbs_;
      use_coefs = orb_type == GBS_orbs;
    }
    else if(shell_map_dfbs_.find(ish) != shell_map_dfbs_.end()) {
      idxbs = wfn_->dfbs_;
      idxshmap = &shell_map_dfbs_;
      use_coefs = orb_type == DFBS_pair;
    }
    else if(shell_map_exbs_.find(ish) != shell_map_exbs_.end()){
      idxbs = exbs_;
      idxshmap = &shell_map_exbs_;
      use_coefs = orb_type == EX_pair;
    }
    else {
      throw ProgrammingError("huh?", __FILE__, __LINE__);
    }
    const int idxish = (*idxshmap)[ish];
    const int bfoff = idxbs->shell_to_function(idxish);
    const int inbf = idxbs->shell(idxish).nfunction();
    const int bfend = idxbs->shell_to_function(idxish) + inbf;
    for(int ibf = bfoff; ibf < bfend; ++ibf) {
      //int ibf_molden;
      //if(inbf > 3) ibf_molden = bfoff + off_molden(ibf-bfoff, inbf);
      //else
      int ibf_molden = ibf;
      out << setw(6) << (fxn_offset+1)
      << setw(16) << (use_coefs ? coefs(ibf_molden) : 0.0)
      << endl;
      ++fxn_offset;
    }
  } // end loop over AOs

}

void
ApproximatePairWriter::compute_pairs_ex()
{
  auto g2exptr = wfn_->ints_to_eigen(
      ShellBlockData<>(exbs_), ShellBlockData<>(exbs_),
      ints2c_ex_, wfn_->metric_oper_type_
  );
  auto& g2ex = *g2exptr;
  CADFCLHF::Decomposition decomp(g2ex);

  for(auto&& pair : pairs_) {
    Eigen::MatrixXd& cmat = excoefs_[pair];
    ShellData ish(pair.first, wfn_->gbs_);
    ShellData jsh(pair.second, wfn_->gbs_);
    cmat.resize(ish.nbf*jsh.nbf, exbs_->nbasis());
    auto g3exptr = wfn_->ints_to_eigen(
      ish, jsh, ShellBlockData<>(exbs_),
      ints3c_ex_, wfn_->metric_oper_type_
    );
    auto& g3ex = *g3exptr;

    for(auto&& mu : function_range(ish)) {
      for(auto&& nu : function_range(jsh)) {
        const int mu_nu = mu.off * jsh.nbf + nu.off;
        cmat.row(mu_nu) = decomp.solve(g3ex.row(mu_nu).transpose());
      }
    }

  }

}

void
ApproximatePairWriter::write_mo_section(std::ostream &out)
{
  out << "[MO]" << endl;

  out << fixed << setprecision(10);
  // Print out the orbitals comprising the pairs
  for(int ipair = 0; ipair < pairs_.size(); ++ipair){
    Eigen::VectorXd coefs(wfn_->gbs_->nbasis());
    coefs.setZero();
    for(int pair_part = 0; pair_part < 2; ++pair_part) {
      const int ish = pair_part == 0 ? pairs_[ipair].first : pairs_[ipair].second;
      const int bfoff = wfn_->gbs_->shell_to_function(ish);
      const int bfend = wfn_->gbs_->shell_to_function(ish) + wfn_->gbs_->shell(ish).nfunction();
      for(int ibf = bfoff; ibf < bfend; ++ibf) {
        std::ostringstream sstr;
        sstr << "AOsh_" << ish << "_bf_" << ibf << "_pair_" << ipair << "_" << pair_part+1;
        coefs[ibf] = 1.0;
        write_mo(out, coefs, GBS_orbs, sstr.str());
        coefs[ibf] = 0.0;
      }
    }
  }

  compute_pairs_ex();

  // Print out the approximate pairs and the "exact" pairs
  for(int ipair = 0; ipair < pairs_.size(); ++ipair){
    ShellData ish(pairs_[ipair].first, wfn_->gbs_, wfn_->dfbs_);
    ShellData jsh(pairs_[ipair].second, wfn_->gbs_, wfn_->dfbs_);
    for(auto&& mu : function_range(ish)) {
      for(auto&& nu : function_range(jsh)) {
        {
          std::ostringstream sstr;
          sstr << "pair_" << ipair << "_bf_" << mu.index << "_" << nu.index << "_approx";
          Eigen::VectorXd cdf(wfn_->dfbs_->nbasis());
          cdf.setZero();
          if(nu <= mu) {
            cdf.segment(ish.atom_dfbfoff, ish.atom_dfnbf) = *(wfn_->coefs_[{mu,nu}].first);
          }
          else {
            cdf.segment(ish.atom_dfbfoff, ish.atom_dfnbf) = *(wfn_->coefs_[{nu,mu}].second);
          }
          if(normalize_pairs_) cdf.array() /= cdf.norm();
          if(ish.center != jsh.center) {
            if(nu <= mu) {
              cdf.segment(jsh.atom_dfbfoff, jsh.atom_dfnbf) = *(wfn_->coefs_[{mu,nu}].second);
            }
            else {
              cdf.segment(jsh.atom_dfbfoff, jsh.atom_dfnbf) = *(wfn_->coefs_[{nu,mu}].first);
            }
          }
          write_mo(out, cdf, DFBS_pair, sstr.str());
        }
        {
          std::ostringstream sstr;
          sstr << "pair_" << ipair << "_bf_" << mu.index << "_" << nu.index << "_exact";
          Eigen::VectorXd cdf(exbs_->nbasis());
          cdf = excoefs_[pairs_[ipair]].row(mu.off*jsh.nbf + nu.off);
          if(normalize_pairs_) cdf.array() /= cdf.norm();
          write_mo(out, cdf, DFBS_pair, sstr.str());
        }
      }
    }
  }

}



