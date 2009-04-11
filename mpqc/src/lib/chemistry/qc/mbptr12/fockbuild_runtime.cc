//
// fockbuild_runtime.cc
//
// Copyright (C) 2009 Edward Valeev
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

#include <sstream>
#include <math/scmat/repl.h>
#include <util/class/scexception.h>
#include <chemistry/qc/mbptr12/fockbuild_runtime.h>
#include <chemistry/qc/mbptr12/registry.timpl.h>
#include <chemistry/qc/mbptr12/orbitalspace.h>
#include <chemistry/qc/mbptr12/fockbuilder.h>

using namespace sc;

static ClassDesc FockBuildRuntime_cd(typeid(FockBuildRuntime),
                                     "FockBuildRuntime", 1,
                                     "virtual public SavableState", 0, 0,
                                     create<FockBuildRuntime> );

FockBuildRuntime::FockBuildRuntime(const Ref<GaussianBasisSet>& refbasis,
                                   const RefSymmSCMatrix& aodensity_alpha,
                                   const RefSymmSCMatrix& aodensity_beta,
                                   const Ref<Integral>& integral,
                                   Ref<MessageGrp> msg,
                                   Ref<ThreadGrp> thr) :
  P_(aodensity_alpha+aodensity_beta), Po_(0), basis_(refbasis), integral_(integral),
  msg_(msg), thr_(thr),
  registry_(FockMatrixRegistry::instance()) {

  RefSymmSCMatrix Po = aodensity_alpha - aodensity_beta;
  spin_polarized_ = Po->maxabs() > DBL_EPSILON;
  if (spin_polarized_)
    Po_ = Po;

}

FockBuildRuntime::FockBuildRuntime(StateIn& si) {
  int spin_polarized;  si.get(spin_polarized);  spin_polarized_ = spin_polarized;
  RefSCDimension pdim;
  pdim << SavableState::restore_state(si);
  Ref<SCMatrixKit> kit = new ReplSCMatrixKit;
  P_ = kit->symmmatrix(pdim);
  P_.restore(si);
  if (spin_polarized_) {
    Po_ = kit->symmmatrix(pdim);
    Po_.restore(si);
  }
  basis_ << SavableState::restore_state(si);
  integral_ << SavableState::restore_state(si);
  registry_ = FockMatrixRegistry::restore_instance(si);

  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();
}

void FockBuildRuntime::save_data_state(StateOut& so) {
  so.put((int)spin_polarized_);
  SavableState::save_state(P_.dim().pointer(), so);
  P_.save(so);
  if (spin_polarized_)
    Po_.save(so);
  SavableState::save_state(basis_.pointer(), so);
  SavableState::save_state(integral_.pointer(), so);
  FockMatrixRegistry::save_instance(registry_, so);
}

void FockBuildRuntime::validate_key(const std::string& key) const {
  try {
    ParsedOneBodyIntKey parsedkey(key);
  } catch (...) {
    std::ostringstream oss;
    oss << "FockBuildRuntime::get() -- key " << key
        << " does not match the format";
    throw ProgrammingError(oss.str().c_str(), __FILE__,__LINE__);
  }
}

namespace {

  std::string
  transposed_key(const std::string& key) {
    ParsedOneBodyIntKey pkey(key);
    std::string tkey = ParsedOneBodyIntKey::key(pkey.ket(), pkey.bra(),
                                                pkey.oper(), pkey.spin());
    return tkey;
  }

  RefSCMatrix
  SymmToRect(const RefSymmSCMatrix& symm) {
    RefSCMatrix rect = symm.kit()->matrix(symm.dim(),symm.dim());
    rect.assign(0.0);
    rect->accumulate(symm.pointer());
    return rect;
  }
}

RefSCMatrix
FockBuildRuntime::get(const std::string& key) {
  validate_key(key);
  if (registry_->key_exists(key)) {
    return registry_->value(key);
  } else { // if not found

    // try transpose first
    const std::string tkey = transposed_key(key);
    if (registry_->key_exists(tkey)) {
      return registry_->value(tkey).t();
    }

    ParsedOneBodyIntKey pkey(key);
    const std::string& bra_key = pkey.bra();
    const std::string& ket_key = pkey.ket();
    const std::string& oper_key = pkey.oper();
    Ref<OrbitalSpaceRegistry> idxreg = OrbitalSpaceRegistry::instance();
    Ref<OrbitalSpace> bra = idxreg->value(bra_key);
    Ref<OrbitalSpace> ket = idxreg->value(ket_key);

    // determine spin
    const SpinCase1 spin = pkey.spin();

    // is the AO matrix available?
    Ref<AOSpaceRegistry> aoidxreg = AOSpaceRegistry::instance();
    Ref<OrbitalSpace> aobra = aoidxreg->value(bra->basis());
    Ref<OrbitalSpace> aoket = aoidxreg->value(ket->basis());
    const std::string& aobra_key = idxreg->key(aobra);
    const std::string& aoket_key = idxreg->key(aoket);
    const std::string aokey = ParsedOneBodyIntKey::key(aobra_key, aoket_key,
                                                       oper_key, spin);
    if (registry_->key_exists(aokey)) {
      RefSCMatrix aofock = registry_->value(aokey);
      RefSCMatrix mofock = bra->coefs().t() * aofock * ket->coefs();
      registry_->add(key, mofock);
      return registry_->value(key);
    }
    const std::string transposed_aokey = transposed_key(aokey);
    if (registry_->key_exists(transposed_aokey)) {
      RefSCMatrix transposed_aofock = registry_->value(transposed_aokey);
      RefSCMatrix mofock = bra->coefs().t() * transposed_aofock.t() * ket->coefs();
      registry_->add(key, mofock);
      return registry_->value(key);
    }

    // AO matrix not found: compute all components of it first, then call itself again
    const Ref<GaussianBasisSet>& bs1 = bra->basis();
    const Ref<GaussianBasisSet>& bs2 = ket->basis();
    const bool bs1_eq_bs2 = bs1->equiv(bs2);
    { // core hamiltonian
      RefSCMatrix H;
      const Ref<GaussianBasisSet>& obs = basis_;
      if (bs1_eq_bs2) {
        Ref<OneBodyFockMatrixBuilder<true> >
            fmb = new OneBodyFockMatrixBuilder<true> (OneBodyFockMatrixBuilder<true>::NonRelativistic,
                bs1, bs2, obs, integral());

        RefSymmSCMatrix Hsymm = fmb->result();
        // convert to H
        H = SymmToRect(Hsymm);
      } else { // result is rectangular already

        Ref<OneBodyFockMatrixBuilder<false> >
            fmb =
                new OneBodyFockMatrixBuilder<false> (OneBodyFockMatrixBuilder<false>::NonRelativistic,
                    bs1, bs2, obs, integral());
        H = fmb->result();
      }
      const std::string hkey = ParsedOneBodyIntKey::key(aobra_key,aoket_key,std::string("H"));
      registry_->add(hkey, H);
    }

    { // J, K, and F
      const Ref<GaussianBasisSet>& obs = basis_;
      bool compute_F = false;
      bool compute_J = (oper_key == "J" || oper_key == "F");
      bool compute_K = (oper_key == "K" || oper_key == "F");

      Ref<TwoBodyFockMatrixDFBuilder> fmb_df;
      if (use_density_fitting() && compute_J)
        fmb_df = new TwoBodyFockMatrixDFBuilder(false,
                                               true,
                                               false,
                                               bs1, bs2, obs,
                                               P_,
                                               Po_,
                                               dfinfo());

      double nints;
      if (bs1_eq_bs2) {
        Ref<TwoBodyFockMatrixBuilder<true> > fmb;
        if ((!use_density_fitting() && compute_J) || compute_K) {
          fmb = new TwoBodyFockMatrixBuilder<true> (compute_F, compute_J,
                                                compute_K, bs1, bs2, obs, P_,
                                                Po_, integral(), msg(), thr());
          nints = fmb->nints();
        }
        {
          RefSCMatrix J;
          if (compute_J) {
            J = use_density_fitting() ? fmb_df->J() : SymmToRect(fmb->J());
            const std::string jkey = ParsedOneBodyIntKey::key(aobra_key,aoket_key,std::string("J"));
            registry_->add(jkey, J);
            if (debug()) {
              J.print(jkey.c_str());
            }
          }

          RefSCMatrix K;
          if (compute_K) {
            K = SymmToRect(fmb->K(spin));
            const std::string kkey = ParsedOneBodyIntKey::key(aobra_key,aoket_key,std::string("K"),spin);
            registry_->add(kkey, K);
            if (debug()) {
              K.print(kkey.c_str());
            }
          }

          RefSCMatrix F;
          if (compute_J && compute_K) {
            F = K.clone(); F.assign(K); F.scale(-1.0); F.accumulate(J);
            const std::string fkey = ParsedOneBodyIntKey::key(aobra_key,aoket_key,std::string("F"),spin);
            registry_->add(fkey, F);
            if (debug()) {
              F.print(fkey.c_str());
            }
          }

        }

      } else { // result is rectangular already

        Ref<TwoBodyFockMatrixBuilder<false> > fmb;
        if ((!use_density_fitting() && compute_J) || compute_K) {
          fmb = new TwoBodyFockMatrixBuilder<false> (compute_F, compute_J,
                                                 compute_K, bs1, bs2, obs, P_,
                                                 Po_, integral(),
                                                 msg(),
                                                 thr());
          nints = fmb->nints();
        }
        {
          RefSCMatrix J;

          if (compute_J) {
            J = use_density_fitting() ? fmb_df->J() : fmb->J();
            const std::string jkey = ParsedOneBodyIntKey::key(aobra_key,aoket_key,std::string("J"));
            registry_->add(jkey, J);
            if (debug()) {
              J.print(jkey.c_str());
            }
          }

          RefSCMatrix K;
          if (compute_K) {
            K = fmb->K(spin);
            const std::string kkey = ParsedOneBodyIntKey::key(aobra_key,aoket_key,std::string("K"),spin);
            registry_->add(kkey, K);
            if (debug()) {
              K.print(kkey.c_str());
            }
          }

          RefSCMatrix F;
          if (compute_J && compute_K) {
            F= K.clone(); F.assign(K); F.scale(-1.0); F.accumulate(J);
            const std::string fkey = ParsedOneBodyIntKey::key(aobra_key,aoket_key,std::string("F"),spin);
            registry_->add(fkey, F);
            if (debug()) {
              F.print(fkey.c_str());
            }
          }
        }
      }
    }

#if 0
    else { // use density-fitting-based Fock builder

      const Ref<GaussianBasisSet>& obs = basis_;
      const bool compute_F = false;
      const bool compute_J = true;
      const bool compute_K = true;

      Ref<TwoBodyFockMatrixDFBuilder> fmb =
          new TwoBodyFockMatrixDFBuilder(compute_F,
                                         compute_J,
                                         compute_K,
                                         bs1, bs2, obs,
                                         P_,
                                         Po_,
                                         dfinfo());
      {
        RefSCMatrix J = fmb->J();
        const std::string jkey = ParsedOneBodyIntKey::key(aobra_key,aoket_key,std::string("J"));
        registry_->add(jkey, J);

        RefSCMatrix K = fmb->K(spin);
        const std::string kkey = ParsedOneBodyIntKey::key(aobra_key,aoket_key,std::string("K"),spin);
        registry_->add(kkey, K);

        RefSCMatrix F = K.clone(); F.assign(K); F.scale(-1.0); F.accumulate(J);
        const std::string fkey = ParsedOneBodyIntKey::key(aobra_key,aoket_key,std::string("F"),spin);
        registry_->add(fkey, F);

        if (debug()) {
          J.print(jkey.c_str());
          K.print(kkey.c_str());
          F.print(fkey.c_str());
        }
      }
    } // end of DF-based computation
#endif

    // now all components are available, call itself again
    return get(key);
  }
  abort(); // unreachable
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

ParsedOneBodyIntKey::ParsedOneBodyIntKey(const std::string& key) :
  key_(key) {
  typedef std::string::const_iterator citer;
  std::string keycopy(key);

  // pop off the leading '<'
  assert(keycopy[0] == '<');
  keycopy.erase(keycopy.begin());
  // get bra
  bra_ = pop_till_token(keycopy, '|');
  // get oper
  oper_ = pop_till_token(keycopy, '|');
  // get ket
  ket_ = pop_till_token(keycopy, '>');
  // get spin
  if (!keycopy.empty()) {
    const std::string tmp = pop_till_token(keycopy, '[');
    const std::string spinkey = pop_till_token(keycopy,']');
    spin_ = to_spincase1(spinkey);
  }
  else {
    spin_ = AnySpinCase1;
  }


#if 0
  ExEnv::out0() << indent << "ParsedOneBodyIntKey::ParsedOneBodyIntKey():" << std::endl << incindent;
  ExEnv::out0() << indent << "key = " << key_ << std::endl;
  ExEnv::out0() << indent << "bra = " << bra_ << std::endl;
  ExEnv::out0() << indent << "ket = " << ket_ << std::endl;
  ExEnv::out0() << indent << "oper = " << oper_ << std::endl;
  ExEnv::out0() << indent << "spin = " << to_string(spin_) << std::endl;
#endif
}

std::string ParsedOneBodyIntKey::key(const std::string& bra,
                                     const std::string& ket,
                                     const std::string& oper,
                                     SpinCase1 spin) {
  std::ostringstream oss;
  oss << "<" << bra << "|" << oper << "|" << ket << ">";
  if (spin != AnySpinCase1)
    oss << "[" << to_string(spin) << "]";
  return oss.str();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
