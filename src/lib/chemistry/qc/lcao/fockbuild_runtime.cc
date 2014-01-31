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

#include <sstream>
#include <cassert>
#include <math/scmat/repl.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/lcao/fockbuild_runtime.h>
#include <util/misc/registry.timpl.h>
#include <chemistry/qc/wfn/orbitalspace.h>
#include <chemistry/qc/lcao/fockbuilder.h>

using namespace sc;

static ClassDesc FockBuildRuntime_cd(typeid(FockBuildRuntime),
                                     "FockBuildRuntime", 1,
                                     "virtual public SavableState", 0, 0,
                                     create<FockBuildRuntime> );

FockBuildRuntime::FockBuildRuntime(const Ref<OrbitalSpaceRegistry>& oreg,
                                   const Ref<AOSpaceRegistry>& aoreg,
                                   const Ref<GaussianBasisSet>& refbasis,
                                   const RefSymmSCMatrix& aodensity_alpha,
                                   const RefSymmSCMatrix& aodensity_beta,
                                   const Ref<Integral>& integral,
                                   const RefSCVector& efield,
                                   Ref<MessageGrp> msg,
                                   Ref<ThreadGrp> thr) :
  oreg_(oreg), aoreg_(aoreg), P_(aodensity_alpha+aodensity_beta), Po_(0),
  basis_(refbasis), integral_(integral),
  efield_(efield), msg_(msg), thr_(thr),
  log2_precision_(-50.0),
  registry_(FockMatrixRegistry::instance()),
  psqrtregistry_(PSqrtRegistry::instance()) {
  RefSymmSCMatrix Po = aodensity_alpha - aodensity_beta;
  spin_polarized_ = Po->maxabs() > DBL_EPSILON;
  if (spin_polarized_)
    Po_ = Po;
}

FockBuildRuntime::FockBuildRuntime(StateIn& si) {
  oreg_ = OrbitalSpaceRegistry::restore_instance(si);
  aoreg_ = AOSpaceRegistry::restore_instance(si);
  si.get(spin_polarized_);
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
  psqrtregistry_ = PSqrtRegistry::restore_instance(si);

  efield_ = SCMatrixKit::default_matrixkit()->vector(new SCDimension(3));
  efield_.restore(si);

  si.get(log2_precision_);

  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();
}

void FockBuildRuntime::save_data_state(StateOut& so) {
  OrbitalSpaceRegistry::save_instance(oreg_, so);
  AOSpaceRegistry::save_instance(aoreg_, so);
  so.put(spin_polarized_);
  SavableState::save_state(P_.dim().pointer(), so);
  P_.save(so);
  if (spin_polarized_)
    Po_.save(so);
  SavableState::save_state(basis_.pointer(), so);
  SavableState::save_state(integral_.pointer(), so);
  FockMatrixRegistry::save_instance(registry_, so);
  PSqrtRegistry::save_instance(psqrtregistry_, so);
  efield_.save(so);
  so.put(log2_precision_);
}

void FockBuildRuntime::set_densities(const RefSymmSCMatrix& aodensity_alpha,
                                     const RefSymmSCMatrix& aodensity_beta) {
  bool new_P = false;
  RefSymmSCMatrix P, Po;

  P = aodensity_alpha.copy(); P.accumulate(aodensity_beta);
  new_P = new_P || (P - P_)->maxabs() > DBL_EPSILON;
  //(P - P_).print("FockBuildRuntime::set_densities() -- \"new - old\" density");
  Po = aodensity_alpha - aodensity_beta;
  RefSymmSCMatrix dPo = Po_ ? (Po - Po_) : Po;
  new_P = new_P || dPo->maxabs() > DBL_EPSILON;
  if (new_P) {
    this->obsolete_density_dependents();
    P_ = P;
    spin_polarized_ = Po->maxabs() > DBL_EPSILON;
    if (spin_polarized_)
      Po_ = Po;
  }
}

void FockBuildRuntime::set_electric_field(const RefSCVector& efield) {
  efield_ = efield;
}

void
FockBuildRuntime::obsolete() {
  registry_->clear();
  psqrtregistry_->clear();
  dfinfo()->obsolete();
}

namespace {
  template <typename T, typename U>
  struct key_equals {
      key_equals(T k) : key(k) {}
      bool operator()(const std::pair<T, U>& i) const {
        return i.first == key;
      }
      T key;
  };
}

void
FockBuildRuntime::obsolete_density_dependents() {
  psqrtregistry_->clear();
  // leave H core intact? too cheap to bother
  registry_->clear();
  // purge density-dependent spaces
  if (dfinfo()) {
    dfinfo()->runtime()->remove_if(std::string("dd"));
    dfinfo()->runtime()->moints_runtime()->runtime_2c()->remove_if(std::string("dd"));
    dfinfo()->runtime()->moints_runtime()->runtime_3c()->remove_if(std::string("dd"));
    key_equals<std::string, RefSymmSCMatrix> pred("dd");
    dfinfo()->runtime()->moints_runtime()->runtime_2c_inv()->remove_if(pred);
  }
}

void
FockBuildRuntime::set_log2_precision(double prec) {
  if (prec < log2_precision_)
    this->obsolete();
  log2_precision_ = prec;
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

  std::string
  alternatespin_key(const std::string& key, SpinCase1 new_spin) {
    ParsedOneBodyIntKey pkey(key);
    std::string askey = ParsedOneBodyIntKey::key(pkey.bra(), pkey.ket(),
                                                 pkey.oper(), new_spin);
    return askey;
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

  // N.B. will ask for field-free matrices first, and add any additional contributions later

  // parse the key
  ParsedOneBodyIntKey pkey(key);
  const std::string& bra_key = pkey.bra();
  const std::string& ket_key = pkey.ket();
  const std::string& oper_key = pkey.oper();
  const SpinCase1 spin = pkey.spin();

  RefSCMatrix result;
  if (registry_->key_exists(key)) {
    result = registry_->value(key);
  }
  else { // if not found

    // try transpose first
    const std::string tkey = transposed_key(key);
    if (registry_->key_exists(tkey)) {
      result = registry_->value(tkey).t();
    }
    else {
      // if spin-nonpolarized and spin != AnySpinCase1, ask for AnySpinCase1
      if (spin_polarized_ != true && spin != AnySpinCase1) {
        return this->get(alternatespin_key(key, AnySpinCase1));
      }
      else {

        Ref<OrbitalSpaceRegistry> idxreg = oreg_;
        Ref<OrbitalSpace> bra = idxreg->value(bra_key);
        Ref<OrbitalSpace> ket = idxreg->value(ket_key);

        // is the AO matrix available?
        Ref<AOSpaceRegistry> aoidxreg = aoreg_;
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
          result = registry_->value(key);
        }
        else {
          const std::string transposed_aokey = transposed_key(aokey);
          if (registry_->key_exists(transposed_aokey)) {
            RefSCMatrix transposed_aofock = registry_->value(transposed_aokey);
            RefSCMatrix mofock = bra->coefs().t() * transposed_aofock.t()
                * ket->coefs();
            registry_->add(key, mofock);
            result = registry_->value(key);
          }
          else {

            // Fock matrices
            if (oper_key == "H" || oper_key == "J" || oper_key == "K"
                || oper_key == "F") {

              // AO matrix not found: compute all components of it first, then call itself again
              const Ref<GaussianBasisSet>& bs1 = bra->basis();
              const Ref<GaussianBasisSet>& bs2 = ket->basis();
              const bool bs1_eq_bs2 = bs1->equiv(bs2);

              const std::string hkey = ParsedOneBodyIntKey::key(
                  aobra_key, aoket_key, std::string("H"));
              const std::string jkey = ParsedOneBodyIntKey::key(
                  aobra_key, aoket_key, std::string("J"));
              const std::string kkey = ParsedOneBodyIntKey::key(
                  aobra_key, aoket_key, std::string("K"), spin);
              const std::string fkey = ParsedOneBodyIntKey::key(
                  aobra_key, aoket_key, std::string("F"), spin);
              const bool have_H = registry_->key_exists(hkey);
              const bool have_J = registry_->key_exists(jkey);
              const bool have_K = registry_->key_exists(kkey);
              const bool have_F = registry_->key_exists(fkey);
              const bool need_H = (oper_key == "H" || oper_key == "F");
              const bool need_J = (oper_key == "J" || oper_key == "F");
              const bool need_K = (oper_key == "K" || oper_key == "F");
              const bool need_F = (oper_key == "F");
              const bool compute_F = false; // tell FockBuilder to not compute F; compute it myself from components
              const bool compute_H = need_H && !have_H;
              const bool compute_J = need_J && !have_J;
              const bool compute_K = need_K && !have_K;

              RefSCMatrix H;
              if (need_H) { // compute core hamiltonian
                if (compute_H) {
                  const Ref<GaussianBasisSet>& obs = basis_;
                  if (bs1_eq_bs2) {
                    Ref<OneBodyFockMatrixBuilder<true> > fmb =
                        new OneBodyFockMatrixBuilder<true>(
                            OneBodyFockMatrixBuilder<true>::NonRelativistic,
                            bs1, bs2, obs, integral(),
                            pow(2.0, log2_precision_));

                    RefSymmSCMatrix Hsymm = fmb->result();
                    // convert to H
                    H = SymmToRect(Hsymm);
                  } else { // result is rectangular already

                    Ref<OneBodyFockMatrixBuilder<false> > fmb =
                        new OneBodyFockMatrixBuilder<false>(
                            OneBodyFockMatrixBuilder<false>::NonRelativistic,
                            bs1, bs2, obs, integral(),
                            pow(2.0, log2_precision_));
                    H = fmb->result();
                  }
                  registry_->add(hkey, H);
                  if (debug()) {
                    H.print(hkey.c_str());
                  }
                } else { // have_H == true
                  H = registry_->value(hkey);
                }
              }

              { // J, K, and F
                const Ref<GaussianBasisSet>& obs = basis_;
                Ref<TwoBodyFockMatrixDFBuilder> fmb_df;
                if (use_density_fitting())
                  fmb_df = new TwoBodyFockMatrixDFBuilder(compute_F, compute_J,
                                                          compute_K, bs1, bs2,
                                                          obs, P_, Po_,
                                                          dfinfo(),
                                                          psqrtregistry_);

                double nints;
                if (bs1_eq_bs2) {
                  Ref<TwoBodyFockMatrixBuilder<true> > fmb;
                  if (!use_density_fitting()) {
                    fmb = new TwoBodyFockMatrixBuilder<true>(
                        compute_F, compute_J, compute_K, bs1, bs2, obs, P_, Po_,
                        integral(), msg(), thr(), pow(2.0, log2_precision_));
                    nints = fmb->nints();
                  }
                  {
                    RefSCMatrix J;
                    if (need_J) {
                      if (compute_J) {
                        J = use_density_fitting() ? fmb_df->J() :
                                                    SymmToRect(fmb->J());
                        registry_->add(jkey, J);
                        if (debug()) {
                          J.print(jkey.c_str());
                        }
                      } else { // have_J == true
                        J = registry_->value(jkey);
                      }
                    }

                    RefSCMatrix K;
                    if (need_K) {
                      if (compute_K) {
                        // since non-DF based Fock builder computes components for exchange of both spins, ask for K matrices for each spin here
                        const int nunique_spins = spin_polarized_ ? 2 : 1;
                        // refer to spin indirectly to properly handle AnySpinCase1
                        std::vector<SpinCase1> spins(nunique_spins);
                        if (spin_polarized_) {
                          spins[0] = Alpha;
                          spins[1] = Beta;
                        } else {
                          spins[0] = AnySpinCase1;
                        }
                        for (int s = 0; s < nunique_spins; ++s) {
                          const SpinCase1 spin1 = spins[s];
                          const std::string kkey = ParsedOneBodyIntKey::key(
                              aobra_key, aoket_key, std::string("K"), spin1);
                          RefSCMatrix KK =
                              use_density_fitting() ? fmb_df->K(spin1) :
                                                      SymmToRect(fmb->K(spin1));
                          registry_->add(kkey, KK);
                          if (debug()) {
                            KK.print(kkey.c_str());
                          }
                          if (spin == spin1)
                            K = KK;
                        }
                      } else { // have_K == true
                        K = registry_->value(kkey);
                      }
                    }

                    RefSCMatrix F;
                    if (need_F && !have_F) {
                      F = K.clone();
                      F.assign(K);
                      F.scale(-1.0);
                      F.accumulate(J);
                      F.accumulate(H);
                      registry_->add(fkey, F);
                      if (debug()) {
                        F.print(fkey.c_str());
                      }
                    }

                  }

                } else { // result is rectangular already

                  Ref<TwoBodyFockMatrixBuilder<false> > fmb;
                  if (!use_density_fitting()) {
                    fmb = new TwoBodyFockMatrixBuilder<false>(
                        compute_F, compute_J, compute_K, bs1, bs2, obs, P_, Po_,
                        integral(), msg(), thr(), pow(2.0, log2_precision_));
                    nints = fmb->nints();
                  }
                  {
                    RefSCMatrix J;
                    if (need_J) {
                      if (compute_J) {
                        J = use_density_fitting() ? fmb_df->J() : fmb->J();
                        registry_->add(jkey, J);
                        if (debug()) {
                          J.print(jkey.c_str());
                        }
                      } else { // have_J == true
                        J = registry_->value(jkey);
                      }
                    }

                    RefSCMatrix K;
                    if (need_K) {
                      if (compute_K) {
                        // since non-DF based Fock builder computes components for exchange of both spins, ask for K matrices for each spin here
                        const int nunique_spins = spin_polarized_ ? 2 : 1;
                        // refer to spin indirectly to properly handle AnySpinCase1
                        std::vector<SpinCase1> spins(nunique_spins);
                        if (spin_polarized_) {
                          spins[0] = Alpha;
                          spins[1] = Beta;
                        } else {
                          spins[0] = AnySpinCase1;
                        }
                        for (int s = 0; s < nunique_spins; ++s) {
                          const SpinCase1 spin1 = spins[s];
                          const std::string kkey = ParsedOneBodyIntKey::key(
                              aobra_key, aoket_key, std::string("K"), spin1);
                          RefSCMatrix KK =
                              use_density_fitting() ? fmb_df->K(spin1) :
                                                      fmb->K(spin1);
                          registry_->add(kkey, KK);
                          if (debug()) {
                            KK.print(kkey.c_str());
                          }
                          if (spin == spin1)
                            K = KK;
                        }
                      } else { // have_K == true
                        K = registry_->value(kkey);
                      }
                    }

                    RefSCMatrix F;
                    if (need_F && !have_F) {
                      F = K.clone();
                      F.assign(K);
                      F.scale(-1.0);
                      F.accumulate(J);
                      F.accumulate(H);
                      const std::string fkey = ParsedOneBodyIntKey::key(
                          aobra_key, aoket_key, std::string("F"), spin);
                      registry_->add(fkey, F);
                      if (debug()) {
                        F.print(fkey.c_str());
                      }
                    }
                  }
                }

              } // J, K, F components
            } // end of Fock matrices
            else if (oper_key.find("mu_") == 0 || // electric dipole moment
                     oper_key.find("q_") == 0  || // electric quadrupole moment
                     oper_key.find("dphi_") == 0 || // electric field
                     oper_key.find("ddphi_") == 0   // electric field gradient
                    ) {

              const unsigned int rank = oper_key.find("mu_") == 0 || oper_key.find("dphi_") == 0 ? 1 : 2;
              const size_t nops = rank == 1 ? 3 : 6;
              std::vector<std::string> operkeys(nops);

              for (int xyz = 0; xyz < nops; ++xyz) {
                const char* xyz_str[] = { "x", "y", "z", "xx", "xy", "xz", "yy", "yz", "zz" };
                std::ostringstream oss;
                if (oper_key.find("mu_") == 0)
                  oss << "mu_" << xyz_str[xyz];
                if (oper_key.find("q_") == 0)
                  oss << "q_" << xyz_str[xyz+3];
                if (oper_key.find("dphi_") == 0)
                  oss << "dphi_" << xyz_str[xyz];
                if (oper_key.find("ddphi_") == 0)
                  oss << "ddphi_" << xyz_str[xyz+3];
                operkeys[xyz] = ParsedOneBodyIntKey::key(aobra_key, aoket_key,
                                                       oss.str());
              }

              std::vector<RefSCMatrix> intmats(nops);
              const bool do_compute = not registry_->key_exists(operkeys[0]);
              if (do_compute) {
                const Ref<GaussianBasisSet>& obs = basis_;
                Ref<IntParamsOrigin> ref;
                if (oper_key.find("mu_") == 0 || oper_key.find("q_") == 0)
                  ref = new IntParamsOrigin();
                else // electric field (gradient) are on the first nucleus
                  ref = new IntParamsOrigin(bra->basis()->molecule()->r(0));
                if (oper_key.find("mu_") == 0)
                  sc::detail::onebodyint_ao<&Integral::dipole>(bra->basis(), ket->basis(),
                                                               integral(), ref, intmats);
                if (oper_key.find("q_") == 0)
                  sc::detail::onebodyint_ao<&Integral::quadrupole>(bra->basis(), ket->basis(),
                                                                   integral(), ref, intmats);
                if (oper_key.find("dphi_") == 0)
                  sc::detail::onebodyint_ao<&Integral::efield>(bra->basis(), ket->basis(),
                                                               integral(), ref, intmats);
                if (oper_key.find("ddphi_") == 0)
                  sc::detail::onebodyint_ao<&Integral::efield_gradient>(bra->basis(), ket->basis(),
                                                                        integral(), ref, intmats);
                for (int xyz = 0; xyz < nops; ++xyz) {
                  RefSCMatrix ao_blk = bra->coefs().kit()->matrix(bra->coefs().rowdim(),ket->coefs().rowdim());
                  ao_blk->convert( intmats[xyz] );
                  registry_->add(operkeys[xyz], ao_blk);
                }
              }
              else { // already computed
                for (int xyz = 0; xyz < nops; ++xyz)
                  intmats[xyz] = registry_->value(operkeys[xyz]);
              }

            }

            // now all components are available, call itself again
            return get(key);
          }
        }
      }
    }
  }

  // add field contribution to the result, if needed
  if (electric_field() && (oper_key == "H" || oper_key == "F") ) {
    RefSCMatrix result_incl_field = result.copy();
    result_incl_field.accumulate( electric_field_contribution(bra_key, ket_key) );
    return result_incl_field;
  }
  else
    return result;
}

RefSCMatrix
FockBuildRuntime::electric_field_contribution(std::string bra_key,
                                              std::string ket_key) {
  MPQC_ASSERT(electric_field());

  // only AO matrices will be cached

  Ref<OrbitalSpaceRegistry> idxreg = oreg_;
  Ref<OrbitalSpace> bra = idxreg->value(bra_key);
  Ref<OrbitalSpace> ket = idxreg->value(ket_key);
  const Ref<GaussianBasisSet>& bs1 = bra->basis();
  const Ref<GaussianBasisSet>& bs2 = ket->basis();
  const bool bs1_eq_bs2 = bs1->equiv(bs2);

  Ref<AOSpaceRegistry> aoidxreg = aoreg_;
  Ref<OrbitalSpace> aobra = aoidxreg->value(bra->basis());
  Ref<OrbitalSpace> aoket = aoidxreg->value(ket->basis());
  const std::string& aobra_key = idxreg->key(aobra);
  const std::string& aoket_key = idxreg->key(aoket);

  RefSCMatrix H_efield_ao;
  std::vector<std::string> mukeys(3);
  for (int xyz = 0; xyz < 3; ++xyz) {
    const char xyz_char[] = { 'x', 'y', 'z' };
    std::ostringstream oss;
    oss << "Mu_" << xyz_char[xyz];
    mukeys[xyz] = ParsedOneBodyIntKey::key(aobra_key, aoket_key,
                                           oss.str());
  }

  std::vector<RefSCMatrix> Mu(3);
  const bool compute_Mu = not registry_->key_exists(mukeys[0]);
  if (compute_Mu) {
    const Ref<GaussianBasisSet>& obs = basis_;
    Ref<IntParamsOrigin> dipole_origin = new IntParamsOrigin();
    sc::detail::onebodyint_ao<&Integral::dipole>(bs1, bs2, integral(), dipole_origin, Mu);
    for (int xyz = 0; xyz < 3; ++xyz) {
      RefSCMatrix mu_ao_blk = bra->coefs().kit()->matrix(bra->coefs().rowdim(),ket->coefs().rowdim());
      mu_ao_blk->convert( Mu[xyz] );
      registry_->add(mukeys[xyz], mu_ao_blk);
    }
  }
  else { // have_Mu == true
    for (int xyz = 0; xyz < 3; ++xyz)
      Mu[xyz] = registry_->value(mukeys[xyz]);
  }

  // electron charge is not included in Mu
  H_efield_ao = efield_.get_element(0) * Mu[0]
              + efield_.get_element(1) * Mu[1]
              + efield_.get_element(2) * Mu[2];

  // convert H_efield_ao from unblocked to blocked dimensions to be able to multiply with coefs
  RefSCMatrix H_efield_ao_blk = bra->coefs().kit()->matrix(bra->coefs().rowdim(),ket->coefs().rowdim());
  H_efield_ao_blk->convert( H_efield_ao );


  // transform to MO basis
  return bra->coefs().t() * H_efield_ao_blk * ket->coefs();
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
  MPQC_ASSERT(keycopy[0] == '<');
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
