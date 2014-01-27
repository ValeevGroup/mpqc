//
// psicc.cc
//
// Copyright (C) 2007 Edward Valeev
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

#include <cassert>
#include <psifiles.h>
#include <ccfiles.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <cstdio>

#include <chemistry/qc/wfn/orbitalspace.h>
#include <math/mmisc/pairiter.impl.h>
#include <util/misc/print.h>
#include <util/misc/consumableresources.h>
#include <math/distarray4/distarray4_node0file.h>
#include <chemistry/qc/psi/psiwfn.h>
#include <chemistry/qc/psi/psicc.h>

using namespace std;
using namespace sc;

namespace {

  /// convert a std::vector<Source> to a C-type array. Container::size() must be defined
  template <typename Result, typename Source>
  Result*
    to_C(const std::vector<Source>& Cpp) {
      const unsigned int size = Cpp.size();
      Result* C = new Result[size];
      std::copy(Cpp.begin(),Cpp.end(),C);
      return C;
    }

  /// compare 2 2-body tensors
  template <sc::fastpairiter::PairSymm BraSymm, sc::fastpairiter::PairSymm KetSymm>
  void
  compare_tbint_tensors(const RefSCMatrix& T2, const RefSCMatrix& T2_ref,
                        unsigned int nbra1, unsigned int nbra2,
                        unsigned int nket1, unsigned int nket2,
                        double zero);


    void _print(SpinCase2 spin,
                const Ref<DistArray4>& mat,
                const char* label);

}

namespace sc {

  //////////////////////////////////////////////////////////////////////////

  static ClassDesc PsiCC_cd(typeid(PsiCC), "PsiCC", 1,
                            "public PsiCorrWavefunction", 0, create<PsiCC>,
                            create<PsiCC>);

  PsiCC::PsiCC(const Ref<KeyVal>&keyval) :
    PsiCorrWavefunction(keyval) {
    maxiter_ = keyval->intvalue("maxiter",KeyValValueint(default_maxiter));
    diis_ = keyval->booleanvalue("diis",KeyValValueboolean(true));
    diis_nvector_ = keyval->intvalue("diis_nvector",KeyValValueint(true));
  }

  PsiCC::~PsiCC() {
  }

  PsiCC::PsiCC(StateIn&s) :
    PsiCorrWavefunction(s) {
    s.get(maxiter_);
  }

  void PsiCC::save_data_state(StateOut&s) {
    PsiCorrWavefunction::save_data_state(s);
    s.put(maxiter_);
  }

  void PsiCC::dpd_start() {
    const int cachelev = 2;
    int* cachefiles = init_int_array(PSIO_MAXUNIT);
    int** cachelist = init_int_matrix(32,32);

    const Ref<OrbitalSpace>& aocc = ((!compute_1rdm_ || nfzc_ == 0)? occ_act_sb(Alpha)
                                    : occ_sb(Alpha));
//    const Ref<OrbitalSpace>& aocc = occ_act_sb(Alpha);
    const Ref<OrbitalSpace>& avir = vir_act_sb(Alpha);

    //const Ref<OrbitalSpace>& aocc = r12eval()->occ(Alpha);

    int* aoccpi = to_C<int>(aocc->block_sizes());
    int* avirpi = to_C<int>(avir->block_sizes());
    int* aocc_sym = to_C<int>(aocc->orbsym());
    int* avir_sym = to_C<int>(avir->orbsym());

    // UHF/ROHF (ROHF always expected in semicanonical orbitals)
    if (reference()->reftype() != PsiSCF::rhf) {
      const Ref<OrbitalSpace>& bocc = ((!compute_1rdm_ || nfzc_ == 0)? occ_act_sb(Beta)
                                      : occ_sb(Beta));
//        const Ref<OrbitalSpace>& bocc = occ_act_sb(Beta);
      const Ref<OrbitalSpace>& bvir = vir_act_sb(Beta);

      int* boccpi = to_C<int>(bocc->block_sizes());
      int* bvirpi = to_C<int>(bvir->block_sizes());
      int* bocc_sym = to_C<int>(bocc->orbsym());
      int* bvir_sym = to_C<int>(bvir->orbsym());

      dpd_init(0, nirrep(), this->memory_, 0, cachefiles, cachelist,
               NULL, 4, aoccpi, aocc_sym, avirpi, avir_sym,
               boccpi, bocc_sym, bvirpi, bvir_sym);
    }
    // RHF
    else {
      dpd_init(0, nirrep(), this->memory_, 0, cachefiles, cachelist,
               NULL, 2, aoccpi, aocc_sym, avirpi, avir_sym);
    }

    psi::PSIO& psio = exenv()->psio();
    for(int i=CC_MIN; i <= CC_MAX; i++) psio.open(i,1);
}

  void PsiCC::dpd_stop() {
    psi::PSIO& psio = exenv()->psio();
    for(int i=CC_MIN; i <= CC_MAX; i++) psio.close(i,1);
    dpd_close(0);
  }

  const RefSCMatrix&PsiCC::T1(SpinCase1 spin1) {
    if (T1_[spin1])
      return T1_[spin1];
    PsiSCF::RefType reftype = reference_->reftype();

    // Grab T matrices
    const char* kwd = (spin1 == Beta && reftype != PsiSCF::rhf) ? "tia" : "tIA";
    // When computing one-electron density, we use non frozen core for CCSD,
    // but only need the frozen core portion of the T1 amplitude
    T1_[spin1] = ((!compute_1rdm_ || nfzc_ == 0)? T1(spin1, kwd) : T1_fzc(spin1, kwd));
//    T1_[spin1] = T1(spin1, kwd);
    if (debug() >= DefaultPrintThresholds::mostN2)
      T1_[spin1].print(prepend_spincase(spin1,"T1 amplitudes").c_str());

    return T1_[spin1];
  }

  const RefSCMatrix& PsiCC::T2(SpinCase2 spin2) {
    if (T2_[spin2])
      return T2_[spin2];
    PsiSCF::RefType reftype = reference_->reftype();

    // If requesting AA or BB T2s, create their forms with unrestricted indices and read those
    if (spin2 != AlphaBeta) {
      dpd_start();

      if (spin2 == AlphaAlpha) {
        dpdbuf4 D;
        dpd_buf4_init(&D, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
        dpd_buf4_copy(&D, CC_TAMPS, "tIJAB (IJ,AB)");
        dpd_buf4_close(&D);
      }
      if (spin2 == BetaBeta) {
        dpdbuf4 D;
        dpd_buf4_init(&D, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
        dpd_buf4_copy(&D, CC_TAMPS, "tijab (ij,ab)");
        dpd_buf4_close(&D);
      }
      dpd_stop();
    }

    // Grab T matrices
    const char* kwd = (spin2 != AlphaBeta && reftype != PsiSCF::rhf) ? (spin2 == AlphaAlpha ? "tIJAB (IJ,AB)" : "tijab (ij,ab)") : "tIjAb";
    T2_[spin2] = T2(spin2, kwd);
    if (debug() >= DefaultPrintThresholds::mostO2N2)
      T2_[spin2].print(prepend_spincase(spin2,"T2 amplitudes").c_str());

    return T2_[spin2];
  }

  Ref<DistArray4> PsiCC::T2_distarray4(SpinCase2 spin2) {
    if (T2_da4_[spin2])
      return T2_da4_[spin2];
    PsiSCF::RefType reftype = reference_->reftype();

    // If requesting AA or BB T2s, create their forms with unrestricted indices and read those
    if (spin2 != AlphaBeta) {
      dpd_start();

      if (spin2 == AlphaAlpha) {
        dpdbuf4 D;
        dpd_buf4_init(&D, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
        dpd_buf4_copy(&D, CC_TAMPS, "tIJAB (IJ,AB)");
        dpd_buf4_close(&D);
      }
      if (spin2 == BetaBeta) {
        dpdbuf4 D;
        dpd_buf4_init(&D, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
        dpd_buf4_copy(&D, CC_TAMPS, "tijab (ij,ab)");
        dpd_buf4_close(&D);
      }
      dpd_stop();
    }

    // Grab T matrices
    const char* kwd = (spin2 != AlphaBeta) ? (spin2 == AlphaAlpha ? "tIJAB (IJ,AB)" : "tijab (ij,ab)") : "tIjAb";
    // When computing one-electron density, we use non frozen core for CCSD,
    // but only need the frozen core portion of the amplitude
    T2_da4_[spin2] = ((!compute_1rdm_ || nfzc_ == 0)? T2_distarray4(spin2, kwd)
                                    : T2_distarray4_fzc(spin2, kwd));
//    T2_da4_[spin2] = T2_distarray4(spin2, kwd);

    return T2_da4_[spin2];
  }

  const RefSCMatrix&PsiCC::Tau2(SpinCase2 spin2) {
    if (Tau2_[spin2])
      return Tau2_[spin2];
    PsiSCF::RefType reftype = reference_->reftype();

    // If requesting AA or BB Tau2s, create their forms with unrestricted indices and read those
    if (spin2 != AlphaBeta) {
      dpd_start();

      if (spin2 == AlphaAlpha) {
        dpdbuf4 D;
        dpd_buf4_init(&D, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tauIJAB");
        dpd_buf4_copy(&D, CC_TAMPS, "tauIJAB (IJ,AB)");
        dpd_buf4_close(&D);
      }
      if (spin2 == BetaBeta) {
        dpdbuf4 D;
        dpd_buf4_init(&D, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tauijab");
        dpd_buf4_copy(&D, CC_TAMPS, "tauijab (ij,ab)");
        dpd_buf4_close(&D);
      }
      dpd_stop();
    }

    const char* kwd = (spin2 != AlphaBeta && reftype != PsiSCF::rhf) ? (spin2 == AlphaAlpha ? "tauIJAB (IJ,AB)" : "tauijab (ij,ab)") : "tauIjAb";
    Tau2_[spin2] = T2(spin2, kwd);
    if (debug() >= DefaultPrintThresholds::mostO2N2)
      Tau2_[spin2].print(prepend_spincase(spin2,"Tau2 amplitudes").c_str());
    return Tau2_[spin2];
  }

  RefSCMatrix PsiCC::T1(SpinCase1 spin, const std::string& dpdlabel) {
    psi::PSIO& psio = exenv()->psio();
    // grab orbital info
    const std::vector<unsigned int>& occpi = reference()->occpi(spin);
    const std::vector<unsigned int>& uoccpi = reference()->uoccpi(spin);
    std::vector<unsigned int> actoccpi(nirrep_);
    std::vector<unsigned int> actuoccpi(nirrep_);
    std::vector<unsigned int> actoccioff(nirrep_);
    std::vector<unsigned int> actuoccioff(nirrep_);
    const std::vector<unsigned int>& frozen_docc = this->frozen_docc();
    const std::vector<unsigned int>& frozen_uocc = this->frozen_uocc();
    unsigned int nocc_act = 0;
    unsigned int nuocc_act = 0;
    for (unsigned int irrep=0; irrep<nirrep_; ++irrep) {
      actoccpi[irrep] = occpi[irrep] - frozen_docc_[irrep];
      actuoccpi[irrep] = uoccpi[irrep] - frozen_uocc_[irrep];
      nocc_act += actoccpi[irrep];
      nuocc_act += actuoccpi[irrep];
    }
    actoccioff[0] = 0;
    actuoccioff[0] = 0;
    for (unsigned int irrep=1; irrep<nirrep_; ++irrep) {
      actoccioff[irrep] = actoccioff[irrep-1] + actoccpi[irrep-1];
      actuoccioff[irrep] = actuoccioff[irrep-1] + actuoccpi[irrep-1];
    }

    RefSCDimension rowdim =  new SCDimension(nocc_act);
    //rowdim->blocks()->set_subdim(0,new SCDimension(rowdim.n()));
    RefSCDimension coldim = new SCDimension(nuocc_act);
    //coldim->blocks()->set_subdim(0,new SCDimension(coldim.n()));
    RefSCMatrix T = matrixkit()->matrix(rowdim, coldim);
    T.assign(0.0);
    // test that T1 is not empty
    if (rowdim.n() && coldim.n()) {
      // read in the i by a matrix in DPD format
      unsigned int nia_dpd = 0;
      for (unsigned int h=0; h<nirrep_; ++h)
      nia_dpd += actoccpi[h] * actuoccpi[h];

      double* T1 = new double[nia_dpd];
      psio.open(CC_OEI, PSIO_OPEN_OLD);
      psio.read_entry(CC_OEI, const_cast<char*>(dpdlabel.c_str()),
                      reinterpret_cast<char*>(T1), nia_dpd*sizeof(double));
      psio.close(CC_OEI, 1);

      // form the full matrix
      unsigned int ia = 0;
      for (unsigned int h=0; h<nirrep_; ++h) {
        const unsigned int i_offset = actoccioff[h];
        const unsigned int a_offset = actuoccioff[h];
        for (int i=0; i<actoccpi[h]; ++i)
          for (int a=0; a<actuoccpi[h]; ++a, ++ia)
            T.set_element(i+i_offset, a+a_offset, T1[ia]);
      }
      delete[] T1;
    }

    return T;
  }

  // When compute_1rdm = true, get the T1 amplitudes (frozen core size)
  // from non frozen core PSI CCSD density calculation
  RefSCMatrix PsiCC::T1_fzc(SpinCase1 spin, const std::string& dpdlabel) {
    psi::PSIO& psio = exenv()->psio();
    // grab orbital info
    const std::vector<unsigned int>& occpi = reference()->occpi(spin);
    const std::vector<unsigned int>& uoccpi = reference()->uoccpi(spin);
    const std::vector<unsigned int>& frozen_docc = this->frozen_docc();
    const std::vector<unsigned int>& frozen_uocc = this->frozen_uocc();

    std::vector<unsigned int> actoccpi(nirrep_);
    std::vector<unsigned int> actuoccpi(nirrep_);
    unsigned int nocc_act = 0;
    unsigned int nuocc_act = 0;
    unsigned int nocc = 0;
    for (unsigned int irrep = 0; irrep < nirrep_; ++irrep) {
      actoccpi[irrep] = occpi[irrep] - frozen_docc_[irrep];
      actuoccpi[irrep] = uoccpi[irrep] - frozen_uocc_[irrep];
      nocc_act += actoccpi[irrep];
      nuocc_act += actuoccpi[irrep];
      nocc += occpi[irrep];
    }

    std::vector<unsigned int> actoccioff(nirrep_);
    std::vector<unsigned int> actuoccioff(nirrep_);
    std::vector<unsigned int> occioff(nirrep_);
    actoccioff[0] = 0;
    actuoccioff[0] = 0;
    occioff[0] = 0;
    for (unsigned int irrep = 1; irrep < nirrep_; ++irrep) {
      actoccioff[irrep] = actoccioff[irrep-1] + actoccpi[irrep-1];
      actuoccioff[irrep] = actuoccioff[irrep-1] + actuoccpi[irrep-1];
      occioff[irrep] = occioff[irrep-1] + occpi[irrep-1];
    }

    RefSCDimension rowdim_occ_act = new SCDimension(nocc_act);
    RefSCDimension rowdim_occ = new SCDimension(nocc);
    RefSCDimension coldim = new SCDimension(nuocc_act);

    // T1 amplitude matrix only exclude the frozen core portion
    RefSCMatrix T = matrixkit()->matrix(rowdim_occ_act, coldim);
    T.assign(0.0);

    // T1 amplitude matrix from the non frozen core PSI CCSD calculation
    RefSCMatrix T_nfzc = matrixkit()->matrix(rowdim_occ, coldim);
    T_nfzc.assign(0.0);

    // test that T1 is not empty
    if (rowdim_occ_act.n() && coldim.n()) {
      // read in the i by a matrix in DPD format
      unsigned int nia_dpd = 0;
      for (unsigned int h=0; h<nirrep_; ++h)
        nia_dpd += occpi[h] * actuoccpi[h];

      double* T1 = new double[nia_dpd];
      psio.open(CC_OEI, PSIO_OPEN_OLD);
      psio.read_entry(CC_OEI, const_cast<char*>(dpdlabel.c_str()),
                      reinterpret_cast<char*>(T1), nia_dpd*sizeof(double));
      psio.close(CC_OEI, 1);

      // form the full matrix
      unsigned int ia = 0;
      for (unsigned int h = 0; h < nirrep_; ++h) {
        const unsigned int i_offset = occioff[h];
        const unsigned int a_offset = actuoccioff[h];

        for (int i = 0; i < occpi[h]; ++i) {
          for (int a = 0; a < actuoccpi[h]; ++a, ++ia) {
            T_nfzc.set_element(i+i_offset, a+a_offset, T1[ia]);
          }
        }
      }
      delete[] T1;

      for (unsigned int h = 0; h < nirrep_; ++h) {
        const unsigned int i_nfzc_offset = occioff[h] + frozen_docc_[h];
        const unsigned int i_offset = actoccioff[h];
        const unsigned int a_offset = actuoccioff[h];

        for (int i = 0; i < actoccpi[h]; ++i) {
          for (int a = 0; a < actuoccpi[h]; ++a) {
            const double Tia = T_nfzc.get_element(i+i_nfzc_offset, a+a_offset);
            T.set_element(i+i_offset, a+a_offset, Tia);
          }
        }
      }
      T_nfzc = 0;

    }

    //T.print(prepend_spincase(spin,"CCSD T1 amplitudes:").c_str());
    return T;
  }

  namespace {

    template <typename T>
    T packed_2index_anti(T i1, T i2) {
      T max_i = std::max(i1, i2);
      T min_i = std::min(i1, i2);
      return max_i * (max_i - 1)/2 + min_i;
    }
  }

  RefSCMatrix PsiCC::T2(SpinCase2 spin12, const std::string& dpdlabel) {
    psi::PSIO& psio = exenv()->psio();
    const SpinCase1 spin1 = case1(spin12);
    const SpinCase1 spin2 = case2(spin12);

    typedef std::vector<unsigned int> uvec;
    // grab orbital info
    const uvec& occpi1 = reference()->occpi(spin1);
    const uvec& occpi2 = reference()->occpi(spin2);
    const uvec& uoccpi1 = reference()->uoccpi(spin1);
    const uvec& uoccpi2 = reference()->uoccpi(spin2);
    uvec actoccpi1(nirrep_);
    uvec actuoccpi1(nirrep_);
    uvec actoccioff1(nirrep_);
    uvec actuoccioff1(nirrep_);
    uvec actoccpi2(nirrep_);
    uvec actuoccpi2(nirrep_);
    uvec actoccioff2(nirrep_);
    uvec actuoccioff2(nirrep_);
    unsigned int nocc1_act = 0;
    unsigned int nuocc1_act = 0;
    unsigned int nocc2_act = 0;
    unsigned int nuocc2_act = 0;
    for (unsigned int irrep=0; irrep<nirrep_; ++irrep) {
      actoccpi1[irrep] = occpi1[irrep] - frozen_docc_[irrep];
      actuoccpi1[irrep] = uoccpi1[irrep] - frozen_uocc_[irrep];
      nocc1_act += actoccpi1[irrep];
      nuocc1_act += actuoccpi1[irrep];
      actoccpi2[irrep] = occpi2[irrep] - frozen_docc_[irrep];
      actuoccpi2[irrep] = uoccpi2[irrep] - frozen_uocc_[irrep];
      nocc2_act += actoccpi2[irrep];
      nuocc2_act += actuoccpi2[irrep];
    }
    actoccioff1[0] = 0;
    actuoccioff1[0] = 0;
    actoccioff2[0] = 0;
    actuoccioff2[0] = 0;
    for (unsigned int irrep=1; irrep<nirrep_; ++irrep) {
      actoccioff1[irrep] = actoccioff1[irrep-1] + actoccpi1[irrep-1];
      actuoccioff1[irrep] = actuoccioff1[irrep-1] + actuoccpi1[irrep-1];
      actoccioff2[irrep] = actoccioff2[irrep-1] + actoccpi2[irrep-1];
      actuoccioff2[irrep] = actuoccioff2[irrep-1] + actuoccpi2[irrep-1];
    }

    // DPD of orbital product spaces
    std::vector<size_t> ijpi(nirrep_);
    std::vector<size_t> abpi(nirrep_);
    size_t nijab_dpd = 0;
    for (unsigned int h=0; h<nirrep_; ++h) {
      size_t nij = 0;
      size_t nab = 0;
      for (unsigned int g=0; g<nirrep_; ++g) {
        nij += (size_t)actoccpi1[g] * actoccpi2[h^g];
        nab += (size_t)actuoccpi1[g] * actuoccpi2[h^g];
      }
      ijpi[h] = nij;
      abpi[h] = nab;
      nijab_dpd += nij*nab;
    }

    const size_t nij = (spin12 == AlphaBeta) ? nocc1_act*nocc2_act : nocc1_act*(nocc1_act-1)/2;
    const size_t nab = (spin12 == AlphaBeta) ? nuocc1_act*nuocc2_act : nuocc1_act*(nuocc1_act-1)/2;
    RefSCDimension rowdim = new SCDimension(nij);
    //rowdim->blocks()->set_subdim(0,new SCDimension(rowdim.n()));
    RefSCDimension coldim = new SCDimension(nab);
    //coldim->blocks()->set_subdim(0,new SCDimension(coldim.n()));
    RefSCMatrix T = matrixkit()->matrix(rowdim, coldim);
    T.assign(0.0);

    // If not empty...
    if (nijab_dpd) {
    // read in T2 in DPD form
    double* T2 = new double[nijab_dpd];
    psio.open(CC_TAMPS, PSIO_OPEN_OLD);
    psio.read_entry(CC_TAMPS, const_cast<char*>(dpdlabel.c_str()),
                    reinterpret_cast<char*>(T2), nijab_dpd*sizeof(double));
    psio.close(CC_TAMPS, 1);

    // convert to the full form
    size_t ijab = 0;
    size_t ij_offset = 0;
    size_t ab_offset = 0;
    for (unsigned int h=0; h<nirrep_; ij_offset+=ijpi[h], ab_offset+=abpi[h],
                                      ++h) {
      for (unsigned int g=0; g<nirrep_; ++g) {
        unsigned int gh = g^h;
        for (int i=0; i<actoccpi1[g]; ++i) {
          const size_t ii = i + actoccioff1[g];

          for (int j=0; j<actoccpi2[gh]; ++j) {
            const size_t jj = j + actoccioff2[gh];
            if (spin12 != AlphaBeta && ii == jj) continue;
            const size_t ij = (spin12 == AlphaBeta) ? ii * nocc2_act + jj
                                                          : packed_2index_anti(ii, jj);

            for (unsigned int f=0; f<nirrep_; ++f) {
              unsigned int fh = f^h;
              for (int a=0; a<actuoccpi1[f]; ++a) {
                const size_t aa = a + actuoccioff1[f];

                for (int b=0; b<actuoccpi2[fh]; ++b, ++ijab) {
                  const size_t bb = b + actuoccioff2[fh];
                  if (spin12 != AlphaBeta && aa == bb) continue;
                  const size_t ab = (spin12 == AlphaBeta) ? aa*nuocc2_act + bb
                                                          : packed_2index_anti(aa, bb);

                  T.set_element(ij, ab, T2[ijab]);
                }
              }
            }
          }
        }
      }
    }
    delete[] T2;
    }

    return T;
  }

  Ref<DistArray4> PsiCC::T2_distarray4(SpinCase2 spin12,
                                       const std::string& dpdlabel) {
    psi::PSIO& psio = exenv()->psio();
    const SpinCase1 spin1 = case1(spin12);
    const SpinCase1 spin2 = case2(spin12);

    typedef std::vector<unsigned int> uvec;
    // grab orbital info
    const uvec& occpi1 = reference()->occpi(spin1);
    const uvec& occpi2 = reference()->occpi(spin2);
    const uvec& uoccpi1 = reference()->uoccpi(spin1);
    const uvec& uoccpi2 = reference()->uoccpi(spin2);
    uvec actoccpi1(nirrep_);
    uvec actuoccpi1(nirrep_);
    uvec actoccioff1(nirrep_);
    uvec actuoccioff1(nirrep_);
    uvec actoccpi2(nirrep_);
    uvec actuoccpi2(nirrep_);
    uvec actoccioff2(nirrep_);
    uvec actuoccioff2(nirrep_);
    unsigned int nocc1_act = 0;
    unsigned int nuocc1_act = 0;
    unsigned int nocc2_act = 0;
    unsigned int nuocc2_act = 0;
    for (unsigned int irrep = 0; irrep < nirrep_; ++irrep) {
      actoccpi1[irrep] = occpi1[irrep] - frozen_docc_[irrep];
      actuoccpi1[irrep] = uoccpi1[irrep] - frozen_uocc_[irrep];
      nocc1_act += actoccpi1[irrep];
      nuocc1_act += actuoccpi1[irrep];
      actoccpi2[irrep] = occpi2[irrep] - frozen_docc_[irrep];
      actuoccpi2[irrep] = uoccpi2[irrep] - frozen_uocc_[irrep];
      nocc2_act += actoccpi2[irrep];
      nuocc2_act += actuoccpi2[irrep];
    }
    actoccioff1[0] = 0;
    actuoccioff1[0] = 0;
    actoccioff2[0] = 0;
    actuoccioff2[0] = 0;
    for (unsigned int irrep = 1; irrep < nirrep_; ++irrep) {
      actoccioff1[irrep] = actoccioff1[irrep - 1] + actoccpi1[irrep - 1];
      actuoccioff1[irrep] = actuoccioff1[irrep - 1] + actuoccpi1[irrep - 1];
      actoccioff2[irrep] = actoccioff2[irrep - 1] + actoccpi2[irrep - 1];
      actuoccioff2[irrep] = actuoccioff2[irrep - 1] + actuoccpi2[irrep - 1];
    }

    // DPD of orbital product spaces
    std::vector<size_t> ijpi(nirrep_);
    std::vector<size_t> abpi(nirrep_);
    size_t nijab_dpd = 0;
    for (unsigned int h = 0; h < nirrep_; ++h) {
      size_t nij = 0;
      size_t nab = 0;
      for (unsigned int g = 0; g < nirrep_; ++g) {
        nij += (size_t)actoccpi1[g] * actoccpi2[h ^ g];
        nab += (size_t)actuoccpi1[g] * actuoccpi2[h ^ g];
      }
      ijpi[h] = nij;
      abpi[h] = nab;
      nijab_dpd += nij * nab;
    }
    const size_t max_nab = *std::max_element(abpi.begin(), abpi.end());

    const size_t nij = nocc1_act * nocc2_act;
    const size_t nab = nuocc1_act * nuocc2_act;
    RefSCDimension rowdim = new SCDimension(nij);
    RefSCDimension coldim = new SCDimension(nab);
    Ref<DistArray4> T;
    // TODO make this work for non-disk-based storage
    {
      const char* spin12_label = (spin12 == AlphaAlpha) ? "aa" : ((spin12 == BetaBeta) ? "bb" : "ab");
      std::string fileext_str(".T2_distarray4_"); fileext_str += spin12_label;
      const std::string default_basename_prefix = SCFormIO::fileext_to_filename(fileext_str.c_str());
      const std::string da4_filename = ConsumableResources::get_default_instance()->disk_location() +
          default_basename_prefix;
      T = new DistArray4_Node0File(da4_filename.c_str(), 1, nocc1_act, nocc2_act,
                                   nuocc1_act, nuocc2_act);
    }

    // If not empty...
    if (nijab_dpd) {

      // read in T2 one ij at a time
      psio.open(CC_TAMPS, PSIO_OPEN_OLD);

      double* t2_ij = new double[max_nab];
      double* T_ij = new double[nab];

      T->activate();

      // convert to the full form
      size_t ij_offset = 0;
      size_t ab_offset = 0;
      psio_address t2_ij_address = PSIO_ZERO;
      for (unsigned int h = 0; h < nirrep_; ij_offset += ijpi[h], ab_offset
          += abpi[h], ++h) {
        for (unsigned int g = 0; g < nirrep_; ++g) {
          unsigned int gh = g ^ h;
          for (int i = 0; i < actoccpi1[g]; ++i) {
            const size_t ii = i + actoccioff1[g];

            for (int j = 0; j < actoccpi2[gh]; ++j) {
              const size_t jj = j + actoccioff2[gh];

              std::fill(T_ij, T_ij+nab, 0.0);

              const size_t nab_h = abpi[h];
              psio.read(CC_TAMPS, const_cast<char*> (dpdlabel.c_str()),
                        reinterpret_cast<char*> (t2_ij), nab_h * sizeof(double),
                        t2_ij_address, &t2_ij_address);

              size_t ab_h = 0;
              for (unsigned int f = 0; f < nirrep_; ++f) {
                unsigned int fh = f ^ h;
                for (int a = 0; a < actuoccpi1[f]; ++a) {
                  const size_t aa = a + actuoccioff1[f];

                  for (int b = 0; b < actuoccpi2[fh]; ++b, ++ab_h) {
                    const size_t bb = b + actuoccioff2[fh];

                    const size_t ab = aa * nuocc2_act + bb;
                    T_ij[ab] = t2_ij[ab_h];
                  }
                }
              }

              T->store_pair_block(ii, jj, 0, T_ij);
            }
          }
        }
      }

      if (T->data_persistent()) T->deactivate();

      delete[] t2_ij;
      delete[] T_ij;

      psio.close(CC_TAMPS, 1);

    }

    return T;
  }

  // When compute_1rdm = true, get the T2 amplitudes (frozen core size)
  // from non frozen core PSI CCSD density calculation
  Ref<DistArray4> PsiCC::T2_distarray4_fzc(SpinCase2 spin12,
                                       const std::string& dpdlabel) {
    psi::PSIO& psio = exenv()->psio();

    const SpinCase1 spin1 = case1(spin12);
    const SpinCase1 spin2 = case2(spin12);

    // grab orbital info
    typedef std::vector<unsigned int> uvec;
    const uvec& occpi1 = reference()->occpi(spin1);
    const uvec& occpi2 = reference()->occpi(spin2);
    const uvec& uoccpi1 = reference()->uoccpi(spin1);
    const uvec& uoccpi2 = reference()->uoccpi(spin2);

    uvec actoccpi1(nirrep_);
    uvec actuoccpi1(nirrep_);
    uvec actoccpi2(nirrep_);
    uvec actuoccpi2(nirrep_);
    unsigned int nocc1_act = 0;
    unsigned int nuocc1_act = 0;
    unsigned int nocc2_act = 0;
    unsigned int nuocc2_act = 0;

    unsigned int nocc1 = 0;
    unsigned int nocc2 = 0;
    for (unsigned int irrep = 0; irrep < nirrep_; ++irrep) {
      actoccpi1[irrep] = occpi1[irrep] - frozen_docc_[irrep];
      actuoccpi1[irrep] = uoccpi1[irrep] - frozen_uocc_[irrep];
      actoccpi2[irrep] = occpi2[irrep] - frozen_docc_[irrep];
      actuoccpi2[irrep] = uoccpi2[irrep] - frozen_uocc_[irrep];
      nocc1_act += actoccpi1[irrep];
      nuocc1_act += actuoccpi1[irrep];
      nocc2_act += actoccpi2[irrep];
      nuocc2_act += actuoccpi2[irrep];

      nocc1 += occpi1[irrep];
      nocc2 += occpi2[irrep];
    }

    uvec actoccioff1(nirrep_);
    uvec actuoccioff1(nirrep_);
    uvec actoccioff2(nirrep_);
    uvec actuoccioff2(nirrep_);

    uvec occioff1(nirrep_);
    uvec occioff2(nirrep_);

    actoccioff1[0] = 0;
    actuoccioff1[0] = 0;
    actoccioff2[0] = 0;
    actuoccioff2[0] = 0;

    occioff1[0] = 0;
    occioff2[0] = 0;
    for (unsigned int irrep = 1; irrep < nirrep_; ++irrep) {
      actoccioff1[irrep] = actoccioff1[irrep - 1] + actoccpi1[irrep - 1];
      actuoccioff1[irrep] = actuoccioff1[irrep - 1] + actuoccpi1[irrep - 1];
      actoccioff2[irrep] = actoccioff2[irrep - 1] + actoccpi2[irrep - 1];
      actuoccioff2[irrep] = actuoccioff2[irrep - 1] + actuoccpi2[irrep - 1];

      occioff1[irrep] = occioff1[irrep - 1] + occpi1[irrep - 1];
      occioff2[irrep] = occioff2[irrep - 1] + occpi2[irrep - 1];
    }

    // DPD of orbital product spaces
    std::vector<size_t> ijpi(nirrep_);
    std::vector<size_t> abpi(nirrep_);
    size_t nijab_dpd = 0;
    for (unsigned int h = 0; h < nirrep_; ++h) {
      size_t nij = 0;
      size_t nab = 0;
      for (unsigned int g = 0; g < nirrep_; ++g) {
        nij += (size_t)occpi1[g] * occpi2[h ^ g];
        nab += (size_t)actuoccpi1[g] * actuoccpi2[h ^ g];
      }
      ijpi[h] = nij;
      abpi[h] = nab;
      nijab_dpd += nij * nab;
    }

    Ref<DistArray4> T;       // T2 amplitude from non frozen core PSI CCSD
    // TODO make this work for non-disk-based storage
    {
      const char* spin12_label = (spin12 == AlphaAlpha) ? "aa" : ((spin12 == BetaBeta) ? "bb" : "ab");
      std::string fileext_str(".T2_distarray4_"); fileext_str += spin12_label;
      const std::string default_basename_prefix = SCFormIO::fileext_to_filename(fileext_str.c_str());
      const std::string da4_filename = ConsumableResources::get_default_instance()->disk_location() +
          default_basename_prefix;
      T = new DistArray4_Node0File(da4_filename.c_str(), 1, nocc1, nocc2,
                                   nuocc1_act, nuocc2_act);
    }

    Ref<DistArray4> T_fzc;   // T2 amplitude excluding the frozen core part
    {
      const char* spin12_label = (spin12 == AlphaAlpha) ? "aa" : ((spin12 == BetaBeta) ? "bb" : "ab");
      std::string fileext_str(".T2_fzc_distarray4_"); fileext_str += spin12_label;
      const std::string default_basename_prefix = SCFormIO::fileext_to_filename(fileext_str.c_str());
      const std::string da4_filename = ConsumableResources::get_default_instance()->disk_location() +
          default_basename_prefix;
      T_fzc = new DistArray4_Node0File(da4_filename.c_str(), 1,
                                       nocc1_act, nocc2_act,
                                       nuocc1_act, nuocc2_act);
    }

    // If not empty...
    if (nijab_dpd) {
      const size_t max_nab = *std::max_element(abpi.begin(), abpi.end());
      const size_t nij = nocc1 * nocc2;
      const size_t nab = nuocc1_act * nuocc2_act;

      // read in T2 one ij at a time
      psio.open(CC_TAMPS, PSIO_OPEN_OLD);

      double* t2_ij = new double[max_nab];
      double* T_ij = new double[nab];

      T->activate();

      // convert to the full form
      size_t ij_offset = 0;
      size_t ab_offset = 0;
      psio_address t2_ij_address = PSIO_ZERO;
      for (unsigned int h = 0; h < nirrep_; ij_offset += ijpi[h], ab_offset
          += abpi[h], ++h) {
        for (unsigned int g = 0; g < nirrep_; ++g) {
          unsigned int gh = g ^ h;
          for (int i = 0; i < occpi1[g]; ++i) {
            const size_t ii = i + occioff1[g];

            for (int j = 0; j < occpi2[gh]; ++j) {
              const size_t jj = j + occioff2[gh];

              std::fill(T_ij, T_ij+nab, 0.0);

              const size_t nab_h = abpi[h];
              psio.read(CC_TAMPS, const_cast<char*> (dpdlabel.c_str()),
                        reinterpret_cast<char*> (t2_ij), nab_h * sizeof(double),
                        t2_ij_address, &t2_ij_address);

              size_t ab_h = 0;
              for (unsigned int f = 0; f < nirrep_; ++f) {
                unsigned int fh = f ^ h;
                for (int a = 0; a < actuoccpi1[f]; ++a) {
                  const size_t aa = a + actuoccioff1[f];

                  for (int b = 0; b < actuoccpi2[fh]; ++b, ++ab_h) {
                    const size_t bb = b + actuoccioff2[fh];

                    const size_t ab = aa * nuocc2_act + bb;
                    T_ij[ab] = t2_ij[ab_h];
                  }
                }
              }

              T->store_pair_block(ii, jj, 0, T_ij);
            }
          }
        }
      }
      psio.close(CC_TAMPS, 1);

      delete[] t2_ij;
      delete[] T_ij;

      T_fzc->activate();

      for (unsigned int h = 0; h < nirrep_; ++h) {
        for (unsigned int g = 0; g < nirrep_; ++g) {
          unsigned int gh = g ^ h;
          const size_t nfzcpi_g = frozen_docc_[g];
          const size_t nfzcpi_gh = frozen_docc_[gh];

          for (int i = 0; i < actoccpi1[g]; ++i) {
            const size_t ii = i + occioff1[g] + nfzcpi_g;
            const size_t ii_fzc = i + actoccioff1[g];

            for (int j = 0; j < actoccpi2[gh]; ++j) {
              const size_t jj = j + occioff2[gh] + nfzcpi_gh;
              const size_t jj_fzc = j + actoccioff2[gh];

              const double* const ij_blk = T->retrieve_pair_block(ii, jj, 0);
              T_fzc->store_pair_block(ii_fzc, jj_fzc, 0, ij_blk);

              T->release_pair_block(ii, jj, 0);
            }
          }
        }
      }

      if (T->data_persistent()) {
        T->deactivate();
        T_fzc->deactivate();
      }

    }
    //_print(spin12, T_fzc, prepend_spincase(spin12,"CCSD T2 amplitudes:").c_str());
    return T_fzc;
  }

  const RefSCMatrix&PsiCC::Lambda1(SpinCase1 spin) {
    if (Lambda1_[spin])
      return Lambda1_[spin];

    throw FeatureNotImplemented("PsiCC::Lambda1() -- cannot read Lambda1 amplitudes yet",__FILE__,__LINE__);
    return Lambda1_[spin];
  }

  const RefSCMatrix&PsiCC::Lambda2(SpinCase2 spin) {
    if (Lambda2_[spin])
      return Lambda2_[spin];

    throw FeatureNotImplemented("PsiCC::Lambda2() -- cannot read Lambda2 amplitudes yet",__FILE__,__LINE__);
    return Lambda2_[spin];
  }

  // import the psi ccsd one-particle density
  RefSCMatrix PsiCC::Onerdm(SpinCase1 spin) {

    // grab orbital info
    psi::PSIO& psio = exenv()->psio();

    // psi ccsd one-particle density does not work for frozen core
    // this code does not work for ROHF
//    MPQC_ASSERT(nfzc_ == 0 && reference()->reftype() != PsiSCF::hsoshf);
    MPQC_ASSERT(reference()->reftype() != PsiSCF::hsoshf);

    // get # of occupied and unoccupied orbitals of spin S per irrep
    const std::vector<unsigned int>& occpi = reference()->occpi(spin);
    const std::vector<unsigned int>& uoccpi = reference()->uoccpi(spin);
//    std::vector<unsigned int> actoccpi(nirrep_);

    // obtain # of orbitals per irrep
    unsigned int nocc = 0;
//    unsigned int nocc_act = 0;
    unsigned int nuocc = 0;
    for (unsigned int irrep = 0; irrep < nirrep_; ++irrep) {
//      actoccpi[irrep] = occpi[irrep] - frozen_docc_[irrep];
//      nocc_act += actoccpi[irrep];
      nocc += occpi[irrep];
      nuocc += uoccpi[irrep];
    }
    unsigned int norb = nocc + nuocc;
    //ExEnv::out0() << indent << spin << " nocc_act:" << nocc_act << " nuocc:"<< nuocc <<endl;

    // Equation for the one-particle density matrix
    // D_ij = -1/2 t_im^ef L^jm_ef - t_i^e L^j_e
    // D_ab = 1/2 L^mn_ae t_mn^be + L^m_a t_m^b
    // D_ia = t_i^a + (t_im^ae - t_i^e t_m^a) L^m_e
    //        - 1/2 L^mn_ef (t_in^ef t_m^a + t_i^e t_mn^af)
    // D_ai = L^i_a
    RefSCDimension rowdim = new SCDimension(norb);
    RefSCDimension coldim = new SCDimension(norb);
    RefSCMatrix D = matrixkit()->matrix(rowdim, coldim);
    D.assign(0.0);

    // test that D is not empty
    if (rowdim.n() && coldim.n()) {
      // read in the p by q matrix in DPD format
      unsigned int npq_dpd = norb * norb;
      double* Dpq = new double[npq_dpd];
      psio.open(PSIF_MO_OPDM, PSIO_OPEN_OLD);

      if (reference()->reftype() == PsiSCF::rhf) {
        // closed-shell: O[I][J] += 2.0 * D.matrix[h][i][j]
        //               O[A][B] += 2.0 * D.matrix[h][a][b]
        //               O[A][I] += 2.0 * D.matrix[h][i][a]
        //               O[I][A] += 2.0 * D.matrix[h][i][a]

        psio.read_entry(PSIF_MO_OPDM, "MO-basis OPDM",
                        reinterpret_cast<char*>(Dpq), npq_dpd * sizeof(double));
        psio.close(PSIF_MO_OPDM, 1);

        unsigned int pq = 0;
        for (int p = 0; p < norb; ++p)
          for (int q = 0; q < norb; ++q, ++pq)
            D.set_element(p, q, 0.5*Dpq[pq]);
      }
      else if (reference()->reftype() == PsiSCF::uhf) {

        const char* Dpq_lbl = (spin == Alpha) ? "MO-basis Alpha OPDM" : "MO-basis Beta OPDM";
        psio.read_entry(PSIF_MO_OPDM, Dpq_lbl,
                        reinterpret_cast<char*>(Dpq), npq_dpd * sizeof(double));
        psio.close(PSIF_MO_OPDM, 1);

        unsigned int pq = 0;
        for (int p = 0; p < norb; ++p)
          for (int q = 0; q < norb; ++q, ++pq)
            D.set_element(p, q, Dpq[pq]);
      }

      delete[] Dpq;
    }
    // end of if (rowdim.n() && coldim.n())

    return D;
  }
  // end of import psi CCSD one-particle density

  // import PSI3 orbital Z-vector for one-electron density
  RefSCMatrix PsiCC::Onerdm_relax_X(SpinCase1 spin) {

    // grab orbital info
    psi::PSIO& psio = exenv()->psio();

    // get # of occupied and unoccupied orbitals of spin S per irrep
    const std::vector<unsigned int>& occpi = reference()->occpi(spin);
    const std::vector<unsigned int>& uoccpi = reference()->uoccpi(spin);

    // obtain # of orbitals per irrep
    std::vector<unsigned int> occioff(nirrep_);
    std::vector<unsigned int> uoccioff(nirrep_);

    occioff[0] = 0;
    uoccioff[0] = 0;
    unsigned int nocc = occpi[0];
    unsigned int nuocc = uoccpi[0];
    for (unsigned int irrep = 1; irrep < nirrep_; ++irrep) {
      occioff[irrep] = occioff[irrep-1] + occpi[irrep-1];
      uoccioff[irrep] = uoccioff[irrep-1] + uoccpi[irrep-1];
      nocc += occpi[irrep];
      nuocc += uoccpi[irrep];
    }

    RefSCDimension rowdim = new SCDimension(nuocc);
    RefSCDimension coldim = new SCDimension(nocc);
    RefSCMatrix X = matrixkit()->matrix(rowdim, coldim);
    X.assign(0.0);

    // test that X is not empty
    if (rowdim.n() && coldim.n()) {
      // read in the i by a matrix in DPD format
      unsigned int nai_dpd = 0;
      for (unsigned int h = 0; h < nirrep_; ++h)
        nai_dpd += uoccpi[h] * occpi[h];

      double* Xai = new double[nai_dpd];

      //ExEnv::out0() << endl << "X nai_dpd: " << nai_dpd << endl;

      psio.open(CC_OEI, PSIO_OPEN_OLD);
      if (reference()->reftype() == PsiSCF::rhf) {
        psio.read_entry(CC_OEI, "XAI",
                       reinterpret_cast<char*>(Xai), nai_dpd*sizeof(double));

      } else if (reference()->reftype() == PsiSCF::uhf) {
        const char* Xai_lbl = (spin == Alpha) ? "XAI" : "Xai";
        psio.read_entry(CC_OEI, Xai_lbl,
                        reinterpret_cast<char*>(Xai), nai_dpd*sizeof(double));
      }
      psio.close(CC_OEI, 1);

      unsigned int ai = 0;
      for (unsigned int h = 0; h < nirrep_; ++h) {
        const unsigned int a_offset = uoccioff[h];
        const unsigned int i_offset = occioff[h];

        for (int a = 0; a < uoccpi[h]; ++a)
          for (int i = 0; i< occpi[h]; ++i, ++ai)
            X.set_element(a+a_offset, i+i_offset,Xai[ai]);
      }
      delete[] Xai;
    }
    // end of if (rowdim.n() && coldim.n())

    return X;
  }
  // end of import PSI3 relaxation for one-particle density

  // import PSI3 relaxation for one-particle density D_orbs
  RefSCMatrix PsiCC::Onerdm_relax_D(SpinCase1 spin) {

    // grab orbital info
    psi::PSIO& psio = exenv()->psio();

    // get # of occupied and unoccupied orbitals of spin S per irrep
    const std::vector<unsigned int>& occpi = reference()->occpi(spin);
    const std::vector<unsigned int>& uoccpi = reference()->uoccpi(spin);

    // obtain # of orbitals per irrep
    std::vector<unsigned int> occioff(nirrep_);
    std::vector<unsigned int> uoccioff(nirrep_);

    occioff[0] = 0;
    uoccioff[0] = 0;
    unsigned int nocc = occpi[0];
    unsigned int nuocc = uoccpi[0];
    for (unsigned int irrep = 1; irrep < nirrep_; ++irrep) {
      occioff[irrep] = occioff[irrep-1] + occpi[irrep-1];
      uoccioff[irrep] = uoccioff[irrep-1] + uoccpi[irrep-1];
      nocc += occpi[irrep];
      nuocc += uoccpi[irrep];
    }

    RefSCDimension rowdim = new SCDimension(nuocc);
    RefSCDimension coldim = new SCDimension(nocc);
    RefSCMatrix D_orbs = matrixkit()->matrix(rowdim, coldim);
    D_orbs.assign(0.0);

    // test that D is not empty
    if (rowdim.n() && coldim.n()) {
      // read in the i by a matrix in DPD format
      unsigned int nai_dpd = 0;
      for (unsigned int h = 0; h < nirrep_; ++h)
        nai_dpd += uoccpi[h] * occpi[h];

      double* Dai = new double[nai_dpd];

      ExEnv::out0() << endl << "D nai_dpd: " << nai_dpd << endl;

      psio.open(CC_OEI, PSIO_OPEN_OLD);
      if (reference()->reftype() == PsiSCF::rhf) {
        psio.read_entry(CC_OEI, "D(orb)(A,I)",
                       reinterpret_cast<char*>(Dai), nai_dpd*sizeof(double));

      } else if (reference()->reftype() == PsiSCF::uhf) {
        const char* Dai_lbl = (spin == Alpha) ? "D(orb)(A,I)" : "D(orb)(A,I)";
        psio.read_entry(CC_OEI, Dai_lbl,
                        reinterpret_cast<char*>(Dai), nai_dpd*sizeof(double));
      }
      psio.close(CC_OEI, 1);

      unsigned int ai = 0;
      for (unsigned int h = 0; h < nirrep_; ++h) {
        const unsigned int a_offset = uoccioff[h];
        const unsigned int i_offset = occioff[h];

        for (int a = 0; a < uoccpi[h]; ++a)
          for (int i = 0; i< occpi[h]; ++i, ++ai)
            D_orbs.set_element(a+a_offset, i+i_offset,Dai[ai]);
      }
      delete[] Dai;
    }
    // end of if (rowdim.n() && coldim.n())

    return D_orbs;
  }
  // end of Onerdm_relax_D

//  // Compute orbital relaxation contribution for
//  // one-electron density
//  void PsiCC::compute_onerdm_relax(RefSCMatrix& Dorbs_alpha,
//                           RefSCMatrix& Dorbs_beta) {
//
//    // grab orbital info
//    psi::PSIO& psio = exenv()->psio();
//
//    const SpinCase1 spin1 = Alpha;
//    // get # of occupied and unoccupied orbitals of spin S per irrep
//    const std::vector<unsigned int>& occ1pi = reference()->occpi(spin1);
//    const std::vector<unsigned int>& uocc1pi = reference()->uoccpi(spin1);
//
//    // obtain # of orbitals per irrep
//    std::vector<unsigned int> occ1pioff(nirrep_);
//    std::vector<unsigned int> uocc1pioff(nirrep_);
//    occ1pioff[0] = 0;
//    uocc1pioff[0] = 0;
//
//    unsigned int nocc1 = occ1pi[0];
//    unsigned int nuocc1 = uocc1pi[0];
//    for (unsigned int irrep = 1; irrep < nirrep_; ++irrep) {
//      occ1pioff[irrep] = occ1pioff[irrep-1] + occ1pi[irrep-1];
//      uocc1pioff[irrep] = uocc1pioff[irrep-1] + uocc1pi[irrep-1];
//      nocc1 += occ1pi[irrep];
//      nuocc1 += uocc1pi[irrep];
//    }
//
//    unsigned int na1i1_dpd = 0;
//    for (unsigned int h = 0; h < nirrep_; ++h)
//      na1i1_dpd += uocc1pi[h] * occ1pi[h];
//
//    // DPD of orbital product spaces
//    std::vector<size_t> a1i1_pi(nirrep_);
//    for (unsigned int h = 0; h < nirrep_; ++h) {
//      size_t nai = 0;
//      for (unsigned int g = 0; g < nirrep_; ++g) {
//        nai += (size_t)uocc1pi[g] * occ1pi[h ^ g];
//      }
//      a1i1_pi[h] = nai;
//    }
//
//    psio.open(CC_OEI, PSIO_OPEN_OLD);
//    psio.open(CC_MISC, PSIO_OPEN_OLD);
//    if (reference()->reftype() == PsiSCF::rhf) {
//
//      double* Xai = new double[na1i1_dpd];
//      ExEnv::out0() << endl << "X nai_dpd: " << na1i1_dpd << endl;
//
//      psio.read_entry(CC_OEI, "XAI",
//                      reinterpret_cast<char*>(Xai), na1i1_dpd*sizeof(double));
//
//      // Grab only irrep 0 of the orbital Hessian
////      unsigned int nai_dpd_h0 = uocc1pi[0] * occ1pi[0];
//      unsigned int nai_dpd_h0 = a1i1_pi[0];
//      ExEnv::out0() << endl << "nai_dpd of irrep 0: " << nai_dpd_h0 << endl;
//      double** A = block_matrix(nai_dpd_h0, nai_dpd_h0);
//
//      psio_address A_address = PSIO_ZERO;
//      for(int ai = 0; ai < nai_dpd_h0; ai++) {
//        psio.read(CC_MISC, "A(EM,AI)",
//                  reinterpret_cast<char*> (A[ai]), nai_dpd_h0*sizeof(double),
//                  A_address, &A_address);
//       }
//
//       FILE* outfile = tmpfile ();
//       pople(A, Xai, nai_dpd_h0, 1, 1e-12, outfile, 0);
//       fclose (outfile);
//
//       psio.close(CC_OEI, 1);
//       psio.close(CC_MISC, 1);
//
//       RefSCDimension rowdim = new SCDimension(nuocc1);
//       RefSCDimension coldim = new SCDimension(nocc1);
//       Dorbs_alpha = matrixkit()->matrix(rowdim, coldim);
//       Dorbs_alpha.assign(0.0);
//
//       unsigned int ai = 0;
//       for (unsigned int h = 0; h < nirrep_; ++h) {
//         const unsigned int a_offset = uocc1pioff[h];
//         const unsigned int i_offset = occ1pioff[h];
//
//         for (int a = 0; a < uocc1pi[h]; ++a)
//           for (int i = 0; i< occ1pi[h]; ++i, ++ai)
//             Dorbs_alpha.set_element(a+a_offset, i+i_offset,Xai[ai]);
//       }
//
//       delete[] Xai;
//
//       Dorbs_beta = Dorbs_alpha;
//       //Dorbs_alpha.print(prepend_spincase(Alpha,"CCSD one-particle density from relaxation:").c_str());
//
//       // test: print the relaxation effect from PSI3
////       RefSCMatrix Dorbs_psi = this->Onerdm_relax_D(Alpha);
////       Dorbs_psi.print(prepend_spincase(Alpha,"PSI3 Dorbs:").c_str());
//
//    } else if (reference()->reftype() == PsiSCF::uhf) {
//
//        const SpinCase1 spin2 = Beta;
//        const std::vector<unsigned int>& occ2pi = reference()->occpi(spin2);
//        const std::vector<unsigned int>& uocc2pi = reference()->uoccpi(spin2);
//
//        std::vector<unsigned int> occ2pioff(nirrep_);
//        std::vector<unsigned int> uocc2pioff(nirrep_);
//        occ2pioff[0] = 0;
//        uocc2pioff[0] = 0;
//
//        unsigned int nocc2 = occ2pi[0];
//        unsigned int nuocc2 = uocc2pi[0];
//        for (unsigned int irrep = 1; irrep < nirrep_; ++irrep) {
//          occ2pioff[irrep] = occ2pioff[irrep-1] + occ2pi[irrep-1];
//          uocc2pioff[irrep] = uocc2pioff[irrep-1] + uocc2pi[irrep-1];
//          nocc2 += occ2pi[irrep];
//          nuocc2 += uocc2pi[irrep];
//        }
//
//        unsigned int na2i2_dpd = 0;
//        for (unsigned int h = 0; h < nirrep_; ++h)
//          na2i2_dpd += uocc2pi[h] * occ2pi[h];
//
//        double* const Xai = new double[na1i1_dpd + na2i2_dpd];
//        ExEnv::out0() << endl << "X na1i1_dpd: " << na1i1_dpd << endl;
//        ExEnv::out0() << endl << "X na2i2_dpd: " << na2i2_dpd << endl;
//
//        double* iter_Xai = Xai;
//        psio.read_entry(CC_OEI, "XAI",
//                        reinterpret_cast<char*>(iter_Xai), na1i1_dpd*sizeof(double));
//
//        iter_Xai += na1i1_dpd;
//        psio.read_entry(CC_OEI, "Xai",
//                        reinterpret_cast<char*>(iter_Xai), na1i1_dpd*sizeof(double));
//
//        std::vector<size_t> a2i2_pi(nirrep_);
//        for (unsigned int h = 0; h < nirrep_; ++h) {
//          size_t nai = 0;
//          for (unsigned int g = 0; g < nirrep_; ++g) {
//            nai += (size_t)uocc2pi[g] * occ2pi[h ^ g];
//          }
//          a2i2_pi[h] = nai;
//        }
//
//        // Grab only irrep 0 of the orbital Hessian
//        unsigned int na1i1_dpd_h0 = a1i1_pi[0];
//        unsigned int na2i2_dpd_h0 = a2i2_pi[0];
//        unsigned int nai_dpd_h0  = na1i1_dpd_h0 + na2i2_dpd_h0;
//        ExEnv::out0() << endl << "na1i1_dpd of irrep 0: " << na1i1_dpd_h0 << endl;
//        ExEnv::out0() << endl << "na2i2_dpd of irrep 0: " << na2i2_dpd_h0 << endl;
//
//        double** A = block_matrix(nai_dpd_h0, nai_dpd_h0);
//
//        psio_address A1_address = PSIO_ZERO;
//        for(int a1i1 = 0; a1i1 < na1i1_dpd_h0; a1i1++) {
//          psio.read(CC_MISC, "A(AI,BJ)",
//                    reinterpret_cast<char*> (A[a1i1]), na1i1_dpd_h0*sizeof(double),
//                    A1_address, &A1_address);
//         }
//
//        psio_address A2_address = PSIO_ZERO;
//        for(int a2i2 = na1i1_dpd_h0; a2i2 < nai_dpd_h0; a2i2++) {
//          double* iter_A = A[a2i2] + na1i1_dpd_h0;
//          psio.read(CC_MISC, "A(ai,bj)",
//                    reinterpret_cast<char*> (iter_A), na2i2_dpd_h0*sizeof(double),
//                    A2_address, &A2_address);
//        }
//
//        psio_address A12_address = PSIO_ZERO;
//        for(int a1i1 = 0; a1i1 < na1i1_dpd_h0; a1i1++) {
//          double* iter_A = A[a1i1] + na1i1_dpd_h0;
//          psio.read(CC_MISC, "A(AI,bj)",
//                    reinterpret_cast<char*> (iter_A), na2i2_dpd_h0*sizeof(double),
//                    A12_address, &A12_address);
//
//          for(int a2i2 = na1i1_dpd_h0; a2i2 < nai_dpd_h0; a2i2++) {
//              A[a2i2][a1i1] = A[a1i1][a2i2];
//          }
//        }
//
//        FILE* outfile = tmpfile ();
//        pople(A, Xai, nai_dpd_h0, 1, 1e-12, outfile, 0);
//        fclose (outfile);
//
//        psio.close(CC_OEI, 1);
//        psio.close(CC_MISC, 1);
//
//        RefSCDimension rowdim1 = new SCDimension(nuocc1);
//        RefSCDimension coldim1 = new SCDimension(nocc1);
//        Dorbs_alpha = matrixkit()->matrix(rowdim1, coldim1);
//        Dorbs_alpha.assign(0.0);
//
//        unsigned int ai = 0;
//        for (unsigned int h = 0; h < nirrep_; ++h) {
//          const unsigned int a_offset = uocc1pioff[h];
//          const unsigned int i_offset = occ1pioff[h];
//
//          for (int a = 0; a < uocc1pi[h]; ++a)
//            for (int i = 0; i< occ1pi[h]; ++i, ++ai)
//              Dorbs_alpha.set_element(a+a_offset, i+i_offset, Xai[ai]);
//        }
//        //Dorbs_alpha.print(prepend_spincase(Alpha,"CCSD one-particle density from relaxation:").c_str());
//
//        RefSCDimension rowdim2 = new SCDimension(nuocc2);
//        RefSCDimension coldim2 = new SCDimension(nocc2);
//        Dorbs_beta = matrixkit()->matrix(rowdim2, coldim2);
//        Dorbs_beta.assign(0.0);
//
//        for (unsigned int h = 0; h < nirrep_; ++h) {
//          const unsigned int a_offset = uocc2pioff[h];
//          const unsigned int i_offset = occ2pioff[h];
//
//          for (int a = 0; a < uocc2pi[h]; ++a)
//            for (int i = 0; i< occ2pi[h]; ++i, ++ai)
//              Dorbs_beta.set_element(a+a_offset, i+i_offset, Xai[ai]);
//        }
//        //Dorbs_beta.print(prepend_spincase(Beta,"CCSD one-particle density from relaxation:").c_str());
//
//        delete[] Xai;
//
//        // test: print the relaxation effect from PSI3
////        RefSCMatrix Dorbs_psi_alpha = this->Onerdm_relax_D(Alpha);
////        Dorbs_psi_alpha.print(prepend_spincase(Alpha,"PSI3 Dorbs:").c_str());
////        RefSCMatrix Dorbs_psi_beta = this->Onerdm_relax_D(Beta);
////        Dorbs_psi_beta.print(prepend_spincase(Beta,"PSI3 Dorbs:").c_str());
//    }
//
//  }
//  // end of one-electron relaxation

#if 0
  RefSCMatrix PsiCC::transform_T1(const SparseMOIndexMap& occ_act_map,
                                  const SparseMOIndexMap& vir_act_map,
                                  const RefSCMatrix& T1,
                                  const Ref<SCMatrixKit>& kit) const {
    RefSCMatrix T1_new;
    T1_new = kit->matrix(T1.rowdim(), T1.coldim());
    T1_new.assign(0.0);

    const unsigned int no1 = occ_act_map.size();
    const unsigned int nv1 = vir_act_map.size();

    // convert T1 to new orbitals
    for (unsigned int i=0; i<no1; ++i) {
      const unsigned int ii = occ_act_map[i].first;
      const double ii_coef = occ_act_map[i].second;

      for (unsigned int a=0; a<nv1; ++a) {
        const unsigned int aa = vir_act_map[a].first;
        const double aa_coef = vir_act_map[a].second;

        const double elem = T1.get_element(ii, aa);
        T1_new.set_element(i, a, elem*ii_coef*aa_coef);
      }
    }

    return T1_new;
  }

  RefSCMatrix PsiCC::transform_T2(const SparseMOIndexMap& occ1_act_map,
                                  const SparseMOIndexMap& occ2_act_map,
                                  const SparseMOIndexMap& vir1_act_map,
                                  const SparseMOIndexMap& vir2_act_map,
                                  const RefSCMatrix& T2,
                                  const Ref<SCMatrixKit>& kit) const {
    RefSCMatrix T2_new;
    T2_new = kit->matrix(T2.rowdim(), T2.coldim());
    T2_new.assign(0.0);

    const unsigned int no1 = occ1_act_map.size();
    const unsigned int nv1 = vir1_act_map.size();
    const unsigned int no2 = occ2_act_map.size();
    const unsigned int nv2 = vir2_act_map.size();

    for (unsigned int i=0; i<no1; ++i) {
      const unsigned int ii = occ1_act_map[i].first;
      const double ii_coef = occ1_act_map[i].second;

      for (unsigned int j=0; j<no2; ++j) {
        const unsigned int jj = occ2_act_map[j].first;
        const double jj_coef = occ2_act_map[j].second;

        const unsigned int ij = i*no2+j;
        const unsigned int iijj = ii*no2+jj;

        const double ij_coef = ii_coef * jj_coef;

        for (unsigned int a=0; a<nv1; ++a) {
          const unsigned int aa = vir1_act_map[a].first;
          const double aa_coef = vir1_act_map[a].second;

          const double ija_coef = ij_coef * aa_coef;

          for (unsigned int b=0; b<nv2; ++b) {
            const unsigned int bb = vir2_act_map[b].first;
            const double bb_coef = vir2_act_map[b].second;

            const unsigned int ab = a*nv2+b;
            const unsigned int aabb = aa*nv2+bb;

            const double t2 = T2.get_element(iijj, aabb);
            const double t2_mpqc = t2 * ija_coef * bb_coef;
            T2_new.set_element(ij, ab, t2_mpqc);
          }
        }
      }
    }
    return T2_new;
  }

  RefSCMatrix PsiCC::transform_T1(const RefSCMatrix& occ_act_tform,
                                  const RefSCMatrix& vir_act_tform,
                                  const RefSCMatrix& T1,
                                  const Ref<SCMatrixKit>& kit) const {
    // convert to raw storage
    double* t1 = new double[T1.rowdim().n() * T1.coldim().n()];
    T1.convert(t1);

    RefSCMatrix T1_ia = kit->matrix(occ_act_tform.coldim(),
                                    vir_act_tform.coldim());
    T1_ia.assign(t1);

    RefSCMatrix T1_Ia = occ_act_tform * T1_ia;
    RefSCMatrix T1_IA = T1_Ia * vir_act_tform.t();

    return T1_IA;
  }

  RefSCMatrix PsiCC::transform_T2(const RefSCMatrix& occ1_act_tform,
                                  const RefSCMatrix& occ2_act_tform,
                                  const RefSCMatrix& vir1_act_tform,
                                  const RefSCMatrix& vir2_act_tform,
                                  const RefSCMatrix& T2,
                                  const Ref<SCMatrixKit>& kit) const {
    MPQC_ASSERT(occ1_act_tform.rowdim().n() == occ1_act_tform.coldim().n());
    MPQC_ASSERT(occ2_act_tform.rowdim().n() == occ2_act_tform.coldim().n());
    MPQC_ASSERT(vir1_act_tform.rowdim().n() == vir1_act_tform.coldim().n());
    MPQC_ASSERT(vir2_act_tform.rowdim().n() == vir2_act_tform.coldim().n());

    // convert to raw storage
    double* t2 = new double[T2.rowdim().n() * T2.coldim().n()];
    T2.convert(t2);

    const unsigned int nij = T2.rowdim().n();
    const unsigned int nab = T2.coldim().n();
    const unsigned int njab = occ2_act_tform.coldim().n() * nab;
    const unsigned int nija = vir1_act_tform.coldim().n() * nij;

    // store as i by jab
    RefSCMatrix T2_i_jab = kit->matrix(occ1_act_tform.coldim(),
                                       new SCDimension(njab));
    T2_i_jab.assign(t2);
    // transform i -> I
    RefSCMatrix T2_I_jab = occ1_act_tform * T2_i_jab;
    T2_i_jab = 0;
    T2_I_jab.convert(t2);

    // for each I store as j by ab
    // transform j -> J
    RefSCMatrix t2_j_ab = kit->matrix(occ2_act_tform.coldim(), T2.coldim());
    RefSCMatrix t2_J_ab = kit->matrix(occ2_act_tform.rowdim(), T2.coldim());
    const unsigned int nI = T2_I_jab.rowdim().n();
    for (unsigned int I=0; I<nI; ++I) {
      t2_j_ab.assign(t2 + I*njab);
      t2_J_ab.assign(0.0);
      t2_J_ab.accumulate_product(occ2_act_tform, t2_j_ab);
      t2_J_ab.convert(t2 + I*njab);
    }

    // for each IJ store as a by b
    // transform a -> A
    RefSCMatrix t2_a_b = kit->matrix(vir1_act_tform.coldim(),
                                     vir2_act_tform.coldim());
    RefSCMatrix t2_A_b = kit->matrix(vir1_act_tform.rowdim(),
                                     vir2_act_tform.coldim());
    const unsigned int nIJ = T2.rowdim().n();
    for (unsigned int IJ=0; IJ<nIJ; ++IJ) {
      t2_a_b.assign(t2 + IJ*nab);
      t2_A_b.assign(0.0);
      t2_A_b.accumulate_product(vir1_act_tform, t2_a_b);
      t2_A_b.convert(t2 + IJ*nab);
    }

    // store as IJA by b
    RefSCMatrix T2_IJA_b = kit->matrix(new SCDimension(nija), vir2_act_tform.coldim());
    T2_IJA_b.assign(t2);
    // transform B -> B
    RefSCMatrix T2_IJA_B = T2_IJA_b * vir2_act_tform.t();
    T2_IJA_b = 0;
    T2_IJA_B.convert(t2);

    // store as IJ by AB
    RefSCMatrix T2_IJ_AB = kit->matrix(T2.rowdim(), T2.coldim());
    T2_IJ_AB.assign(t2);
    return T2_IJ_AB;
  }
#endif

  namespace {
    bool gtzero(double a) {
      return a > 0.0;
    }
  }

  void PsiCC::compare_T2(const RefSCMatrix& T2, const RefSCMatrix& T2_ref,
                         SpinCase2 spin12, unsigned int no1, unsigned int no2, unsigned int nv1,
                         unsigned int nv2, double zero) const {
    using namespace sc::fastpairiter;
    if (spin12 != AlphaBeta)
      compare_tbint_tensors<AntiSymm,AntiSymm>(T2,T2_ref,no1,no2,nv1,nv2,zero);
    else
      compare_tbint_tensors<ASymm,ASymm>(T2,T2_ref,no1,no2,nv1,nv2,zero);
  }

  const Ref<OrbitalSpace>&
  PsiCC::occ_act_sb(SpinCase1 spin) {
    if (occ_act_sb_[spin])
      return occ_act_sb_[spin];
    if (reference_->reftype() == PsiSCF::rhf && spin == Beta)
      return occ_act_sb(Alpha);

    const unsigned int nmo = reference_->nmo();
    const std::vector<unsigned int>& mopi = reference_->mopi();
    const std::vector<unsigned int>& occpi = reference_->occpi(spin);
    const int nirreps = occpi.size();
    std::vector<bool> occ_mask(nmo, false);
    for (int irrep = 0, irrep_offset = 0; irrep < nirreps; ++irrep) {
      for (int i = 0; i < occpi[irrep]; ++i) {
        occ_mask[i + irrep_offset] = true;
      }
      irrep_offset += mopi[irrep];
    }

    const std::string id(spin == Alpha ? "I" : "i");
    Ref<OrbitalSpace> orbs_sb = new OrbitalSpace("", "",
                                                 reference_->coefs(spin),
                                                 basis(), integral(),
                                                 reference_->evals(spin), 0, 0,
                                                 OrbitalSpace::symmetry);
    Ref<OrbitalSpace> occ_sb =
        new MaskedOrbitalSpace("", "", orbs_sb, occ_mask);
    occ_act_sb_[spin]
        = new OrbitalSpace(id, prepend_spincase(spin,
                                                "active occupied MOs (Psi3)"),
                           occ_sb->coefs(), occ_sb->basis(),
                           occ_sb->integral(), occ_sb->evals(), nfzc_, 0,
                           OrbitalSpace::symmetry);

    return occ_act_sb_[spin];
  }

  const Ref<OrbitalSpace>&
  PsiCC::vir_act_sb(SpinCase1 spin) {
    if (vir_act_sb_[spin])
      return vir_act_sb_[spin];
    if (reference_->reftype() == PsiSCF::rhf && spin == Beta)
      return vir_act_sb(Alpha);

    const unsigned int nmo = reference_->nmo();
    const std::vector<unsigned int>& mopi = reference_->mopi();
    const std::vector<unsigned int>& uoccpi = reference_->uoccpi(spin);
    const int nirreps = uoccpi.size();
    std::vector<bool> uocc_mask(nmo, false);
    for (int irrep = 0, irrep_offset = 0; irrep < nirreps; ++irrep) {
      const unsigned int nocc = mopi[irrep] - uoccpi[irrep];
      for (int i = 0; i < uoccpi[irrep]; ++i) {
        uocc_mask[i + irrep_offset + nocc] = true;
      }
      irrep_offset += mopi[irrep];
    }

    const std::string id(spin == Alpha ? "A" : "a");
    Ref<OrbitalSpace> orbs_sb = new OrbitalSpace("", "",
                                                 reference_->coefs(spin),
                                                 basis(), integral(),
                                                 reference_->evals(spin), 0, 0,
                                                 OrbitalSpace::symmetry);
    Ref<OrbitalSpace> uocc_sb = new MaskedOrbitalSpace("", "", orbs_sb,
                                                       uocc_mask);
    vir_act_sb_[spin]
        = new OrbitalSpace(id, prepend_spincase(spin,
                                                "active virtual MOs (Psi3)"),
                           uocc_sb->coefs(), uocc_sb->basis(),
                           uocc_sb->integral(), uocc_sb->evals(), 0, nfzv_,
                           OrbitalSpace::symmetry);

    return vir_act_sb_[spin];
  }

  const Ref<OrbitalSpace>&
  PsiCC::occ_sb(SpinCase1 spin) {
    if (occ_sb_[spin])
      return occ_sb_[spin];
    if (reference_->reftype() == PsiSCF::rhf && spin == Beta)
      return occ_sb(Alpha);

    const unsigned int nmo = reference_->nmo();
    const std::vector<unsigned int>& mopi = reference_->mopi();
    const std::vector<unsigned int>& occpi = reference_->occpi(spin);
    const int nirreps = occpi.size();
    std::vector<bool> occ_mask(nmo, false);
    for (int irrep = 0, irrep_offset = 0; irrep < nirreps; ++irrep) {
      for (int i = 0; i < occpi[irrep]; ++i) {
        occ_mask[i + irrep_offset] = true;
      }
      irrep_offset += mopi[irrep];
    }

    Ref<OrbitalSpace> orbs_sb = new OrbitalSpace("", "",
                                                 reference_->coefs(spin),
                                                 basis(), integral(),
                                                 reference_->evals(spin), 0, 0,
                                                 OrbitalSpace::symmetry);
    occ_sb_[spin] = new MaskedOrbitalSpace("", "", orbs_sb, occ_mask);

    return occ_sb_[spin];
  }

  void
  PsiCC::obsolete() {
    occ_act_sb_[Alpha] = occ_act_sb_[Beta] = 0;
    vir_act_sb_[Alpha] = vir_act_sb_[Beta] = 0;
    occ_sb_[Alpha] = occ_sb_[Beta] = 0;
    T1_[Alpha] = T1_[Beta] = 0;
    T2_[Alpha] = T2_[Beta] = 0;
    T2_da4_[Alpha] = T2_da4_[Beta] = 0;
    Tau2_[Alpha] = Tau2_[Beta] = 0;
    Lambda1_[Alpha] = Lambda1_[Beta] = 0;
    Lambda2_[Alpha] = Lambda2_[Beta] = 0;
    PsiCorrWavefunction::obsolete();
  }

  //////////////////////////////////////////////////////////////////////////

  static ClassDesc PsiCCSD_cd(typeid(PsiCCSD), "PsiCCSD", 1, "public PsiCC", 0,
                              create<PsiCCSD>, create<PsiCCSD>);

  PsiCCSD::PsiCCSD(const Ref<KeyVal>&keyval) :
    PsiCC(keyval) {
    pccsd_alpha_ = keyval->doublevalue("pccsd_alpha", KeyValValuedouble(1.0));
    pccsd_beta_ = keyval->doublevalue("pccsd_beta", KeyValValuedouble(1.0));
    pccsd_gamma_ = keyval->doublevalue("pccsd_gamma", KeyValValuedouble(1.0));
  }

  PsiCCSD::~PsiCCSD() {
  }

  PsiCCSD::PsiCCSD(StateIn&s) :
    PsiCC(s) {
  }

  bool PsiCCSD::analytic_gradient_implemented() const {
    bool impl = false;
    PsiSCF::RefType reftype = reference_->reftype();
    if (reftype == PsiSCF::rhf || reftype == PsiSCF::hsoshf)
      impl = true;
    return impl;
  }

  void PsiCCSD::save_data_state(StateOut&s) {
    PsiCC::save_data_state(s);
  }

  void PsiCCSD::write_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->open();
    PsiCorrWavefunction::write_input(convergence);
    input->write_keyword("psi:wfn", "ccsd");
    input->write_keyword("ccenergy:convergence", convergence);
    input->write_keyword("ccenergy:maxiter", maxiter_);
    input->write_keyword("ccenergy:diis", diis_);
    input->write_keyword("ccenergy:diis_nvector", diis_nvector_);
    input->write_keyword("ccenergy:pccsd_alpha", pccsd_alpha_);
    input->write_keyword("ccenergy:pccsd_beta", pccsd_beta_);
    input->write_keyword("ccenergy:pccsd_gamma", pccsd_gamma_);
    input->close();
  }

  //////////////////////////////////////////////////////////////////////////

  static ClassDesc PsiCCSD_T_cd(typeid(PsiCCSD_T), "PsiCCSD_T", 1,
                                "public PsiCC", 0, create<PsiCCSD_T>, create<
                                    PsiCCSD_T>);

  PsiCCSD_T::PsiCCSD_T(const Ref<KeyVal>&keyval) :
    PsiCC(keyval) {
  }

  PsiCCSD_T::~PsiCCSD_T() {
  }

  PsiCCSD_T::PsiCCSD_T(StateIn&s) :
    PsiCC(s) {
  }

  void PsiCCSD_T::save_data_state(StateOut&s) {
    PsiCC::save_data_state(s);
    SavableState::save_state(reference_.pointer(), s);
  }

  void PsiCCSD_T::write_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->open();
    PsiCorrWavefunction::write_input(convergence);
    input->write_keyword("psi:wfn", "ccsd_t");
    input->write_keyword("ccenergy:convergence", convergence);
    input->write_keyword("ccenergy:maxiter", maxiter_);
    input->close();
  }

  //////////////////////////////////////////////////////////////////////////

  static ClassDesc PsiCC2_cd(typeid(PsiCC2), "PsiCC2", 1,
                             "public PsiCC", 0, create<PsiCC2>, create<
                             PsiCC2>);

  PsiCC2::PsiCC2(const Ref<KeyVal>&keyval) :
    PsiCC(keyval) {
  }

  PsiCC2::~PsiCC2() {
  }

  PsiCC2::PsiCC2(StateIn&s) :
    PsiCC(s) {
  }

  void PsiCC2::save_data_state(StateOut&s) {
    PsiCC::save_data_state(s);
    SavableState::save_state(reference_.pointer(), s);
  }

  void PsiCC2::write_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->open();
    PsiCorrWavefunction::write_input(convergence);
    input->write_keyword("psi:wfn", "cc2");
    input->write_keyword("ccenergy:convergence", convergence);
    input->write_keyword("ccenergy:maxiter", maxiter_);
    input->close();
  }

  //////////////////////////////////////////////////////////////////////////

  static ClassDesc PsiCC3_cd(typeid(PsiCC3), "PsiCC3", 1,
                             "public PsiCC", 0, create<PsiCC3>, create<
                             PsiCC3>);

  PsiCC3::PsiCC3(const Ref<KeyVal>&keyval) :
    PsiCC(keyval) {
  }

  PsiCC3::~PsiCC3() {
  }

  PsiCC3::PsiCC3(StateIn&s) :
    PsiCC(s) {
  }

  void PsiCC3::save_data_state(StateOut&s) {
    PsiCC::save_data_state(s);
    SavableState::save_state(reference_.pointer(), s);
  }

  void PsiCC3::write_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->open();
    PsiCorrWavefunction::write_input(convergence);
    input->write_keyword("psi:wfn", "cc3");
    input->write_keyword("ccenergy:convergence", convergence);
    input->write_keyword("ccenergy:maxiter", maxiter_);
    input->close();
  }

} // namespace

//////////////////

namespace {
  using namespace sc::fastpairiter;
  template <PairSymm BraSymm, PairSymm KetSymm>
  void
  compare_tbint_tensors(const RefSCMatrix& T2, const RefSCMatrix& T2_ref,
                        unsigned int nb1, unsigned int nb2, unsigned int nk1,
                        unsigned int nk2, double zero)
  {
    T2_ref.print("compare_tbint_tensors() -- T2(ref)");
    T2.print("compare_tbint_tensors() -- T2");
    sc::fastpairiter::MOPairIter<BraSymm> biter(nb1,nb2);
    sc::fastpairiter::MOPairIter<KetSymm> kiter(nk1,nk2);
    for(biter.start(); biter; biter.next()) {
      const int b12 = biter.ij();
      for(kiter.start(); kiter; kiter.next()) {
        const int k12 = kiter.ij();

        const double t2 = T2.get_element(b12, k12);
        const double t2_ref = T2_ref.get_element(b12, k12);

        if (fabs(t2_ref - t2) > zero) {
          T2_ref.print("compare_tbint_tensors() -- T2(ref)");
          T2.print("compare_tbint_tensors() -- T2");
          throw ProgrammingError("2-body tensors do not match",__FILE__,__LINE__);
        }
      }
    }
  }

    void _print(SpinCase2 spin,
               const Ref<DistArray4>& mat,
               const char* label) {
      if (mat->msg()->me() == 0) {
        const size_t nij = (spin != AlphaBeta && mat->ni() == mat->nj()) ? mat->ni() * (mat->ni()-1) / 2 : mat->ni() * mat->nj();
        const size_t nxy = (spin != AlphaBeta && mat->nx() == mat->ny()) ? mat->nx() * (mat->nx()-1) / 2 : mat->nx() * mat->ny();
        RefSCMatrix scmat = SCMatrixKit::default_matrixkit()->matrix(new SCDimension(nij), new SCDimension(nxy));
        scmat << mat;
        scmat.print(label);
      }
    }

} // anonymous namespace

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
