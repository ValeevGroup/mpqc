//
// contract_tbint_tensors_to_obtensor.h
//
// Copyright (C) 2008 Martin Torheyden
//
// Author: Martin Torheyden <mtorhey@vt.edu>
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


#ifndef _chemistry_qc_mbptr12_contract_tbint_tensors_to_obtensor_h
#define _chemistry_qc_mbptr12_contract_tbint_tensors_to_obtensor_h

#include <cmath>
#include <math/mmisc/pairiter.h>
#include <chemistry/qc/lcao/utils.h>
#include <chemistry/qc/lcao/utils.impl.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <util/misc/print.h>


namespace sc {

  template <typename DataProcess_Bra,
            typename DataProcess_Ket,
            bool CorrFactorInBra,
            bool CorrFactorInKet,
            bool CorrFactorInInt>
    void R12IntEval::contract_tbint_tensors_to_obtensor(
                                            RefSCMatrix& T,
                                            SpinCase2 pairspin,
                                            TwoBodyTensorInfo tbtensor_type_bra,
                                            TwoBodyTensorInfo tbtensor_type_ket,
                                            const Ref<OrbitalSpace>& space1_bra,
                                            const Ref<OrbitalSpace>& space1_intb,
                                            const Ref<OrbitalSpace>& space2_intb,
                                            const Ref<OrbitalSpace>& space3_intb,
                                            const Ref<OrbitalSpace>& space1_ket,
                                            const Ref<OrbitalSpace>& space1_intk,
                                            const Ref<OrbitalSpace>& space2_intk,
                                            const Ref<OrbitalSpace>& space3_intk,
                                            const Ref<mbptr12::TwoParticleContraction>& tpcontract,
                                            const std::vector<std::string>& tformkeys_bra,
                                            const std::vector<std::string>& tformkeys_ket) {

    bool r12coeffs_bra = (tbtensor_type_bra.twobodytensor_type()==TwoBodyTensorInfo::geminalcoeff) ? true : false;
    bool r12coeffs_ket = (tbtensor_type_ket.twobodytensor_type()==TwoBodyTensorInfo::geminalcoeff) ? true : false;

    TwoBodyOper::type tbint_type_bra = (r12coeffs_bra) ? corrfactor()->tbint_type_eri()       // some dummy integral type
                                                       : tbtensor_type_bra.twobodyint_type();
    TwoBodyOper::type tbint_type_ket = (r12coeffs_ket) ? corrfactor()->tbint_type_eri()       // some dummy integral type
                                                       : tbtensor_type_ket.twobodyint_type();

    // are external spaces of particles 1 and 2 equivalent?
    const bool part1_strong_equiv_part2 = (space1_bra==space2_intb &&
                                           space1_ket==space2_intk);

    // can external spaces of particles 1 and 2 be equivalent?
    const bool part1_weak_equiv_part2 = false;

    bool correct_semantics = ( (space1_intb->rank() == space1_intk->rank()) &&
                               (space2_intb->rank() == space2_intk->rank()) &&
                               (space3_intb->rank() == space3_intk->rank()) );

    // also dimensions of tpcontract must match those of space2_intb and space3_intb
    correct_semantics = ( correct_semantics &&
                          (tpcontract->nrow() == space2_intb->rank()) &&
                          (tpcontract->ncol() == space3_intb->rank()) );

    if (!correct_semantics)
      throw ProgrammingError("R12IntEval::contract_tbint_tensor_() -- incorrect call semantics",
                             __FILE__,__LINE__);

    // are internal spaces of particles 1 and 2 equivalent?
    const bool part1_intequiv_part2 = (space2_intb==space3_intb &&
                                       space2_intk==space3_intk);

    const bool alphabeta = !(part1_strong_equiv_part2 &&
                             part1_intequiv_part2);
    //const bool alphabeta = (pairspin==AlphaBeta) ? true : false;
    // NOTE! Even if computing in AlphaBeta, internal sums can be over AlphaAlpha!!!

    // Using spinorbital iterators means I don't take into account perm symmetry
    //const SpinCase2 S = (alphabeta ? AlphaBeta : AlphaAlpha);
    const SpinCase2 S = pairspin;

    const bool antisymmetrize = !alphabeta;

    const bool CorrFactorInBraInt = CorrFactorInBra && CorrFactorInInt;
    const bool CorrFactorInKetInt = CorrFactorInKet && CorrFactorInInt;

    const unsigned int nbrasets = (CorrFactorInBra ? corrfactor()->nfunctions() : 1);
    const unsigned int nketsets = (CorrFactorInKet ? corrfactor()->nfunctions() : 1);
    const unsigned int nintsets = (CorrFactorInInt ? corrfactor()->nfunctions() : 1);
    const unsigned int nbraintsets = (CorrFactorInBraInt ? UINT_MAX : nbrasets*nintsets);
    const unsigned int nketintsets = (CorrFactorInKetInt ? UINT_MAX : nketsets*nintsets);


    //
    // create transforms, if needed
    //
    typedef std::vector<Ref<TwoBodyMOIntsTransform> > tformvec;

    // bra transforms
    const size_t num_tforms_bra = tformkeys_bra.size();
    tformvec transforms_bra(num_tforms_bra);
    for (unsigned int t = 0; t < num_tforms_bra; ++t) {
      transforms_bra[t] = moints_runtime4()->get(tformkeys_bra[t]);
    }

    // ket transforms
    const size_t num_tforms_ket = tformkeys_ket.size();
    tformvec transforms_ket(num_tforms_ket);
    for (unsigned int t = 0; t < num_tforms_ket; ++t) {
      transforms_ket[t] = moints_runtime4()->get(tformkeys_ket[t]);
    }

    //
    // Generate contract label
    //
    Timer timer("Generic tensor contract");
    std::string label;
    {
      std::ostringstream oss_bra;
      oss_bra << "<" << space1_bra->id() << " " << space1_intb->id() << (antisymmetrize ? "||" : "|")
              << space2_intb->id() << " " << space3_intb->id() << ">";
      const std::string label_bra = oss_bra.str();
      std::ostringstream oss_ket;
      oss_ket << "<" << space1_ket->id() << " " << space1_intk->id() << (antisymmetrize ? "||" : "|")
              << space2_intk->id() << " " << space3_intk->id() << ">";
      const std::string label_ket = oss_ket.str();
      std::ostringstream oss;
      oss << "<" << space1_bra->id() << (antisymmetrize ? "||" : "|")
              << space1_ket->id() << "> = "
              << label_bra << " . " << label_ket << "^T";
      label = oss.str();
    }
    ExEnv::out0() << std::endl << indent
                  << "Entered generic contraction (" << label << ")" << std::endl;
    ExEnv::out0() << incindent;

    //
    // Construct maps
    //
    // WARNING: Assuming all transforms are over same spaces!!!
    //
    Ref<OrbitalSpace> tspace1_bra = transforms_bra[0]->space1();
    Ref<OrbitalSpace> tspace1_intb = transforms_bra[0]->space3();
    Ref<OrbitalSpace> tspace2_intb = transforms_bra[0]->space2();
    Ref<OrbitalSpace> tspace3_intb = transforms_bra[0]->space4();
    Ref<OrbitalSpace> tspace1_ket = transforms_ket[0]->space1();
    Ref<OrbitalSpace> tspace1_intk = transforms_ket[0]->space3();
    Ref<OrbitalSpace> tspace2_intk = transforms_ket[0]->space2();
    Ref<OrbitalSpace> tspace3_intk = transforms_ket[0]->space4();
    // maps spaceX to spaceX of the transform
    std::vector<unsigned int> map1_bra, map1_ket,
                              map1_intb, map2_intb, map3_intb, map1_intk, map2_intk, map3_intk;
    // maps space3_intb to space2_intb of transform
    std::vector<unsigned int> map23_intb;
    // maps space2_intb to space3_intb of transform
    std::vector<unsigned int> map32_intb;
    // maps space3_intk to space2_intk of transform
    std::vector<unsigned int> map23_intk;
    // maps space2_intk to space3_intk of transform
    std::vector<unsigned int> map32_intk;

    { // bra maps
      map1_bra = *tspace1_bra<<*space1_bra;
      map1_intb = *tspace1_intb<<*space1_intb;
      map2_intb = *tspace2_intb<<*space2_intb;
      map3_intb = *tspace3_intb<<*space3_intb;
      // Will antisymmetrize the integrals? Then need ijkl AND ijlk
      if(S!=AlphaBeta) {
        if (tspace2_intb == tspace3_intb) {
          map23_intb = map2_intb;
          map32_intb = map3_intb;
        }
        else {
          map23_intb = *tspace2_intb<<*space3_intb;
          map32_intb = *tspace3_intb<<*space2_intb;
        }
      }
    }

    { // ket maps
      map1_ket = *tspace1_ket<<*space1_ket;
      map1_intk = *tspace1_intk<<*space1_intk;
      map2_intk = *tspace2_intk<<*space2_intk;
      map3_intk = *tspace3_intk<<*space3_intk;
      // Will antisymmetrize the integrals? Then need ijkl AND ijlk
      if(S!=AlphaBeta) {
        if (tspace2_intk == tspace3_intk) {
          map23_intk = map2_intk;
          map32_intk = map3_intk;
        }
        else {
          map23_intk = *tspace2_intk<<*space3_intk;
          map32_intk = *tspace3_intk<<*space2_intk;
        }
      }
    }

    const unsigned int bratform_block_ncols = tspace3_intb->rank();
    const unsigned int kettform_block_ncols = tspace3_intk->rank();
    const RefDiagSCMatrix evals1_bra = space1_bra->evals();
    const RefDiagSCMatrix evals1_ket = space1_ket->evals();
    const RefDiagSCMatrix evals1_intb = space1_intb->evals();
    const RefDiagSCMatrix evals2_intb = space2_intb->evals();
    const RefDiagSCMatrix evals3_intb = space3_intb->evals();
    const RefDiagSCMatrix evals1_intk = space1_intk->evals();
    const RefDiagSCMatrix evals2_intk = space2_intk->evals();
    const RefDiagSCMatrix evals3_intk = space3_intk->evals();

    SpinMOPairIter iterint(space2_intb->rank(),(alphabeta ? space3_intb : space2_intb)->rank(),S);
    // size of one block of |space1_int space2_int>
    const unsigned int nint = iterint.nij();
    RefSCDimension space1_bra_dim = space1_bra->dim();
    RefSCDimension space1_ket_dim = space1_ket->dim();
    const unsigned int nbra = space1_bra_dim->n();
    const unsigned int nket = space1_ket_dim->n();

    // size of integral blocks to contract
    const unsigned int blksize_int = space2_intb->rank() * space3_intb->rank();
    double* T_alphap = new double[blksize_int];
    double* T_qp = new double[blksize_int];


    int nmo_space1_bra = space1_bra->rank();
    int nmo_space1_intb = space1_intb->rank();
    int nmo_space1_ket = space1_ket->rank();
    int nmo_space1_intk = space1_intk->rank();

    // outermost loop over contraction blocks to minimize number of activate()/deactivate() calls
    for(unsigned int fint=0; fint<nintsets; ++fint) {
      unsigned int fbraoffset = 0;
      for(unsigned int fbra=0; fbra<nbrasets; ++fbra,fbraoffset+=nbra) {
        const unsigned int fbraint = fbra*nintsets+fint;
        Ref<TwoBodyMOIntsTransform> tformb = transforms_bra[fbraint];
        const Ref<TwoBodyIntDescr>& intdescrb = tformb->intdescr();
        const unsigned int intsetidx_bra = intdescrb->intset(tbint_type_bra);

        Ref<DistArray4> accumb = tformb->ints_distarray4();
        // if transforms have not been computed yet, compute
        if (accumb.null()) {
          tformb->compute();
          accumb = tformb->ints_distarray4();
        }
        accumb->activate();

        unsigned int fketoffset = 0;
        for(unsigned int fket=0; fket<nketsets; ++fket,fketoffset+=nket) {
          const unsigned int fketint = fket*nintsets+fint;
          Ref<TwoBodyMOIntsTransform> tformk = transforms_ket[fketint];
          const Ref<TwoBodyIntDescr>& intdescrk = tformk->intdescr();
          const unsigned int intsetidx_ket = intdescrk->intset(tbint_type_ket);

          Ref<DistArray4> accumk = tformk->ints_distarray4();
          if (accumk.null()) {
            tformk->compute();
            accumk = tformk->ints_distarray4();
          }
          accumk->activate();

          if (debug_ >= DefaultPrintThresholds::diagnostics) {
            ExEnv::out0() << indent << "Using transforms "
                                    << tformb->name() << " and "
                                    << tformk->name() << std::endl;
          }

          // split work over tasks which have access to integrals
          // WARNING: assuming same accessibility for both bra and ket transforms
          std::vector<int> proc_with_ints;
          const int nproc_with_ints = accumb->tasks_with_access(proc_with_ints);
          const int me = r12world()->world()->msg()->me();
          const bool nket_ge_nevals = (nket >= nproc_with_ints);

          Timer tim_mo_ints_retrieve;
          tim_mo_ints_retrieve.set_default("MO ints retrieve");
          if (accumb->has_access(me)) {
            unsigned int alphapq = 0;
            for(unsigned int alpha=0; alpha<nmo_space1_bra; alpha++) {  // external bra index
              for(unsigned int p=0; p<nmo_space1_intb; p++) {  // internal bra index
                unsigned int alphap = alpha*nmo_space1_intb + p;
                bool fetch_alphap_block = false;
                if (nket_ge_nevals) {
                  fetch_alphap_block = true;
                }
                else {
                  const int first_alphap_task = alphapq%nproc_with_ints;
                  const int last_alphap_task = (alphapq+nket-1)%nproc_with_ints;
                  // figure out if for this alphap there is alphapq to be handled by me
                  if ( !(first_alphap_task > proc_with_ints[me] && proc_with_ints[me] > last_alphap_task) )
                    fetch_alphap_block = true;
                }
                if(!fetch_alphap_block)
                  continue;

                const unsigned int alphaalpha = map1_bra[alpha];
                const unsigned int pp = map1_intb[p];

                const double *alphap_buf;
                if(!r12coeffs_bra) {
                  tim_mo_ints_retrieve.enter_default();
                  alphap_buf = accumb->retrieve_pair_block(alphaalpha,pp,intsetidx_bra);
                  tim_mo_ints_retrieve.exit_default();
                }

                if (debug_ >= DefaultPrintThresholds::mostO4)
                  ExEnv::outn() << indent << "task " << me << ": obtained alphap blocks" << std::endl;

                for(unsigned int q=0; q<nmo_space1_ket; q++) {  // external ket index
                  unsigned int pq = p*nmo_space1_ket + q;  // pq is a rectangular index pair
                  const int alphapq_proc = alphapq%nproc_with_ints;
                  if(pairspin==AlphaBeta) {
                    if (alphapq_proc != proc_with_ints[me])
                      continue;
                    const unsigned int qq = map1_ket[q];

                    if (debug_ >= DefaultPrintThresholds::mostO4)
                      ExEnv::outn() << indent << "task " << me << ": working on <alpha,p | q,p> = <"
                                    << alpha << "," << p << " | " << q << "," << p << ">" << std::endl;

                    const double *qp_buf;
                    if(!r12coeffs_ket) {
                      tim_mo_ints_retrieve.enter_default();
                      qp_buf = accumk->retrieve_pair_block(qq,pp,intsetidx_ket);
                      tim_mo_ints_retrieve.exit_default();
                    }

                    if (debug_ >= DefaultPrintThresholds::mostO4)
                      ExEnv::outn() << indent << "task " << me << ": obtained qp blocks" << std::endl;

                    // zero out intblocks
                    memset(T_alphap, 0, blksize_int*sizeof(double));
                    memset(T_qp, 0, blksize_int*sizeof(double));

                    if (debug_ >= DefaultPrintThresholds::mostO4) {
                      ExEnv::out0() << indent << "alpha = " << alpha << " p = " << p << " q = " << q
                                    << incindent << std::endl;
                    }

                    for(iterint.start(); iterint; iterint.next()) {  // internal pair index. Add \f$ \bar{r}_{p_{\gamma f\gamma} \alpha}^{rs} t_{rs}^{pq} \f$. rs is a rectangular index pair.
                      unsigned int r = iterint.i();
                      unsigned int s = iterint.j();
                      unsigned int rs = iterint.ij();

                      int RS_bra,RS_ket;
                      {
                        const unsigned int rr = map2_intb[r];
                        const unsigned int ss = map3_intb[s];
                        RS_bra = rr*bratform_block_ncols + ss;
                      }
                      {
                        const unsigned int rr = map2_intk[r];
                        const unsigned int ss = map3_intk[s];
                        RS_ket = rr*kettform_block_ncols + ss;
                      }
                      double I_alphaprs;
                      if(!r12coeffs_bra) {
                        I_alphaprs = alphap_buf[RS_bra];
                      }
                      else {
                        I_alphaprs = C_CuspConsistent(alpha,p,r,s,S);
                      }
                      double I_qprs;
                      if(!r12coeffs_ket) {
                        I_qprs = qp_buf[RS_ket];
                      }
                      else {
                        I_qprs = C_CuspConsistent(q,p,r,s,S);
                      }

                      if (debug_ >= DefaultPrintThresholds::mostO6) {
                        ExEnv::out0() << indent << " r = " << r << " s = " << s << std::endl;
                      }

                      if (debug_ >= DefaultPrintThresholds::mostO6) {
                        ExEnv::out0() << indent << " <alphap|rs> = " << I_alphaprs << std::endl
                                      << indent << "     <qp|rs> = " << I_qprs << std::endl;
                      }

                      double T_alphaprs;
                      if(!r12coeffs_bra) {
                        T_alphaprs = DataProcess_Bra::I2T(I_alphaprs,alpha,p,r,s,
                                                      evals1_bra,evals2_intb,evals1_intb,evals3_intb);
                      }
                      else {
                        T_alphaprs = I_alphaprs;
                      }

                      double T_qprs;
                      if(!r12coeffs_ket) {
                        T_qprs = DataProcess_Ket::I2T(I_qprs,q,p,r,s,
                                                      evals1_ket,evals2_intk,evals1_intk,evals3_intk);
                      }
                      else {
                        T_qprs = I_qprs;
                      }

                      if (debug_ >= DefaultPrintThresholds::mostO6) {
                        ExEnv::out0() << indent << " <alphap|T|rs> = " << T_alphaprs << std::endl
                                      << indent << "     <qp|T|rs> = " << T_qprs << std::endl;
                      }

                      T_alphap[rs] = T_alphaprs;
                      T_qp[rs] = T_qprs;

                    }  // loop over iterint, i.e. rs

                    // contract matrices
                    double T_alphapqp = tpcontract->contract(T_alphap,T_qp);
                    T_alphapqp *= 2.0;

                    if (debug_ >= DefaultPrintThresholds::mostO4) {
                      ExEnv::out0() << decindent << indent
                                    << " <alphap|qp> = " << T_alphapqp << std::endl;
                    }

                    T.accumulate_element(alpha,q,T_alphapqp);

                    if(!r12coeffs_ket) {
                      accumk->release_pair_block(qq,pp,intsetidx_ket);
                    }

                    alphapq += 1;
                  }  // if(pairspin==AlphaBeta)
                  else {  // if(pairspin!=AlphaBeta)
                    if(q==p) continue;  // q==p is not allowed for pairspin!=AlphaBeta -> skip.
                    if (alphapq_proc != proc_with_ints[me])
                      continue;
                    const unsigned int qq = map1_ket[q];

                    if (debug_ >= DefaultPrintThresholds::mostO4)
                      ExEnv::outn() << indent << "task " << me << ": working on <alpha,p | q,p> = <"
                                    << alpha << "," << p << " | " << q << "," << p << ">" << std::endl;

                    const double *qp_buf;
                    const double *pq_buf;
                    if(!r12coeffs_ket) {
                      tim_mo_ints_retrieve.enter_default();
                      if(q>p){
                        qp_buf = accumk->retrieve_pair_block(qq,pp,intsetidx_ket);
                      }
                      else {
                        pq_buf = accumk->retrieve_pair_block(pp,qq,intsetidx_ket);
                      }
                      tim_mo_ints_retrieve.exit_default();
                    }

                    if (debug_ >= DefaultPrintThresholds::mostO4)
                      ExEnv::outn() << indent << "task " << me << ": obtained qp blocks" << std::endl;

                    // zero out intblocks
                    memset(T_alphap, 0, blksize_int*sizeof(double));
                    memset(T_qp, 0, blksize_int*sizeof(double));

                    if (debug_ >= DefaultPrintThresholds::mostO4) {
                      ExEnv::out0() << indent << "alpha = " << alpha << " p = " << p << " q = " << q
                                    << incindent << std::endl;
                    }

                    for(iterint.start(); iterint; iterint.next()) {  // internal pair index
                      unsigned int r = iterint.i();
                      unsigned int s = iterint.j();
                      unsigned int rs = iterint.ij();

                      int RS_bra,RS_ket;
                      {
                        const unsigned int rr = map2_intb[r];
                        const unsigned int ss = map3_intb[s];
                        RS_bra = rr*bratform_block_ncols + ss;
                      }
                      {
                        const unsigned int rr = map2_intk[r];
                        const unsigned int ss = map3_intk[s];
                        RS_ket = rr*kettform_block_ncols + ss;
                      }
                      double I_alphaprs;
                      if(!r12coeffs_bra) {
                        I_alphaprs = alphap_buf[RS_bra];
                      }
                      else {
                        I_alphaprs = C_CuspConsistent(alpha,p,r,s,S);
                      }
                      // getting integrals with exchanged indices
                      int SR_bra, SR_ket;
                      {
                        const unsigned int rr = map32_intb[r];
                        const unsigned int ss = map23_intb[s];
                        SR_bra = ss*bratform_block_ncols + rr;
                      }
                      {
                        const unsigned int rr = map32_intk[r];
                        const unsigned int ss = map23_intk[s];
                        SR_ket = ss*kettform_block_ncols + rr;
                      }

                      double I_alphapsr;
                      if(!r12coeffs_bra) {
                        I_alphapsr = alphap_buf[SR_bra];
                      }
                      else {
                        I_alphapsr = C_CuspConsistent(alpha,p,s,r,S);
                      }

                      if (debug_ >= DefaultPrintThresholds::mostO6) {
                        ExEnv::out0() << indent << " r = " << r << " s = " << s << std::endl;
                      }

                      if(q>p) {  /// Add \f$ - \bar{r}_{p_{\gamma f\gamma} \alpha}^{rs} t_{rs}^{qp} \f$. qp is a triangular index.
                        double I_qprs;
                        if(!r12coeffs_ket) {
                          I_qprs = qp_buf[RS_ket];
                        }
                        else {
                          I_qprs = C_CuspConsistent(q,p,r,s,S);
                        }

                        if (debug_ >= DefaultPrintThresholds::mostO6) {
                          ExEnv::out0() << indent << " <alphap|rs> = " << I_alphaprs << std::endl
                                        << indent << "     <qp|rs> = " << I_qprs << std::endl;
                        }

                        double I_qpsr;
                        if(!r12coeffs_ket) {
                          I_qpsr = qp_buf[SR_ket];
                        }
                        else {
                          I_qpsr = C_CuspConsistent(q,p,s,r,S);
                        }

                        if (debug_ >= DefaultPrintThresholds::mostO6) {
                          ExEnv::out0() << " <alphap|rs> = " << I_alphaprs << std::endl
                                        << " <alphap|sr> = " << I_alphapsr << std::endl
                                        << "     <qp|rs> = " << I_qprs << std::endl
                                        << "     <qp|sr> = " << I_qpsr << std::endl;
                        }

                        double T_alphaprs;
                        if(!r12coeffs_bra) {
                          T_alphaprs = DataProcess_Bra::I2T(I_alphaprs-I_alphapsr,alpha,p,r,s,
                                                        evals1_bra,evals2_intb,evals1_intb,evals3_intb);
                        }
                        else {
                          T_alphaprs = I_alphaprs;
                        }

                        double T_qprs;
                        if(!r12coeffs_ket) {
                          T_qprs = DataProcess_Ket::I2T(I_qprs-I_qpsr,q,p,r,s,
                                                        evals1_ket,evals2_intk,evals1_intk,evals3_intk);
                        }
                        else {
                          T_qprs = I_qprs;
                        }

                        T_alphap[rs] = T_alphaprs;
                        T_qp[rs] = T_qprs;

                      }  // if(q>p)
                      else {  // if(p>q). Add \f$ \bar{r}_{p_{\gamma f\gamma} \alpha}^{rs} t_{rs}^{pq} \f$. pq is a triangular index.
                        double I_pqrs;
                        if(!r12coeffs_ket) {
                          I_pqrs = pq_buf[RS_ket];
                        }
                        else {
                          I_pqrs = C_CuspConsistent(p,q,r,s,S);
                        }

                        double I_pqsr;
                        if(!r12coeffs_ket) {
                          I_pqsr = pq_buf[SR_ket];
                        }
                        else {
                          I_pqsr = C_CuspConsistent(p,q,s,r,S);
                        }

                        if (debug_ >= DefaultPrintThresholds::mostO6) {
                          ExEnv::out0() << " <alphap|rs> = " << I_alphaprs << std::endl
                                        << " <alphap|rs> = " << I_alphapsr << std::endl
                                        << "     <pq|rs> = " << I_pqrs << std::endl
                                        << "     <pq|sr> = " << I_pqsr << std::endl;
                        }

                        double T_alphaprs;
                        if(!r12coeffs_bra) {
                          T_alphaprs = DataProcess_Bra::I2T(I_alphaprs-I_alphapsr,alpha,p,r,s,
                                                        evals1_bra,evals2_intb,evals1_intb,evals3_intb);
                        }
                        else {
                          T_alphaprs = I_alphaprs;
                        }

                        double T_qprs;
                        if(!r12coeffs_ket) {
                          T_qprs = DataProcess_Ket::I2T(I_pqsr-I_pqrs,q,p,r,s,
                                                        evals1_ket,evals2_intk,evals1_intk,evals3_intk);
                        }
                        else {
                          T_qprs = -I_pqrs;
                        }

                        T_alphap[rs] = T_alphaprs;
                        T_qp[rs] = T_qprs;

                      }  // if(p>q)
                    }  // loop over iterint, i.e. rs

                    // contract matrices
                    double T_alphapqp = tpcontract->contract(T_alphap,T_qp);

                    if (debug_ >= DefaultPrintThresholds::mostO4) {
                      ExEnv::out0() << decindent << indent
                                    << " <alphap|qp> = " << T_alphapqp << std::endl;
                    }

                    T.accumulate_element(alpha,q,T_alphapqp);

                    if(!r12coeffs_ket) {
                      if(q>p) {
                        accumk->release_pair_block(qq,pp,intsetidx_ket);
                      }
                      else {
                        accumk->release_pair_block(pp,qq,intsetidx_ket);
                      }
                    }

                    alphapq += 1;
                  }  // if(pairspin!=AlphaBeta)
                }  // loop over q

                if(fetch_alphap_block && (!r12coeffs_bra)) {
                  accumb->release_pair_block(alphaalpha,pp,intsetidx_bra);
                }

              }  // loop over p
            }  // loop over alpha

          }

          if (accumb != accumk) {
            if (accumk->data_persistent()) accumk->deactivate();
          }

        }  // fket loop
        if (accumb->data_persistent()) accumb->deactivate();
      }  // fbra loop
    }  // fint loop

    delete[] T_alphap;
    delete[] T_qp;

    ExEnv::out0() << decindent;
    ExEnv::out0() << indent << "Exited generic contraction (" << label << ")" << std::endl;

    timer.exit();
  }

}

#endif /* _chemistry_qc_mbptr12_contract_tbint_tensors_to_obtensor_h */
