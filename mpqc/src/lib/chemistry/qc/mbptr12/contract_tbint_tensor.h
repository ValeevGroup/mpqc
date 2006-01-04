//
// contract_tbint_tensor.h
//
// Copyright (C) 2005 Edward Valeev
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
#pragma interface
#endif

#ifndef _chemistry_qc_mbptr12_contracttbinttensor_h
#define _chemistry_qc_mbptr12_contracttbinttensor_h

#include <cmath>
#include <util/misc/timer.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>

namespace sc {
  
  template <typename DataProcess_Bra,
            typename DataProcess_Ket,
            typename DataProcess_BraKet,
            bool CorrFactorInBra,
            bool CorrFactorInKet,
            bool CorrFactorInInt>
    void
    R12IntEval::contract_tbint_tensor(
      RefSCMatrix& T,
      unsigned int tbint_type_bra,
      unsigned int tbint_type_ket,
      const Ref<MOIndexSpace>& space1_bra,
      const Ref<MOIndexSpace>& space2_bra,
      const Ref<MOIndexSpace>& space1_int,
      const Ref<MOIndexSpace>& space2_int,
      const Ref<MOIndexSpace>& space1_ket,
      const Ref<MOIndexSpace>& space2_ket,
      const Ref<LinearR12::TwoParticleContraction>& tpcontract,
      bool antisymmetrize,
      const std::vector< Ref<TwoBodyMOIntsTransform> >& tforms_bra,
      const std::vector< Ref<TwoBodyMOIntsTransform> >& tforms_ket,
      const std::vector< Ref<TwoBodyIntDescr> >& intdescrs_bra,
      const std::vector< Ref<TwoBodyIntDescr> >& intdescrs_ket
    )
    {
      // if antisymmetrize==yes not sure about correct semantics of the function call -- disable for now
      if (antisymmetrize)
        throw FeatureNotImplemented("R12IntEval::contract_tbint_tensor_() -- antisymmetrization not supported yet",
                                    __FILE__,__LINE__);
      
      // are particles 1 and 2 equivalent?
      const bool part1_equiv_part2 = (space1_bra==space2_bra &&
                                      space1_ket==space2_ket &&
                                      space1_int==space2_int);
      // Check correct semantics of this call : if antisymmetrize then particles must be equivalent
      bool correct_semantics = (antisymmetrize && part1_equiv_part2) ||
                               !antisymmetrize;
      // also dimensions of tpcontract must match those of space1_int and space2_int
      correct_semantics = ( correct_semantics &&
                            (tpcontract->nrow() == space1_int->rank()) &&
                            (tpcontract->ncol() == space2_int->rank()) );
      if (!correct_semantics)
        throw ProgrammingError("R12IntEval::contract_tbint_tensor_() -- incorrect call semantics",
                               __FILE__,__LINE__);
      
      const bool CorrFactorInBraInt = CorrFactorInBra && CorrFactorInInt;
      const bool CorrFactorInKetInt = CorrFactorInKet && CorrFactorInInt;
      
      const unsigned int nbrasets = (CorrFactorInBra ? corrfactor()->nfunctions() : 1);
      const unsigned int nketsets = (CorrFactorInKet ? corrfactor()->nfunctions() : 1);
      const unsigned int nintsets = (CorrFactorInInt ? corrfactor()->nfunctions() : 1);
      const unsigned int nbraintsets = (CorrFactorInBraInt ? -1 : nbrasets*nintsets);
      const unsigned int nketintsets = (CorrFactorInKetInt ? -1 : nketsets*nintsets);
      
      
      //
      // create transforms, if needed
      //
      typedef std::vector< Ref<TwoBodyMOIntsTransform> > tformvec;
      
      // bra transforms
      tformvec transforms_bra = tforms_bra;
      if (transforms_bra.empty()) {
        if (CorrFactorInBraInt) {
          unsigned int fbraint = 0;
          for(unsigned int fbra=0; fbra<nbrasets; ++fbra) {
            for(unsigned int fint=0; fint<nintsets; ++fint, ++fbraint) {
              std::string tlabel(transform_label(space1_bra,space1_int,space2_bra,space2_int,fbra,fint));
              try {
                transforms_bra.push_back(get_tform_(tlabel));
              }
              catch (TransformNotFound& a){
                Ref<MOIntsTransformFactory> tfactory = r12info()->tfactory();
                tfactory->set_spaces(space1_bra,space1_int,space2_bra,space2_int);
                Ref<TwoBodyMOIntsTransform> tform = tfactory->twobody_transform(
                                                      MOIntsTransformFactory::StorageType_13,
                                                      tlabel,
                                                      intdescrs_bra[fbraint]
                                                    );
                transforms_bra.push_back(tform);
              }
            }
          }
        }
        else {
          for(int f=0; f<nbraintsets; f++) {
            std::string tlabel(transform_label(space1_bra,space1_int,space2_bra,space2_int,f));
            try {
              transforms_bra.push_back(get_tform_(tlabel));
            }
            catch (TransformNotFound& a){
              Ref<MOIntsTransformFactory> tfactory = r12info()->tfactory();
              tfactory->set_spaces(space1_bra,space1_int,space2_bra,space2_int);
              Ref<TwoBodyMOIntsTransform> tform = tfactory->twobody_transform(
                                                    MOIntsTransformFactory::StorageType_13,
                                                    tlabel,
                                                    intdescrs_bra[f]
                                                  );
              transforms_bra.push_back(tform);
            }
          }
        }
      }
      
      // ket transforms
      tformvec transforms_ket = tforms_ket;
      if (transforms_ket.empty()) {
        if (CorrFactorInKetInt) {
          unsigned int fketint = 0;
          for(unsigned int fket=0; fket<nketsets; ++fket) {
            for(unsigned int fint=0; fint<nintsets; ++fint, ++fketint) {
              std::string tlabel(transform_label(space1_ket,space1_int,space2_ket,space2_int,fket,fint));
              try {
                transforms_ket.push_back(get_tform_(tlabel));
              }
              catch (TransformNotFound& a){
                Ref<MOIntsTransformFactory> tfactory = r12info()->tfactory();
                tfactory->set_spaces(space1_ket,space1_int,space2_ket,space2_int);
                Ref<TwoBodyMOIntsTransform> tform = tfactory->twobody_transform(
                                                      MOIntsTransformFactory::StorageType_13,
                                                      tlabel,
                                                      intdescrs_ket[fketint]
                                                    );
                transforms_ket.push_back(tform);
              }
            }
          }
        }
        else {
          for(int f=0; f<nketintsets; f++) {
            std::string tlabel(transform_label(space1_ket,space1_int,space2_ket,space2_int,f));
            try {
              transforms_ket.push_back(get_tform_(tlabel));
            }
            catch (TransformNotFound& a){
              Ref<MOIntsTransformFactory> tfactory = r12info()->tfactory();
              tfactory->set_spaces(space1_ket,space1_int,space2_ket,space2_int);
              Ref<TwoBodyMOIntsTransform> tform = tfactory->twobody_transform(
                                                    MOIntsTransformFactory::StorageType_13,
                                                    tlabel,
                                                    intdescrs_ket[f]
                                                  );
              transforms_ket.push_back(tform);
            }
          }
        }
      }
      
      
      //
      // Generate contract label
      //
      tim_enter("Generic tensor contract");
      std::string label;
      {
        std::ostringstream oss_bra;
        oss_bra << "<" << space1_bra->id() << " " << space2_bra->id() << (antisymmetrize ? "||" : "|")
                << space1_int->id() << " " << space2_int->id() << ">";
        const std::string label_bra = oss_bra.str();
        std::ostringstream oss_ket;
        oss_ket << "<" << space1_ket->id() << " " << space2_ket->id() << (antisymmetrize ? "||" : "|")
                << space1_int->id() << " " << space2_int->id() << ">";
        const std::string label_ket = oss_ket.str();
        std::ostringstream oss;
        oss << "<" << space1_ket->id() << " " << space2_ket->id() << (antisymmetrize ? "||" : "|")
                << space1_int->id() << " " << space2_int->id() << "> = "
                << label_bra << " . " << label_ket << "^T";
        label = oss.str();
      }
      ExEnv::out0() << endl << indent
                    << "Entered generic contraction (" << label << ")" << endl;
      ExEnv::out0() << incindent;
      
      //
      // Construct maps
      //
      // WARNING: Assuming all transforms are over same spaces!!!
      //
      Ref<MOIndexSpace> tspace1_bra = transforms_bra[0]->space1();
      Ref<MOIndexSpace> tspace2_bra = transforms_bra[0]->space3();
      Ref<MOIndexSpace> tspace1_intb = transforms_bra[0]->space2();
      Ref<MOIndexSpace> tspace2_intb = transforms_bra[0]->space4();
      Ref<MOIndexSpace> tspace1_ket = transforms_ket[0]->space1();
      Ref<MOIndexSpace> tspace2_ket = transforms_ket[0]->space3();
      Ref<MOIndexSpace> tspace1_intk = transforms_ket[0]->space2();
      Ref<MOIndexSpace> tspace2_intk = transforms_ket[0]->space4();
      // maps spaceX to spaceX of the transform
      std::vector<unsigned int> map1_bra, map2_bra, map1_ket, map2_ket,
                                map1_intb, map2_intb, map1_intk, map2_intk;
      // maps space2_int to space1_intb of transform
      std::vector<unsigned int> map12_intb;
      // maps space1_int to space2_intb of transform
      std::vector<unsigned int> map21_intb;
      // maps space2_int to space1_intk of transform
      std::vector<unsigned int> map12_intk;
      // maps space1_int to space2_intk of transform
      std::vector<unsigned int> map21_intk;
      
      { // bra maps
        map1_bra = *tspace1_bra<<*space1_bra;
        map2_bra = *tspace2_bra<<*space2_bra;
        map1_intb = *tspace1_intb<<*space1_int;
        map2_intb = *tspace2_intb<<*space2_int;
        if (antisymmetrize) {
          if (tspace1_intb == tspace2_intb) {
            map12_intb = map1_intb;
            map21_intb = map2_intb;
          }
          else {
            map12_intb = *tspace1_intb<<*space2_int;
            map21_intb = *tspace2_intb<<*space1_int;
          }
        }
      }
      { // ket maps
        map1_ket = *tspace1_ket<<*space1_ket;
        map2_ket = *tspace2_ket<<*space2_ket;
        map1_intk = *tspace1_intk<<*space1_int;
        map2_intk = *tspace2_intk<<*space2_int;
        if (antisymmetrize) {
          if (tspace1_intk == tspace2_intk) {
            map12_intk = map1_intk;
            map21_intk = map2_intk;
          }
          else {
            map12_intk = *tspace1_intk<<*space2_int;
            map21_intk = *tspace2_intk<<*space1_int;
          }
        }
      }
      
      const unsigned int bratform_block_ncols = tspace2_intb->rank();
      const unsigned int kettform_block_ncols = tspace2_intk->rank();
      const RefDiagSCMatrix evals1_bra = space1_bra->evals();
      const RefDiagSCMatrix evals2_bra = space2_bra->evals();
      const RefDiagSCMatrix evals1_ket = space1_ket->evals();
      const RefDiagSCMatrix evals2_ket = space2_ket->evals();
      const RefDiagSCMatrix evals1_int = space1_int->evals();
      const RefDiagSCMatrix evals2_int = space2_int->evals();
      
      // Using spinorbital iterators means I don't take into account perm symmetry
      // More efficient algorithm will require generic code
      const SpinCase2 S = (antisymmetrize ? AlphaAlpha : AlphaBeta);
      SpinMOPairIter iterbra(space1_bra,space2_bra,S);
      SpinMOPairIter iterket(space1_ket,space2_ket,S);
      SpinMOPairIter iterint(space1_int,space2_int,S);
      // size of one block of <space1_bra space2_bra|
      const unsigned int nbra = iterbra.nij();
      // size of one block of <space1_ket space2_ket|
      const unsigned int nket = iterket.nij();
      // size of one block of |space1_int space2_int>
      const unsigned int nint = iterint.nij();

      // size of integral blocks to contract
      const unsigned int blksize_int = space1_int->rank() * space2_int->rank();
      double* T_ij = new double[blksize_int];
      double* T_kl = new double[blksize_int];
      
      //
      // Contraction loops
      //
      unsigned int fbraint = 0;
      unsigned int fbraoffset = 0;
      for(unsigned int fbra=0; fbra<nbrasets; ++fbra,fbraoffset+=nbra) {
        
        unsigned int fketint = 0;
        unsigned int fketoffset = 0;
        for(unsigned int fket=0; fket<nketsets; ++fket,fketoffset+=nket) {
          
          for(unsigned int fint=0; fint<nintsets; ++fint, ++fbraint, ++fketint) {
            
            Ref<TwoBodyMOIntsTransform> tformb = transforms_bra[fbraint];
            Ref<TwoBodyMOIntsTransform> tformk = transforms_ket[fketint];
            
            if (debug_ > 0) {
              ExEnv::out0() << indent << "Using transforms "
                                      << tformb->name() << " and "
                                      << tformk->name() << std::endl;
            }
            
            Ref<R12IntsAcc> accumb = tformb->ints_acc();
            Ref<R12IntsAcc> accumk = tformk->ints_acc();
            // if transforms have not been computed yet, compute
            if (accumb.null() || !accumb->is_committed()) {
              tformb->compute();
            }
            if (!accumb->is_active())
              accumb->activate();
            if (accumk.null() || !accumk->is_committed()) {
              tformk->compute();
            }
            if (!accumk->is_active())
              accumk->activate();
            
            // split work over tasks which have access to integrals
            // WARNING: assuming same accessibility for both bra and ket transforms
            vector<int> proc_with_ints;
            const int nproc_with_ints = tasks_with_ints_(accumb,proc_with_ints);
            const int me = r12info()->msg()->me();
            const bool nket_ge_nevals = (nket >= nproc_with_ints);
            
            if (accumb->has_access(me)) {
              
              unsigned int ijkl = 0;
              for(iterbra.start(); iterbra; iterbra.next()) {
                const int ij = iterbra.ij();
                
                bool fetch_ij_block = false;
                if (nket_ge_nevals) {
                  fetch_ij_block = true;
                }
                else {
                  const int first_ij_task = ijkl%nproc_with_ints;
                  const int last_ij_task = (ijkl+nket-1)%nproc_with_ints;
                  // figure out if for this ij there is ijkl to be handled by me
                  if ( !(first_ij_task > proc_with_ints[me] && proc_with_ints[me] > last_ij_task) )
                    fetch_ij_block = true;
                }
                if (!fetch_ij_block)
                  continue;
                
                const unsigned int i = iterbra.i();
                const unsigned int j = iterbra.j();
                const unsigned int ii = map1_bra[i];
                const unsigned int jj = map2_bra[j];

                tim_enter("MO ints retrieve");
                const double *ij_buf = accumb->retrieve_pair_block(ii,jj,tbint_type_bra);
                tim_exit("MO ints retrieve");
                if (debug_)
                  ExEnv::outn() << indent << "task " << me << ": obtained ij blocks" << endl;
                
                for(iterket.start(); iterket; iterket.next(),ijkl++) {
                  const int kl = iterket.ij();
                
                  const int ijkl_proc = ijkl%nproc_with_ints;
                  if (ijkl_proc != proc_with_ints[me])
                    continue;
                
                  const unsigned int k = iterket.i();
                  const unsigned int l = iterket.j();
                  const unsigned int kk = map1_ket[k];
                  const unsigned int ll = map2_ket[l];
                
                  if (debug_)
                    ExEnv::outn() << indent << "task " << me << ": working on (i,j | k,l) = ("
                                  << i << "," << j << " | " << k << "," << l << ")" << endl;
                  tim_enter("MO ints retrieve");
                  const double *kl_buf = accumk->retrieve_pair_block(kk,ll,tbint_type_ket);
                  tim_exit("MO ints retrieve");
                  if (debug_)
                    ExEnv::outn() << indent << "task " << me << ": obtained kl blocks" << endl;
                  
                  // zero out intblocks
                  memset(T_ij, 0, blksize_int*sizeof(double));
                  memset(T_kl, 0, blksize_int*sizeof(double));
                  
                  if (debug_ > 3) {
                    ExEnv::out0() << indent << "i = " << i << " j = " << j << " k = " << k << " l = " << l
                                  << incindent << endl;
                  }
                  
                  for(iterint.start(); iterint; iterint.next()) {
                    const unsigned int m = iterint.i();
                    const unsigned int n = iterint.j();
                    const int mn = iterint.ij();
                    int MN_bra, MN_ket;
                    {
                      const unsigned int mm = map1_intb[m];
                      const unsigned int nn = map2_intb[n];
                      MN_bra = mm*bratform_block_ncols+nn;
                    }
                    {
                      const unsigned int mm = map1_intk[m];
                      const unsigned int nn = map2_intk[n];
                      MN_ket = mm*kettform_block_ncols+nn;
                    }
                    
                    const double I_ijmn = ij_buf[MN_bra];
                    const double I_klmn = kl_buf[MN_ket];

                    if (debug_ > 3) {
                      ExEnv::out0() << indent << " m = " << m << " n = " << n << endl;
                    }
                    
                    if (!antisymmetrize) {
                      if (debug_ > 3) {
                        ExEnv::out0() << indent << " <ij|mn> = " << I_ijmn << endl
                                      << indent << " <kl|mn> = " << I_klmn << endl;
                      }
                      const double T_ijmn = DataProcess_Bra::I2T(I_ijmn,i,j,m,n,
                                                                 evals1_bra,evals1_int,evals2_bra,evals2_int);
                      const double T_klmn = DataProcess_Ket::I2T(I_klmn,k,l,m,n,
                                                                 evals1_ket,evals1_int,evals2_ket,evals2_int);
                      if (debug_ > 3) {
                        ExEnv::out0() << indent << " <ij|T|mn> = " << T_ijmn << endl
                                      << indent << " <kl|T|mn> = " << T_klmn << endl;
                      }
                      
                      T_ij[mn] = T_ijmn;
                      T_kl[mn] = T_klmn;
                    }
                    else {

                      int NM_bra, NM_ket;
                      {
                        const unsigned int mm = map21_intb[m];
                        const unsigned int nn = map12_intb[n];
                        NM_bra = nn*bratform_block_ncols+mm;
                      }
                      {
                        const unsigned int mm = map21_intk[m];
                        const unsigned int nn = map12_intk[n];
                        NM_ket = nn*kettform_block_ncols+mm;
                      }
                      
                      const double I_ijnm = ij_buf[NM_bra];
                      const double I_klnm = kl_buf[NM_ket];
                      
                      if (debug_ > 3) {
                        ExEnv::out0() << " <ij|mn> = " << I_ijmn << endl
                                      << " <ij|nm> = " << I_ijnm << endl
                                      << " <kl|mn> = " << I_klmn << endl
                                      << " <kl|nm> = " << I_klnm << endl;
                      }
                      
                      const double T_ijmn = DataProcess_Bra::I2T(I_ijmn-I_ijnm,i,j,m,n,
                                                                 evals1_bra,evals1_int,evals2_bra,evals2_int);
                      const double T_klmn = DataProcess_Ket::I2T(I_klmn-I_klnm,k,l,m,n,
                                                                 evals1_ket,evals1_int,evals2_ket,evals2_int);
                      T_ij[mn] = T_ijmn;
                      T_kl[mn] = T_klmn;
                    }
                    
                  } // int loop
                  
                  // contract matrices
                  double T_ijkl = tpcontract->contract(T_ij,T_kl);
                  if (debug_ > 3) {
                    ExEnv::out0() << decindent << indent
                                  << " <ij|kl> = " << T_ijkl << endl;
                  }
                  T_ijkl = DataProcess_BraKet::I2T(T_ijkl,i,j,k,l,
                                                   evals1_bra,evals1_ket,evals2_bra,evals2_ket);
                  if (debug_ > 3) {
                    ExEnv::out0() << indent << " <ij|T|kl> = " << T_ijkl << endl;
                  }
                  T.accumulate_element(ij+fbraoffset,kl+fketoffset,T_ijkl);
                  
                  accumk->release_pair_block(kk,ll,tbint_type_ket);
                  
                } // ket loop
                
                if (fetch_ij_block)
                  accumb->release_pair_block(ii,jj,tbint_type_bra);

              } // bra loop
            } // loop over tasks with access
          } // int blocks
        } // ket blocks
      } // bra blocks
      
      delete[] T_ij;
      delete[] T_kl;
      
      ExEnv::out0() << decindent;
      ExEnv::out0() << indent << "Exited generic contraction (" << label << ")" << endl;
      tim_exit("Generic tensor contract");
    }
              
}         

#endif

