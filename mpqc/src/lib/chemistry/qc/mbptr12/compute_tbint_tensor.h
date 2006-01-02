//
// compute_tbint_tensor.h
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

#include <util/misc/timer.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>

#ifndef _chemistry_qc_mbptr12_computetbinttensor_h
#define _chemistry_qc_mbptr12_computetbinttensor_h

namespace sc {
  
  template <typename DataProcess,
            bool CorrFactorInBra,
            bool CorrFactorInKet>
    void
    R12IntEval::compute_tbint_tensor(RefSCMatrix& T,
                                     int tbint_type,
                                     const Ref<MOIndexSpace>& space1,
                                     const Ref<MOIndexSpace>& space2,
                                     const Ref<MOIndexSpace>& space3,
                                     const Ref<MOIndexSpace>& space4,
                                     bool antisymmetrize,
                                     const std::vector< Ref<TwoBodyMOIntsTransform> >& transvec,
                                     const std::vector< Ref<TwoBodyIntDescr> >& intdescrs)
    {
      // are particles 1 and 2 equivalent?
      const bool part1_equiv_part2 = (space1==space3 && space2==space4);
      // Check correct semantics of this call : if antisymmetrize then particles must be equivalent
      const bool correct_semantics = (antisymmetrize && part1_equiv_part2) ||
                                     !antisymmetrize;
      if (!correct_semantics)
        throw ProgrammingError("R12IntEval::compute_tbint_tensor_() -- incorrect call semantics",
                               __FILE__,__LINE__);
      
      const bool CorrFactorInBraKet = CorrFactorInBra && CorrFactorInKet;
      
      const unsigned int num13sets = (CorrFactorInBra ? corrfactor()->nfunctions() : 1);
      const unsigned int num24sets = (CorrFactorInKet ? corrfactor()->nfunctions() : 1);
      const unsigned int nsets = (CorrFactorInBraKet ? -1 : num13sets*num24sets);
      
      // create transforms, if needed
      typedef std::vector< Ref<TwoBodyMOIntsTransform> > tformvec;
      tformvec transforms = transvec;
      if (transforms.empty()) {
        if (CorrFactorInBraKet) {
          for(unsigned int f13=0; f13<num13sets; f13++) {
            for(unsigned int f24=0; f24<num24sets; f24++) {
              transforms.push_back(get_tform_(transform_label(space1,space2,space3,space4,f13,f24)));
            }
          }
        }
        else {
          for(int f=0; f<nsets; f++)
            transforms.push_back(get_tform_(transform_label(space1,space2,space3,space4,f)));
        }
      }
      
      tim_enter("Generic tensor");
      std::ostringstream oss;
      oss << "<" << space1->id() << " " << space3->id() << (antisymmetrize ? "||" : "|")
      << space2->id() << " " << space4->id() << ">";
      const std::string label = oss.str();
      ExEnv::out0() << endl << indent
      << "Entered generic tensor (" << label << ") evaluator" << endl;
      ExEnv::out0() << incindent;
      
      //
      // WARNING: Assuming all transforms are over same spaces!!!
      //
      Ref<MOIndexSpace> tspace1 = transforms[0]->space1();
      Ref<MOIndexSpace> tspace2 = transforms[0]->space2();
      Ref<MOIndexSpace> tspace3 = transforms[0]->space3();
      Ref<MOIndexSpace> tspace4 = transforms[0]->space4();
      
      // maps spaceX to spaceX of the transform
      std::vector<unsigned int> map1, map2, map3, map4;
      // maps space4 to space2 of transform
      std::vector<unsigned int> map24;
      // maps space2 to space4 of transform
      std::vector<unsigned int> map42;
      map1 = *tspace1<<*space1;
      map2 = *tspace2<<*space2;
      map3 = *tspace3<<*space3;
      map4 = *tspace4<<*space4;
      if (antisymmetrize) {
        if (tspace2 == tspace4) {
          map24 = map2;
          map42 = map4;
        }
        else {
          map24 = *tspace2<<*space4;
          map42 = *tspace4<<*space2;
        }
      }
      
      const unsigned int trank4 = tspace4->rank();
      const RefDiagSCMatrix evals1 = space1->evals();
      const RefDiagSCMatrix evals2 = space2->evals();
      const RefDiagSCMatrix evals3 = space3->evals();
      const RefDiagSCMatrix evals4 = space4->evals();
      
      // Using spinorbital iterators means I don't take into account perm symmetry
      // More efficient algorithm will require generic code
      const SpinCase2 S = (antisymmetrize ? AlphaAlpha : AlphaBeta);
      SpinMOPairIter iter13(space1,space3,S);
      SpinMOPairIter iter24(space2,space4,S);
      // size of one block of <space1 space3|
      const unsigned int n13 = iter13.nij();
      // size of one block of |space3 space4>
      const unsigned int n24 = iter24.nij();
      
      unsigned int f13f24 = 0;
      unsigned int f13offset = 0;
      for(unsigned int f13=0; f13<num13sets; ++f13,f13offset+=n13) {
        
        unsigned int f24offset = 0;
        for(unsigned int f24=0; f24<num24sets; ++f24,f24offset+=n24,++f13f24) {
          
          Ref<TwoBodyMOIntsTransform> tform = transforms[f13f24];
          
          if (debug_ > 0)
            ExEnv::out0() << indent << "Using transform " << tform->name() << std::endl;
          
          Ref<R12IntsAcc> accum = tform->ints_acc();
          // if transforms have not been computed yet, compute
          if (accum.null() || !accum->is_committed()) {
            tform->compute(intdescrs[f13f24]);
          }
          if (!accum->is_active())
            accum->activate();
          
          // split work over tasks which have access to integrals
          vector<int> proc_with_ints;
          const int nproc_with_ints = tasks_with_ints_(accum,proc_with_ints);
          const int me = r12info()->msg()->me();
          
          if (accum->has_access(me)) {
            for(iter13.start(); iter13; iter13.next()) {
              const int ij = iter13.ij();
              
              const int ij_proc = ij%nproc_with_ints;
              if (ij_proc != proc_with_ints[me])
                continue;
              
              const unsigned int i = iter13.i();
              const unsigned int j = iter13.j();
              const unsigned int ii = map1[i];
              const unsigned int jj = map3[j];
              
              if (debug_)
                ExEnv::outn() << indent << "task " << me << ": working on (i,j) = " << i << "," << j << " " << endl;
              tim_enter("MO ints retrieve");
              const double *ij_buf = accum->retrieve_pair_block(ii,jj,tbint_type);
              tim_exit("MO ints retrieve");
              if (debug_)
                ExEnv::outn() << indent << "task " << me << ": obtained ij blocks" << endl;
              
              for(iter24.start(); iter24; iter24.next()) {
                const unsigned int a = iter24.i();
                const unsigned int b = iter24.j();
                const unsigned int aa = map2[a];
                const unsigned int bb = map4[b];
                const int AB = aa*trank4+bb;
                const int ab = iter24.ij();
                
                const double I_ijab = ij_buf[AB];
                
                if (!antisymmetrize) {
                  if (debug_ > 2) {
                    ExEnv::out0() << "i = " << i << " j = " << j << " a = " << a << " b = " << b
                    << " <ij|ab> = " << I_ijab << endl;
                  }
                  const double T_ijab = DataProcess::I2T(I_ijab,i,j,a,b,evals1,evals2,evals3,evals4);
                  T.accumulate_element(ij+f13offset,ab+f24offset,T_ijab);
                }
                else {
                  const int aa = map42[a];
                  const int bb = map24[b];
                  const int BA = bb*trank4+aa;
                  const double I_ijba = ij_buf[BA];
                  if (debug_ > 2) {
                    ExEnv::out0() << "i = " << i << " j = " << j << " a = " << a << " b = " << b
                    << " <ij|ab> = " << I_ijab << " <ij|ba> = " << I_ijba << endl;
                  }
                  const double I_anti = I_ijab - I_ijba;
                  const double T_ijab = DataProcess::I2T(I_anti,i,j,a,b,evals1,evals2,evals3,evals4);
                  T.accumulate_element(ij+f13offset,ab+f24offset,T_ijab);
                }
                
              } // 24 loop
              accum->release_pair_block(ii,jj,tbint_type);
              
            } // 13 loop
          } // loop over tasks with access
          
        }
      }
      
      ExEnv::out0() << decindent;
      ExEnv::out0() << indent << "Exited generic tensor (" << label << ") evaluator" << endl;
      tim_exit("Generic tensor");
    }
  
#if 1
  /// Contains classes used to compute many-body tensors
  namespace ManyBodyTensors {
    
    /// Tensor elements are <pq||rs>
    class I_to_T {
    public:
      static double I2T(double I, int i1, int i3, int i2, int i4,
      const RefDiagSCMatrix& evals1,
      const RefDiagSCMatrix& evals2,
      const RefDiagSCMatrix& evals3,
      const RefDiagSCMatrix& evals4)
      {
        return I;
      }
    };
    
    /// MP2 T2 tensor elements are <ij||ab> /(e_i + e_j - e_a - e_b)
    class ERI_to_T2 {
    public:
      static double I2T(double I, int i1, int i3, int i2, int i4,
      const RefDiagSCMatrix& evals1,
      const RefDiagSCMatrix& evals2,
      const RefDiagSCMatrix& evals3,
      const RefDiagSCMatrix& evals4)
      {
        const double denom = 1.0/(evals1(i1) + evals3(i3) - evals2(i2) - evals4(i4));
        return I*denom;
      }
    };
    
    /// MP2 pseudo-T2 (S2) tensor elements are <ij||ab> /sqrt(|e_i + e_j - e_a - e_b|) such
    /// that MP2 pair energies are the diagonal elements of S2 * S2.t()
    class ERI_to_S2 {
    public:
      static double I2T(double I, int i1, int i3, int i2, int i4,
                 const RefDiagSCMatrix& evals1,
                 const RefDiagSCMatrix& evals2,
                 const RefDiagSCMatrix& evals3,
                 const RefDiagSCMatrix& evals4)
      {
        const double denom = 1.0/sqrt(fabs(evals2(i2) + evals4(i4) - evals1(i1) - evals3(i3)));
        return I*denom;
      }
    };
    
  }
#endif

}

#endif

