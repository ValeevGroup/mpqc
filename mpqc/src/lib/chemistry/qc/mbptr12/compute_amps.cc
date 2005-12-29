//
// compute_amps.cc
//
// Copyright (C) 2004 Edward Valeev
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

#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <scconfig.h>
#include <util/misc/formio.h>
#include <util/misc/timer.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/matrix.h>
#include <math/scmat/local.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/mbptr12/r12ia.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>

using namespace std;
using namespace sc;

void
R12IntEval::compute_T2_(RefSCMatrix& T2,
                        const Ref<MOIndexSpace>& space1,
                        const Ref<MOIndexSpace>& space2,
                        const Ref<MOIndexSpace>& space3,
                        const Ref<MOIndexSpace>& space4,
                        const Ref<TwoBodyMOIntsTransform>& transform)
{
  T2.assign(0.0);
  // maps spaceX to spaceX of the transform
  std::vector<unsigned int> map1, map2, map3, map4;
  Ref<TwoBodyMOIntsTransform> tform = transform;
  if (tform.null()) {
    // Only need 1/r12 integrals, hence doesn't matter which f12 to use
    tform = get_tform_(transform_label(space1,space2,space3,space4,0));
  }
  Ref<MOIndexSpace> tspace1 = tform->space1();
  Ref<MOIndexSpace> tspace2 = tform->space2();
  Ref<MOIndexSpace> tspace3 = tform->space3();
  Ref<MOIndexSpace> tspace4 = tform->space4();
  map1 = *tspace1<<*space1;
  map2 = *tspace2<<*space2;
  map3 = *tspace3<<*space3;
  map4 = *tspace4<<*space4;
  
  const unsigned int rank1 = space1->rank();
  const unsigned int rank2 = space2->rank();
  const unsigned int rank3 = space3->rank();
  const unsigned int rank4 = space4->rank();
  const unsigned int trank4 = tspace4->rank();
  const RefDiagSCMatrix evals1 = space1->evals();
  const RefDiagSCMatrix evals2 = space2->evals();
  const RefDiagSCMatrix evals3 = space3->evals();
  const RefDiagSCMatrix evals4 = space4->evals();
  
  // Using spinorbital iterators means I don't take into account perm symmetry
  // More efficient algorithm will require generic code
  SpinMOPairIter iter13(space1,space3,AlphaBeta);
  SpinMOPairIter iter24(space2,space4,AlphaBeta);
  
  Ref<R12IntsAcc> accum = tform->ints_acc();
  if (accum.null() || !accum->is_committed()) {
    // only need ERIs
    Ref<TwoBodyIntDescr> tbintdescr = new TwoBodyIntDescrERI(r12info()->integral());
    tform->compute(tbintdescr);
  }
  if (!accum->is_active())
    accum->activate();
  
  tim_enter("T2 amplitudes");
  std::ostringstream oss;
  oss << "<" << space1->id() << " " << space3->id() << "|T2|"
      << space2->id() << " " << space4->id() << ">";
  const std::string label = oss.str();
  ExEnv::out0() << endl << indent
	       << "Entered MP2 T2 amplitude (" << label << ") evaluator" << endl;
  ExEnv::out0() << incindent;
  if (debug_ > 0)
    ExEnv::out0() << indent << "Using transform " << tform->name() << std::endl;
  
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
      const double *ij_buf_eri = accum->retrieve_pair_block(ii,jj,corrfactor()->tbint_type_eri());
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
        
        const double ERI_iajb = ij_buf_eri[AB];
#if 1
        const double denom = 1.0/(evals1(i) + evals3(j) - evals2(a) - evals4(b));
#else
        // use this to test T2
        const double denom = 1.0/sqrt(fabs(evals1(i) + evals3(j) - evals2(a) - evals4(b)));
#endif
        
        if (debug_ > 2) {
          ExEnv::out0() << "i = " << i << " j = " << j << " a = " << a << " b = " << b
          << " <ij|ab> = " << ERI_iajb
          << " denom = " << denom << endl;
        }
        
        T2.set_element(ij,ab,ERI_iajb*denom);
      }
      accum->release_pair_block(ii,jj,corrfactor()->tbint_type_eri());
      
    }
  }
  
  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited MP2 T2 amplitude (" << label << ") evaluator" << endl;
  tim_exit("T2 amplitudes");
}


void
R12IntEval::compute_F12_(RefSCMatrix& F12,
                        const Ref<MOIndexSpace>& space1,
                        const Ref<MOIndexSpace>& space2,
                        const Ref<MOIndexSpace>& space3,
                        const Ref<MOIndexSpace>& space4,
                        const std::vector< Ref<TwoBodyMOIntsTransform> >& transvec)
{
  F12.assign(0.0);
  const unsigned int nf12 = corrfactor()->nfunctions();
  // create transforms, if needed
  std::vector< Ref<TwoBodyMOIntsTransform> > transforms = transvec;
  if (transforms.empty()) {
    for(unsigned int f12=0; f12<nf12; f12++) {
      transforms.push_back(get_tform_(transform_label(space1,space2,space3,space4,f12)));
    }
  }
  
  tim_enter("F12 amplitudes");
  std::ostringstream oss;
  oss << "<" << space1->id() << " " << space3->id() << "|F12|"
      << space2->id() << " " << space4->id() << ">";
  const std::string label = oss.str();
  ExEnv::out0() << endl << indent
                << "Entered F12 amplitude (" << label << ") evaluator" << endl
                << incindent;

  // Using spinorbital iterators means I don't take into account perm symmetry
  // More efficient algorithm will require generic code
  SpinMOPairIter iter13(space1,space3,AlphaBeta);
  SpinMOPairIter iter24(space2,space4,AlphaBeta);
  const unsigned int rank1 = space1->rank();
  const unsigned int rank2 = space2->rank();
  const unsigned int rank3 = space3->rank();
  const unsigned int rank4 = space4->rank();
  // size of one block of <space1 space3|
  const unsigned int n13 = rank1*rank3;
  
  // counts how many rows of F12 have written
  unsigned int f12offset = 0;
  for(unsigned int f12=0; f12<nf12; f12++,f12offset+=n13) {
    Ref<TwoBodyMOIntsTransform> tform = transforms[f12];
    if (debug_ > 0)
      ExEnv::out0() << indent << "Using transform " << tform->name() << std::endl;

    // maps spaceX to spaceX of the transform
    std::vector<unsigned int> map1, map2, map3, map4;
    Ref<MOIndexSpace> tspace1 = tform->space1();
    Ref<MOIndexSpace> tspace2 = tform->space2();
    Ref<MOIndexSpace> tspace3 = tform->space3();
    Ref<MOIndexSpace> tspace4 = tform->space4();
    map1 = *tspace1<<*space1;
    map2 = *tspace2<<*space2;
    map3 = *tspace3<<*space3;
    map4 = *tspace4<<*space4;
    const unsigned int trank4 = tspace4->rank();
    
    Ref<R12IntsAcc> accum = tform->ints_acc();
    if (accum.null() || !accum->is_committed()) {
      // need F12 integral
      Ref<TwoBodyIntDescr> tbintdescr = corrfactor()->tbintdescr(r12info()->integral(),f12);
      tform->compute(tbintdescr);
    }
    if (!accum->is_active())
      accum->activate();
    
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
        const double *ij_buf_f12 = accum->retrieve_pair_block(ii,jj,corrfactor()->tbint_type_f12());
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
          
          const double F12_iajb = ij_buf_f12[AB];
          
          if (debug_ > 2) {
            ExEnv::out0() << "i = " << i << " j = " << j << " a = " << a << " b = " << b
            << " <ij| f12[" << f12 << "]| ab> = " << F12_iajb << endl;
          }
          
          F12.set_element(ij+f12offset,ab,F12_iajb);
        }
        accum->release_pair_block(ii,jj,corrfactor()->tbint_type_f12());
        
      }
    }
    
  }
  ExEnv::out0() << decindent << indent << "Exited F12 amplitude (" << label << ") evaluator" << endl;
  tim_exit("F12 amplitudes");
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
