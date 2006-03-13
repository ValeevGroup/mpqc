//
// ebc_contribs.cc
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
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/mbptr12/r12ia.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/r12_amps.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>
#include <chemistry/qc/mbptr12/print_scmat_norms.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>
#include <chemistry/qc/mbptr12/utils.h>
#include <chemistry/qc/mbptr12/utils.impl.h>

using namespace std;
using namespace sc;

#define TEST_T2 0
#define TEST_A 0
// if set to 1 then use f+k rather than f to compute A
#define A_DIRECT_EXCLUDE_K 0

//
// these are for testing purposes only
//
#define ACOMM_INCLUDE_TR_ONLY 0
#define ACOMM_INCLUDE_R_ONLY 0


void
R12IntEval::compute_A_direct_(RefSCMatrix& A,
                              const Ref<MOIndexSpace>& space1,
                              const Ref<MOIndexSpace>& space2,
                              const Ref<MOIndexSpace>& space3,
                              const Ref<MOIndexSpace>& space4,
                              const Ref<MOIndexSpace>& fspace2,
                              const Ref<MOIndexSpace>& fspace4,
                              bool antisymmetrize)
{
  // are particles 1 and 2 equivalent?
  const bool part1_equiv_part2 = (space1==space3 && space2 == space4);
  
  const unsigned int nf12 = corrfactor()->nfunctions();
  // create transforms, if needed
  std::vector< Ref<TwoBodyIntDescr> > descrs; // get 1 3 |F12| 2 4_f
  TwoBodyIntDescrCreator descr_creator(corrfactor(),
                                       r12info()->integral(),
                                       true,false);
  fill_container(descr_creator,descrs);
  
  tim_enter("A intermediate (direct)");
  std::ostringstream oss;
  oss << "<" << space1->id() << " " << space3->id() << "|A|"
      << space2->id() << " " << space4->id() << ">";
  const std::string label = oss.str();
  ExEnv::out0() << endl << indent
                << "Entered \"direct\" A intermediate (" << label << ") evaluator" << endl
                << incindent;
  //
  // ij|A|kl = ij|f12|kl_f, symmetrized if part1_equiv_part2
  //
  std::vector< Ref<TwoBodyMOIntsTransform> > tforms4f; // get 1 3 |F12| 2 4_f
  compute_F12_(A,space1,space2,space3,fspace4,antisymmetrize,tforms4f,descrs);
  if (part1_equiv_part2) {
    // no need to symmetrize if computing antisymmetric matrix -- compute_tbint_tensor takes care of that
    if (!antisymmetrize)
      symmetrize<false>(A,A,space1,space2);
    A.scale(2.0);
  }
  else {
    std::vector< Ref<TwoBodyMOIntsTransform> > tforms2f;
    compute_F12_(A,space1,fspace2,space3,space4,antisymmetrize,tforms2f,descrs);
  }

  ExEnv::out0() << decindent << indent << "Exited \"direct\" A intermediate (" << label << ") evaluator" << endl;
  tim_exit("A intermediate (direct)");
}

void
R12IntEval::compute_A_viacomm_(RefSCMatrix& A,
                               const Ref<MOIndexSpace>& space1,
                               const Ref<MOIndexSpace>& space2,
                               const Ref<MOIndexSpace>& space3,
                               const Ref<MOIndexSpace>& space4,
                               bool antisymmetrize,
                               const std::vector< Ref<TwoBodyMOIntsTransform> >& tforms)
{
  tim_enter("A intermediate (via commutator)");
  std::ostringstream oss;
  oss << "<" << space1->id() << " " << space3->id() << "|A|"
      << space2->id() << " " << space4->id() << ">";
  const std::string label = oss.str();
  ExEnv::out0() << endl << indent
                << "Entered \"commutator\" A intermediate (" << label << ") evaluator" << endl
                << incindent;
  
  // create descriptors, if needed
  std::vector< Ref<TwoBodyIntDescr> > descrs; // get 1 3 |F12| 2 4
  if (tforms.empty()) {
    TwoBodyIntDescrCreator descr_creator(corrfactor(),
                                         r12info()->integral(),
                                         true,false);
    fill_container(descr_creator,descrs);
  }
  
  //
  // ij|A|kl = - ij|[t1+t2,f12]|kl - (e_k+e_l-e_i-e_j) * ij|f12|kl
  //
  using ManyBodyTensors::Plus;
  using ManyBodyTensors::Minus;
#if !ACOMM_INCLUDE_TR_ONLY
  compute_tbint_tensor<ManyBodyTensors::Apply_H0minusE0<Minus>,true,false>(
    A, corrfactor()->tbint_type_f12(),
    space1, space2, space3, space4,
    antisymmetrize, tforms, descrs
  );
#endif
#if !ACOMM_INCLUDE_R_ONLY
  compute_tbint_tensor<ManyBodyTensors::Apply_Identity<Plus>,true,false>(
    A, corrfactor()->tbint_type_t1f12(),
    space1, space2, space3, space4,
    antisymmetrize, tforms, descrs
  );
  compute_tbint_tensor<ManyBodyTensors::Apply_Identity<Plus>,true,false>(
    A, corrfactor()->tbint_type_t2f12(),
    space1, space2, space3, space4,
    antisymmetrize, tforms, descrs
  );
#endif

  ExEnv::out0() << decindent << indent << "Exited \"commutator\" A intermediate (" << label << ") evaluator" << endl;
  tim_exit("A intermediate (via commutator)");
}

void
R12IntEval::AT2_contrib_to_V_()
{
  if (evaluated_)
    return;
  if (r12info_->msg()->me() == 0) {
    for(unsigned int s=0; s<nspincases2(); s++) {
      SpinCase2 spin = static_cast<SpinCase2>(s);
      
      // Use normal or commutator form of A, depending on the approach
      RefSCMatrix A;
      if (follow_ks_ebcfree_)
        A = Ac_[s];
      else
        A = A_[s];
      RefSCMatrix V = A * amps()->T2(spin).t();
      
      if (debug_ > 0) {
        std::string label = prepend_spincase(spin,"AT2 contribution to V");
        print_scmat_norms(V,label.c_str());
      }
      V_[s].accumulate(V);
    }
  }
  globally_sum_intermeds_();
}

void
R12IntEval::AF12_contrib_to_B_()
{
  if (evaluated_)
    return;
  if (r12info_->msg()->me() == 0) {
    for(unsigned int s=0; s<nspincases2(); s++) {
      SpinCase2 spin = static_cast<SpinCase2>(s);
      
      // Use normal or commutator form of A, depending on the approach
      RefSCMatrix A;
      if (follow_ks_ebcfree_)
        A = Ac_[s];
      else
        A = A_[s];
      RefSCMatrix AF = A * amps()->Fvv(spin).t();

      // B^{EBC} implies summation over all ab, not just unique ones, hence a factor of 2
      const double spin_pfac = 2.0;
      // minus 1/2 for Symmetrize AND -1/2 factor in B^{EBC} expression: B^{EBC} = -0.5 A . F12^t
      const double scale = -0.25 * spin_pfac;
      RefSCMatrix B = B_[s].clone();  B.assign(0.0);
      AF.scale(scale); B.accumulate(AF);
      RefSCMatrix AFt = AF.t();
      B.accumulate(AFt);
      
      const std::string label = prepend_spincase(spin,"B^{EBC} contribution");
      if (debug_ > 1) {
        B.print(label.c_str());
      }
      B_[s].accumulate(B);
      if (debug_ > 0) {
        print_scmat_norms(B,label.c_str());
      }
    }
  }
  globally_sum_intermeds_();
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
