//
// exchange.cc
//
// Copyright (C) 2004 Edward Valeev
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

#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <scconfig.h>
#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <math/scmat/blas.h>
#include <chemistry/qc/mbptr12/distarray4.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/print.h>

using namespace std;
using namespace sc;

RefSCMatrix R12IntEval::exchange_(const SpinCase1 &spin, const Ref<OrbitalSpace>& bra_space,
                      const Ref<OrbitalSpace>& ket_space) {

  Ref<MessageGrp> msg = r12world()->world()->msg();

  int me = msg->me();
  int nproc = msg->n();

  Ref<OrbitalSpace> contr_space;
  if(r12world()->sdref()) {
    contr_space = occ(spin);
    return(exchange_(contr_space,bra_space,ket_space));
  }
  else {
    contr_space = orbs(spin);
  }
  RefSymmSCMatrix opdm = this->ordm(spin);

  Timer tim_exchange("general exchange");
  ExEnv::out0() << endl << indent
                << "Entered general exchange matrix evaluator" << endl;
  ExEnv::out0() << incindent;
  const int ncontr = contr_space->rank();
  const int nbra = bra_space->rank();
  const int nket = ket_space->rank();
  const int nbraket = nbra*nket;
  if (debug_ >= DefaultPrintThresholds::fine) {
    ExEnv::out0() << indent << "nbra = " << nbra << endl;
    ExEnv::out0() << indent << "nket = " << nket << endl;
    ExEnv::out0() << indent << "ncontr = " << ncontr << endl;
  }

  // Do the AO->MO transform
  Ref<MOIntsTransformFactory> tfactory = this->tfactory();
  // Gaussians are real, hence contr_space and bra_space can be swapped
  tfactory->set_spaces(contr_space,bra_space,
                       contr_space,ket_space);
  // Only need 1/r12 integrals
  Ref<TwoBodyMOIntsTransform> rpsq_tform = tfactory->twobody_transform_13("(pr|qs)");
  rpsq_tform->compute();
  Ref<DistArray4> rpsq_acc = rpsq_tform->ints_distarray4();

  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  vector<int> proc_with_ints;
  const int nproc_with_ints = rpsq_acc->tasks_with_access(proc_with_ints);

  Timer tim_mo_ints_retrieve;
  tim_mo_ints_retrieve.set_default("MO ints retrieve");

  double* K_pq = new double[nbraket];
  memset(K_pq,0,nbraket*sizeof(double));
  if (rpsq_acc->has_access(me)) {
    for(int r=0; r<ncontr; r++) {
      for(int s=0; s<ncontr; s++) {
        int rs = r*ncontr+s;
        const int rs_proc = rs%nproc_with_ints;
        if (rs_proc != proc_with_ints[me])
          continue;

        if (debug_ >= DefaultPrintThresholds::fine)
          ExEnv::outn() << indent << "task " << me << ": working on (r,s) = " << r << "," << s << " " << endl;

        const double onepdm = opdm.get_element(s,r);

        // Get (|1/r12|) integrals
        tim_mo_ints_retrieve.enter_default();
        const double *rspq_buf_eri = rpsq_acc->retrieve_pair_block(r,s,corrfactor()->tbint_type_eri());
        tim_mo_ints_retrieve.exit_default();

        if (debug_ >= DefaultPrintThresholds::fine)
          ExEnv::outn() << indent << "task " << me << ": obtained mm block" << endl;

        //for(int p=0; p<nbra; p++) {
        //  for(int q=0; q<nket; q++) {
        //    int pq = p*nket+q;
        //    K_pq[pq] += rspq_buf_eri[pq]*onepdm;
        //  }
        //}
        const int unit_stride = 1;
        F77_DAXPY(&nbraket,&onepdm,rspq_buf_eri,&unit_stride,K_pq,&unit_stride);
        rpsq_acc->release_pair_block(r,s,corrfactor()->tbint_type_eri());
      }
    }
  }

  ExEnv::out0() << indent << "End of computation of exchange matrix" << endl;
  if (rpsq_acc->data_persistent()) rpsq_acc->deactivate();

  msg->sum(K_pq,nbraket);

  RefSCMatrix K(bra_space->coefs()->coldim(), ket_space->coefs()->coldim(), bra_space->coefs()->kit());
  K.assign(K_pq);
  delete [] K_pq;

  if (debug_ >= DefaultPrintThresholds::allN2) {
    K.print("Exchange matrix");
  }

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited general exchange matrix evaluator" << endl;

  tim_exchange.exit();

  return(K);
}

RefSCMatrix
R12IntEval::exchange_(const Ref<OrbitalSpace>& occ_space, const Ref<OrbitalSpace>& bra_space,
                      const Ref<OrbitalSpace>& ket_space)
{
  Ref<MessageGrp> msg = r12world()->world()->msg();

  Timer tim_exchange("exchange");

  int me = msg->me();
  int nproc = msg->n();
  ExEnv::out0() << endl << indent
	       << "Entered exchange matrix evaluator" << endl;
  ExEnv::out0() << incindent;

  // Do the AO->MO transform
  Ref<MOIntsTransformFactory> tfactory = this->tfactory();
  // Gaussians are real, hence occ_space and bra_space can be swapped
  tfactory->set_spaces(occ_space,bra_space,
                       occ_space,ket_space);
  // Only need 1/r12 integrals
  Ref<TwoBodyMOIntsTransform> mxny_tform = tfactory->twobody_transform_13("(mx|ny)");
  mxny_tform->compute();
  Ref<DistArray4> mnxy_acc = mxny_tform->ints_distarray4();
  mnxy_acc->activate();

  const int nocc = occ_space->rank();
  const int nbra = bra_space->rank();
  const int nket = ket_space->rank();
  const int nbraket = nbra*nket;

  ExEnv::out0() << indent << "Begin computation of exchange matrix" << endl;
  if (debug_ >= DefaultPrintThresholds::fine) {
    ExEnv::out0() << indent << "nbra = " << nbra << endl;
    ExEnv::out0() << indent << "nket = " << nket << endl;
    ExEnv::out0() << indent << "nocc = " << nocc << endl;
  }

  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  vector<int> proc_with_ints;
  int nproc_with_ints = mnxy_acc->tasks_with_access(proc_with_ints);

  //////////////////////////////////////////////////////////////
  //
  // Evaluation of the exchange matrix proceeds as follows:
  //
  //    loop over batches of mm, 0<=m<nocc
  //      load (mmxy)=(xm|my) into memory
  //
  //      loop over xy, 0<=x<nbra, 0<=y<nket
  //        compute K[x][y] += (mmxy)
  //      end xy loop
  //    end mm loop
  //
  /////////////////////////////////////////////////////////////////////////////////

  double* K_xy = new double[nbraket];
  memset(K_xy,0,nbraket*sizeof(double));
  Timer tim_mo_ints_retrieve;
  tim_mo_ints_retrieve.set_default("MO ints retrieve");
  if (mnxy_acc->has_access(me)) {

    for(int m=0; m<nocc; m++) {

      const int mm = m*nocc+m;
      const int mm_proc = mm%nproc_with_ints;
      if (mm_proc != proc_with_ints[me])
        continue;

      if (debug_ >= DefaultPrintThresholds::fine)
        ExEnv::outn() << indent << "task " << me << ": working on (m) = " << m << " " << endl;

      // Get (|1/r12|) integrals
      tim_mo_ints_retrieve.enter_default();
      const double *mmxy_buf_eri = mnxy_acc->retrieve_pair_block(m,m,corrfactor()->tbint_type_eri());
      tim_mo_ints_retrieve.exit_default();

      if (debug_ >= DefaultPrintThresholds::fine)
        ExEnv::outn() << indent << "task " << me << ": obtained mm block" << endl;

      const double one = 1.0;
      const int unit_stride = 1;
      F77_DAXPY(&nbraket,&one,mmxy_buf_eri,&unit_stride,K_xy,&unit_stride);

      mnxy_acc->release_pair_block(m,m,corrfactor()->tbint_type_eri());
    }
  }

  ExEnv::out0() << indent << "End of computation of exchange matrix" << endl;

  msg->sum(K_xy,nbraket);

  RefSCMatrix K(bra_space->coefs()->coldim(), ket_space->coefs()->coldim(), bra_space->coefs()->kit());
  K.assign(K_xy);
  delete[] K_xy;

  if (debug_ >= DefaultPrintThresholds::allN2) {
    K.print("Exchange matrix");
  }

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited exchange matrix evaluator" << endl;
  tim_exchange.exit();

  return K;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
