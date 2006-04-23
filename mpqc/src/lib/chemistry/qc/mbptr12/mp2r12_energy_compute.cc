//
// mp2r12_energy_compute.cc
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
#pragma implementation
#endif

#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#include <chemistry/qc/mbptr12/svd.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/utils.h>
#include <chemistry/qc/mbptr12/mp2r12_energy_util.h>

using namespace std;
using namespace sc;

#define USE_INVERT 0

void
MP2R12Energy::compute()
{
  if (evaluated_)
    return;
  
  Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  Ref<MessageGrp> msg = r12info->msg();
  int me = msg->me();
  int ntasks = msg->n();
  
  const bool ebc = r12eval()->ebc();
  /// KS approach to EBC-free method differs from mine if std approx != B || C
  const bool ks_ebcfree = r12eval()->ks_ebcfree() &&
                          (stdapprox_ != LinearR12::StdApprox_B &&
                           stdapprox_ != LinearR12::StdApprox_C);
  // Diagonal ansatz?
  const bool diag = r12info->ansatz()->diag();
  // WARNING only RHF and UHF are considered
  const int num_unique_spincases2 = (r12eval()->spin_polarized() ? 3 : 2);
  
  //
  // Evaluate pair energies:
  // distribute workload among nodes by pair index
  //
  for(int spin=0; spin<num_unique_spincases2; spin++) {
    
    const SpinCase2 spincase2 = static_cast<SpinCase2>(spin);
    
    Ref<LinearR12::NullCorrelationFactor> nullcorrptr; nullcorrptr << r12info->corrfactor();
    // if no explicit correlation -- just get MP2 energies
    if (nullcorrptr.nonnull()) {
      ef12_[spin].assign(0.0);
      RefSCVector emp2 = r12eval()->emp2(spincase2);
      double* buf = new double[emp2.dim().n()];
      emp2.convert(buf);
      emp2f12_[spin].assign(buf);
      delete[] buf;
    }
    else {
      
      const Ref<MOIndexSpace>& occ1_act = r12eval()->occ_act(case1(spincase2));
      const Ref<MOIndexSpace>& vir1_act = r12eval()->vir_act(case1(spincase2));
      const Ref<MOIndexSpace>& occ2_act = r12eval()->occ_act(case2(spincase2));
      const Ref<MOIndexSpace>& vir2_act = r12eval()->vir_act(case2(spincase2));
      int nocc1_act = occ1_act->rank();
      int nvir1_act = vir1_act->rank();
      int nocc2_act = occ2_act->rank();
      int nvir2_act = vir2_act->rank();
      const std::vector<double> evals_act_occ1 = convert(occ1_act->evals());
      const std::vector<double> evals_act_vir1 = convert(vir1_act->evals());
      const std::vector<double> evals_act_occ2 = convert(occ2_act->evals());
      const std::vector<double> evals_act_vir2 = convert(vir2_act->evals());
      
      // Get the intermediates V and X
      RefSCMatrix V = r12eval()->V(spincase2);
      RefSCMatrix X = r12eval()->X(spincase2);
      // The choice of B depends intricately on approximations
      RefSymmSCMatrix B;
      if (stdapprox_ == LinearR12::StdApprox_C) {
        B = r12eval()->BC(spincase2);
      }
      else {
        B = r12eval()->B(spincase2);
        // in standard approximation B, add up [K1+K2,F12] term
        if (stdapprox_ == LinearR12::StdApprox_B) {
          RefSymmSCMatrix BB = r12eval()->BB(spincase2);
          B.accumulate(BB);
        }
      }
      
      // In Klopper-Samson method, two kinds of A are used: app B (I replace with my A) and app XX (my Ac)
      RefSCMatrix A, Ac;
      if (ebc == false) {
        A = r12eval()->A(spincase2);
        if (ks_ebcfree) {
          Ac = r12eval()->Ac(spincase2);
        }
      }
      
      // Prepare total and R12 pairs
      const RefSCDimension dim_oo = V.coldim();
      const RefSCDimension dim_xc = V.rowdim();
      const int noo = dim_oo.n();
      if (noo == 0)
        continue;
      const int nxc = dim_xc.n();
      const int num_f12 = r12info->corrfactor()->nfunctions();
      
      Ref<MP2R12EnergyUtil_base> util;
      if (diag)
        util = new MP2R12EnergyUtil<true>(dim_oo, dim_xc);
      else
        util = new MP2R12EnergyUtil<false>(dim_oo, dim_xc);
      
      double* ef12_vec = new double[noo];
      memset(ef12_vec,0,sizeof(double)*noo);
      
      if (debug_ > 1) {
        util->print(prepend_spincase(spincase2,"V matrix").c_str(),V);
        util->print(prepend_spincase(spincase2,"MP2-F12 B matrix").c_str(),B);
        // A is too big to print out
      }
      
      // Prepare the B matrix:
      // 1) the B matrix is the same for all pairs in approximation A (EBC assumed) or if
      //    the ansatz is diagonal
      // 2) in approximations A', B, and C the B matrix is pair-specific
      const bool same_B_for_all_pairs = ( (stdapprox_ == LinearR12::StdApprox_A &&
                                           ebc == true) ||
                                          diag);
      
      if (same_B_for_all_pairs) {
        RefSymmSCMatrix B_ij = B.clone();
	B_ij.assign(B);
        SpinMOPairIter ij_iter(occ1_act, occ2_act, spincase2);
        for(ij_iter.start(); int(ij_iter); ij_iter.next()) {
          const int ij = ij_iter.ij();
          const int i = ij_iter.i();
          const int j = ij_iter.j();
          
          for(int f=0; f<num_f12; f++) {
            const int f_off = f*noo;
            const int ijf = f_off + ij;
            
            for(int g=0; g<=f; g++) {
              const int g_off = g*noo;
              const int ijg = g_off + ij;
              
              // contribution from X still appears in approximations C and A''
              if (stdapprox_ != LinearR12::StdApprox_A) {
                if (stdapprox_ == LinearR12::StdApprox_C ||
		    stdapprox_ == LinearR12::StdApprox_App) {
                  double fx = - (evals_act_occ1[i] + evals_act_occ2[j]) * X.get_element(ijf,ijg);
                  B_ij.accumulate_element(ijf,ijg,fx);
		}
              }
                  
              // If EBC is not assumed add 2.0*Akl,cd*Acd,ow/(ec+ed-ei-ej)
              if (ebc == false) {
                double fy = 0.0;
                SpinMOPairIter cd_iter(vir1_act, vir2_act, spincase2);
                if (ks_ebcfree) {
                  for(cd_iter.start(); cd_iter; cd_iter.next()) {
                    const int cd = cd_iter.ij();
                    const int c = cd_iter.i();
                    const int d = cd_iter.j();
                    
                    fy -= 0.5*( A.get_element(ijf,cd)*Ac.get_element(ijg,cd) + Ac.get_element(ijf,cd)*A.get_element(ijg,cd) )/
                    (evals_act_vir1[c] + evals_act_vir2[d]
                    - evals_act_occ1[i] - evals_act_occ2[j]);
                  }
                }
                else {
                  for(cd_iter.start(); cd_iter; cd_iter.next()) {
                    const int cd = cd_iter.ij();
                    const int c = cd_iter.i();
                    const int d = cd_iter.j();
                    
                    fy -= A.get_element(ijf,cd)*A.get_element(ijg,cd)/(evals_act_vir1[c] + evals_act_vir2[d]
                    - evals_act_occ1[i] - evals_act_occ2[j]);
                  }
                }
                
                B_ij.accumulate_element(ijf,ijg,fy);
              }
            }
          }
        }
        if (debug_ > 1)
          util->print(prepend_spincase(spincase2,"MP2-F12 B matrix").c_str(),B_ij);
#if USE_INVERT
        util->invert(B_ij);
        if (debug_ > 1)
          util->print("Inverse MP2-F12/A B matrix",B_ij);
        RefSCMatrix C = C_[spin].clone();
        util->times(B_ij,V,C);
        C_[spin].assign(C);  C = 0;
        C_[spin].scale(-1.0);
#else
        // solve B * C = V
        RefSCMatrix C = C_[spin].clone();
        util->solve_linear_system(B_ij, C, V);
        C_[spin].assign(C);  C = 0;
        C_[spin].scale(-1.0);
#endif
      }
      // B is pair specific
      //
      else {
        RefSymmSCMatrix B_ij = B.clone();
        SpinMOPairIter ij_iter(occ1_act, occ2_act, spincase2);
        for(ij_iter.start(); int(ij_iter); ij_iter.next()) {
          const int ij = ij_iter.ij();
          if (ij%ntasks != me)
            continue;
          const int i = ij_iter.i();
          const int j = ij_iter.j();
          
          // 
          
          // In approximations A', B, or C matrices B are pair-specific:
          // app A' or B: form B(ij)kl,ow = Bkl,ow + 1/2(ek + el + eo + ew - 2ei - 2ej)Xkl,ow
          // app C:       form B(ij)kl,ow = Bkl,ow - (ei + ej)Xkl,ow
          B_ij.assign(B);
          
          for(int f=0; f<num_f12; f++) {
            const int f_off = f*noo;
            
            SpinMOPairIter kl_iter(occ1_act, occ2_act, spincase2);
            for(kl_iter.start(); kl_iter; kl_iter.next()) {
              const int kl = kl_iter.ij() + f_off;
              const int k = kl_iter.i();
              const int l = kl_iter.j();
              
              for(int g=0; g<=f; g++) {
                const int g_off = g*noo;
                
                SpinMOPairIter ow_iter(occ1_act, occ2_act, spincase2);
                for(ow_iter.start(); ow_iter; ow_iter.next()) {
                  const int ow = ow_iter.ij() + g_off;
                  const int o = ow_iter.i();
                  const int w = ow_iter.j();
                  
                  // This is a symmetric matrix
                  if (ow > kl)
                    continue;

                  // contribution from X vanishes in approximation A
                  if (stdapprox_ != LinearR12::StdApprox_A) {
                    double fx;
                    if (stdapprox_ != LinearR12::StdApprox_C &&
                      stdapprox_ != LinearR12::StdApprox_App)
                    fx = 0.5 * (evals_act_occ1[k] + evals_act_occ2[l] + evals_act_occ1[o] + evals_act_occ2[w]
                    - 2.0*evals_act_occ1[i] - 2.0*evals_act_occ2[j])
                    * X.get_element(kl,ow);
                    else
                      fx = - (evals_act_occ1[i] + evals_act_occ2[j]) * X.get_element(kl,ow);
                    
                    B_ij.accumulate_element(kl,ow,fx);
                  }
                  
                  // If EBC is not assumed add 2.0*Akl,cd*Acd,ow/(ec+ed-ei-ej)
                  if (ebc == false) {
                    double fy = 0.0;
                    SpinMOPairIter cd_iter(vir1_act, vir2_act, spincase2);
                    if (ks_ebcfree) {
                      for(cd_iter.start(); cd_iter; cd_iter.next()) {
                        const int cd = cd_iter.ij();
                        const int c = cd_iter.i();
                        const int d = cd_iter.j();
                        
                        fy -= 0.5*( A.get_element(kl,cd)*Ac.get_element(ow,cd) + Ac.get_element(kl,cd)*A.get_element(ow,cd) )/
                          (evals_act_vir1[c] + evals_act_vir2[d]
                           - evals_act_occ1[i] - evals_act_occ2[j]);
                      }
                    }
                    else {
                      for(cd_iter.start(); cd_iter; cd_iter.next()) {
                        const int cd = cd_iter.ij();
                        const int c = cd_iter.i();
                        const int d = cd_iter.j();
                        
                        fy -= A.get_element(kl,cd)*A.get_element(ow,cd)/(evals_act_vir1[c] + evals_act_vir2[d]
                                                                         - evals_act_occ1[i] - evals_act_occ2[j]);
                      }
                    }
                    
                    B_ij.accumulate_element(kl,ow,fy);
                  }
                }
              }
            }
          }
          if (debug_ > 1)
            util->print(prepend_spincase(spincase2,"MP2-F12 B matrix").c_str(),B_ij);

          /// Block in which I compute ef12
          {
            RefSCVector V_ij = V.get_column(ij);
#if USE_INVERT
            // The r12 amplitudes B^-1 * V
            util->invert(B_ij);
            if (debug_ > 1)
              util->print("Inverse MP2-F12 B matrix",B_ij);
            RefSCVector Cij = -1.0*(B_ij * V_ij);
            for(int kl=0; kl<nxc; kl++)
              C_[spin].set_element(kl,ij,Cij.get_element(kl));
#else
            RefSCVector Cij = V_ij.clone();
            // solve B * C = V
            Cij = V_ij.clone();
            sc::exp::lapack_linsolv_symmnondef(B_ij, Cij, V_ij);
            Cij.scale(-1.0);
            for(int kl=0; kl<nxc; kl++)
              C_[spin].set_element(kl,ij,Cij.get_element(kl));
#endif
          }
        }
      } // pair-specific B block
      
      // Compute pair energies
      RefSCVector ef12 = util->dot_product(C_[spin],V);
      for(unsigned int ij=0; ij<noo; ij++)
        if (ij%ntasks == me)
          ef12_vec[ij] = ef12(ij);
      msg->sum(ef12_vec,noo,0,-1);
      ef12_[spin]->assign(ef12_vec);

      RefSCVector emp2 = r12eval()->emp2(spincase2);
      double* buf = new double[emp2.dim().n()];
      emp2.convert(buf);
      emp2f12_[spin].assign(buf);
      delete[] buf;
      emp2f12_[spin]->accumulate(ef12_[spin]);

      delete[] ef12_vec;
    }
    
  } // end of spincase loop
  
  // Set beta-beta energies to alpha-alpha for closed-shell
  if (!r12info->refinfo()->ref()->spin_polarized()) {
    emp2f12_[BetaBeta] = emp2f12_[AlphaAlpha];
    ef12_[BetaBeta] = ef12_[AlphaAlpha];
  }

  evaluated_ = true;
  
  return;
}

namespace {

  // Extracts the diagonal d of a symmetric matrix A
  void diagonal(const RefSymmSCMatrix& A, RefSCVector& d)
  {
    const int n = A.dim().n();
    // allocate d if necessary
    if (d.null())
      d = A.kit()->vector(A.dim());
    else {
      if (d.dim().n() != A.dim().n())
        throw ProgrammingError("anonymous::diagonal() -- dimensions of d and A don't match",__FILE__,__LINE__);
    }
    for(int i=0; i<n; i++)
      d[i] = A.get_element(i,i);
  }

  // Extracts the "diagonal" d of a matrix A -- the diagonal is defined as concatenation of diagonals of square blocks
  void diagonal(const RefSCMatrix& A, RefSCVector& d)
  {
    const int nrow = A.rowdim().n();
    const int ncol = A.coldim().n();
    const int nmin = std::min(nrow,ncol);
    if ( nrow%nmin && ncol%nmin)
      throw ProgrammingError("anonymous::diagonal() -- dimensions of A are not perfect divisors of each other");
    const int nmax = std::max(nrow,ncol);
    
    // allocate d if necessary
    if (d.null())
      d = A.kit()->vector(new SCDimension(nmax));
    for(int i=0; i<nmax; i++)
      d[i] = A.get_element(i%nrow,i%ncol);
  }

  // Puts back the "diagonal" d into matrix A -- the diagonal is defined as concatenation of diagonals of square blocks
  void diagonal(const RefSCVector& d, RefSCMatrix& A)
  {
    const int nrow = A.rowdim().n();
    const int ncol = A.coldim().n();
    const int nmin = std::min(nrow,ncol);
    if ( nrow%nmin && ncol%nmin)
      throw ProgrammingError("anonymous::diagonal() -- dimensions of A are not perfect divisors of each other");
    const int nmax = std::max(nrow,ncol);
    if (nmax != d.dim().n())
      throw ProgrammingError("anonymous::diagonal() -- dimension of d and A do not match",__FILE__,__LINE__);

    A.assign(0.0);
    for(int i=0; i<nmax; i++)
      A.set_element(i%nrow,i%ncol,d[i]);
  }

  // Prints out matrix A
  template <typename Matrix> void print(bool diag_only, const char* label, Matrix& A) {
    if (!diag_only)
      A.print(label);
    else {
      RefSCVector d;
      diagonal(A,d);
      d.print(label);
    }
  }
  
  template <typename Matrix> void invert(bool diag_only, Matrix& A) {
    if (!diag_only)
      A->gen_invert_this();
    else {
      RefSCVector d;
      diagonal(A,d);
      const int nd = d.dim().d();
      for(int i=0; i<nd; i++) {
        const double oovalue = 1.0/d.get_element(i);
        d.set_element(oovalue);
      }
      diagonal(d,A);
    }
  }

  // Solves A*X = B
  void solve_linear_system(bool diag_only, const RefSymmSCMatrix& A, RefSCMatrix& X, const RefSCMatrix& B)
  {
    if (!diag_only)
      sc::exp::lapack_linsolv_symmnondef(A, X, B);
    else {
      RefSCVector Ad; diagonal(A,Ad);
      RefSCVector Xd; diagonal(X,Xd);
      RefSCVector Bd; diagonal(B,Bd);
      const int n = Ad.dim().n();
      if (n != Xd.dim().n() || n != Bd.dim().n())
        throw ProgrammingError("anonymous::solve_linear_system() -- dimensions don't match");
      for(int i=0; i<n; i++) {
        const double a = Ad.get_element(i);
        const double b = Bd.get_element(i);
        const double tol = 1.0e-10;
        if (a < tol)
          throw ToleranceExceeded("ananymous::solve_linear_system() -- a is below tolerance or negative",__FILE__,__LINE__,tol,a);
        Xd.set_element(i,b/a);
      }
      diagonal(Xd,X);
    }
  }

  // computes y = A x
  void times(bool diag_only, const RefSymmSCMatrix& A, const RefSCMatrix& x, RefSCMatrix& y)
  {
    if (!diag_only)
      y = A*x;
    else {
      RefSCVector Ad; diagonal(A,Ad);
      RefSCVector xd; diagonal(A,xd);
      RefSCVector yd = xd.kit()->vector(xd.dim());
      const int n = Ad.dim().n();
      if (n != xd.dim().n())
        throw ProgrammingError("anonymous::times() -- dimensions don't match");
      for(int i=0; i<n; i++) {
        yd(i) = Ad(i) * xd(i);
      }
      diagonal(yd,y);
    }
  }
  
}
