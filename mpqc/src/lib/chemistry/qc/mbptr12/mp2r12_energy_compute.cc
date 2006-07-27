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
#include <chemistry/qc/mbptr12/print.h>

using namespace std;
using namespace sc;

#define USE_INVERT 0

namespace {
  /// Eigenvalues statistics
  struct EvalStats {
    EvalStats(const RefDiagSCMatrix& eigenvalues, double thresh, unsigned int nevals_below, double emin, double emax) :
      evals(eigenvalues), threshold(thresh), num_below_threshold(nevals_below), min(emin), max(emax) {}
    RefDiagSCMatrix evals;
    double threshold;
    unsigned int num_below_threshold;
    double min;
    double max;
  };
  /// Compute EvalsStats for eigenvalues evals, given threshold (default -- 0.0)
  EvalStats eigenvalue_stats(const RefDiagSCMatrix& evals, double threshold = 0.0);
  /// Print eigenstats
  void print_eigenstats(const EvalStats& stats,
			const char* label,
			const Ref<MP2R12EnergyUtil_base>& util,
			int debug,
			bool warn,
			std::ostream& os = ExEnv::out0());
};

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
  // 1) the B matrix is the same for all pairs in approximation A (EBC assumed) or if
  //    the ansatz is diagonal
  // 2) in approximations A', B, and C the B matrix is pair-specific
  const bool same_B_for_all_pairs = ( (stdapprox_ == LinearR12::StdApprox_A && 
                                       ebc == true) ||
                                      diag);
  // check positive definiteness of B? -- cannot yet do for diagonal ansatz
  const bool check_posdef_B = (r12info->safety_check() && !diag);
  // make B positive definite? -- cannot yet do for diagonal ansatz
  const bool posdef_B = ((r12info->posdef_B() != LinearR12::PositiveDefiniteB_no) && !diag && same_B_for_all_pairs);
  // make ~B(ij) positive definite? -- cannot yet do for diagonal ansatz
  const bool posdef_Bij = ((r12info->posdef_B() == LinearR12::PositiveDefiniteB_yes) && !diag && !same_B_for_all_pairs);
  if (r12info->safety_check()) {
    if (diag)
      ExEnv::out0() << indent << "WARNING: cannot yet check positive definitess of the B matrix when using the diagonal ansatz." << endl;
    if ( (r12info->posdef_B() != LinearR12::PositiveDefiniteB_no) && diag)
      ExEnv::out0() << indent << "WARNING: cannot yet ensure positive definite B matrix when using the diagonal ansatz." << endl;
  }

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
      RefSymmSCMatrix X = r12eval()->X(spincase2);
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

      // util class treats abstractly on dense and "diagonal" matrices used in orb-invariant and non-orb-invariant ansatze
      Ref<MP2R12EnergyUtil_base> util;
      if (diag)
        util = new MP2R12EnergyUtil<true>(dim_oo, dim_xc);
      else
        util = new MP2R12EnergyUtil<false>(dim_oo, dim_xc);
      
      double* ef12_vec = new double[noo];
      memset(ef12_vec,0,sizeof(double)*noo);
      
      if (debug_ >= DefaultPrintThresholds::mostO4) {
        util->print(prepend_spincase(spincase2,"V matrix").c_str(),V);
        util->print(prepend_spincase(spincase2,"X matrix").c_str(),X);
        util->print(prepend_spincase(spincase2,"B matrix (excludes X and coupling)").c_str(),B);
        // A is too big to print out
      }

      //
      // SAFETY FIRST! Check that X is nonsigular and X and B are positive definite. Rectify, if needed
      //
      unsigned int nlindep_g12 = 0;
      bool lindep_g12 = false;
      unsigned int nnegevals_B = 0;
      bool negevals_B = false;
      RefSCMatrix UX;                // orthonormalizes the geminal space
      if (r12info->safety_check()) {

	//
	// Check if XC pair basis is linearly dependent and compute the orthonormalizer, if possible
	//
	if (!diag) {
	  ExEnv::out0() << indent << "Computing orthonormal " << prepend_spincase(spincase2,"geminal space.",ToLowerCase) << endl << incindent;
	  Ref<SCF> ref = r12info->refinfo()->ref();
	  OverlapOrthog orthog(OverlapOrthog::Canonical,
			       X,
			       X.kit(),
			       ref->lindep_tol());
	  nlindep_g12 = orthog.nlindep();
	  UX = orthog.basis_to_orthog_basis();
	  ExEnv::out0() << decindent;
	}
	else {
	  // For diagonal ansatz cannot yet get rid of linear dependencies, only check for the presence
	  RefDiagSCMatrix Xevals = util->eigenvalues(X);
	  const double thresh = 1.0e-12;
	  EvalStats stats = eigenvalue_stats(Xevals,thresh);
	  print_eigenstats(stats,prepend_spincase(spincase2,"X matrix").c_str(),
			   util,debug_,true);
	  if (stats.num_below_threshold)
	    ExEnv::out0() << indent << "WARNING: Cannot yet get rid of linear dependencies for the nonorbital invariant ansatz yet." << endl;
	  nlindep_g12 = 0;
	}
	lindep_g12 = (nlindep_g12 > 0);

	// Check if B matrix is positive definite -- check eigenvalues of B
	if (check_posdef_B) {
	  RefSymmSCMatrix Borth = B.kit()->symmmatrix(UX.rowdim());
	  Borth.assign(0.0);
	  Borth.accumulate_transform(UX,B);

	  RefDiagSCMatrix Bevals;
	  RefSCMatrix Bevecs;
	  if (posdef_B) {
	    util->diagonalize(Borth,Bevals,Bevecs);
	  }
	  else {
	    Bevals = util->eigenvalues(Borth);
	  }
	  Borth = 0;
	  EvalStats stats = eigenvalue_stats(Bevals,0.0);
	  nnegevals_B = stats.num_below_threshold;
	  // only warn if !posdef_B
	  print_eigenstats(stats,prepend_spincase(spincase2,"B matrix in orthonormal basis").c_str(),
			   util,debug_,!posdef_B);
	  negevals_B = (nnegevals_B > 0);

	  // If found negative eigenvalues, throw away the vectors corresponding to negative eigenvalues
	  if (posdef_B && negevals_B) {
	    RefSCDimension bevdim = new SCDimension(Bevals.dim().n() - nnegevals_B);
	    RefSCMatrix V = B.kit()->matrix(Bevecs.rowdim(),bevdim);
	    const unsigned int nbevdim = bevdim.n();
	    unsigned int ivec = nnegevals_B;
	    for(unsigned int i=0; i<nbevdim; ++i, ++ivec) {
	      RefSCVector col = Bevecs.get_column(ivec);
	      V.assign_column(col,i);
	      col = 0;
	    }
	    
	    RefSCMatrix Vt = V.t(); V=0;
	    UX = Vt * UX;
	    
	    // Test new transform matrix
	    if (debug_ >= DefaultPrintThresholds::allO4) {
	      RefSymmSCMatrix Borth = B.kit()->symmmatrix(UX.rowdim());
	      Borth.assign(0.0);
	      Borth.accumulate_transform(UX,B);

	      RefDiagSCMatrix Bevals = util->eigenvalues(Borth);
	      Borth = 0;
	      Bevals.print("Eigenvalues of orthonormalized B matrix");

	      RefSymmSCMatrix Xorth = X.kit()->symmmatrix(UX.rowdim());
	      Xorth.assign(0.0);
	      Xorth.accumulate_transform(UX,X);
	      Xorth.print("Should be 1: UX * X * UX.t()");
	    }

	    ExEnv::out0() << indent << "WARNING: eliminated " << nnegevals_B << " eigenvalues from B matrix" << endl;

	  } // getting rid of negative eigenvalues of B
	} // check of positive definiteness of B

      } // safety checks of X and B matrices
      
      // Prepare the B matrix:
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
        if (debug_ >= DefaultPrintThresholds::mostO4)
          util->print(prepend_spincase(spincase2,"~B matrix").c_str(),B_ij);

	// Compute first-order wave function
	const bool need_to_transform_f12dim = (lindep_g12 || (posdef_B && negevals_B));
#if USE_INVERT
	if (need_to_transform_f12dim && r12info->safety_check())
	  throw ProgrammingError("MP2R12Energy::compute() -- safe evaluation of the MP2-R12 energy using inversion is not implemented yet");
        util->invert(B_ij);
        if (debug_ >= DefaultPrintThresholds::mostO4)
          util->print("Inverse MP2-F12/A B matrix",B_ij);
        RefSCMatrix C = C_[spin].clone();
        util->times(B_ij,V,C);
        C_[spin].assign(C);  C = 0;
        C_[spin].scale(-1.0);
#else
        // solve B * C = V
	if (!need_to_transform_f12dim) {
	  RefSCMatrix C = C_[spin].clone();
	  util->solve_linear_system(B_ij, C, V);
	  C_[spin].assign(C);  C = 0;
	}
	else {
	  RefSCMatrix Vo = UX * V;
	  RefSymmSCMatrix Bo = B_ij.kit()->symmmatrix(UX.rowdim());
	  Bo.assign(0.0);
	  Bo.accumulate_transform(UX,B_ij);
	  RefSCMatrix Co = B_ij.kit()->matrix(UX.rowdim(),C_[spin].coldim());
	  util->solve_linear_system(Bo, Co, Vo);
	  RefSCMatrix C = UX.t() * Co;
	  C_[spin].assign(C);
	}
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
	  std::string B_ij_label;
	  {
	    ostringstream oss; oss << "~B(" << ij << ") matrix";
	    B_ij_label = oss.str();
	  }
          if (debug_ >= DefaultPrintThresholds::mostO4) {
            util->print(prepend_spincase(spincase2,B_ij_label).c_str(),B_ij);
	  }

          /// Block in which I compute ef12
          {
            RefSCVector V_ij = V.get_column(ij);

	    //
	    // Check if ~B(ij) is positive definite
	    //
	    RefSCMatrix UXB = UX;
	    unsigned int nnegevals_B_ij = 0;
	    bool negevals_B_ij = false;
	    if (r12info->safety_check() && check_posdef_B) {
	      RefSymmSCMatrix Borth = B.kit()->symmmatrix(UX.rowdim());
	      Borth.assign(0.0);
	      Borth.accumulate_transform(UX,B_ij);
	      
	      RefDiagSCMatrix Bevals;
	      RefSCMatrix Bevecs;
	      if (posdef_Bij) {
		util->diagonalize(Borth,Bevals,Bevecs);
	      }
	      else {
		Bevals = util->eigenvalues(Borth);
	      }
	      Borth = 0;
	      EvalStats stats = eigenvalue_stats(Bevals,0.0);
	      nnegevals_B_ij = stats.num_below_threshold;
	      negevals_B_ij = (nnegevals_B_ij > 0);
	      {
		ostringstream oss;  oss << B_ij_label << " in orthonormal basis";
		// only warn if !posdef_Bij
		print_eigenstats(stats,prepend_spincase(spincase2,oss.str()).c_str(),
				 util,debug_,!posdef_Bij);
	      }

	      // If found negative eigenvalues, throw away the vectors corresponding to negative eigenvalues
	      if (posdef_Bij && negevals_B_ij) {
		RefSCDimension bevdim = new SCDimension(Bevals.dim().n() - nnegevals_B_ij);
		RefSCMatrix V = B.kit()->matrix(Bevecs.rowdim(),bevdim);
		const unsigned int nbevdim = bevdim.n();
		unsigned int ivec = nnegevals_B_ij;
		for(unsigned int i=0; i<nbevdim; ++i, ++ivec) {
		  RefSCVector col = Bevecs.get_column(ivec);
		  V.assign_column(col,i);
		  col = 0;
		}
		
		RefSCMatrix Vt = V.t(); V=0;
		UXB = Vt * UX;

		ExEnv::out0() << indent << "WARNING: eliminated " << nnegevals_B_ij << " eigenvalues from " << B_ij_label << endl;

	      } // getting rid of negative eigenvalues of B
	    } // check of positive definiteness of B

	    const bool need_to_transform_f12dim = (lindep_g12 || (posdef_Bij && negevals_B_ij));
#if USE_INVERT
	    if (need_to_transform_f12dim && r12info->safety_check())
	      throw ProgrammingError("MP2R12Energy::compute() -- safe evaluation of the MP2-R12 energy using inversion is not implemented yet");
            // The r12 amplitudes B^-1 * V
            util->invert(B_ij);
            if (debug_ >= DefaultPrintThresholds::allO4) {
	      osringstream oss; oss << "Inverse MP2-F12 matrix B(" << ij << ")";
              util->print(oss.str(),B_ij);
	    }
            RefSCVector Cij = -1.0*(B_ij * V_ij);
            for(int kl=0; kl<nxc; kl++)
              C_[spin].set_element(kl,ij,Cij.get_element(kl));
#else
	    // solve B * C = V
	    if (!need_to_transform_f12dim) {
	      RefSCVector Cij = V_ij.clone();
	      sc::exp::lapack_linsolv_symmnondef(B_ij, Cij, V_ij);
	      Cij.scale(-1.0);
	      for(int kl=0; kl<nxc; kl++)
		C_[spin].set_element(kl,ij,Cij.get_element(kl));
	    }
	    else {
	      RefSCVector Vo = UXB * V_ij;
	      RefSymmSCMatrix Bo = B_ij.kit()->symmmatrix(UXB.rowdim());
	      Bo.assign(0.0);
	      Bo.accumulate_transform(UXB,B_ij);
	      RefSCVector Co = Vo.clone();
	      sc::exp::lapack_linsolv_symmnondef(Bo, Co, Vo);
	      RefSCVector Cij = UXB.t() * Co;
	      Cij.scale(-1.0);
	      for(int kl=0; kl<nxc; kl++)
		C_[spin].set_element(kl,ij,Cij.get_element(kl));
	    }
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

  EvalStats eigenvalue_stats(const RefDiagSCMatrix& evals, double threshold)
  {
    const unsigned int n = evals.dim().n();
    unsigned int nneg = 0;
    double eval_max = -1.0e100;
    double eval_min =  1.0e100;
    for(unsigned int i=0; i<n; i++) {
      const double eval = evals.get_element(i);
      eval_max = std::max(eval_max,eval);
      eval_min = std::min(eval_min,eval);
      nneg += ( eval < threshold ? 1 : 0);
    }
    return EvalStats(evals,threshold,nneg,eval_min,eval_max);
  }

  void print_eigenstats(const EvalStats& stats,
			const char* label,
			const Ref<MP2R12EnergyUtil_base>& util,
			int debug,
			bool warn,
			std::ostream& os)
  {
    bool below_thresh = (stats.num_below_threshold > 0);
    bool zero_thresh = (stats.threshold == 0.0);
    
    if (debug >= DefaultPrintThresholds::mostN) {
      ostringstream oss;  oss << "Eigenvalues of " << label;
      util->print(oss.str().c_str(),stats.evals);
      os << indent << ( below_thresh && warn ? "WARNING: " : "")
	 << "min = " << stats.min
	 << ",  max = " << stats.max;
      if (below_thresh) {
	os << ",  " << stats.num_below_threshold;
	if (zero_thresh)
	  os << " negative";
	else
	  os << " below " << stats.threshold;
      }
      os << endl;
    }
    else {
      if (below_thresh && warn) {
	os << indent << "WARNING: " << label << " has " << stats.num_below_threshold;
	if (zero_thresh)
	  os << " negative eigenvalues";
	else
	  os << " eigenvalues below " << stats.threshold;
	os << endl;
      }
    }
  }
  
}
