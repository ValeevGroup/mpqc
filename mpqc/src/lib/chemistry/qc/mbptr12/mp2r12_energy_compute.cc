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
#include <chemistry/qc/mbptr12/utils.impl.h>
#include <chemistry/qc/mbptr12/mp2r12_energy_util.h>
#include <chemistry/qc/mbptr12/print.h>
#include <math/scmat/local.h>

using namespace std;
using namespace sc;

#define USE_INVERT 0
#define USE_INVERT_B 0

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
MP2R12Energy_SpinOrbital::compute()
{
  if (evaluated_)
    return;
  
  Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  Ref<MessageGrp> msg = r12info->msg();
  int me = msg->me();
  int ntasks = msg->n();
  
  const bool obs_eq_vbs = r12info->basis_vir()->equiv(r12info->basis());
  const bool obs_eq_ribs = r12info->basis_ri()->equiv(r12info->basis());
  const bool cabs_empty = obs_eq_vbs && obs_eq_ribs;
  bool ebc = r12eval()->ebc();
  if (r12eval()->vir(Alpha)->rank()==0 ||
      r12eval()->vir(Beta)->rank()==0 ||
      cabs_empty)
    ebc = true;
  /// KS approach to EBC-free method differs from mine if std approx != B || C
  const bool ks_ebcfree = r12eval()->ks_ebcfree() &&
                          (stdapprox_ != LinearR12::StdApprox_B &&
                           stdapprox_ != LinearR12::StdApprox_C);
  const bool diag=r12info->r12tech()->ansatz()->diag();
  if (diag)
    throw ProgrammingError("MP2R12Energy_SpinOrbital::comute() -- diag=true is not supported. Use MP2R12Energy_SpinOrbital_new.");

  // WARNING only RHF and UHF are considered
  const int num_unique_spincases2 = (r12eval()->spin_polarized() ? 3 : 2);
  // 1) the B matrix is the same if the ansatz is diagonal
  // 2) in approximations A', A'', B, and C the B matrix is pair-specific
  const bool same_B_for_all_pairs = diag;
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
    if (r12eval()->dim_oo(spincase2).n() == 0)
      continue;
    
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
      const Ref<MOIndexSpace>& xspace1  = r12eval()->xspace(case1(spincase2));
      const Ref<MOIndexSpace>& occ2_act = r12eval()->occ_act(case2(spincase2));
      const Ref<MOIndexSpace>& vir2_act = r12eval()->vir_act(case2(spincase2));
      const Ref<MOIndexSpace>& xspace2  = r12eval()->xspace(case2(spincase2));
      int nocc1_act = occ1_act->rank();
      int nvir1_act = vir1_act->rank();
      int nx1 = xspace1->rank();
      int nocc2_act = occ2_act->rank();
      int nvir2_act = vir2_act->rank();
      int nx2 = xspace2->rank();
      const std::vector<double> evals_act_occ1 = convert(occ1_act->evals());
      const std::vector<double> evals_act_vir1 = convert(vir1_act->evals());
      const std::vector<double> evals_xspace1  = convert(xspace1->evals());
      const std::vector<double> evals_act_occ2 = convert(occ2_act->evals());
      const std::vector<double> evals_act_vir2 = convert(vir2_act->evals());
      const std::vector<double> evals_xspace2  = convert(xspace2->evals());
      
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
      const int nxy = nxc / num_f12;
      RefSCDimension dim_xy = new SCDimension(nxy);

      // util class treats abstractly on dense and "diagonal" matrices used in orb-invariant and non-orb-invariant ansatze
      Ref<MP2R12EnergyUtil_base> util;
      util = generate_MP2R12EnergyUtil(spincase2,dim_oo, dim_xy, dim_xc, nocc1_act, diag);
      //if (diag_)
      //  util = new MP2R12EnergyUtil_Diag(dim_oo, dim_oo, dim_xc,nocc1_act);
      //else
      //  util = new MP2R12EnergyUtil_Nondiag(dim_oo, dim_xy, dim_xc,nocc1_act);
      
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
	  OverlapOrthog orthog(OverlapOrthog::Canonical,
			       X,
			       X.kit(),
			       r12info->refinfo()->ref()->lindep_tol());
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
        SpinMOPairIter xy_iter(xspace1, xspace2, spincase2);
        for(xy_iter.start(); int(xy_iter); xy_iter.next()) {
          const int xy = xy_iter.ij();
          const int x = xy_iter.i();
          const int y = xy_iter.j();
          
          for(int f=0; f<num_f12; f++) {
            const int f_off = f*nxy;
            const int xyf = f_off + xy;
            
            for(int g=0; g<=f; g++) {
              const int g_off = g*nxy;
              const int xyg = g_off + xy;
              
              const double fx = - (evals_xspace1[x] + evals_xspace2[y]) * X.get_element(xyf,xyg);
              B_ij.accumulate_element(xyf,xyg,fx);
                  
              // If EBC is not assumed add 2.0*Akl,cd*Acd,ow/(ec+ed-ex-ey)
              if (ebc == false) {
                double fy = 0.0;
                SpinMOPairIter cd_iter(vir1_act, vir2_act, spincase2);
                if (ks_ebcfree) {
                  for(cd_iter.start(); cd_iter; cd_iter.next()) {
                    const int cd = cd_iter.ij();
                    const int c = cd_iter.i();
                    const int d = cd_iter.j();
                    
                    fy -= 0.5*( A.get_element(xyf,cd)*Ac.get_element(xyg,cd) + Ac.get_element(xyf,cd)*A.get_element(xyg,cd) )/
                    (evals_act_vir1[c] + evals_act_vir2[d]
                    - evals_xspace1[x] - evals_xspace2[y]);
                  }
                }
                else {
                  for(cd_iter.start(); cd_iter; cd_iter.next()) {
                    const int cd = cd_iter.ij();
                    const int c = cd_iter.i();
                    const int d = cd_iter.j();
                    
                    fy -= A.get_element(xyf,cd)*A.get_element(xyg,cd)/(evals_act_vir1[c] + evals_act_vir2[d]
                    - evals_act_occ1[x] - evals_act_occ2[y]);
                  }
                }
                
                B_ij.accumulate_element(xyf,xyg,fy);
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
          // app A' or B:  form B(ij)kl,ow = Bkl,ow + 1/2(ek + el + eo + ew - 2ei - 2ej)Xkl,ow
          // app A'' or C: form B(ij)kl,ow = Bkl,ow - (ei + ej)Xkl,ow
          B_ij.assign(B);
          
          for(int f=0; f<num_f12; f++) {
            const int f_off = f*nxy;
            
            SpinMOPairIter kl_iter(xspace1, xspace2, spincase2);
            for(kl_iter.start(); kl_iter; kl_iter.next()) {
              const int kl = kl_iter.ij() + f_off;
              const int k = kl_iter.i();
              const int l = kl_iter.j();
              
              for(int g=0; g<=f; g++) {
                const int g_off = g*nxy;
                
                SpinMOPairIter ow_iter(xspace1, xspace2, spincase2);
                for(ow_iter.start(); ow_iter; ow_iter.next()) {
                  const int ow = ow_iter.ij() + g_off;
                  const int o = ow_iter.i();
                  const int w = ow_iter.j();
                  
                  // This is a symmetric matrix
                  if (ow > kl)
                    continue;

                  const double fx = - (evals_act_occ1[i] + evals_act_occ2[j])
                                      * X.get_element(kl, ow);
                  
                  B_ij.accumulate_element(kl,ow,fx);
                  
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
            util->print(prepend_spincase(spincase2,B_ij_label).c_str(),B_ij,ExEnv::outn());
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
  if (!r12info->refinfo()->spin_polarized()) {
    emp2f12_[BetaBeta] = emp2f12_[AlphaAlpha];
    ef12_[BetaBeta] = ef12_[AlphaAlpha];
  }

  evaluated_ = true;
  
  return;
}

RefSymmSCMatrix MP2R12Energy_SpinOrbital_new::compute_B_non_pairspecific(const RefSymmSCMatrix &B,
                                                                         const RefSymmSCMatrix &X,
                                                                         const RefSCMatrix &V,
                                                                         const RefSCMatrix &A,
                                                                         const SpinCase2 &spincase2) {
  const Ref<MOIndexSpace> &xspace1  = r12eval()->xspace(case1(spincase2));
  const Ref<MOIndexSpace> &xspace2  = r12eval()->xspace(case2(spincase2));
  const Ref<MOIndexSpace> &vir1_act = r12eval()->vir_act(case1(spincase2));
  const Ref<MOIndexSpace> &vir2_act = r12eval()->vir_act(case2(spincase2));
  const Ref<MOIndexSpace> &occ1_act = r12eval()->occ_act(case1(spincase2));
  const Ref<MOIndexSpace> &occ2_act = r12eval()->occ_act(case2(spincase2));
  const std::vector<double> evals_xspace1  = convert(xspace1->evals());
  const std::vector<double> evals_xspace2  = convert(xspace2->evals());
  const std::vector<double> evals_act_occ1 = convert(occ1_act->evals());
  const std::vector<double> evals_act_vir1 = convert(vir1_act->evals());
  const std::vector<double> evals_act_occ2 = convert(occ2_act->evals());
  const std::vector<double> evals_act_vir2 = convert(vir2_act->evals());
  const int nocc1_act=occ1_act->rank();
  int nx2 = xspace2->rank();
  const RefSCDimension dim_oo = V.coldim();
  const RefSCDimension dim_xc = V.rowdim();
  const int noo = dim_oo.n();
  const int nxc = dim_xc.n();
  const int num_f12 = r12eval_->r12info()->corrfactor()->nfunctions();
  const int nxy = nxc / num_f12;
  RefSCDimension dim_xy = new SCDimension(nxy);
  
  const bool obs_eq_vbs = r12eval_->r12info()->basis_vir()->equiv(r12eval_->r12info()->basis());
  const bool obs_eq_ribs = r12eval_->r12info()->basis_ri()->equiv(r12eval_->r12info()->basis());
  const bool cabs_empty = obs_eq_vbs && obs_eq_ribs;
  bool ebc = r12eval()->ebc();
  if (r12eval()->vir(Alpha)->rank()==0 ||
      r12eval()->vir(Beta)->rank()==0 ||
      cabs_empty)
    ebc = true;
  SpinMOPairIter xy_iter(xspace1, xspace2, spincase2);
  
  Ref<MP2R12EnergyUtil_Diag> util = generate_MP2R12EnergyUtil_Diag(spincase2,dim_oo,dim_xy,dim_xc,nocc1_act);

  RefSymmSCMatrix B_ij = B.clone();
  B_ij.assign(B);
  for(xy_iter.start(); int(xy_iter); xy_iter.next()) {
    const int xy = xy_iter.ij();
    const int x = xy_iter.i();
    const int y = xy_iter.j();
        
    for(int f=0; f<num_f12; f++) {
      const int f_off = f*nxy;
      const int xyf = f_off + xy;
          
      for(int g=0; g<=f; g++) {
        const int g_off = g*nxy;
        const int xyg = g_off + xy;

        const double fx = - (evals_xspace1[x] + evals_xspace2[y]) * X.get_element(xyf,xyg);
        B_ij.accumulate_element(xyf,xyg,fx);
        if((spincase2==AlphaBeta) && (x!=y)){
      const int yx=y*nx2+x;
      const int yxg = g_off + yx;
      if (xyf > yxg) {
            const double fx = - (evals_xspace1[x] + evals_xspace2[y]) * X.get_element(xyf,yxg);
            B_ij.accumulate_element(xyf,yxg,fx);
          }
        }

        // If EBC is not assumed add 2.0*Akl,cd*Acd,ow/(ec+ed-ex-ey)
        if (ebc == false) {
          double fy = 0.0;
          SpinMOPairIter cd_iter(vir1_act, vir2_act, spincase2);
          for(cd_iter.start(); cd_iter; cd_iter.next()) {
            const int cd = cd_iter.ij();
            const int c = cd_iter.i();
            const int d = cd_iter.j();
              
            fy -= A.get_element(xyf,cd)*A.get_element(xyg,cd)/(evals_act_vir1[c] + evals_act_vir2[d]
            - evals_act_occ1[x] - evals_act_occ2[y]);
          }  // loop over cd_iter

          B_ij.accumulate_element(xyf,xyg,fy);
        }  // ebc == false

      }  // loop over geminal index g
    }  // loop over geminal index f
  }  // loop over xy_iter
  //B_ij.print("B_ij:");
  
  return(B_ij);
}

RefSymmSCMatrix MP2R12Energy_SpinOrbital_new::compute_B_pairspecific(const SpinMOPairIter &ij_iter,
                                                                     const RefSymmSCMatrix &B,
                                                                     const RefSymmSCMatrix &X,
                                                                     const RefSCMatrix &V,
                                                                     const RefSCMatrix &A,
                                                                     const SpinCase2 &spincase2) {
  const Ref<MOIndexSpace> &xspace1  = r12eval()->xspace(case1(spincase2));
  const Ref<MOIndexSpace> &xspace2  = r12eval()->xspace(case2(spincase2));
  const Ref<MOIndexSpace> &vir1_act = r12eval()->vir_act(case1(spincase2));
  const Ref<MOIndexSpace> &vir2_act = r12eval()->vir_act(case2(spincase2));
  const Ref<MOIndexSpace> &occ1_act = r12eval()->occ_act(case1(spincase2));
  const Ref<MOIndexSpace> &occ2_act = r12eval()->occ_act(case2(spincase2));
  const std::vector<double> evals_xspace1  = convert(xspace1->evals());
  const std::vector<double> evals_xspace2  = convert(xspace2->evals());
  const std::vector<double> evals_act_occ1 = convert(occ1_act->evals());
  const std::vector<double> evals_act_vir1 = convert(vir1_act->evals());
  const std::vector<double> evals_act_occ2 = convert(occ2_act->evals());
  const std::vector<double> evals_act_vir2 = convert(vir2_act->evals());
  const int nocc1_act=occ1_act->rank();
  int nx2 = xspace2->rank();
  const RefSCDimension dim_oo = V.coldim();
  const RefSCDimension dim_xc = V.rowdim();
  const int noo = dim_oo.n();
  const int nxc = dim_xc.n();
  const int num_f12 = r12eval_->r12info()->corrfactor()->nfunctions();
  const int nxy = nxc / num_f12;
  RefSCDimension dim_xy = new SCDimension(nxy);
  
  const bool obs_eq_vbs = r12eval_->r12info()->basis_vir()->equiv(r12eval_->r12info()->basis());
  const bool obs_eq_ribs = r12eval_->r12info()->basis_ri()->equiv(r12eval_->r12info()->basis());
  const bool cabs_empty = obs_eq_vbs && obs_eq_ribs;
  bool ebc = r12eval()->ebc();
  if (r12eval()->vir(Alpha)->rank()==0 ||
      r12eval()->vir(Beta)->rank()==0 ||
      cabs_empty)
    ebc = true;
  SpinMOPairIter xy_iter(xspace1, xspace2, spincase2);
  
  RefSymmSCMatrix B_ij = B.clone();
  
  const int i = ij_iter.i();
  const int j = ij_iter.j();


  // In approximations A', B, or C matrices B are pair-specific:
  // app A' or B:  form B(ij)kl,ow = Bkl,ow + 1/2(ek + el + eo + ew - 2ei - 2ej)Xkl,ow
  // app A'' or C: form B(ij)kl,ow = Bkl,ow - (ei + ej)Xkl,ow
  B_ij.assign(B);

  for(int f=0; f<num_f12; f++) {
    const int f_off = f*nxy;

    SpinMOPairIter kl_iter(xspace1, xspace2, spincase2);
    for(kl_iter.start(); kl_iter; kl_iter.next()) {
      const int kl = kl_iter.ij() + f_off;
      const int k = kl_iter.i();
      const int l = kl_iter.j();
  
      for(int g=0; g<=f; g++) {
        const int g_off = g*nxy;
    
        SpinMOPairIter ow_iter(xspace1, xspace2, spincase2);
        for(ow_iter.start(); ow_iter; ow_iter.next()) {
          const int ow = ow_iter.ij() + g_off;
          const int o = ow_iter.i();
          const int w = ow_iter.j();
                
          // This is a symmetric matrix
          if (ow > kl)
            continue;

          const double fx = - (evals_act_occ1[i] + evals_act_occ2[j])
                              * X.get_element(kl, ow);
              
          B_ij.accumulate_element(kl,ow,fx);
              
          // If EBC is not assumed add 2.0*Akl,cd*Acd,ow/(ec+ed-ei-ej)
          if (ebc == false) {
            double fy = 0.0;
            SpinMOPairIter cd_iter(vir1_act, vir2_act, spincase2);
            for(cd_iter.start(); cd_iter; cd_iter.next()) {
              const int cd = cd_iter.ij();
              const int c = cd_iter.i();
              const int d = cd_iter.j();
                
              fy -= A.get_element(kl,cd)*A.get_element(ow,cd)/(evals_act_vir1[c] + evals_act_vir2[d]
                                                               - evals_act_occ1[i] - evals_act_occ2[j]);
            }
            
            B_ij.accumulate_element(kl,ow,fy);
          }  // ebc == false

        }  // loop over ow_iter
      }  // loop over geminal index g
    }  // loop over kl_iter
  }  // loop over geminal index f
  
  return(B_ij);
}

void MP2R12Energy_SpinOrbital_new::compute_MP2R12_old(SpinCase2 &spincase2) {
  Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  const Ref<MOIndexSpace> &occ1_act = r12eval()->occ_act(case1(spincase2));
  const Ref<MOIndexSpace> &vir1_act = r12eval()->vir_act(case1(spincase2));
  const Ref<MOIndexSpace> &xspace1  = r12eval()->xspace(case1(spincase2));
  const Ref<MOIndexSpace> &occ2_act = r12eval()->occ_act(case2(spincase2));
  const Ref<MOIndexSpace> &vir2_act = r12eval()->vir_act(case2(spincase2));
  const Ref<MOIndexSpace> &xspace2  = r12eval()->xspace(case2(spincase2));
  int nocc1_act = occ1_act->rank();
  int nvir1_act = vir1_act->rank();
  int nx1 = xspace1->rank();
  int nocc2_act = occ2_act->rank();
  int nvir2_act = vir2_act->rank();
  int nx2 = xspace2->rank();
  const std::vector<double> evals_act_occ1 = convert(occ1_act->evals());
  const std::vector<double> evals_act_vir1 = convert(vir1_act->evals());
  const std::vector<double> evals_xspace1  = convert(xspace1->evals());
  const std::vector<double> evals_act_occ2 = convert(occ2_act->evals());
  const std::vector<double> evals_act_vir2 = convert(vir2_act->evals());
  const std::vector<double> evals_xspace2  = convert(xspace2->evals());

  const bool obs_eq_vbs = r12info->basis_vir()->equiv(r12info->basis());
  const bool obs_eq_ribs = r12info->basis_ri()->equiv(r12info->basis());
  const bool cabs_empty = obs_eq_vbs && obs_eq_ribs;
  bool ebc = r12eval()->ebc();
  if (r12eval()->vir(Alpha)->rank()==0 ||
      r12eval()->vir(Beta)->rank()==0 ||
      cabs_empty)
    ebc = true;

  // WARNING only RHF and UHF are considered
  const int num_unique_spincases2 = (r12eval()->spin_polarized() ? 3 : 2);
  ExEnv::out0() << "num_unique_spincases2 = " << num_unique_spincases2 << endl;
  // 1) the B matrix is the same if the ansatz is diagonal
  // 2) in approximations A', A'', B, and C the B matrix is pair-specific
  const bool same_B_for_all_pairs = diag_;

  // Get the intermediates V and X
  RefSCMatrix V = r12eval()->V(spincase2);
  RefSymmSCMatrix X = r12eval()->X(spincase2);
  RefSymmSCMatrix B;
  if(stdapprox_==LinearR12::StdApprox_B){
    B = r12eval()->BB(spincase2);
  }
  else if(stdapprox_==LinearR12::StdApprox_C) {
    B = r12eval()->BC(spincase2);
  }
  else {
    B = r12eval()->B(spincase2);
  }

  RefSCMatrix A, Ac;
  if (ebc == false) {
    A = r12eval()->A(spincase2);
  }

  // Prepare total and R12 pairs
  const RefSCDimension dim_oo = V.coldim();
  const RefSCDimension dim_xc = V.rowdim();
  const int noo = dim_oo.n();
  const int nxc = dim_xc.n();
  const int num_f12 = r12eval_->r12info()->corrfactor()->nfunctions();
  const int nxy = nxc / num_f12;
  RefSCDimension dim_xy = new SCDimension(nxy);

  // Prepare the B matrix:
  ExEnv::out0() << "ebc = " << ebc << endl;
  ExEnv::out0() << "diag_ = " << diag_ << endl;
  if(diag_) {  // and same_B_for_all_pairs
    Ref<MP2R12EnergyUtil_Diag> util = generate_MP2R12EnergyUtil_Diag(spincase2,dim_oo,dim_xy,dim_xc,nocc1_act);

    RefSymmSCMatrix B_ij = B.clone();
    B_ij.assign(B);
    SpinMOPairIter xy_iter(xspace1, xspace2, spincase2);
    for(xy_iter.start(); int(xy_iter); xy_iter.next()) {
      const int xy = xy_iter.ij();
      const int x = xy_iter.i();
      const int y = xy_iter.j();
          
      for(int f=0; f<num_f12; f++) {
        const int f_off = f*nxy;
        const int xyf = f_off + xy;
            
        for(int g=0; g<=f; g++) {
          const int g_off = g*nxy;
          const int xyg = g_off + xy;

          const double fx = - (evals_xspace1[x] + evals_xspace2[y]) * X.get_element(xyf,xyg);
          B_ij.accumulate_element(xyf,xyg,fx);
          if((spincase2==AlphaBeta) && (x!=y)){
        const int yx=y*nx2+x;
        const int yxg = g_off + yx;
        if (xyf > yxg) {
              const double fx = - (evals_xspace1[x] + evals_xspace2[y]) * X.get_element(xyf,yxg);
          B_ij.accumulate_element(xyf,yxg,fx);
            }
          }

          // If EBC is not assumed add 2.0*Akl,cd*Acd,ow/(ec+ed-ex-ey)
          if (ebc == false) {
            double fy = 0.0;
            SpinMOPairIter cd_iter(vir1_act, vir2_act, spincase2);
            for(cd_iter.start(); cd_iter; cd_iter.next()) {
              const int cd = cd_iter.ij();
              const int c = cd_iter.i();
              const int d = cd_iter.j();
                
              fy -= A.get_element(xyf,cd)*A.get_element(xyg,cd)/(evals_act_vir1[c] + evals_act_vir2[d]
              - evals_act_occ1[x] - evals_act_occ2[y]);
            }  // loop over cd_iter

            B_ij.accumulate_element(xyf,xyg,fy);
          }  // ebc == false

        }  // loop over geminal index g
      }  // loop over geminal index f
    }  // loop over xy_iter
    //B_ij.print("B_ij:");
    //V.print("V");

#if USE_INVERT_B
    util->invert(B_ij);
    util->times(B_ij,V,C_[spincase2]);
    C_[spincase2].scale(-1.0);
#else   /* not USE_INVERT_B */
    util->solve_linear_system(B_ij, C_[spincase2], V);
    C_[spincase2].scale(-1.0);
#endif  /* not USE_INVERT_B */
    ef12_[spincase2]->assign(util->dot_product(C_[spincase2],V));
  }  // diag_ and same_B_for_all_pairs
  else {  // not diag_ and not same_B_for_all_pairs
    RefSymmSCMatrix B_ij = B.clone();
    SpinMOPairIter ij_iter(occ1_act, occ2_act, spincase2);
    for(ij_iter.start(); int(ij_iter); ij_iter.next()) {
      const int ij = ij_iter.ij();
      const int i = ij_iter.i();
      const int j = ij_iter.j();

      // In approximations A', B, or C matrices B are pair-specific:
      // app A' or B:  form B(ij)kl,ow = Bkl,ow + 1/2(ek + el + eo + ew - 2ei - 2ej)Xkl,ow
      // app A'' or C: form B(ij)kl,ow = Bkl,ow - (ei + ej)Xkl,ow
      B_ij.assign(B);

      for(int f=0; f<num_f12; f++) {
        const int f_off = f*nxy;

    SpinMOPairIter kl_iter(xspace1, xspace2, spincase2);
    for(kl_iter.start(); kl_iter; kl_iter.next()) {
          const int kl = kl_iter.ij() + f_off;
      const int k = kl_iter.i();
      const int l = kl_iter.j();

      for(int g=0; g<=f; g++) {
            const int g_off = g*nxy;

        SpinMOPairIter ow_iter(xspace1, xspace2, spincase2);
        for(ow_iter.start(); ow_iter; ow_iter.next()) {
              const int ow = ow_iter.ij() + g_off;
          const int o = ow_iter.i();
          const int w = ow_iter.j();
                  
              // This is a symmetric matrix
              if (ow > kl)
                continue;

              const double fx = - (evals_act_occ1[i] + evals_act_occ2[j])
                                  * X.get_element(kl, ow);
                  
              B_ij.accumulate_element(kl,ow,fx);
                  
              // If EBC is not assumed add 2.0*Akl,cd*Acd,ow/(ec+ed-ei-ej)
              if (ebc == false) {
                double fy = 0.0;
                SpinMOPairIter cd_iter(vir1_act, vir2_act, spincase2);
                for(cd_iter.start(); cd_iter; cd_iter.next()) {
                  const int cd = cd_iter.ij();
                  const int c = cd_iter.i();
                  const int d = cd_iter.j();
                    
                  fy -= A.get_element(kl,cd)*A.get_element(ow,cd)/(evals_act_vir1[c] + evals_act_vir2[d]
                                                                   - evals_act_occ1[i] - evals_act_occ2[j]);
                }
                
                B_ij.accumulate_element(kl,ow,fy);
              }  // ebc == false

            }  // loop over ow_iter
          }  // loop over geminal index g
        }  // loop over kl_iter
      }  // loop over geminal index f


      // computation of ef12
      RefSCVector V_ij = V.get_column(ij);

      //B_ij.print("B_ij:");
      //V_ij.print("V_ij:");
#if USE_INVERT_B
      // The r12 amplitudes B^-1 * V
      B_ij->gen_invert_this();
      RefSCVector Cij = -1.0*(B_ij * V_ij);
      for(int kl=0; kl<nxc; kl++) {
        C_[spincase2].set_element(kl,ij,Cij.get_element(kl));
      }
#else  /* not USE_INVERT_B */
      // solve B * C = V
      RefSCVector Cij = V_ij.clone();
      sc::exp::lapack_linsolv_symmnondef(B_ij, Cij, V_ij);
      Cij.scale(-1.0);
      for(int kl=0; kl<nxc; kl++) {
        C_[spincase2].set_element(kl,ij,Cij.get_element(kl));
      }
#endif  /* not USE_INVERT_B */
    }  // loop over ij_iter
    // compute ef12_
    RefSCMatrix product = C_[spincase2].t()*V;
    for(int ij=0; ij<noo; ij++){
      ef12_[spincase2].set_element(ij,product.get_element(ij,ij));
    }
  }  // not diag_ and not same_B_for_all_pairs
  emp2f12_[spincase2].assign(r12eval()->emp2(spincase2));
  emp2f12_[spincase2]->accumulate(ef12_[spincase2]);
}

void MP2R12Energy_SpinOrbital_new::determine_C_non_pairspecific(const RefSymmSCMatrix &B_ij,
                                                                const RefSCMatrix &V,
                                                                const SpinCase2 &spincase2,
                                                                const Ref<MP2R12EnergyUtil_Diag> &util) {
  util->solve_linear_system(B_ij, C_[spincase2], V);
  C_[spincase2].scale(-1.0);
}

void MP2R12Energy_SpinOrbital_new::determine_C_fixed_non_pairspecific(const SpinCase2 &spincase2) {
  RefSCDimension oodim = r12eval()->dim_oo(spincase2);
  RefSCDimension xydim = r12eval()->dim_xy(spincase2);
  RefSCDimension f12dim = r12eval()->dim_f12(spincase2);
  
  const unsigned int nij = oodim.n();
  const unsigned int nxy = xydim.n();
  const unsigned int nf12 = f12dim.n();
  const unsigned int ngem = nf12/nxy;  // number of geminals
  
  RefSCDimension gdim=new SCDimension(ngem);
  RefSCDimension two_gdim=new SCDimension(2*gdim.n());
  
  const Ref<MOIndexSpace> &occ1_act = r12eval()->occ_act(case1(spincase2));
  const Ref<MOIndexSpace> &occ2_act = r12eval()->occ_act(case2(spincase2));
  SpinMOPairIter ij_iter(occ1_act, occ2_act, spincase2);
  
  int nocc1_act = occ1_act->rank();
  Ref<MP2R12EnergyUtil_Diag> util = generate_MP2R12EnergyUtil_Diag(spincase2,oodim,xydim,f12dim,nocc1_act);
  
  RefSCVector C_pair;
  if(!STG(r12eval_->r12info()->r12tech()->corrfactor()->geminaldescriptor())){
    throw ProgrammingError("MP2R12Energy_SpinOrbital_new::determine_C_fixed_non_pairspecific -- fixed coefficients can be only determined for Slater type geminals.",__FILE__,__LINE__);
  }
  double gamma=single_slater_exponent(r12eval_->r12info()->r12tech()->corrfactor()->geminaldescriptor());
  double Cp_ij_ij=-1.0/(2.0*gamma);  // singlett
  double Cm_ij_ij=-1.0/(4.0*gamma);  // triplett
  
  C_[spincase2].assign(0.0);
  
  for(ij_iter.start(); int(ij_iter); ij_iter.next()) {
    int ij = ij_iter.ij();
    int i = ij_iter.i();
    int j = ij_iter.j();
    
    if((spincase2==AlphaBeta) && (i!=j)){
      C_pair = C_[spincase2].kit()->vector(two_gdim);
      for(int g=0; g<ngem; g++){
        int goffset=2*g;
        C_pair.set_element(goffset  ,0.5*(Cp_ij_ij+Cm_ij_ij));
        C_pair.set_element(goffset+1,0.5*(Cp_ij_ij-Cm_ij_ij));
      }
      util->put(ij,C_[spincase2],C_pair);
    }  // (spincase2==AlphaBeta) && (i!=j)
    else {  // (spincase2!=AlphaBeta) || (i==j)
      C_pair = C_[spincase2].kit()->vector(gdim);
      if(spincase2==AlphaBeta){
        for(int g=0; g<ngem; g++){
          C_pair.set_element(g,Cp_ij_ij);
        }
        util->put(ij,C_[spincase2],C_pair);
      }  // spincase2==AlphaBeta
      else {  // spincase2!=AlphaBeta
        for(int g=0; g<ngem; g++){
          C_pair.set_element(g,Cm_ij_ij);
        }
        util->put(ij,C_[spincase2],C_pair);
      }  // spincase2!=AlphaBeta
    }  // (spincase2!=AlphaBeta) || (i==j)
  }
}

void MP2R12Energy_SpinOrbital_new::determine_ef12_fixedcoeff_hylleraas(const RefSymmSCMatrix &B_ij,
                                                                       const RefSCMatrix &V,
                                                                       const SpinCase2 &spincase2,
                                                                       const Ref<MP2R12EnergyUtil_Diag> &util) {
  RefSCVector ef12_tmp = util->dot_product(C_[spincase2],V);
  RefSCMatrix product=C_[spincase2].clone();
  util->times(B_ij,C_[spincase2],product);
  ef12_[spincase2]->assign(ef12_tmp*2.0+util->dot_product(C_[spincase2],product));
}

void MP2R12Energy_SpinOrbital_new::determine_C_pairspecific(const int ij,
                                                            const RefSymmSCMatrix &B_ij,
                                                            const RefSCMatrix &V,
                                                            const SpinCase2 &spincase2) {
  const RefSCDimension dim_xc = V.rowdim();
  const int nxc=dim_xc.n();
  
  RefSCVector V_ij = V.get_column(ij);

  //B_ij.print("B_ij:");
  //V_ij.print("V_ij:");

  // solve B * C = V
  RefSCVector Cij = V_ij.clone();
  sc::exp::lapack_linsolv_symmnondef(B_ij, Cij, V_ij);
  Cij.scale(-1.0);
  for(int kl=0; kl<nxc; kl++) {
    C_[spincase2].set_element(kl,ij,Cij.get_element(kl));
  }
}

void MP2R12Energy_SpinOrbital_new::compute_MP2R12_nondiag(){
  const int num_unique_spincases2 = (r12eval()->spin_polarized() ? 3 : 2);
  for(int spin2=0; spin2<num_unique_spincases2; spin2++){
    SpinCase2 spincase2=static_cast<SpinCase2>(spin2);
    Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
    const Ref<MOIndexSpace> &occ1_act = r12eval()->occ_act(case1(spincase2));
    const Ref<MOIndexSpace> &occ2_act = r12eval()->occ_act(case2(spincase2));
    int nocc1_act = occ1_act->rank();
    int nocc2_act = occ2_act->rank();
  
    const bool obs_eq_vbs = r12info->basis_vir()->equiv(r12info->basis());
    const bool obs_eq_ribs = r12info->basis_ri()->equiv(r12info->basis());
    const bool cabs_empty = obs_eq_vbs && obs_eq_ribs;
    bool ebc = r12eval()->ebc();
    if (r12eval()->vir(Alpha)->rank()==0 ||
        r12eval()->vir(Beta)->rank()==0 ||
        cabs_empty)
      ebc = true;
    
    // WARNING only RHF and UHF are considered
    const int num_unique_spincases2 = (r12eval()->spin_polarized() ? 3 : 2);
    ExEnv::out0() << "num_unique_spincases2 = " << num_unique_spincases2 << endl;
    // 1) the B matrix is the same if the ansatz is diagonal
    // 2) in approximations A', A'', B, and C the B matrix is pair-specific
    const bool same_B_for_all_pairs = diag_;
  
    // Get the intermediates V and X
    RefSCMatrix V = r12eval()->V(spincase2);
    RefSymmSCMatrix X = r12eval()->X(spincase2);
    RefSymmSCMatrix B;
    if(stdapprox_==LinearR12::StdApprox_B){
      B = r12eval()->BB(spincase2);
    }
    else if(stdapprox_==LinearR12::StdApprox_C) {
      B = r12eval()->BC(spincase2);
    }
    else {
      B = r12eval()->B(spincase2);
    }
  
    RefSCMatrix A, Ac;
    if (ebc == false) {
      A = r12eval()->A(spincase2);
    }
  
    // Prepare total and R12 pairs
    const RefSCDimension dim_oo = V.coldim();
    const RefSCDimension dim_xc = V.rowdim();
    const int noo = dim_oo.n();
    const int nxc = dim_xc.n();
    const int num_f12 = r12eval_->r12info()->corrfactor()->nfunctions();
    const int nxy = nxc / num_f12;
    RefSCDimension dim_xy = new SCDimension(nxy);
    
    SpinMOPairIter ij_iter(occ1_act, occ2_act, spincase2);
    for(ij_iter.start(); int(ij_iter); ij_iter.next()) {
      const int ij = ij_iter.ij();
      const int i = ij_iter.i();
      const int j = ij_iter.j();
  
      RefSymmSCMatrix B_ij = compute_B_pairspecific(ij_iter,B,X,V,A,spincase2);
  
      determine_C_pairspecific(ij,B_ij,V,spincase2);
    }  // loop over ij_iter
    // compute ef12_
    RefSCMatrix product = C_[spincase2].t()*V;
    for(int ij=0; ij<noo; ij++){
      ef12_[spincase2].set_element(ij,product.get_element(ij,ij));
    }
    
    emp2f12_[spincase2].assign(r12eval()->emp2(spincase2));
    emp2f12_[spincase2]->accumulate(ef12_[spincase2]);
  }
}

void MP2R12Energy_SpinOrbital_new::compute_MP2R12_diag_nonfixed() {
  const int num_unique_spincases2 = (r12eval()->spin_polarized() ? 3 : 2);
  for(int spin2=0; spin2<num_unique_spincases2; spin2++){
    SpinCase2 spincase2=static_cast<SpinCase2>(spin2);
    
    Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
    const Ref<MOIndexSpace> &occ1_act = r12eval()->occ_act(case1(spincase2));
    const Ref<MOIndexSpace> &occ2_act = r12eval()->occ_act(case2(spincase2));
    int nocc1_act = occ1_act->rank();
    int nocc2_act = occ2_act->rank();

    const bool obs_eq_vbs = r12info->basis_vir()->equiv(r12info->basis());
    const bool obs_eq_ribs = r12info->basis_ri()->equiv(r12info->basis());
    const bool cabs_empty = obs_eq_vbs && obs_eq_ribs;
    bool ebc = r12eval()->ebc();
    if (r12eval()->vir(Alpha)->rank()==0 ||
        r12eval()->vir(Beta)->rank()==0 ||
        cabs_empty)
      ebc = true;
    
    // WARNING only RHF and UHF are considered
    const int num_unique_spincases2 = (r12eval()->spin_polarized() ? 3 : 2);
    ExEnv::out0() << "num_unique_spincases2 = " << num_unique_spincases2 << endl;
    // 1) the B matrix is the same if the ansatz is diagonal
    // 2) in approximations A', A'', B, and C the B matrix is pair-specific
    const bool same_B_for_all_pairs = diag_;

    // Get the intermediates V and X
    RefSCMatrix V = r12eval()->V(spincase2);
    RefSymmSCMatrix X = r12eval()->X(spincase2);
    RefSymmSCMatrix B;
    if(stdapprox_==LinearR12::StdApprox_B){
      B = r12eval()->BB(spincase2);
    }
    else if(stdapprox_==LinearR12::StdApprox_C) {
      B = r12eval()->BC(spincase2);
    }
    else {
      B = r12eval()->B(spincase2);
    }

    RefSCMatrix A, Ac;
    if (ebc == false) {
      A = r12eval()->A(spincase2);
    }

    // Prepare total and R12 pairs
    const RefSCDimension dim_oo = V.coldim();
    const RefSCDimension dim_xc = V.rowdim();
    const int noo = dim_oo.n();
    const int nxc = dim_xc.n();
    const int num_f12 = r12eval_->r12info()->corrfactor()->nfunctions();
    const int nxy = nxc / num_f12;
    RefSCDimension dim_xy = new SCDimension(nxy);
    
    RefSymmSCMatrix B_ij = compute_B_non_pairspecific(B,X,V,A,spincase2);

    //B_ij.print("B_ij:");
    //V.print("V");
    
    Ref<MP2R12EnergyUtil_Diag> util = generate_MP2R12EnergyUtil_Diag(spincase2,dim_oo,dim_xy,dim_xc,nocc1_act);

    // compute C_
    determine_C_non_pairspecific(B_ij,V,spincase2,util);

    ef12_[spincase2]->assign(util->dot_product(C_[spincase2],V));
    
    emp2f12_[spincase2].assign(r12eval()->emp2(spincase2));
    emp2f12_[spincase2]->accumulate(ef12_[spincase2]);
  }
}

void MP2R12Energy_SpinOrbital_new::compute_MP2R12_diag_fixed_hylleraas() {
  const int num_unique_spincases2 = (r12eval()->spin_polarized() ? 3 : 2);
  for(int spin2=0; spin2<num_unique_spincases2; spin2++){
    SpinCase2 spincase2=static_cast<SpinCase2>(spin2);
    
    Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
    const Ref<MOIndexSpace> &occ1_act = r12eval()->occ_act(case1(spincase2));
    const Ref<MOIndexSpace> &occ2_act = r12eval()->occ_act(case2(spincase2));
    int nocc1_act = occ1_act->rank();
    int nocc2_act = occ2_act->rank();

    const bool obs_eq_vbs = r12info->basis_vir()->equiv(r12info->basis());
    const bool obs_eq_ribs = r12info->basis_ri()->equiv(r12info->basis());
    const bool cabs_empty = obs_eq_vbs && obs_eq_ribs;
    bool ebc = r12eval()->ebc();
    if (r12eval()->vir(Alpha)->rank()==0 ||
        r12eval()->vir(Beta)->rank()==0 ||
        cabs_empty)
      ebc = true;
    
    // WARNING only RHF and UHF are considered
    const int num_unique_spincases2 = (r12eval()->spin_polarized() ? 3 : 2);
    ExEnv::out0() << "num_unique_spincases2 = " << num_unique_spincases2 << endl;
    // 1) the B matrix is the same if the ansatz is diagonal
    // 2) in approximations A', A'', B, and C the B matrix is pair-specific
    const bool same_B_for_all_pairs = diag_;

    // Get the intermediates V and X
    RefSCMatrix V = r12eval()->V(spincase2);
    RefSymmSCMatrix X = r12eval()->X(spincase2);
    RefSymmSCMatrix B;
    if(stdapprox_==LinearR12::StdApprox_B){
      B = r12eval()->BB(spincase2);
    }
    else if(stdapprox_==LinearR12::StdApprox_C) {
      B = r12eval()->BC(spincase2);
    }
    else {
      B = r12eval()->B(spincase2);
    }

    RefSCMatrix A, Ac;
    if (ebc == false) {
      A = r12eval()->A(spincase2);
    }

    // Prepare total and R12 pairs
    const RefSCDimension dim_oo = V.coldim();
    const RefSCDimension dim_xc = V.rowdim();
    const int noo = dim_oo.n();
    const int nxc = dim_xc.n();
    const int num_f12 = r12eval_->r12info()->corrfactor()->nfunctions();
    const int nxy = nxc / num_f12;
    RefSCDimension dim_xy = new SCDimension(nxy);
    
    Ref<MP2R12EnergyUtil_Diag> util = generate_MP2R12EnergyUtil_Diag(spincase2,dim_oo,dim_xy,dim_xc,nocc1_act);
    RefSymmSCMatrix B_ij = compute_B_non_pairspecific(B,X,V,A,spincase2);
    determine_C_fixed_non_pairspecific(spincase2);
    determine_ef12_fixedcoeff_hylleraas(B_ij,V,spincase2,util);
    
    emp2f12_[spincase2].assign(r12eval()->emp2(spincase2));
    emp2f12_[spincase2]->accumulate(ef12_[spincase2]);
  }
}

void MP2R12Energy_SpinOrbital_new::compute_MP2R12_diag_fixed_nonhylleraas() {
  const int num_unique_spincases2 = (r12eval()->spin_polarized() ? 3 : 2);
  for(int spin2=0; spin2<num_unique_spincases2; spin2++){
    SpinCase2 spincase2=static_cast<SpinCase2>(spin2);
    
    Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
    const Ref<MOIndexSpace> &occ1_act = r12eval()->occ_act(case1(spincase2));
    const Ref<MOIndexSpace> &occ2_act = r12eval()->occ_act(case2(spincase2));
    int nocc1_act = occ1_act->rank();
    int nocc2_act = occ2_act->rank();

    const bool obs_eq_vbs = r12info->basis_vir()->equiv(r12info->basis());
    const bool obs_eq_ribs = r12info->basis_ri()->equiv(r12info->basis());
    const bool cabs_empty = obs_eq_vbs && obs_eq_ribs;
    bool ebc = r12eval()->ebc();
    if (r12eval()->vir(Alpha)->rank()==0 ||
        r12eval()->vir(Beta)->rank()==0 ||
        cabs_empty)
      ebc = true;
    
    // WARNING only RHF and UHF are considered
    const int num_unique_spincases2 = (r12eval()->spin_polarized() ? 3 : 2);
    ExEnv::out0() << "num_unique_spincases2 = " << num_unique_spincases2 << endl;
    // 1) the B matrix is the same if the ansatz is diagonal
    // 2) in approximations A', A'', B, and C the B matrix is pair-specific
    const bool same_B_for_all_pairs = diag_;

    // Get the intermediates V and X
    RefSCMatrix V = r12eval()->V(spincase2);
    RefSymmSCMatrix X = r12eval()->X(spincase2);
    //RefSymmSCMatrix B;
    //if(stdapprox_==LinearR12::StdApprox_B){
    //  B = r12eval()->BB(spincase2);
    //}
    //else if(stdapprox_==LinearR12::StdApprox_C) {
    //  B = r12eval()->BC(spincase2);
    //}
    //else {
    // B = r12eval()->B(spincase2);
    //}

    //RefSCMatrix A, Ac;
    //if (ebc == false) {
    //  A = r12eval()->A(spincase2);
    //}

    // Prepare total and R12 pairs
    const RefSCDimension dim_oo = V.coldim();
    const RefSCDimension dim_xc = V.rowdim();
    const int noo = dim_oo.n();
    const int nxc = dim_xc.n();
    const int num_f12 = r12eval_->r12info()->corrfactor()->nfunctions();
    const int nxy = nxc / num_f12;
    RefSCDimension dim_xy = new SCDimension(nxy);
    
    Ref<MP2R12EnergyUtil_Diag> util = generate_MP2R12EnergyUtil_Diag(spincase2,dim_oo,dim_xy,dim_xc,nocc1_act);
    determine_C_fixed_non_pairspecific(spincase2);
    ef12_[spincase2]->assign(util->dot_product(C_[spincase2],V));
    
    emp2f12_[spincase2].assign(r12eval()->emp2(spincase2));
    emp2f12_[spincase2]->accumulate(ef12_[spincase2]);
  }
}

void MP2R12Energy_SpinOrbital_new::compute_MP2R12(SpinCase2 &spincase2) {
  Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  const Ref<MOIndexSpace> &occ1_act = r12eval()->occ_act(case1(spincase2));
  const Ref<MOIndexSpace> &occ2_act = r12eval()->occ_act(case2(spincase2));
  int nocc1_act = occ1_act->rank();
  int nocc2_act = occ2_act->rank();

  const bool obs_eq_vbs = r12info->basis_vir()->equiv(r12info->basis());
  const bool obs_eq_ribs = r12info->basis_ri()->equiv(r12info->basis());
  const bool cabs_empty = obs_eq_vbs && obs_eq_ribs;
  bool ebc = r12eval()->ebc();
  if (r12eval()->vir(Alpha)->rank()==0 ||
      r12eval()->vir(Beta)->rank()==0 ||
      cabs_empty)
    ebc = true;
  
  // WARNING only RHF and UHF are considered
  const int num_unique_spincases2 = (r12eval()->spin_polarized() ? 3 : 2);
  ExEnv::out0() << "num_unique_spincases2 = " << num_unique_spincases2 << endl;
  // 1) the B matrix is the same if the ansatz is diagonal
  // 2) in approximations A', A'', B, and C the B matrix is pair-specific
  const bool same_B_for_all_pairs = diag_;

  // Get the intermediates V and X
  RefSCMatrix V = r12eval()->V(spincase2);
  RefSymmSCMatrix X = r12eval()->X(spincase2);
  RefSymmSCMatrix B;
  if(stdapprox_==LinearR12::StdApprox_B){
    B = r12eval()->BB(spincase2);
  }
  else if(stdapprox_==LinearR12::StdApprox_C) {
    B = r12eval()->BC(spincase2);
  }
  else {
    B = r12eval()->B(spincase2);
  }

  RefSCMatrix A, Ac;
  if (ebc == false) {
    A = r12eval()->A(spincase2);
  }

  // Prepare total and R12 pairs
  const RefSCDimension dim_oo = V.coldim();
  const RefSCDimension dim_xc = V.rowdim();
  const int noo = dim_oo.n();
  const int nxc = dim_xc.n();
  const int num_f12 = r12eval_->r12info()->corrfactor()->nfunctions();
  const int nxy = nxc / num_f12;
  RefSCDimension dim_xy = new SCDimension(nxy);

  // Prepare the B matrix:
  //ExEnv::out0() << "ebc = " << ebc << endl;
  //ExEnv::out0() << "diag_ = " << diag_ << endl;
  if(diag_) {  // and same_B_for_all_pairs
    RefSymmSCMatrix B_ij = compute_B_non_pairspecific(B,X,V,A,spincase2);

    //B_ij.print("B_ij:");
    //V.print("V");
    
    Ref<MP2R12EnergyUtil_Diag> util = generate_MP2R12EnergyUtil_Diag(spincase2,dim_oo,dim_xy,dim_xc,nocc1_act);

    // compute C_
    determine_C_non_pairspecific(B_ij,V,spincase2,util);

    ef12_[spincase2]->assign(util->dot_product(C_[spincase2],V));
    
    if(fixedcoeff()){
      determine_C_fixed_non_pairspecific(spincase2);
      //determine_ef12_fixedcoeff(B_ij,V,spincase2,util);
      ef12_[spincase2]->assign(util->dot_product(C_[spincase2],V));
      
      emp2f12_[spincase2].assign(r12eval()->emp2(spincase2));
      emp2f12_[spincase2]->accumulate(ef12_[spincase2]);
      emp2f12_[spincase2].assign(r12eval()->emp2(spincase2));
      emp2f12_[spincase2]->accumulate(ef12_[spincase2]);
    }
  }  // diag_ and same_B_for_all_pairs
  else {  // not diag_ and not same_B_for_all_pairs
    SpinMOPairIter ij_iter(occ1_act, occ2_act, spincase2);
    for(ij_iter.start(); int(ij_iter); ij_iter.next()) {
      const int ij = ij_iter.ij();
      const int i = ij_iter.i();
      const int j = ij_iter.j();

      RefSymmSCMatrix B_ij = compute_B_pairspecific(ij_iter,B,X,V,A,spincase2);

      determine_C_pairspecific(ij,B_ij,V,spincase2);
    }  // loop over ij_iter
    // compute ef12_
    RefSCMatrix product = C_[spincase2].t()*V;
    for(int ij=0; ij<noo; ij++){
      ef12_[spincase2].set_element(ij,product.get_element(ij,ij));
    }
  }  // not diag_ and not same_B_for_all_pairs
  emp2f12_[spincase2].assign(r12eval()->emp2(spincase2));
  emp2f12_[spincase2]->accumulate(ef12_[spincase2]);
}

void MP2R12Energy_SpinOrbital_new::compute() {
  if(evaluated_==true) {
    return;
  }

  if(diag_){
    if(fixed_coeff_){
      if(hylleraas_) {
        ExEnv::out0() << "Calling MP2R12Energy_SpinOrbital_new::compute_MP2R12_diag_fixed_hylleraas()" << endl;
        compute_MP2R12_diag_fixed_hylleraas();
      }  // hylleraas_
      else {  // not hylleraas_
        ExEnv::out0() << "Calling MP2R12Energy_SpinOrbital_new::compute_MP2R12_diag_fixed_nonhylleraas()" << endl;
        compute_MP2R12_diag_fixed_nonhylleraas();
      }  // not hylleraas_
    }   // fixed_coeff_
    else {  // not fixed_coeff_
      ExEnv::out0() << "Calling MP2R12Energy_SpinOrbital_new::compute_MP2R12_diag_nonfixed();" << endl;
      compute_MP2R12_diag_nonfixed();
    }  // not fixed_coeff_
  }  // diag_
  else {  // not diag_
    ExEnv::out0() << "Calling MP2R12Energy_SpinOrbital_new::compute_MP2R12_nondiag()" << endl;
    compute_MP2R12_nondiag();
  }  // not diag_
                          
  C_[BetaBeta] = C_[AlphaAlpha];
  
  ef12_[BetaBeta] = ef12_[AlphaAlpha];
  
  emp2f12_[BetaBeta] = emp2f12_[AlphaAlpha];
  
  evaluated_=true;
  
  return;
}

void MP2R12Energy_SpinAdapted::compute() {
  if (evaluated_)
    return;
  
  int spin=0;   /// only matrices of spincase AlphaBeta necessary in closed shell case
  RefSCDimension dim_oo = r12eval()->dim_oo(static_cast<SpinCase2>(spin));
  RefSCDimension dim_f12 = r12eval()->dim_f12(static_cast<SpinCase2>(spin));
  Ref<SCMatrixKit> kit = new LocalSCMatrixKit;
  RefSCMatrix C_alphabeta = kit->matrix(dim_f12,dim_oo);
  
  Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  Ref<MessageGrp> msg = r12info->msg();
  int me = msg->me();
  int ntasks = msg->n();
  
  const bool obs_eq_vbs = r12info->basis_vir()->equiv(r12info->basis());
  const bool obs_eq_ribs = r12info->basis_ri()->equiv(r12info->basis());
  const bool cabs_empty = obs_eq_vbs && obs_eq_ribs;
  bool ebc = r12eval()->ebc();
  if (r12eval()->vir(Alpha)->rank()==0 ||
      r12eval()->vir(Beta)->rank()==0 ||
      cabs_empty)
    ebc = true;
  /// KS approach to EBC-free method differs from mine if std approx != B || C
  const bool ks_ebcfree = r12eval()->ks_ebcfree() &&
                          (stdapprox_ != LinearR12::StdApprox_B &&
                           stdapprox_ != LinearR12::StdApprox_C);
  // Diagonal ansatz?
  const bool diag = false;   /// In the intermediate matrices, all elements are needed. Beause of
                             /// this, diag is set to false even if diagonal spinadapted version is
                             /// used.
  // WARNING only RHF and UHF are considered
  const int num_unique_spincases2 = (r12eval()->spin_polarized() ? 3 : 2);
  
  // 1) the B matrix is the same for all pairs in approximation A (EBC assumed) or if
  //    the ansatz is diagonal
  // 2) in approximations A', A'', B, and C the B matrix is pair-specific
  const bool same_B_for_all_pairs = diag;
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
  
  const SpinCase2 spincase2 = static_cast<SpinCase2>(spin);
  
  RefSCDimension dim_xc;
  
  // Prepare total and R12 pairs
  const int noo = dim_oo.n();
  if (noo == 0)
    throw ProgrammingError("MP2R12Energy_SpinAdapted::compute(): noo is equal zero",__FILE__,__LINE__);
  const int nxc = dim_xc.n();
  const int num_f12 = r12info->corrfactor()->nfunctions();
  const int nxy = nxc / num_f12;
  RefSCDimension dim_xy = new SCDimension(nxy);
  const Ref<MOIndexSpace>& occ1_act = r12eval()->occ_act(case1(spincase2));
  unsigned int nocc1_act = occ1_act->rank();
  
  // util class treats abstractly on dense and "diagonal" matrices used in orb-invariant and non-orb-invariant ansatze
  Ref<MP2R12EnergyUtil_base> util;
  util = generate_MP2R12EnergyUtil(spincase2,dim_oo, dim_xy, dim_xc, nocc1_act, diag);
  //if (diag)
  //  util = new MP2R12EnergyUtil_Diag(dim_oo, dim_oo, dim_xc,nocc1_act);
  //else
  //  util = new MP2R12EnergyUtil_Nondiag(dim_oo, dim_xy, dim_xc,nocc1_act);
  
  double* ef12_vec = new double[noo];
  memset(ef12_vec,0,sizeof(double)*noo);
  
  // Get the intermediates V and X
  RefSCMatrix V = r12eval()->V(spincase2);
  RefSymmSCMatrix X = r12eval()->X(spincase2);
  
  Ref<LinearR12::NullCorrelationFactor> nullcorrptr; nullcorrptr << r12info->corrfactor();
  // if no explicit correlation -- just get MP2 energies
  if (nullcorrptr.nonnull()) {
    throw ProgrammingError("MP2R12Energy_SpinAdapted::compute(): spinadapted case for standard MP2 not implemented",__FILE__,__LINE__);
  }
  else {
    
    const Ref<MOIndexSpace>& occ1_act = r12eval()->occ_act(case1(spincase2));
    const Ref<MOIndexSpace>& vir1_act = r12eval()->vir_act(case1(spincase2));
    const Ref<MOIndexSpace>& xspace1  = r12eval()->xspace(case1(spincase2));
    const Ref<MOIndexSpace>& occ2_act = r12eval()->occ_act(case2(spincase2));
    const Ref<MOIndexSpace>& vir2_act = r12eval()->vir_act(case2(spincase2));
    const Ref<MOIndexSpace>& xspace2  = r12eval()->xspace(case2(spincase2));
    int nocc1_act = occ1_act->rank();
    int nvir1_act = vir1_act->rank();
    int nx1 = xspace1->rank();
    int nocc2_act = occ2_act->rank();
    int nvir2_act = vir2_act->rank();
    int nx2 = xspace2->rank();
    const std::vector<double> evals_act_occ1 = convert(occ1_act->evals());
    const std::vector<double> evals_act_vir1 = convert(vir1_act->evals());
    const std::vector<double> evals_xspace1  = convert(xspace1->evals());
    const std::vector<double> evals_act_occ2 = convert(occ2_act->evals());
    const std::vector<double> evals_act_vir2 = convert(vir2_act->evals());
    const std::vector<double> evals_xspace2  = convert(xspace2->evals());
    
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
    dim_oo = V.coldim();
    dim_xc = V.rowdim();
    
    // In Klopper-Samson method, two kinds of A are used: app B (I replace with my A) and app XX (my Ac)
    RefSCMatrix A, Ac;
    if (ebc == false) {
      A = r12eval()->A(spincase2);
      if (ks_ebcfree) {
        Ac = r12eval()->Ac(spincase2);
      }
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
      SpinMOPairIter xy_iter(xspace1, xspace2, spincase2);
      for(xy_iter.start(); int(xy_iter); xy_iter.next()) {
        const int xy = xy_iter.ij();
        const int x = xy_iter.i();
        const int y = xy_iter.j();
        
        for(int f=0; f<num_f12; f++) {
          const int f_off = f*nxy;
          const int xyf = f_off + xy;
          
          for(int g=0; g<=f; g++) {
            const int g_off = g*nxy;
            const int xyg = g_off + xy;
            
            double fx = - (evals_xspace1[x] + evals_xspace2[y]) * X.get_element(xyf,xyg);
            B_ij.accumulate_element(xyf,xyg,fx);
                
            // If EBC is not assumed add 2.0*Akl,cd*Acd,ow/(ec+ed-ex-ey)
            if (ebc == false) {
              double fy = 0.0;
              SpinMOPairIter cd_iter(vir1_act, vir2_act, spincase2);
              if (ks_ebcfree) {
                for(cd_iter.start(); cd_iter; cd_iter.next()) {
                  const int cd = cd_iter.ij();
                  const int c = cd_iter.i();
                  const int d = cd_iter.j();
                  
                  fy -= 0.5*( A.get_element(xyf,cd)*Ac.get_element(xyg,cd) + Ac.get_element(xyf,cd)*A.get_element(xyg,cd) )/
                  (evals_act_vir1[c] + evals_act_vir2[d]
                  - evals_xspace1[x] - evals_xspace2[y]);
                }
              }
              else {
                for(cd_iter.start(); cd_iter; cd_iter.next()) {
                  const int cd = cd_iter.ij();
                  const int c = cd_iter.i();
                  const int d = cd_iter.j();
                  
                  fy -= A.get_element(xyf,cd)*A.get_element(xyg,cd)/(evals_act_vir1[c] + evals_act_vir2[d]
                  - evals_act_occ1[x] - evals_act_occ2[y]);
                }
              }
              
              B_ij.accumulate_element(xyf,xyg,fy);
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
            RefSCMatrix C = C_alphabeta.clone();
            util->times(B_ij,V,C);
            C_alphabeta.assign(C);  C = 0;
            C_alphabeta.scale(-1.0);
    #else
            // solve B * C = V
        if (!need_to_transform_f12dim) {
          RefSCMatrix C = C_alphabeta.clone();
          util->solve_linear_system(B_ij, C, V);
          C_alphabeta.assign(C);  C = 0;
        }
        else {
          RefSCMatrix Vo = UX * V;
          RefSymmSCMatrix Bo = B_ij.kit()->symmmatrix(UX.rowdim());
          Bo.assign(0.0);
          Bo.accumulate_transform(UX,B_ij);
          RefSCMatrix Co = B_ij.kit()->matrix(UX.rowdim(),C_alphabeta.coldim());
          util->solve_linear_system(Bo, Co, Vo);
          RefSCMatrix C = UX.t() * Co;
          C_alphabeta.assign(C);
        }
            C_alphabeta.scale(-1.0);
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
              // app A' or B:  form B(ij)kl,ow = Bkl,ow + 1/2(ek + el + eo + ew - 2ei - 2ej)Xkl,ow
              // app A'' or C: form B(ij)kl,ow = Bkl,ow - (ei + ej)Xkl,ow
              B_ij.assign(B);
              
              for(int f=0; f<num_f12; f++) {
                const int f_off = f*nxy;
                
                SpinMOPairIter kl_iter(xspace1, xspace2, spincase2);
                for(kl_iter.start(); kl_iter; kl_iter.next()) {
                  const int kl = kl_iter.ij() + f_off;
                  const int k = kl_iter.i();
                  const int l = kl_iter.j();
                  
                  for(int g=0; g<=f; g++) {
                    const int g_off = g*nxy;
                    
                    SpinMOPairIter ow_iter(xspace1, xspace2, spincase2);
                    for(ow_iter.start(); ow_iter; ow_iter.next()) {
                      const int ow = ow_iter.ij() + g_off;
                      const int o = ow_iter.i();
                      const int w = ow_iter.j();
                      
                      // This is a symmetric matrix
                      if (ow > kl)
                        continue;

                      const double fx = - (evals_act_occ1[i] + evals_act_occ2[j]) * X.get_element(kl,ow);
                      B_ij.accumulate_element(kl,ow,fx);
                      
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
                util->print(prepend_spincase(spincase2,B_ij_label).c_str(),B_ij,ExEnv::outn());
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
                  C_alphabeta.set_element(kl,ij,Cij.get_element(kl));
    #else
            // solve B * C = V
            if (!need_to_transform_f12dim) {
              RefSCVector Cij = V_ij.clone();
              sc::exp::lapack_linsolv_symmnondef(B_ij, Cij, V_ij);
              Cij.scale(-1.0);
              for(int kl=0; kl<nxc; kl++)
            C_alphabeta.set_element(kl,ij,Cij.get_element(kl));
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
              C_alphabeta.set_element(kl,ij,Cij.get_element(kl));
            }
    #endif 
          }
        }
      }  // pair-specific B block
    
  }
  
  // Compute pure spin energies
  for(int purespin2=0; purespin2<NPureSpinCases2; purespin2++){
    const PureSpinCase2 purespincase2 = static_cast<PureSpinCase2>(purespin2);
    
    const Ref<MOIndexSpace>& occ1_act = r12eval()->occ_act(case1(spincase2));
    const Ref<MOIndexSpace>& vir1_act = r12eval()->vir_act(case1(spincase2));
    const Ref<MOIndexSpace>& xspace1  = r12eval()->xspace(case1(spincase2));
    const Ref<MOIndexSpace>& occ2_act = r12eval()->occ_act(case2(spincase2));
    const Ref<MOIndexSpace>& vir2_act = r12eval()->vir_act(case2(spincase2));
    const Ref<MOIndexSpace>& xspace2  = r12eval()->xspace(case2(spincase2));
    
    RefSCMatrix V_spinadapt;
    if(purespin2==0){
      C_[purespin2]=sc::spinadapt<Singlet>(C_alphabeta,occ1_act,xspace1);
      V_spinadapt=spinadapt<Singlet>(V,occ1_act,xspace1);
    }
    else if (purespin2==1){
      C_[purespin2]=spinadapt<Triplet>(C_alphabeta,occ1_act,xspace1);
      V_spinadapt=spinadapt<Triplet>(V,occ1_act,xspace1);      
    }
    
    delete util;
    util = generate_MP2R12EnergyUtil(spincase2,dim_oo, dim_xy, dim_xc, nocc1_act, true);
    //util = new MP2R12EnergyUtil_Diag(dim_oo, dim_oo, dim_xc,nocc1_act);  //?
    
    // Compute pair energies
    RefSCVector ef12 = util->dot_product(C_[purespin2],V_spinadapt);
    for(unsigned int ij=0; ij<noo; ij++)
      if (ij%ntasks == me)
        ef12_vec[ij] = ef12(ij);
    msg->sum(ef12_vec,noo,0,-1);
    ef12_[purespin2]->assign(ef12_vec);
  }
  
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
