//
// ri_basis.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdexcept>
#include <sstream>

#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/symmint.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbptr12/linearr12.h>
#include <chemistry/qc/mbptr12/svd.h>
#include <chemistry/qc/mbptr12/transform_factory.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>

using namespace sc;
using namespace std;

void
R12IntEvalInfo::construct_ri_basis_(bool safe)
{
  if (bs_aux_->equiv(bs_)) {
    bs_ri_ = bs_;
    if (abs_method_ == LinearR12::ABS_CABS ||
	abs_method_ == LinearR12::ABS_CABSPlus)
      throw std::runtime_error("R12IntEvalInfo::construct_ri_basis_ -- ABS methods CABS and CABS+ can only be used when ABS != OBS");
  }
  else {
    switch(abs_method_) {
      case LinearR12::ABS_ABS:
	construct_ri_basis_ks_(safe);
	break;
      case LinearR12::ABS_ABSPlus:
	construct_ri_basis_ksplus_(safe);
	break;
      case LinearR12::ABS_CABS:
	construct_ri_basis_ev_(safe);
	break;
      case LinearR12::ABS_CABSPlus:
	construct_ri_basis_evplus_(safe);
	break;
      default:
	throw std::runtime_error("R12IntEvalInfo::construct_ri_basis_ -- invalid ABS method");
    }
  }
}

void
R12IntEvalInfo::construct_ri_basis_ks_(bool safe)
{
  bs_ri_ = bs_aux_;
  if (!abs_spans_obs_()) {
    ExEnv::out0() << endl << indent << "WARNING: the auxiliary basis is not safe to use with the given orbital basis" << endl << endl;
    if (safe)
      throw std::runtime_error("R12IntEvalInfo::construct_ri_basis_ks_ -- auxiliary basis is not safe to use with the given orbital basis");
  }
}

void
R12IntEvalInfo::construct_ri_basis_ksplus_(bool safe)
{
  GaussianBasisSet& abs = *(bs_aux_.pointer());
  bs_ri_ = abs + bs_;
  construct_orthog_ri_();
}

void
R12IntEvalInfo::construct_ri_basis_ev_(bool safe)
{
  bs_ri_ = bs_aux_;
  if (!abs_spans_obs_()) {
    ExEnv::out0() << endl << indent << "WARNING: the auxiliary basis is not safe to use with the given orbital basis" << endl << endl;
    if (safe)
      throw std::runtime_error("R12IntEvalInfo::construct_ri_basis_ev_ -- auxiliary basis is not safe to use with the given orbital basis");
  }
  construct_ortho_comp_svd_();
}

void
R12IntEvalInfo::construct_ri_basis_evplus_(bool safe)
{
  GaussianBasisSet& abs = *(bs_aux_.pointer());
  bs_ri_ = abs + bs_;
  construct_ortho_comp_svd_();
}

void
R12IntEvalInfo::construct_orthog_aux_()
{
  if (abs_space_.nonnull())
    return;

  abs_space_ = orthogonalize("ABS", bs_aux_, ref_->orthog_method(), ref_->lindep_tol(), nlindep_aux_);

  if (bs_aux_ == bs_ri_)
    ribs_space_ = abs_space_;
}

void
R12IntEvalInfo::construct_orthog_vir_()
{
  if (vir_space_.nonnull())
    return;

  if (bs_ == bs_vir_) {
    // If virtuals are from the same space as occupieds, then everything is easy
    vir_space_ = new MOIndexSpace("unoccupied MOs sorted by energy", mo_space_->coefs(),
                                  mo_space_->basis(), mo_space_->evals(), ndocc(), 0);
    // If virtuals are from the same space as occupieds, then everything is easy
    vir_space_symblk_ = new MOIndexSpace("unoccupied MOs symmetry-blocked", mo_space_->coefs(),
                                         mo_space_->basis(), mo_space_->evals(), ndocc(), 0, MOIndexSpace::symmetry);

    if (nfzv() == 0)
      act_vir_space_ = vir_space_;
    else
      act_vir_space_ = new MOIndexSpace("active unoccupied MOs sorted by energy", mo_space_->coefs(),
                                  mo_space_->basis(), mo_space_->evals(), ndocc(), nfzv());
    nlindep_vir_ = 0;
  }
  else {
    // This is a set of orthonormal functions that span VBS
    Ref<MOIndexSpace> vir_space = orthogonalize("VBS", bs_vir_, ref_->orthog_method(), ref_->lindep_tol(), nlindep_vir_);
    // Now project out occupied MOs
    vir_space_symblk_ = orthog_comp(occ_space_symblk_, vir_space, "VBS", ref_->lindep_tol());

    // Design flaw!!! Need to compute Fock matrix right here but can't since Fock is built into R12IntEval
    // Need to move all relevant code outside of MBPT2-R12 code
    if (nfzv_ != 0)
      throw std::runtime_error("R12IntEvalInfo::construct_orthog_vir_() -- nfzv_ != 0 is not allowed yet");
    vir_space_ = vir_space_symblk_;
    act_vir_space_ = vir_space_symblk_;
  }
}

void
R12IntEvalInfo::construct_orthog_ri_()
{
  if (bs_ri_.null())
    throw std::runtime_error("R12IntEvalInfo::construct_orthog_ri_ -- RI basis has not been set yet");
  if (bs_aux_ == bs_ri_)
    construct_orthog_aux_();
  if (ribs_space_.nonnull())
    return;

  ribs_space_ = orthogonalize("RI-BS", bs_ri_, ref_->orthog_method(), ref_->lindep_tol(), nlindep_ri_);
}

bool
R12IntEvalInfo::abs_spans_obs_()
{
  construct_orthog_aux_();

  // Compute the bumber of linear dependencies in OBS+ABS
  GaussianBasisSet& abs = *(bs_aux_.pointer());
  Ref<GaussianBasisSet> ri_basis = abs + bs_;
  int nlindep_ri = 0;
  if (bs_ri_.nonnull() && ri_basis->equiv(bs_ri_)) {
    construct_orthog_ri_();
    nlindep_ri = nlindep_ri_;
  }
  else {
    Ref<MOIndexSpace> ribs_space = orthogonalize("OBS+ABS", ri_basis, ref_->orthog_method(), ref_->lindep_tol(), nlindep_ri);
  }

  if (nlindep_ri - nlindep_aux_ - mo_space_->rank() == 0)
    return true;
  else
    return false;
}

/////////////////////////////////////////////////////////////////////////////

void
R12IntEvalInfo::construct_ortho_comp_svd_()
{
   construct_orthog_aux_();
   construct_orthog_vir_();
   construct_orthog_ri_();

   if (debug_ > 1) {
     occ_space_symblk_->coefs().print("Occupied MOs (symblocked)");
     vir_space_symblk_->coefs().print("Virtual MOs (symblocked)");
     obs_space_->coefs().print("All MOs");
     act_occ_space_->coefs().print("Active occupied MOs");
     act_vir_space_->coefs().print("Active virtual MOs");
     ribs_space_->coefs().print("Orthogonal RI-BS");
   }

   ribs_space_ = orthog_comp(occ_space_symblk_, ribs_space_, "RI-BS", ref_->lindep_tol());
   ribs_space_ = orthog_comp(vir_space_symblk_, ribs_space_, "RI-BS", ref_->lindep_tol());
}

Ref<MOIndexSpace>
R12IntEvalInfo::orthogonalize(const std::string& name, const Ref<GaussianBasisSet>& bs,
                              OverlapOrthog::OrthogMethod orthog_method, double lindep_tol,
                              int& nlindep)
{
  // Make an Integral and initialize with bs_aux
  Ref<Integral> integral = Integral::get_default_integral()->clone();
  integral->set_basis(bs);
  Ref<PetiteList> plist = integral->petite_list();
  Ref<OneBodyInt> ov_engine = integral->overlap();

  // form skeleton s matrix
  Ref<SCMatrixKit> matrixkit = bs->matrixkit();
  RefSymmSCMatrix s(bs->basisdim(), matrixkit);
  Ref<SCElementOp> ov =
    new OneBodyIntOp(new SymmOneBodyIntIter(ov_engine, plist));

  s.assign(0.0);
  s.element_op(ov);
  ov=0;
  //if (debug_ > 1) {
  //  std::string s_label = "AO skeleton overlap (" + name + "/" + name + ")";
  //  s.print(s_label.c_str());
  //}

  // then symmetrize it
  RefSCDimension sodim = plist->SO_basisdim();
  Ref<SCMatrixKit> so_matrixkit = bs->so_matrixkit();
  RefSymmSCMatrix overlap(sodim, so_matrixkit);
  plist->symmetrize(s,overlap);

  // and clean up a bit
  ov_engine = 0;
  s = 0;

  //
  // Compute orthogonalizer for bs
  //
  ExEnv::out0() << indent << "Orthogonalizing basis for space " << name << ":" << endl << incindent;
  OverlapOrthog orthog(orthog_method,
                       overlap,
                       so_matrixkit,
                       lindep_tol,
                       0);
  RefSCMatrix orthog_so = orthog.basis_to_orthog_basis();
  orthog_so = orthog_so.t();
  RefSCMatrix orthog_ao = plist->evecs_to_AO_basis(orthog_so);
  orthog_so = 0;
  ExEnv::out0() << decindent;

  nlindep = orthog.nlindep();
  Ref<MOIndexSpace> space = new MOIndexSpace(name,orthog_ao,bs);

  return space;
}


Ref<MOIndexSpace>
R12IntEvalInfo::orthog_comp(const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2,
                            const std::string& name, double lindep_tol)
{
  // Both spaces must be ordered in the same way
  if (space1->moorder() != space2->moorder())
    throw std::runtime_error("R12IntEvalInfo::orthog_comp() -- space1 and space2 are ordered differently ");

  ExEnv::out0() << indent
                << "SVD-projecting out " << space1->name() << " out of " << space2->name()
                << " to obtain space " << name << endl << incindent;

  // C12 = C1 * S12 * C2
  RefSCMatrix C12;
  compute_overlap_ints(space1,space2,C12);

  //
  // SVDecompose C12 = U Sigma V and throw out columns of V
  //
  Ref<SCMatrixKit> ao_matrixkit = space1->basis()->matrixkit();
  Ref<SCMatrixKit> so_matrixkit = space1->basis()->so_matrixkit();
  int nblocks = C12.nblock();
  const double toler = lindep_tol;
  double min_sigma = 1.0;
  double max_sigma = 0.0;
  int* nvec_per_block = new int[nblocks];
  // basis for orthogonal complement is a vector of nvecs by nbasis2
  // we don't know nvecs yet, so use rank2
  RefSCMatrix orthog2 = space2->coefs();
  int rank2 = orthog2.coldim().n();
  int nbasis2 = orthog2.rowdim().n();
  double* vecs = new double[rank2 * nbasis2];
  int nlindep = 0;

  int v_offset = 0;
  for(int b=0; b<nblocks; b++) {

    RefSCDimension rowd = C12.rowdim()->blocks()->subdim(b);
    RefSCDimension cold = C12.coldim()->blocks()->subdim(b);
    int nrow = rowd.n();
    int ncol = cold.n();
    if (nrow && ncol) {

      RefSCMatrix C12_b = C12.block(b);
      RefSCDimension sigd = nrow < ncol ? rowd : cold;
      int nsigmas = sigd.n();

      RefSCMatrix U(rowd, rowd, ao_matrixkit);
      RefSCMatrix V(cold, cold, ao_matrixkit);
      RefDiagSCMatrix Sigma(sigd, ao_matrixkit);

      // C12_b.svd(U,Sigma,V);
      exp::lapack_svd(C12_b,U,Sigma,V);

      // Transform V into AO basis. Vectors are in rows
      RefSCMatrix orthog2_b = orthog2.block(b);
      V = V * orthog2_b.t();

      // Figure out how many sigmas are too small, i.e. how many vectors from space2 overlap
      // only weakly with space1.
      // NOTE: Sigma values returned by svd() are in descending order
      int nzeros = 0;
      for(int s=0; s<nsigmas; s++) {
        double sigma = Sigma(s);
        if (sigma < toler)
          nzeros++;
        if (sigma < min_sigma)
          min_sigma = sigma;
        if (sigma > max_sigma)
          max_sigma = sigma;
      }

      // number of vectors that span the orthogonal space
      nvec_per_block[b] = nzeros + ncol - nsigmas;
      nlindep += nsigmas - nzeros;

      if (nvec_per_block[b]) {
        int v_first = nsigmas - nzeros;
        int v_last = ncol - 1;
        double* v_ptr = vecs + v_offset*nbasis2;
        RefSCMatrix vtmp = V.get_subblock(v_first,v_last,0,nbasis2-1);
        vtmp.convert(v_ptr);
      }
    }
    else {
      nvec_per_block[b] = ncol;

      if (nvec_per_block[b]) {
        RefSCMatrix orthog2_b = orthog2.block(b);
        orthog2_b = orthog2_b.t();
        double* v_ptr = vecs + v_offset*nbasis2;
        orthog2_b.convert(v_ptr);
      }
    }

    v_offset += nvec_per_block[b];
  }

  // Modify error message
  if (v_offset == 0) {
    const std::string errmsg = "R12IntEvalInfo::orthog_comp() -- " + space2->name()
    + " has null projection on orthogonal complement to " + space2->name()
    + "Modify/increase basis for " + space2->name() + ".";
    throw std::runtime_error(errmsg.c_str());
  }

  // convert vecs into orthog2
  // modify for the dimension
  RefSCDimension orthog_dim = new SCDimension(v_offset, nblocks, nvec_per_block, "");
  for(int b=0; b<nblocks; b++)
    orthog_dim->blocks()->set_subdim(b, new SCDimension(nvec_per_block[b]));
  RefSCMatrix orthog_vecs(orthog_dim,orthog2.rowdim(),so_matrixkit);
  orthog_vecs.assign(vecs);
  orthog2 = orthog_vecs.t();

  ExEnv::out0() << indent
    << nlindep << " basis function"
    << (nlindep>1?"s":"")
    << " projected out of " << space2->name() << "."
    << endl;
  ExEnv::out0() << indent
    << "n(basis):        ";
  for (int i=0; i<orthog_dim->blocks()->nblock(); i++) {
    ExEnv::out0() << scprintf(" %5d", orthog_dim->blocks()->size(i));
  }
  ExEnv::out0() << endl;
  ExEnv::out0() << indent
    << "Maximum singular value = "
    << max_sigma << endl
    << indent
    << "Minimum singular value = "
    << min_sigma << endl;
  ExEnv::out0() << decindent;

  delete[] vecs;
  delete[] nvec_per_block;

  Ref<MOIndexSpace> orthog_comp_space = new MOIndexSpace(name,orthog2,space2->basis());
  
  return orthog_comp_space;
}


Ref<MOIndexSpace>
R12IntEvalInfo::gen_project(const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2,
                            const std::string& name, double lindep_tol)
{
  //
  // Projection works as follows:
  // 1) Compute overlap matrix between orthonormal spaces 1 and 2: C12 = C1 * S12 * C2
  // 2) SVDecompose C12 = U Sigma V^t, throw out (near)singular triplets
  // 3) Projected vectors (in AO basis) are X2 = C2 * V * Sigma^{-1} * U^t, where Sigma^{-1} is the generalized inverse
  //


  
  // Both spaces must be ordered in the same way
  if (space1->moorder() != space2->moorder())
    throw std::runtime_error("R12IntEvalInfo::orthog_comp() -- space1 and space2 are ordered differently ");

  ExEnv::out0() << indent
                << "Projecting " << space1->name() << " onto " << space2->name()
                << " exactly to obtain space " << name << endl << incindent;

  // C12 = C1 * S12 * C2
  RefSCMatrix C12;
  compute_overlap_ints(space1,space2,C12);
  C12.print("C12 matrix");

  // Check dimensions of C12 to make sure that projection makes sense
  
  
  Ref<SCMatrixKit> ao_matrixkit = space1->basis()->matrixkit();
  Ref<SCMatrixKit> so_matrixkit = space1->basis()->so_matrixkit();
  int nblocks = C12.nblock();
  const double toler = lindep_tol;
  double min_sigma = 1.0;
  double max_sigma = 0.0;
  int* nvec_per_block = new int[nblocks];

  // projected vectors are a matrix of nvecs by nbasis2
  // we don't know nvecs yet, so use rank1
  RefSCMatrix C1 = space1->coefs();
  RefSCMatrix C2 = space2->coefs();
  int rank1 = space1->coefs()->ncol();
  int nbasis2 = C2->nrow();
  double* vecs = new double[rank1 * nbasis2];
  int nweakovlp = 0;

  int v_offset = 0;
  for(int b=0; b<nblocks; b++) {

    RefSCDimension rowd = C12.rowdim()->blocks()->subdim(b);
    RefSCDimension cold = C12.coldim()->blocks()->subdim(b);
    int nrow = rowd.n();
    int ncol = cold.n();
    
    // Cannot project if rank of the target space is smaller than the rank of the source space
    if (nrow > ncol)
      throw std::runtime_error("R12IntEvalInfo::svd_project() -- rank of the target space is smaller than the rank of the source space");
    
    if (nrow && ncol) {

      RefSCMatrix C12_b = C12.block(b);
      RefSCDimension sigd = rowd;
      int nsigmas = sigd.n();

      RefSCMatrix U(rowd, rowd, ao_matrixkit);
      RefSCMatrix V(cold, cold, ao_matrixkit);
      RefDiagSCMatrix Sigma(sigd, ao_matrixkit);
      
      //
      // Compute C12 = U * Sigma * V
      //
      /* C12_b.svd(U,Sigma,V); */
      exp::lapack_svd(C12_b,U,Sigma,V);

      // Figure out how many sigmas are too small, i.e. how many vectors from space2 overlap
      // only weakly with space1.
      // NOTE: Sigma values returned by svd() are in descending order
      int nzeros = 0;
      for(int s=0; s<nsigmas; s++) {
        double sigma = Sigma(s);
        if (sigma < toler)
          nzeros++;
        if (sigma < min_sigma)
          min_sigma = sigma;
        if (sigma > max_sigma)
          max_sigma = sigma;
      }

      // number of vectors that span the projected space
      nvec_per_block[b] = nsigmas - nzeros;
      if (nvec_per_block[b] < nrow)
        throw std::runtime_error("R12IntEvalInfo::gen_project() -- space 1 is not fully spanned by space 2");
      nweakovlp += nzeros + ncol - nrow;

      if (nvec_per_block[b]) {
        int s_first = 0;
        int s_last = nvec_per_block[b]-1;
        RefSCMatrix vtmp = V.get_subblock(s_first,s_last,0,ncol-1);
        RefSCDimension rowdim = vtmp.rowdim();
        RefDiagSCMatrix stmp = vtmp.kit()->diagmatrix(rowdim);
        for(int i=0; i<nvec_per_block[b]; i++)
          stmp(i) = 1.0/(Sigma(i));
        RefSCMatrix utmp = U.get_subblock(0,nrow-1,s_first,s_last);
        RefSCMatrix C12_inv_t = (utmp * stmp) * vtmp;
        
        (C12_b * C12_inv_t.t()).print("C12 * C12^{-1}");
        (C12_inv_t * C12_b.t()).print("C12^{-1} * C12");
        
        // Transform V into AO basis and transpose so that vectors are in rows
        RefSCMatrix C2_b = C2.block(b);
        RefSCMatrix X2_t = C12_inv_t * C2_b.t();
        double* x2t_ptr = vecs + v_offset*nbasis2;
        X2_t.convert(x2t_ptr);
      }
    }
    else {
      nvec_per_block[b] = 0;
    }

    
    v_offset += nvec_per_block[b];
  }

  // convert vecs into proj
  RefSCMatrix proj(C1.coldim(),C2.rowdim(),so_matrixkit);
  proj.assign(vecs);
  proj = proj.t();

  ExEnv::out0() << indent
    << nweakovlp << " basis function"
    << (nweakovlp>1?"s":"")
    << " in " << space2->name() << " did not overlap significantly with "
    << space1->name() << "." << endl;
  ExEnv::out0() << indent
    << "n(basis):        ";
  for (int i=0; i<proj.coldim()->blocks()->nblock(); i++) {
    ExEnv::out0() << scprintf(" %5d", proj.coldim()->blocks()->size(i));
  }
  ExEnv::out0() << endl;
  ExEnv::out0() << indent
    << "Maximum singular value = "
    << max_sigma << endl
    << indent
    << "Minimum singular value = "
    << min_sigma << endl;
  ExEnv::out0() << decindent;

  delete[] vecs;
  delete[] nvec_per_block;

  Ref<MOIndexSpace> proj_space = new MOIndexSpace(name,proj,space2->basis());

  RefSCMatrix S12;  compute_overlap_ints(space1,proj_space,S12);
  S12.print("Check: overlap between space1 and projected space");
  
  return proj_space;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
