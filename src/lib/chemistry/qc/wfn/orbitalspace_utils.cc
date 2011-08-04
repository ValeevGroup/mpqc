//
// orbitalspace_utils.cc
//
// Copyright (C) 2009 Edward Valeev
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

#include <chemistry/qc/wfn/orbitalspace_utils.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/obintfactory.h>
#include <math/scmat/svd.h>

using namespace std;

namespace sc {

  Ref<OrbitalSpace>
  orthogonalize(const std::string& id, const std::string& name, const Ref<GaussianBasisSet>& bs,
                const Ref<Integral>& ints,
                OverlapOrthog::OrthogMethod orthog_method, double lindep_tol,
                int& nlindep)
  {
    // Make an Integral and initialize with bs_aux
    Ref<Integral> integral = ints->clone();
    integral->set_basis(bs);
    Ref<PetiteList> plist = integral->petite_list();
    RefSymmSCMatrix overlap = compute_onebody_matrix<&Integral::overlap>(plist);

    // if nlindep is non-negative AND orthog_method is symmetric/canonical
    // tell OverlapOrthog to make the number of linear dependencies be exactly nlindep
    if (nlindep >= 0 && (orthog_method == OverlapOrthog::Symmetric || orthog_method == OverlapOrthog::Canonical))
      lindep_tol = -nlindep;

    //
    // Compute orthogonalizer for bs
    //
    ExEnv::out0() << indent << "Orthogonalizing basis for space " << name << ":" << endl << incindent;
    OverlapOrthog orthog(orthog_method,
                         overlap,
                         bs->so_matrixkit(),
                         lindep_tol,
                         0);
    RefSCMatrix orthog_so = orthog.basis_to_orthog_basis();
    orthog_so = orthog_so.t();
    RefSCMatrix orthog_ao = plist->evecs_to_AO_basis(orthog_so);
    orthog_so = 0;
    ExEnv::out0() << decindent;

    nlindep = orthog.nlindep();
    Ref<OrbitalSpace> space = new OrbitalSpace(id,name,orthog_ao,bs,integral);

    return space;
  }


  Ref<OrbitalSpace>
  orthog_comp(const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2,
              const std::string& id, const std::string& name, double lindep_tol)
  {
    if (!space1->integral()->equiv(space2->integral()))
      throw ProgrammingError("Two OrbitalSpaces use incompatible Integral factories");
    // Both spaces must have same blocking
    if (space1->nblocks() != space2->nblocks())
      throw std::runtime_error("sc::orthog_comp() -- space1 and space2 have incompatible blocking");

    ExEnv::out0() << indent
                  << "SVD-projecting out " << space1->name() << " out of " << space2->name()
                  << " to obtain space " << name << endl << incindent;

    // If space1 is void, return a copy of the original space
    if (space1->rank() == 0)
        return new OrbitalSpace(id,name,space2->coefs(),space2->basis(),space2->integral());

    // C12 = C1 * S12 * C2
    RefSCMatrix C12 = compute_overlap_ints(space1,space2);

    //
    // SVDecompose C12 = U Sigma V and throw out columns of V
    //
    Ref<SCMatrixKit> ao_matrixkit = space1->basis()->matrixkit();
    Ref<SCMatrixKit> so_matrixkit = space1->basis()->so_matrixkit();
    int nblocks = C12.nblock();
    const double toler = lindep_tol;
    double min_sigma = 1.0;
    double max_sigma = 0.0;
    vector<int> nvec_per_block(nblocks);
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
        lapack_svd(C12_b,U,Sigma,V);

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
      const std::string errmsg = "R12WavefunctionWorld::orthog_comp() -- " + space2->name()
      + " has null projection on orthogonal complement to " + space2->name()
      + "Modify/increase basis for " + space2->name() + ".";
      throw std::runtime_error(errmsg.c_str());
    }

    // convert vecs into orthog2
    // modify for the dimension
    RefSCDimension orthog_dim = new SCDimension(v_offset, nblocks, &nvec_per_block[0], "");
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

    Ref<OrbitalSpace> orthog_comp_space = new OrbitalSpace(id,name,orthog2,space2->basis(),space2->integral());

    return orthog_comp_space;
  }


  Ref<OrbitalSpace>
  gen_project(const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2,
              const std::string& id, const std::string& name, double lindep_tol)
  {
    //
    // Projection works as follows:
    // 1) Compute overlap matrix between orthonormal spaces 1 and 2: C12 = C1 * S12 * C2
    // 2) SVDecompose C12 = U Sigma V^t, throw out (near)singular triplets
    // 3) Projected vectors (in AO basis) are X2 = C2 * V * Sigma^{-1} * U^t, where Sigma^{-1} is the generalized inverse
    //

    // Check integral factories
    if (!space1->integral()->equiv(space2->integral()))
      throw ProgrammingError("Two OrbitalSpaces use incompatible Integral factories");
    // Both spaces must have same blocking
    if (space1->nblocks() != space2->nblocks())
      throw std::runtime_error("R12WavefunctionWorld::orthog_comp() -- space1 and space2 have incompatible blocking");

    ExEnv::out0() << indent
                  << "Projecting " << space1->name() << " onto " << space2->name()
                  << " exactly to obtain space " << name << endl << incindent;

    // C12 = C1 * S12 * C2
    RefSCMatrix C12 = compute_overlap_ints(space1,space2);

    // Check dimensions of C12 to make sure that projection makes sense


    Ref<SCMatrixKit> ao_matrixkit = space1->basis()->matrixkit();
    Ref<SCMatrixKit> so_matrixkit = space1->basis()->so_matrixkit();
    int nblocks = C12.nblock();
    const double toler = lindep_tol;
    double min_sigma = 1.0;
    double max_sigma = 0.0;
    vector<int> nvec_per_block(nblocks);

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
        throw std::runtime_error("R12WavefunctionWorld::svd_project() -- rank of the target space is smaller than the rank of the source space");

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
        lapack_svd(C12_b,U,Sigma,V);

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
          throw std::runtime_error("R12WavefunctionWorld::gen_project() -- space 1 is not fully spanned by space 2");
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

    Ref<OrbitalSpace> proj_space = new OrbitalSpace(id,name,proj,space2->basis(),space2->integral());

    return proj_space;
  }

  RefSCMatrix
  compute_overlap_ints(const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2) {
    if (!space1->integral()->equiv(space2->integral()))
      throw ProgrammingError("two OrbitalSpaces use incompatible Integral factories");
    const Ref<GaussianBasisSet> bs1 = space1->basis();
    const Ref<GaussianBasisSet> bs2 = space2->basis();
    const bool bs1_eq_bs2 = (bs1 == bs2);
    int nshell1 = bs1->nshell();
    int nshell2 = bs2->nshell();

    RefSCMatrix vec1t = space1->coefs().t();
    RefSCMatrix vec2 = space2->coefs();

    Ref<Integral> localints = space1->integral()->clone();
    localints->set_basis(bs1,bs2);

    Ref<OneBodyInt> ov_ints = localints->overlap();

    // form AO moment matrices
    RefSCDimension aodim1 = vec1t.coldim();
    RefSCDimension aodim2 = vec2.rowdim();
    Ref<SCMatrixKit> aokit = bs1->so_matrixkit();
    RefSCMatrix s(aodim1, aodim2, aokit);
    s.assign(0.0);

    for(int sh1=0; sh1<nshell1; sh1++) {
      int bf1_offset = bs1->shell_to_function(sh1);
      int nbf1 = bs1->shell(sh1).nfunction();

      int sh2max;
      if (bs1_eq_bs2)
        sh2max = sh1;
      else
        sh2max = nshell2-1;

      for(int sh2=0; sh2<=sh2max; sh2++) {
        int bf2_offset = bs2->shell_to_function(sh2);
        int nbf2 = bs2->shell(sh2).nfunction();

        ov_ints->compute_shell(sh1,sh2);
        const double *ovintsptr = ov_ints->buffer();

        int bf1_index = bf1_offset;
        for(int bf1=0; bf1<nbf1; bf1++, bf1_index++, ovintsptr+=nbf2) {
          int bf2_index = bf2_offset;
          const double *ptr = ovintsptr;
          int bf2max;
          if (bs1_eq_bs2 && sh1 == sh2)
            bf2max = bf1;
          else
            bf2max = nbf2-1;
          for(int bf2=0; bf2<=bf2max; bf2++, bf2_index++) {

            s.set_element(bf1_index, bf2_index, *(ptr++));

          }
        }
      }
    }

    // and clean up a bit
    ov_ints = 0;

    // Symmetrize matrices, if necessary
    if (bs1_eq_bs2) {

      const int nbasis = bs1->nbasis();

      for(int bf1=0; bf1<nbasis; bf1++)
        for(int bf2=0; bf2<=bf1; bf2++) {
          s(bf2,bf1) = s(bf1,bf2);
        }

    }


      // finally, transform
      RefSCMatrix S = vec1t * s * vec2;

      // and clean up a bit
      s = 0;

      //if (debug_ > 1) {
      //  S.print("Overlap");
      //}
      return S;
  }

  void
  compute_multipole_ints(const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2,
                         RefSCMatrix& MX, RefSCMatrix& MY, RefSCMatrix& MZ,
                         RefSCMatrix& MXX, RefSCMatrix& MYY, RefSCMatrix& MZZ,
                         RefSCMatrix& MXY, RefSCMatrix& MXZ, RefSCMatrix& MYZ)
  {
    if (!space1->integral()->equiv(space2->integral()))
      throw ProgrammingError("two OrbitalSpaces use incompatible Integral factories");
    const Ref<GaussianBasisSet> bs1 = space1->basis();
    const Ref<GaussianBasisSet> bs2 = space2->basis();
    const bool bs1_eq_bs2 = (bs1 == bs2);
    int nshell1 = bs1->nshell();
    int nshell2 = bs2->nshell();

    RefSCMatrix vec1t = space1->coefs().t();
    RefSCMatrix vec2 = space2->coefs();

    Ref<Integral> localints = space1->integral()->clone();
    localints->set_basis(bs1,bs2);

    Ref<OneBodyInt> m1_ints = localints->dipole(0);
    Ref<OneBodyInt> m2_ints = localints->quadrupole(0);

    // form AO moment matrices
    RefSCDimension aodim1 = vec1t.coldim();
    RefSCDimension aodim2 = vec2.rowdim();
    Ref<SCMatrixKit> aokit = bs1->so_matrixkit();
    RefSCMatrix mx(aodim1, aodim2, aokit);
    RefSCMatrix my(aodim1, aodim2, aokit);
    RefSCMatrix mz(aodim1, aodim2, aokit);
    RefSCMatrix mxx(aodim1, aodim2, aokit);
    RefSCMatrix myy(aodim1, aodim2, aokit);
    RefSCMatrix mzz(aodim1, aodim2, aokit);
    mx.assign(0.0);
    my.assign(0.0);
    mz.assign(0.0);
    mxx.assign(0.0);
    myy.assign(0.0);
    mzz.assign(0.0);

    for(int sh1=0; sh1<nshell1; sh1++) {
      int bf1_offset = bs1->shell_to_function(sh1);
      int nbf1 = bs1->shell(sh1).nfunction();

      int sh2max;
      if (bs1_eq_bs2)
        sh2max = sh1;
      else
        sh2max = nshell2-1;

      for(int sh2=0; sh2<=sh2max; sh2++) {
        int bf2_offset = bs2->shell_to_function(sh2);
        int nbf2 = bs2->shell(sh2).nfunction();

        m1_ints->compute_shell(sh1,sh2);
        const double *m1intsptr = m1_ints->buffer();

        m2_ints->compute_shell(sh1,sh2);
        const double *m2intsptr = m2_ints->buffer();

        int bf1_index = bf1_offset;
        for(int bf1=0; bf1<nbf1; bf1++, bf1_index++, m1intsptr+=3*nbf2, m2intsptr+=6*nbf2) {
      int bf2_index = bf2_offset;
      const double *ptr1 = m1intsptr;
          const double *ptr2 = m2intsptr;
      int bf2max;
          if (bs1_eq_bs2 && sh1 == sh2)
            bf2max = bf1;
          else
        bf2max = nbf2-1;
      for(int bf2=0; bf2<=bf2max; bf2++, bf2_index++) {

        mx.set_element(bf1_index, bf2_index, *(ptr1++));
        my.set_element(bf1_index, bf2_index, *(ptr1++));
        mz.set_element(bf1_index, bf2_index, *(ptr1++));

            mxx.set_element(bf1_index, bf2_index, *(ptr2++));
            ptr2 += 2;
            myy.set_element(bf1_index, bf2_index, *(ptr2++));
            ptr2++;
            mzz.set_element(bf1_index, bf2_index, *(ptr2++));

          }
        }
      }
    }

    // and clean up a bit
    m1_ints = 0;
    m2_ints = 0;

    // Symmetrize matrices, if necessary
    if (bs1_eq_bs2) {

      const int nbasis = bs1->nbasis();

      for(int bf1=0; bf1<nbasis; bf1++)
        for(int bf2=0; bf2<=bf1; bf2++) {
          mx(bf2,bf1) = mx(bf1,bf2);
          my(bf2,bf1) = my(bf1,bf2);
          mz(bf2,bf1) = mz(bf1,bf2);
          mxx(bf2,bf1) = mxx(bf1,bf2);
          myy(bf2,bf1) = myy(bf1,bf2);
          mzz(bf2,bf1) = mzz(bf1,bf2);
        }

    }


    // finally, transform
    MX = vec1t * mx * vec2;
    MY = vec1t * my * vec2;
    MZ = vec1t * mz * vec2;
    MXX = vec1t * mxx * vec2;
    MYY = vec1t * myy * vec2;
    MZZ = vec1t * mzz * vec2;

    // and clean up a bit
    mx = 0;
    my = 0;
    mz = 0;
    mxx = 0;
    myy = 0;
    mzz = 0;

    //if (debug_ > 1) {
    //  MX.print("mu(X)");
    //  MY.print("mu(Y)");
    //  MZ.print("mu(Z)");
    //  MXX.print("mu(XX)");
    //  MYY.print("mu(YY)");
    //  MZZ.print("mu(ZZ)");
    //}
  }

};

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
