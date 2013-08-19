//
                                                                                                                                // fockbuilder.cc
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

#define EIGEN_NO_AUTOMATIC_RESIZING 0
#define TIMER_DEPTH 1
#define timer_change(str, depth) \
  if(TIMER_DEPTH >= depth) tim.change(str);
#define timer_enter(str, depth) \
  if(TIMER_DEPTH >= depth) tim.enter(str);
#define timer_exit(depth) \
  if(TIMER_DEPTH >= depth) tim.exit();
#define DO_EIGEN_OUT 0
#if DO_EIGEN_OUT
#  define eigenout(label, var) cout << "@@@@@ " label << " @@@@@" << endl << setprecision(20) << var << endl << "%%%%%%%%%%" << endl
#else
#  define eigenout(label, var)
#endif

#include <util/misc/sharedptr.h>
#include <cfloat>
#include<chemistry/qc/basis/symmint.h>
#include<chemistry/qc/basis/orthog.h>
#include<math/scmat/blas.h>
#include<math/scmat/svd.h>
#include<chemistry/qc/lcao/fockbuilder.h>
#include<util/misc/consumableresources.h>
#include<Eigen/Dense>

typedef Eigen::Map<Eigen::MatrixXd> EigenMatrixMap;
typedef Eigen::Map<const Eigen::MatrixXd> ConstEigenMatrixMap;
typedef Eigen::Map<Eigen::VectorXd> VectorMap;
typedef Eigen::MatrixXd EigenMatrix;
typedef std::pair<int, int> IntPair;
typedef Eigen::HouseholderQR<EigenMatrix> Decomposition;
typedef std::map<IntPair, Decomposition*> DecompositionMap;

typedef DecompositionMap::iterator DMap_iter;

using namespace std;
using namespace sc;

namespace sc {
  namespace detail {

    RefSymmSCMatrix overlap(const Ref<GaussianBasisSet>& bas,
                            const Ref<Integral>& integral) {

      Ref<Integral> localints = integral->clone();
      localints->set_basis(bas);
      Ref<PetiteList> pl = localints->petite_list();

      // form skeleton overlap in AO basis
      RefSymmSCMatrix sao(bas->basisdim(), bas->matrixkit());
      sao.assign(0.0);
      Ref<SCElementOp> sc =
          new OneBodyIntOp(new SymmOneBodyIntIter(localints->overlap(), pl));
      sao.element_op(sc);
      sc = 0;

      // now symmetrize Sso
      RefSymmSCMatrix s(pl->SO_basisdim(), bas->so_matrixkit());
      pl->symmetrize(sao, s);

      return s;
    }

    RefSCMatrix overlap(const Ref<GaussianBasisSet>& brabas,
                        const Ref<GaussianBasisSet>& ketbas,
                        const Ref<Integral>& integral) {

      Ref<Integral> localints = integral->clone();
      Ref<GPetiteList2> pl12 = GPetiteListFactory::plist2(brabas,ketbas);
      localints->set_basis(brabas,ketbas);

      // form overlap in AO basis
      RefSCMatrix sao(brabas->basisdim(), ketbas->basisdim(), brabas->matrixkit());
      sao.assign(0.0);

      abort();
    }

    void edipole(const Ref<GaussianBasisSet>& bas,
                 const Ref<Integral>& integral,
                 RefSymmSCMatrix& mu_x,
                 RefSymmSCMatrix& mu_y,
                 RefSymmSCMatrix& mu_z) {

      Ref<Integral> localints = integral->clone();
      localints->set_basis(bas);
      Ref<PetiteList> pl = localints->petite_list();
      Ref<OneBodyInt> m1_ints = localints->dipole(0);

      // form skeleton mu_i in AO basis
      RefSymmSCMatrix mu_x_ao(bas->basisdim(), bas->matrixkit()); mu_x_ao.assign(0.0);
      RefSymmSCMatrix mu_y_ao(bas->basisdim(), bas->matrixkit()); mu_y_ao.assign(0.0);
      RefSymmSCMatrix mu_z_ao(bas->basisdim(), bas->matrixkit()); mu_z_ao.assign(0.0);

      const int nshell = bas->nshell();
      for(int sh1=0; sh1<nshell; sh1++) {
        int bf1_offset = bas->shell_to_function(sh1);
        int nbf1 = bas->shell(sh1).nfunction();

        int sh2max = sh1;
        for(int sh2=0; sh2<=sh2max; sh2++) {
          int bf2_offset = bas->shell_to_function(sh2);
          int nbf2 = bas->shell(sh2).nfunction();

          m1_ints->compute_shell(sh1,sh2);
          const double *m1intsptr = m1_ints->buffer();

          int bf1_index = bf1_offset;
          for(int bf1=0; bf1<nbf1; bf1++, bf1_index++, m1intsptr+=3*nbf2) {
            int bf2_index = bf2_offset;
            const double *ptr1 = m1intsptr;
            int bf2max;
            if (sh1 == sh2)
              bf2max = bf1;
            else
              bf2max = nbf2-1;
            for(int bf2=0; bf2<=bf2max; bf2++, bf2_index++) {

              // the negative charge of the electron is not included
              mu_x_ao.set_element(bf1_index, bf2_index, ptr1[0]);
              mu_y_ao.set_element(bf1_index, bf2_index, ptr1[1]);
              mu_z_ao.set_element(bf1_index, bf2_index, ptr1[2]);
              ptr1 += 3;

            }
          }
        }
      }
      m1_ints = 0;

      const int nbasis = bas->nbasis();
      for(int bf1=0; bf1<nbasis; bf1++)
        for(int bf2=0; bf2<=bf1; bf2++) {
          mu_x_ao(bf2,bf1) = mu_x_ao(bf1,bf2);
          mu_y_ao(bf2,bf1) = mu_y_ao(bf1,bf2);
          mu_z_ao(bf2,bf1) = mu_z_ao(bf1,bf2);
        }


      // now symmetrize
      mu_x = bas->so_matrixkit()->symmmatrix(pl->SO_basisdim());
      mu_y = bas->so_matrixkit()->symmmatrix(pl->SO_basisdim());
      mu_z = bas->so_matrixkit()->symmmatrix(pl->SO_basisdim());
      pl->symmetrize(mu_x_ao, mu_x);
      pl->symmetrize(mu_y_ao, mu_y);
      pl->symmetrize(mu_z_ao, mu_z);
    }

    RefSymmSCMatrix nonrelativistic(const Ref<GaussianBasisSet>& bas,
                                    const Ref<Integral>& integral) {

      Ref<Integral> localints = integral->clone();
      localints->set_basis(bas);
      Ref<PetiteList> pl = localints->petite_list();

      // form skeleton Hcore in AO basis
      RefSymmSCMatrix hao(bas->basisdim(), bas->matrixkit());
      hao.assign(0.0);

      // kinetic energy
      Ref<SCElementOp> hc =
          new OneBodyIntOp(new SymmOneBodyIntIter(localints->kinetic(), pl));
      hao.element_op(hc);
      hc = 0;

      // molecular potential
      Ref<OneBodyInt> nuc = localints->nuclear();
      nuc->reinitialize();
      hc = new OneBodyIntOp(new SymmOneBodyIntIter(nuc, pl));
      hao.element_op(hc);
      hc = 0;

      // now symmetrize Hso
      RefSymmSCMatrix h(pl->SO_basisdim(), bas->so_matrixkit());
      pl->symmetrize(hao, h);

      return h;
    }

    RefSymmSCMatrix twobody_twocenter_coulomb(const Ref<GaussianBasisSet>& bas,
                                              const Ref<Integral>& integral) {

      Ref<Integral> localints = integral->clone();
      localints->set_basis(bas);
      Ref<PetiteList> pl = localints->petite_list();

      // form skeleton 2-center Coulomb in AO basis
      RefSymmSCMatrix ao(bas->basisdim(), bas->matrixkit());
      ao.assign(0.0);
      Ref<SCElementOp> sc =
          new TwoBodyTwoCenterIntOp(new SymmTwoBodyTwoCenterIntIter(localints->electron_repulsion2(),
                                                                    pl));
      ao.element_op(sc);
      sc = 0;

      // now symmetrize
      RefSymmSCMatrix result(pl->SO_basisdim(), bas->so_matrixkit());
      pl->symmetrize(ao, result);

      return result;
    }

    RefSymmSCMatrix pauli(const Ref<GaussianBasisSet>& bas,
                          const Ref<Integral>& integral) {
      const double c = 137.0359895; // speed of light in a vacuum in a.u.

      Ref<Integral> localints = integral->clone();
      localints->set_basis(bas);
      Ref<PetiteList> pl = localints->petite_list();

      // form skeleton Hcore in AO basis
      RefSymmSCMatrix hao(bas->basisdim(), bas->matrixkit());
      hao.assign(0.0);

      // mass-velocity
      Ref<OneBodyInt> mv = localints->p4();
      mv->reinitialize();
      Ref<SCElementOp> hc = new OneBodyIntOp(new SymmOneBodyIntIter(mv, pl));
      RefSymmSCMatrix mv_ao(bas->basisdim(), bas->matrixkit());
      mv_ao.assign(0.0);
      mv_ao.element_op(hc);
      mv_ao.scale((-1.0) / (8.0 * c * c));
      hao.accumulate(mv_ao);
      mv_ao = 0;
      hc = 0;

      // Darwin term
      RefSymmSCMatrix Darwin(bas->basisdim(), bas->matrixkit());
      Darwin.assign(0.0);
      {
        const double darwin_prefac = M_PI / (c * c * 2.0);

        GaussianBasisSet::ValueData* vdata1 =
            new GaussianBasisSet::ValueData(bas, localints);
        const int nbasis1 = bas->nbasis();
        double* values1 = new double[nbasis1];
        const Ref<Molecule>& molecule = bas->molecule();
        const int natom = molecule->natom();

        for (int iatom = 0; iatom < natom; iatom++) {

          SCVector3 R_iatom = molecule->r(iatom);
          // this puts values of basis functions evaluated at R_iatom into values1
          bas->values(R_iatom, vdata1, values1);
          const double prefac = darwin_prefac * molecule->charge(iatom);

          for (int ibasis = 0; ibasis < bas->nbasis(); ibasis++) {
            for (int jbasis = 0; jbasis <= ibasis; jbasis++) {

              const double d = values1[ibasis] * values1[jbasis] * prefac;
              Darwin.accumulate_element(ibasis, jbasis, d);
            }
          }
        }
        delete[] values1;
        delete vdata1;
      }
      hao.accumulate(Darwin);

      // now symmetrize Hso
      RefSymmSCMatrix h(pl->SO_basisdim(), bas->so_matrixkit());
      pl->symmetrize(hao, h);

      // add the nonrelativistic contribution
      RefSymmSCMatrix hnr = nonrelativistic(bas, integral);
      h->accumulate(hnr);

      return h;
    }

    namespace {

      void dk2_contrib(const RefSymmSCMatrix &h_pbas,
                       const RefDiagSCMatrix &E,
                       const RefDiagSCMatrix &K,
                       const RefDiagSCMatrix &p2,
                       const RefDiagSCMatrix &p2K2,
                       const RefDiagSCMatrix &p2K2_inv,
                       const RefSymmSCMatrix &AVA_pbas,
                       const RefSymmSCMatrix &BpVpB_pbas) {
        RefSCDimension pdim = AVA_pbas.dim();

        RefSymmSCMatrix AVA_prime = AVA_pbas.clone();
        RefSymmSCMatrix BpVpB_prime = BpVpB_pbas.clone();
        int npbas = pdim.n();
        for (int i=0; i<npbas; i++) {
          double Ei = E(i);
          for (int j=0; j<=i; j++) {
            double Ej = E(j);
            AVA_prime(i,j) = AVA_pbas(i,j)/(Ei+Ej);
            BpVpB_prime(i,j) = BpVpB_pbas(i,j)/(Ei+Ej);
          }
        }

        RefSCMatrix h_contrib;
        h_contrib
          =
          - 1.0 * BpVpB_prime * E * AVA_prime
          - 0.5 * BpVpB_prime * AVA_prime * E
          - 0.5 * AVA_prime * BpVpB_prime * E
          + 0.5 * AVA_prime * p2K2 * AVA_prime * E
          + 0.5 * BpVpB_prime * p2K2_inv * BpVpB_prime * E
          + 0.5 * AVA_prime * (p2K2 * E) * AVA_prime
          + 0.5 * BpVpB_prime * (p2K2_inv * E) * BpVpB_prime
          ;

        h_pbas.accumulate_symmetric_sum(h_contrib);
      }

      void dk3_contrib(const RefSymmSCMatrix &h_pbas,
                       const RefDiagSCMatrix &E,
                       const RefDiagSCMatrix &B,
                       const RefDiagSCMatrix &p2K2_inv,
                       const RefSCMatrix &so_to_p,
                       const RefSymmSCMatrix &pxVp) {
        RefSCDimension p_oso_dim = so_to_p.rowdim();
        Ref<SCMatrixKit> p_so_kit = so_to_p.kit();
        int noso = p_oso_dim.n();
        RefSymmSCMatrix pxVp_pbas(p_oso_dim, p_so_kit);
        pxVp_pbas.assign(0.0);
        pxVp_pbas.accumulate_transform(so_to_p, pxVp);

        RefSymmSCMatrix BpxVpB_prime(p_oso_dim, p_so_kit);
        for (int i=0; i<noso; i++) {
          double Ei = E(i);
          for (int j=0; j<=i; j++) {
            double Ej = E(j);
            BpxVpB_prime(i,j) = pxVp_pbas(i,j)*B(i)*B(j)/(Ei+Ej);
          }
        }

        RefSCDimension pdim = E.dim();

        RefSCMatrix h_contrib;
        h_contrib
          =
          - 0.5 * BpxVpB_prime * E * p2K2_inv * BpxVpB_prime
          - 0.5 * BpxVpB_prime * p2K2_inv * BpxVpB_prime * E
          ;

        h_pbas.accumulate_symmetric_sum(h_contrib);
      }

    } // anonymous namespace

    RefSymmSCMatrix dk(int dklev, const Ref<GaussianBasisSet>& bas,
                       const Ref<GaussianBasisSet>& p_bas,
                       const Ref<Integral>& integral) {

#define DK_DEBUG 0

      ExEnv::out0() << indent << "Including DK" << dklev
          << (dklev == 1 ? " (free particle projection)" : "")
          << (dklev == 2 ? " (Douglas-Kroll-Hess)" : "")
          << (dklev == 3 ? " (complete spin-free Douglas-Kroll)" : "")
          << " terms in the one body Hamiltonian." << std::endl;

      if (dklev > 2) {
        throw FeatureNotImplemented("dklev must be 0, 1, or 2", __FILE__, __LINE__);
      }

      Ref<Integral> localints = integral->clone();
      // The one electron integrals will be computed in the momentum basis.
      localints->set_basis(p_bas);

      Ref<PetiteList> p_pl = localints->petite_list();

      RefSCDimension p_so_dim = p_pl->SO_basisdim();
      RefSCDimension p_ao_dim = p_pl->AO_basisdim();
      Ref<SCMatrixKit> p_kit = p_bas->matrixkit();
      Ref<SCMatrixKit> p_so_kit = p_bas->so_matrixkit();

      // Compute the overlap in the momentum basis.
      RefSymmSCMatrix S_skel(p_ao_dim, p_kit);
      S_skel.assign(0.0);
      Ref<SCElementOp> hc =
          new OneBodyIntOp(new SymmOneBodyIntIter(localints->overlap(), p_pl));
      S_skel.element_op(hc);
      hc = 0;
      RefSymmSCMatrix S(p_so_dim, p_so_kit);
      p_pl->symmetrize(S_skel, S);

      ExEnv::out0() << indent << "The momentum basis is:" << std::endl;
      ExEnv::out0() << incindent;
      p_bas->print_brief(ExEnv::out0());
      ExEnv::out0() << decindent;

      ExEnv::out0() << indent << "Orthogonalizing the momentum basis"
          << std::endl;
      const int debug = 0;
      Ref<OverlapOrthog> p_orthog = new OverlapOrthog(OverlapOrthog::default_orthog_method(), S,
                                                      p_so_kit, OverlapOrthog::default_lindep_tol(),
                                                      debug);

      RefSCDimension p_oso_dim = p_orthog->orthog_dim();

      // form skeleton Hcore in the momentum basis
      RefSymmSCMatrix T_skel(p_ao_dim, p_kit);
      T_skel.assign(0.0);

      hc = new OneBodyIntOp(new SymmOneBodyIntIter(localints->kinetic(), p_pl));
      T_skel.element_op(hc);
      hc = 0;

      // finish constructing the kinetic energy integrals,
      // for which the skeleton is in hao
      RefSymmSCMatrix T(p_so_dim, p_so_kit);
      p_pl->symmetrize(T_skel, T);
      T_skel = 0;

      // Transform T into an orthogonal basis
      RefSymmSCMatrix T_oso(p_oso_dim, p_so_kit);
      T_oso.assign(0.0);
      T_oso.accumulate_transform(p_orthog->basis_to_orthog_basis(), T);

      // diagonalize the T integrals to get a momentum basis
      RefDiagSCMatrix Tval(p_oso_dim, p_so_kit);
      RefSCMatrix Tvec(p_oso_dim, p_oso_dim, p_so_kit);
      // Tvec * Tval * Tvec.t() = T_oso
      T_oso.diagonalize(Tval, Tvec);

      T_oso = 0;

#if DK_DEBUG
      T.print("T");
      Tval.print("Tval");
#endif

      // Compute the kinematic factors
      RefDiagSCMatrix A(p_oso_dim, p_so_kit);
      RefDiagSCMatrix B(p_oso_dim, p_so_kit);
      RefDiagSCMatrix E(p_oso_dim, p_so_kit);
      RefDiagSCMatrix K(p_oso_dim, p_so_kit);
      RefDiagSCMatrix p2(p_oso_dim, p_so_kit);
      RefDiagSCMatrix Emc2(p_oso_dim, p_so_kit);
      const double c = 137.0359895; // speed of light in a vacuum in a.u.
      int noso = p_oso_dim.n();
      for (int i = 0; i < noso; i++) {
        double T_val = Tval(i);
        // momentum basis sets with near linear dependencies may
        // have T_val approximately equal to zero.  These can be
        // negative, which will cause a SIGFPE in sqrt.
        if (T_val < DBL_EPSILON)
          T_val = 0.0;
        double p = sqrt(2.0 * T_val);
        double E_val = c * sqrt(p * p + c * c);
        double A_val = sqrt((E_val + c * c) / (2.0 * E_val));
        double K_val = c / (E_val + c * c);
        double B_val = A_val * K_val;
        double Emc2_val = c * c * p * p / (E_val + c * c); // = E - mc^2
        A( i) = A_val;
        B( i) = B_val;
        E( i) = E_val;
        K( i) = K_val;
        p2( i) = p * p;
        Emc2( i) = Emc2_val;
      }

#if DK_DEBUG
      A.print("A");
      B.print("B");
      E.print("E");
      K.print("K");
      Emc2.print("Emc2");
#endif

      // Construct the transform from the coordinate to the momentum
      // representation in the momentum basis
      RefSCMatrix so_to_p = Tvec.t() * p_orthog->basis_to_orthog_basis();

#if DK_DEBUG
      so_to_p.print("so_to_p");
#endif

      // compute the V integrals
      Ref<OneBodyInt> V_obi = localints->nuclear();
      V_obi->reinitialize();
      hc = new OneBodyIntOp(new SymmOneBodyIntIter(V_obi, p_pl));
      RefSymmSCMatrix V_skel(p_ao_dim, p_kit);
      V_skel.assign(0.0);
      V_skel.element_op(hc);
      V_obi = 0;
      hc = 0;
      RefSymmSCMatrix V(p_so_dim, p_so_kit);
      p_pl->symmetrize(V_skel, V);
      V_skel = 0;

#if DK_DEBUG
      V.print("V");
#endif

      // transform V to the momentum basis
      RefSymmSCMatrix V_pbas(p_oso_dim, p_so_kit);
      V_pbas.assign(0.0);
      V_pbas.accumulate_transform(so_to_p, V);

      // compute the p.Vp integrals
      Ref<OneBodyInt> pVp_obi = localints->p_dot_nuclear_p();
      hc = new OneBodyIntOp(new SymmOneBodyIntIter(pVp_obi, p_pl));
      RefSymmSCMatrix pVp_skel(p_ao_dim, p_kit);
      pVp_skel.assign(0.0);
      pVp_skel.element_op(hc);
#if DK_DEBUG
      const double *buf = pVp_obi->buffer();
      for (int I=0,Ii=0; I<p_bas->nshell(); I++) {
        for (int i=0; i<p_bas->shell(I).nfunction(); i++,Ii++) {
          for (int J=0,Jj=0; J<p_bas->nshell(); J++) {
            pVp_obi->compute_shell(I,J);
            int ij = i*p_bas->shell(J).nfunction();
            for (int j=0; j<p_bas->shell(J).nfunction(); j++,ij++,Jj++) {
              std::cout << "pVp["<<Ii<<"]["<<Jj<<"][0]= " << buf[ij]
              << std::endl;
            }
          }
        }
      }
#endif
      pVp_obi = 0;
      hc = 0;
      RefSymmSCMatrix pVp(p_so_dim, p_so_kit);
      p_pl->symmetrize(pVp_skel, pVp);
      pVp_skel = 0;

#if DK_DEBUG
      pVp.print("pVp");
      (-2.0*T).print("-2*T");
#endif

      // transform p.Vp to the momentum basis
      RefSymmSCMatrix pVp_pbas(p_oso_dim, p_so_kit);
      pVp_pbas.assign(0.0);
      pVp_pbas.accumulate_transform(so_to_p, pVp);

      RefSymmSCMatrix AVA_pbas(p_oso_dim, p_so_kit);
      RefSymmSCMatrix BpVpB_pbas(p_oso_dim, p_so_kit);
      for (int i = 0; i < noso; i++) {
        for (int j = 0; j <= i; j++) {
          AVA_pbas(i, j) = V_pbas(i, j) * A(i) * A(j);
          BpVpB_pbas(i, j) = pVp_pbas(i, j) * B(i) * B(j);
        }
      }

      V_pbas = 0;
      pVp_pbas = 0;

      // form the momentum basis hamiltonian
      RefSymmSCMatrix h_pbas(p_oso_dim, p_so_kit);
      h_pbas = AVA_pbas + BpVpB_pbas;

      // Add the kinetic energy
      for (int i = 0; i < noso; i++) {
        h_pbas(i, i) = h_pbas(i, i) + Emc2(i);
      }

      if (dklev > 1) {
        RefDiagSCMatrix p2K2 = p2 * K * K;
        RefDiagSCMatrix p2K2_inv = p2K2->clone();

        for (int i = 0; i < noso; i++) {
          double p2K2_val = p2K2(i);
          if (fabs(p2K2_val) > DBL_EPSILON)
            p2K2_inv( i) = 1.0 / p2K2_val;
          else
            p2K2_inv( i) = 0.0;
        }

        dk2_contrib(h_pbas, E, K, p2, p2K2, p2K2_inv, AVA_pbas, BpVpB_pbas);

        if (dklev > 2) {
          Ref<OneBodyInt> pxVp_obi = localints->p_cross_nuclear_p();
          Ref<SCElementOp3> hc3;
          hc3 = new OneBody3IntOp(new SymmOneBodyIntIter(pxVp_obi, p_pl));
          RefSymmSCMatrix pxVp_x_skel(p_ao_dim, p_kit);
          RefSymmSCMatrix pxVp_y_skel(p_ao_dim, p_kit);
          RefSymmSCMatrix pxVp_z_skel(p_ao_dim, p_kit);
          pxVp_x_skel.assign(0.0);
          pxVp_y_skel.assign(0.0);
          pxVp_z_skel.assign(0.0);
          pxVp_x_skel.element_op(hc3, pxVp_y_skel, pxVp_z_skel);
          RefSymmSCMatrix pxVp_x(p_so_dim, p_so_kit);
          RefSymmSCMatrix pxVp_y(p_so_dim, p_so_kit);
          RefSymmSCMatrix pxVp_z(p_so_dim, p_so_kit);
          p_pl->symmetrize(pxVp_x_skel, pxVp_x);
          p_pl->symmetrize(pxVp_y_skel, pxVp_y);
          p_pl->symmetrize(pxVp_z_skel, pxVp_z);
          pxVp_x_skel = 0;
          pxVp_y_skel = 0;
          pxVp_z_skel = 0;

          dk3_contrib(h_pbas, E, B, p2K2_inv, so_to_p, pxVp_x);
          dk3_contrib(h_pbas, E, B, p2K2_inv, so_to_p, pxVp_y);
          dk3_contrib(h_pbas, E, B, p2K2_inv, so_to_p, pxVp_z);
        }

      }

#if DK_DEBUG
      h_pbas.print("h_pbas");
#endif

      AVA_pbas = 0;
      BpVpB_pbas = 0;
      A = 0;
      B = 0;
      E = 0;
      K = 0;
      Emc2 = 0;

      // Construct the transform from the momentum representation to the
      // coordinate representation in the momentum basis
      RefSCMatrix p_to_so = p_orthog->basis_to_orthog_basis_inverse() * Tvec;

      // Construct the transform from the momentum basis to the
      // coordinate basis.
      localints->set_basis(bas, p_bas);
      Ref<PetiteList> pl = localints->petite_list();
      RefSCMatrix S_ao_p(pl->AO_basisdim(), p_ao_dim, p_kit);
      S_ao_p.assign(0.0);
      hc = new OneBodyIntOp(localints->overlap());
      S_ao_p.element_op(hc);
      hc = 0;
      // convert s_ao_p into the so ao and so p basis
      RefSCMatrix blocked_S_ao_p(pl->AO_basisdim(), p_pl->AO_basisdim(),
                                 p_so_kit);
      blocked_S_ao_p->convert(S_ao_p);
      RefSCMatrix S_ao_p_so_l = pl->sotoao() * blocked_S_ao_p;
      RefSCMatrix S_ao_p_so = S_ao_p_so_l * p_pl->aotoso();
      S_ao_p_so_l = 0;

      // transform h_pbas back to the so basis
      RefSymmSCMatrix h_dk_so(pl->SO_basisdim(), bas->so_matrixkit());
      h_dk_so.assign(0.0);
      h_dk_so.accumulate_transform(S_ao_p_so * p_orthog->overlap_inverse()
          * p_to_so, h_pbas);

      // Compute the overlap in bas
      localints->set_basis(bas);
      RefSymmSCMatrix S_bas;
      {
        Ref<SCMatrixKit> kit = bas->matrixkit();
        Ref<SCMatrixKit> so_kit = bas->so_matrixkit();
        RefSCDimension so_dim = pl->SO_basisdim();
        RefSCDimension ao_dim = pl->AO_basisdim();
        RefSymmSCMatrix S_skel(ao_dim, kit);
        S_skel.assign(0.0);
        hc
            = new OneBodyIntOp(
                               new SymmOneBodyIntIter(localints->overlap(), pl));
        S_skel.element_op(hc);
        hc = 0;
        S_bas = so_kit->symmmatrix(so_dim);
        pl->symmetrize(S_skel, S_bas);
      }

      localints->set_basis(bas);

#if DK_DEBUG
      {
        RefSCMatrix tmp = S_ao_p_so * p_orthog->overlap_inverse() * S_ao_p_so.t();
        tmp.print("S(OBS,pbasis) * S(pbasis)^-1 * S(pbasis,OBS)");
        ExEnv::out0() << indent << " trace = " << tmp.trace() << endl;
        S_bas.print("S(OBS)");

        ExEnv::out0() << indent << " trace = " << S_bas.trace() << endl;

        ExEnv::out0() << indent << "nso = " << pl->SO_basisdim()->n() << endl;
      }
#endif

      // Check to see if the momentum basis spans the coordinate basis.  The
      // following approach seems reasonable, but a more careful mathematical
      // analysis would be desirable.
      const double S_ao_projected_trace = (S_ao_p_so
          * p_orthog->overlap_inverse() * S_ao_p_so.t()).trace()
          / pl->SO_basisdim()->n();
      const double S_ao_trace = S_bas.trace() / pl->SO_basisdim()->n();
      const double completeness_diagnostic = S_ao_projected_trace / S_ao_trace;
      ExEnv::out0() << indent
          << "Tr(basis overlap projected into momentum basis)/Tr(basis overlap) = "
          << completeness_diagnostic << std::endl;
      if (fabs(1.0 - completeness_diagnostic) > OverlapOrthog::default_lindep_tol()) {
        ExEnv::out0() << indent
            << "WARNING: the momentum basis does not span the orbital basis"
            << std::endl;
      }

#if DK_DEBUG
      S_ao_p_so.print("S_ao_p_so");
      p_to_so.print("p_to_so");
      //(p_to_so*so_to_p).print("p_to_so*so_to_p");
      (S_ao_p_so*S.gi()*p_to_so).print("S_ao_p_so*S.gi()*p_to_so");
#endif

#if DK_DEBUG
      (T+V).print("T+V");
      h_dk_so.print("h_dk_so");
#endif

      return h_dk_so;

    }

    RefSCMatrix coulomb(const Ref<TwoBodyFourCenterMOIntsRuntime>& ints_rtime,
                        const Ref<OrbitalSpace>& occ_space,
                        const Ref<OrbitalSpace>& bra_space,
                        const Ref<OrbitalSpace>& ket_space) {

      Ref<MessageGrp> msg = MessageGrp::get_default_messagegrp();

      Timer tim_coulomb("coulomb");

      int me = msg->me();
      int nproc = msg->n();
      //ExEnv::out0() << endl << indent
      //         << "Entered Coulomb matrix evaluator" << endl;
      //ExEnv::out0() << incindent;

      // Only need 1/r12 integrals. In principle, almost any Descr will do. For now ask for ERI
      const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(occ_space->id(),bra_space->id(),
                                                             occ_space->id(),ket_space->id(),
                                                             std::string("ERI"),
                                                             std::string(TwoBodyIntLayout::b1k1_b2k2));
      Ref<TwoBodyMOIntsTransform> mnxy_tform = ints_rtime->get(tform_key);
      Ref<DistArray4> mnxy_acc = mnxy_tform->ints_distarray4();

      const int nocc = occ_space->rank();
      const int nbra = bra_space->rank();
      const int nket = ket_space->rank();
      const blasint nbraket = nbra*nket;

      ExEnv::out0() << indent << "Begin computation of Coulomb matrix" << endl;

      // Compute the number of tasks that have full access to the integrals
      // and split the work among them
      vector<int> proc_with_ints;
      int nproc_with_ints = mnxy_acc->tasks_with_access(proc_with_ints);

      //////////////////////////////////////////////////////////////
      //
      // Evaluation of the coulomb matrix proceeds as follows:
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

      double* J_xy = new double[nbraket];
      memset(J_xy,0,nbraket*sizeof(double));
      Timer tim_mo_ints_retrieve;
      tim_mo_ints_retrieve.set_default("MO ints retrieve");
      if (mnxy_acc->has_access(me)) {

        for(int m=0; m<nocc; m++) {

          const int mm = m*nocc+m;
          const int mm_proc = mm%nproc_with_ints;
          if (mm_proc != proc_with_ints[me])
            continue;

          // Get (|1/r12|) integrals
          tim_mo_ints_retrieve.enter_default();
          const double *mmxy_buf_eri = mnxy_acc->retrieve_pair_block(m,m,TwoBodyOper::eri);
          tim_mo_ints_retrieve.exit_default();

          const double one = 1.0;
          const blasint unit_stride = 1;
          F77_DAXPY(&nbraket,&one,mmxy_buf_eri,&unit_stride,J_xy,&unit_stride);

          mnxy_acc->release_pair_block(m,m,TwoBodyOper::eri);
        }
      }

      ExEnv::out0() << indent << "End of computation of Coulomb matrix" << endl;

      msg->sum(J_xy,nbraket);

      RefSCMatrix J(bra_space->coefs()->coldim(), ket_space->coefs()->coldim(), bra_space->coefs()->kit());
      J.assign(J_xy);
      delete[] J_xy;

      //ExEnv::out0() << decindent;
      //ExEnv::out0() << indent << "Exited Coulomb matrix evaluator" << endl;
      tim_coulomb.exit();

      return J;
    }

    RefSCMatrix exchange(const Ref<TwoBodyFourCenterMOIntsRuntime>& ints_rtime,
                         const Ref<OrbitalSpace>& occ_space,
                         const Ref<OrbitalSpace>& bra_space,
                         const Ref<OrbitalSpace>& ket_space) {

      Ref<MessageGrp> msg = MessageGrp::get_default_messagegrp();

      Timer tim_exchange("exchange");

      int me = msg->me();
      int nproc = msg->n();
      //ExEnv::out0() << endl << indent
      //         << "Entered exchange matrix evaluator" << endl;
      //ExEnv::out0() << incindent;

      // Only need 1/r12 integrals. In principle, almost any Descr will do. For now ask for ERI
      const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(occ_space->id(),occ_space->id(),
                                                             bra_space->id(),ket_space->id(),
                                                             std::string("ERI"),
                                                             std::string(TwoBodyIntLayout::b1b2_k1k2));
      Ref<TwoBodyMOIntsTransform> mxny_tform = ints_rtime->get(tform_key);
      Ref<DistArray4> mnxy_acc = mxny_tform->ints_distarray4();

      const int nocc = occ_space->rank();
      const int nbra = bra_space->rank();
      const int nket = ket_space->rank();
      const blasint nbraket = nbra*nket;

      ExEnv::out0() << indent << "Begin computation of exchange matrix" << endl;

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

          // Get (|1/r12|) integrals
          tim_mo_ints_retrieve.enter_default();
          const double *mmxy_buf_eri = mnxy_acc->retrieve_pair_block(m,m,TwoBodyOper::eri);
          tim_mo_ints_retrieve.exit_default();

          const double one = 1.0;
          const blasint unit_stride = 1;
          F77_DAXPY(&nbraket,&one,mmxy_buf_eri,&unit_stride,K_xy,&unit_stride);

          mnxy_acc->release_pair_block(m,m,TwoBodyOper::eri);
        }
      }

      ExEnv::out0() << indent << "End of computation of exchange matrix" << endl;

      msg->sum(K_xy,nbraket);

      RefSCMatrix K(bra_space->coefs()->coldim(), ket_space->coefs()->coldim(), bra_space->coefs()->kit());
      K.assign(K_xy);
      delete[] K_xy;

      //ExEnv::out0() << decindent;
      //ExEnv::out0() << indent << "Exited exchange matrix evaluator" << endl;
      tim_exchange.exit();

      return K;
    }

    RefSCMatrix coulomb_df(const Ref<DensityFittingInfo>& df_info,
                           const RefSymmSCMatrix& P,
                           const Ref<GaussianBasisSet>& brabs,
                           const Ref<GaussianBasisSet>& ketbs,
                           const Ref<GaussianBasisSet>& obs) {

      Ref<MessageGrp> msg = MessageGrp::get_default_messagegrp();

      Timer tim("coulomb(DF)");
      const double wall_time_start = tim.wall_time("coulomb(DF)");

      int me = msg->me();
      int nproc = msg->n();
      //ExEnv::out0() << endl << indent
      //         << "Entered Coulomb(DF) matrix evaluator" << endl;
      //ExEnv::out0() << incindent;

      //////////////////////////////////////////////////////////////
      //
      // Evaluation of the coulomb matrix proceeds as follows:
      //
      // density-fit (mu nu| . P_{mu nu} -> R^Mu
      // contract cC_ab^Mu with R^Mu to produce J_ab
      //
      /////////////////////////////////////////////////////////////////////////////////

      const Ref<AOSpaceRegistry> ao_registry = df_info->runtime()->moints_runtime()->factory()->ao_registry();
      const Ref<OrbitalSpace>& braspace = ao_registry->value(brabs);
      const Ref<OrbitalSpace>& ketspace = ao_registry->value(ketbs);
      const Ref<OrbitalSpace>& dfspace = ao_registry->value(df_info->params()->basis());
      const Ref<OrbitalSpace>& obs_space = ao_registry->value(obs);
      const blasint ndf = dfspace->rank();

      const Ref<DensityFittingRuntime>& df_rtime = df_info->runtime();
      const Ref<TwoBodyThreeCenterMOIntsRuntime>& int3c_rtime = df_rtime->moints_runtime()->runtime_3c();

      std::string kernel_key = df_info->params()->kernel_key();
      const bool noncoulomb_kernel = (kernel_key.find("exp") != std::string::npos);
      std::string params_key = df_info->params()->intparams_key();
      std::string operset_key = noncoulomb_kernel ? "G12'" : "ERI";
      TwoBodyOper::type kernel_oper = noncoulomb_kernel ? TwoBodyOper::r12_0_g12 : TwoBodyOper::eri;
      unsigned int ints_type_idx = TwoBodyOperSetDescr::instance(noncoulomb_kernel ? TwoBodyOperSet::G12NC : TwoBodyOperSet::ERI)->opertype(kernel_oper);

      std::vector<double> R(ndf,0.0);
      {
        ///////////////////////////////////////
        // contract 3-center fitting kernel_key ints with density
        ///////////////////////////////////////
        const std::string C_key = ParsedTwoBodyThreeCenterIntKey::key(obs_space->id(),
                                                                      dfspace->id(),
                                                                      obs_space->id(),
                                                                      operset_key, params_key);
        tim.enter(C_key);
        const Ref<TwoBodyThreeCenterMOIntsTransform>& C_tform = int3c_rtime->get(C_key);
        C_tform->compute();
        Ref<DistArray4> C = C_tform->ints_acc();  C->activate();
        tim.exit();

        tim.enter("contract density");
        const blasint nobs = obs_space->rank();
        const blasint ndf = dfspace->rank();
        std::vector<double> Q(ndf,0.0);
        {
          // Compute the number of tasks that have full access to the integrals
          // and split the work among them
          vector<int> proc_with_ints;
          int nproc_with_ints = C->tasks_with_access(proc_with_ints);

          if (C->has_access(me)) {

            double** P_blk = new double*[nobs];
            P_blk[0] = new double[nobs*nobs];
            for(int i=1; i<nobs; ++i) {
              P_blk[i] = P_blk[i-1] + nobs;
            }
            for(int i=0; i<nobs; ++i)
              for(int j=0; j<=i; ++j)
                P_blk[i][j] = P_blk[j][i] = P.get_element(i,j);

            for (int p = 0; p < nobs; ++p) {
              if (not C->is_local(0,p))
                continue;

              const double* C_p_qR_buf = C->retrieve_pair_block(0, p, ints_type_idx);
              const double* P_p_q = P_blk[p];

              const char notrans = 'n';
              const double one = 1.0;
              const blasint unit_stride = 1;
              F77_DGEMV(&notrans, &ndf, &nobs, &one, C_p_qR_buf, &ndf, P_p_q,
                        &unit_stride, &one, &(Q[0]), &unit_stride);

              C->release_pair_block(0, p, ints_type_idx);

            }

            // cleanup
            delete[] P_blk[0];
            delete[] P_blk;

          }

          // sum all contributions
          msg->sum(&(Q[0]), ndf);
        }
        tim.exit();

        // intermediate cleanup
        if (C->data_persistent()) C->deactivate();
        C = 0;

        ///////
        // compute 2-center fitting kernel_key
        ///////
        RefSymmSCMatrix kernel = dfspace->coefs_nb()->kit()->symmmatrix(dfspace->dim());
        {
          const std::string kernel_key = ParsedTwoBodyTwoCenterIntKey::key(dfspace->id(),
                                                                           dfspace->id(),
                                                                           operset_key, params_key);
          tim.enter(kernel_key);
          RefSCMatrix kernel_rect = df_rtime->moints_runtime()->runtime_2c()->get(kernel_key);
          kernel.assign_subblock(kernel_rect, 0, ndf-1, 0, ndf-1);
          tim.exit();
        }

        //////////////////////////////////////////////
        // density fit Q using Cholesky inverse
        //////////////////////////////////////////////
        tim.enter("density fit");
        {
          // check if kernel_key fit already exists
          RefSymmSCMatrix kernel_i_mat;
          if ( not df_rtime->moints_runtime()->runtime_2c_inv()->key_exists(dfspace->id()) ) {
            kernel_i_mat = kernel.copy();
            lapack_invert_symmposdef(kernel_i_mat, 1e10);
            df_rtime->moints_runtime()->runtime_2c_inv()->add(dfspace->id(), kernel_i_mat);
          }
          kernel_i_mat = df_rtime->moints_runtime()->runtime_2c_inv()->value(dfspace->id());
          double* kernel_i_rect = allocate<double>(ndf * ndf);
          for(int r=0, rc=0; r<ndf; ++r) {
            for(int c=r; c<ndf; ++c, ++rc) {
              double value = kernel_i_mat(r,c);
              kernel_i_rect[r*ndf + c] = value;
              kernel_i_rect[c*ndf + r] = value;
            }
          }
          {
            char notransp = 'n';
            blasint n = ndf;
            double one = 1.0;
            blasint ione = 1;
            double zero = 0.0;
            F77_DGEMV(&notransp, &n, &n, &one, kernel_i_rect, &n, &Q[0], &ione, &zero, &R[0], &ione);
          }
          deallocate(kernel_i_rect);
#if 0
          {
          std::vector<double> kernel_packed;  // only needed for factorized methods
          std::vector<double> kernel_factorized;
          // convert kernel_ to a packed upper-triangle form
          kernel_packed.resize(ndf * (ndf + 1) / 2);
          kernel_key->convert(&(kernel_packed[0]));
          // factorize kernel_ using diagonal pivoting from LAPACK's DSPTRF
          kernel_factorized.resize(ndf * (ndf + 1) / 2);
          sc::lapack_cholesky_symmposdef(kernel_key,
                                         &(kernel_factorized[0]),
                                         1e10);

          const bool refine_solution = false;
          sc::lapack_linsolv_cholesky_symmposdef(&(kernel_packed[0]), ndf,
                                                 &(kernel_factorized[0]),
                                                 &(R[0]), &(Q[0]), 1,
                                                 refine_solution);
          }

          for(int i=0; i<ndf; ++i) {
            std::cout << R[i] << " " << Rnew[i] << std::endl;
          }
#endif

        }
        tim.exit();

      }

      // finish out by multiplying with the 3-center fitting kernel_key integrals
      const std::string cC_key = ParsedTwoBodyThreeCenterIntKey::key(braspace->id(),
                                                                     dfspace->id(),
                                                                     ketspace->id(),
                                                                     "ERI","");
      tim.enter(cC_key);
      const Ref<TwoBodyThreeCenterMOIntsTransform>& cC_tform = int3c_rtime->get(cC_key);
      cC_tform->compute();
      Ref<DistArray4> cC = cC_tform->ints_acc();  cC->activate();
      tim.exit();

      tim.enter("assemble J");
      const int nbra = braspace->rank();
      const blasint nket = ketspace->rank();
      const int nbraket = nbra * nket;
      std::vector<double> J(nbraket, 0.0);
      {
        // Compute the number of tasks that have full access to the integrals
        // and split the work among them
        vector<int> proc_with_ints;
        int nproc_with_ints = cC->tasks_with_access(proc_with_ints);

        if (cC->has_access(me)) {

          for (int a = 0; a < nbra; ++a) {
            if (not cC->is_local(0, a))
              continue;

            const double* cC_a_bR_buf = cC->retrieve_pair_block(0, a, 0);
            double* J_a = &(J[a*nket]);

            const char trans = 't';
            const double one = 1.0;
            const blasint unit_stride = 1;
            F77_DGEMV(&trans, &ndf, &nket, &one, cC_a_bR_buf, &ndf, &(R[0]),
                      &unit_stride, &one, J_a, &unit_stride);

            cC->release_pair_block(0, a, 0);

          }

        }

        // sum all contributions
        msg->sum(&(J[0]), nbraket);
      }
      tim.exit();

      if (cC->data_persistent()) cC->deactivate();  cC = 0;

      Ref<Integral> localints = int3c_rtime->factory()->integral()->clone();
      localints->set_basis(brabs);
      Ref<PetiteList> brapl = localints->petite_list();
      localints->set_basis(ketbs);
      Ref<PetiteList> ketpl = localints->petite_list();
      RefSCDimension bradim = brapl->AO_basisdim();
      RefSCDimension ketdim = ketpl->AO_basisdim();
      RefSCMatrix result(bradim,
                         ketdim,
                         brabs->so_matrixkit());

      result.assign(&(J[0]));

      tim.exit();
      //ExEnv::out0() << decindent;
      //ExEnv::out0() << indent << "Exited Coulomb(DF) matrix evaluator ("
      //              << (tim.wall_time("coulomb(DF)") - wall_time_start) << " sec)" << endl;

      return result;
    }

    RefSCMatrix coulomb_df_local(const Ref<DensityFittingInfo>& df_info,
                                 const RefSymmSCMatrix& P,
                                 const Ref<GaussianBasisSet>& brabs,
                                 const Ref<GaussianBasisSet>& ketbs,
                                 const Ref<GaussianBasisSet>& obs)
    {
      /*=======================================================================================*/
      /* Setup and stuff                                       		                        {{{1 */ #if 1 // begin fold
      //----------------------------------------//
      Timer tim("coulomb(DF local)");
      timer_enter("01 - setup", 1);
      //----------------------------------------//
      Ref<MessageGrp> msg = MessageGrp::get_default_messagegrp();
      int me = msg->me();
      int nproc = msg->n();
      //----------------------------------------//
      const Ref<DensityFittingRuntime>& df_rtime = df_info->runtime();
      const Ref<AOSpaceRegistry> ao_registry = df_rtime->moints_runtime()->factory()->ao_registry();
      const Ref<OrbitalSpace>& braspace = ao_registry->value(brabs);
      const Ref<OrbitalSpace>& ketspace = ao_registry->value(ketbs);
      const Ref<GaussianBasisSet>& dfbs = df_info->params()->basis();
      assert(dfbs.nonnull());
      const Ref<OrbitalSpace>& dfspace = ao_registry->value(dfbs);
      const Ref<OrbitalSpace>& obs_space = ao_registry->value(obs);
      //----------------------------------------//
      const blasint dfnbf = dfspace->rank();
      const int branbf = brabs->nbasis();
      const int ketnbf = ketbs->nbasis();
      const int obsnbf = obs->nbasis();
      //----------------------------------------//
      const Ref<TwoBodyThreeCenterMOIntsRuntime>& int3c_rtime = df_rtime->moints_runtime()->runtime_3c();
      const Ref<TwoBodyTwoCenterMOIntsRuntime>& int2c_rtime = df_rtime->moints_runtime()->runtime_2c();
      //----------------------------------------//
      std::string metric_key = df_info->params()->kernel_key();
      const bool noncoulomb_kernel = (metric_key.find("exp") != std::string::npos);
      std::string params_key = df_info->params()->intparams_key();
      std::string operset_key = noncoulomb_kernel ? "G12'" : "ERI";
      TwoBodyOper::type metric_oper =
          noncoulomb_kernel ? TwoBodyOper::r12_0_g12 : TwoBodyOper::eri;
      unsigned int ints_type_idx = TwoBodyOperSetDescr::instance(
          noncoulomb_kernel ? TwoBodyOperSet::G12NC : TwoBodyOperSet::ERI
      )->opertype(metric_oper);
      /*****************************************************************************************/ #endif //1}}}
      /*=======================================================================================*/
      /* The old way                                           		                        {{{1 */ #if 0 // begin fold
      /*=======================================================================================*/
      /* Get the coefficients                                  		                        {{{1 */
      //========================================//
      timer_change("02 - get coefficients", 1);
      //----------------------------------------//
      // Setup
      timer_enter("01 - (mu nu | M | X) setup", 2);
      //----------------------------------------//
      // Go ahead and compute them all for now.  It will be faster to
      //   compute only the ones we need and save them for later, but
      //   that requires significant infrastructure changes.
      const std::string munu_M_X_key = ParsedTwoBodyThreeCenterIntKey::key(obs_space->id(),
                                                                           dfspace->id(),
                                                                           obs_space->id(),
                                                                           operset_key, params_key);
      timer_change("02 - (mu nu | M | X) compute", 2);
      const Ref<TwoBodyThreeCenterMOIntsTransform>& munu_M_X_tform = int3c_rtime->get(munu_M_X_key);
      munu_M_X_tform->compute();
      Ref<DistArray4> munu_M_X = munu_M_X_tform->ints_acc(); munu_M_X->activate();
      //----------------------------------------//
      // < X Y | metric | >
      // in chemists notation: ( X | metric | Y )
      timer_change("03 - compute (X | M | Y)", 2);
      const std::string metric2c_key = ParsedTwoBodyTwoCenterIntKey::key(
          dfspace->id(),
          dfspace->id(),
          operset_key, params_key
      );
      RefSCMatrix metric_2c_ints = int2c_rtime->get(metric2c_key);
      //----------------------------------------//
      // loop over atom pairs...
      timer_change("04 - loop over basis function pairs", 2);
      EigenMatrix C(dfnbf, branbf * ketnbf);
      C = EigenMatrix::Zero(dfnbf, branbf * ketnbf);
      if(munu_M_X->has_access(me)){
        DecompositionMap decomps;
        for(int ibf = 0; ibf < branbf; ++ibf){
          if(not munu_M_X->is_local(0, ibf))
            continue;
          //----------------------------------------//
          // Convenience variables for atom A
          const int ishA = brabs->function_to_shell(ibf);
          const int atomA = brabs->shell_to_center(ishA);
          const int nshA = brabs->nshell_on_center(atomA);
          const int nbfA = brabs->nbasis_on_center(atomA);
          const int shoffA = brabs->shell_on_center(atomA, 0);
          const int bfoffA = brabs->shell_to_function(shoffA);
          const int dfnshA = dfbs->nshell_on_center(atomA);
          const int dfnbfA = dfbs->nbasis_on_center(atomA);
          const int dfshoffA = dfbs->shell_on_center(atomA, 0);
          const int dfbfoffA = dfbs->shell_to_function(dfshoffA);
          //----------------------------------------//
          // TODO Permutational symmetry
          for(int jbf = 0; jbf < ketnbf; ++jbf){
            //----------------------------------------//
            // Convenience variables for atom B
            const int jshB = ketbs->function_to_shell(jbf);
            const int atomB = ketbs->shell_to_center(jshB);
            const int nshB = ketbs->nshell_on_center(atomB);
            const int nbfB = ketbs->nbasis_on_center(atomB);
            const int shoffB = ketbs->shell_on_center(atomB, 0);
            const int bfoffB = ketbs->shell_to_function(shoffB);
            const int dfnshB = dfbs->nshell_on_center(atomB);
            const int dfnbfB = dfbs->nbasis_on_center(atomB);
            const int dfshoffB = dfbs->shell_on_center(atomB, 0);
            const int dfbfoffB = dfbs->shell_to_function(dfshoffB);
            int dfpart_size = dfnbfA;
            if(atomA != atomB) dfpart_size += dfnbfB;
            /* old code */ #if 0
            IntPair atomAB(atomA, atomB);
            /*---------------------------------------------------------------*/
            /*  Get and Decompose ( X_(ab) | M | Y_(ab) )               {{{2 */ #if 2 // begin fold
            //===============================================================//
            timer_enter("01 - get ( X_(ab) | M | Y_(ab) ) block", 3);
            if(decomps.count(atomAB) == 0){
              double* metric_2c_part_ptr = allocate<double>(dfpart_size*dfpart_size);
              double** blockptr = allocate<double*>(dfpart_size);
              //---------------------------------------------------------------//
              // Get the block corresponding to (A|m|A)
              for(int idfpart = 0; idfpart < dfnbfA; ++idfpart) {
                blockptr[idfpart] = &(metric_2c_part_ptr[idfpart * dfpart_size]);
              }
              {
                RefSCMatrix block = metric_2c_ints.get_subblock(
                    dfbfoffA, dfbfoffA+dfnbfA-1,
                    dfbfoffA, dfbfoffA+dfnbfA-1
                );
                block.convert(blockptr);
              }
              // if atomA == atomB, we're done.  Otherwise,
              if(atomA != atomB){
                // Get the block corresponding to (A|m|B)
                for(int idfpart = 0; idfpart < dfnbfA; ++idfpart) {
                  blockptr[idfpart] = &(metric_2c_part_ptr[idfpart * dfpart_size + dfnbfA]);
                }
                {
                  RefSCMatrix block = metric_2c_ints.get_subblock(
                      dfbfoffA, dfbfoffA+dfnbfA-1,
                      dfbfoffB, dfbfoffB+dfnbfB-1
                  );
                  block.convert(blockptr);
                }
                //---------------------------------------------------------------//
                // Get the block corresponding to (B|m|A)
                for(int idfpart = 0; idfpart < dfnbfB; ++idfpart) {
                  blockptr[idfpart] = &(metric_2c_part_ptr[(idfpart+dfnbfA) * dfpart_size]);
                }
                {
                  RefSCMatrix block = metric_2c_ints.get_subblock(
                      dfbfoffB, dfbfoffB+dfnbfB-1,
                      dfbfoffA, dfbfoffA+dfnbfA-1
                  );
                  block.convert(blockptr);
                }
                //---------------------------------------------------------------//
                // Get the block corresponding to (B|m|B)
                for(int idfpart = 0; idfpart < dfnbfB; ++idfpart) {
                  blockptr[idfpart] = &(metric_2c_part_ptr[(idfpart+dfnbfA) * dfpart_size + dfnbfA]);
                }
                {
                  RefSCMatrix block = metric_2c_ints.get_subblock(
                      dfbfoffB, dfbfoffB+dfnbfB-1,
                      dfbfoffB, dfbfoffB+dfnbfB-1
                  );
                  block.convert(blockptr);
                }
              }
              //-----------------------------------------------------------------//
              EigenMatrixMap X_m_Y(metric_2c_part_ptr, dfpart_size, dfpart_size);
              timer_change("02 - decompose ( X(ab) | M | Y(ab) )", 3);
              decomps[atomAB] = new Decomposition(X_m_Y);
              deallocate(blockptr);
              deallocate(metric_2c_part_ptr);
            }
            /*****************************************************************/ #endif //2}}}
            /*---------------------------------------------------------------*/
            timer_change("03 - solve equations", 3);
            //----------------------------------------//
            // Allocate (ibf jbf | M | X(AB) )
            Eigen::VectorXd ibfjbf_M_X(dfpart_size);
            //----------------------------------------//
            // get (ibf jbf | M | X(A) )
            double* ibfjbf_M_X_buffer_A = allocate<double>(dfnbfA);
            VectorMap ibfjbf_M_X_A(ibfjbf_M_X_buffer_A, dfnbfA);
            munu_M_X->retrieve_pair_subblock(
                0, ibf,     // index in unit basis, index in mu
                ints_type_idx,
                jbf,        jbf+1,             // nu start, nu fence
                dfbfoffA,   dfbfoffA + dfnbfA, // X start,  X fence
                ibfjbf_M_X_buffer_A
            );
            ibfjbf_M_X.head(dfnbfA) = ibfjbf_M_X_A;
            deallocate(ibfjbf_M_X_buffer_A);
            //----------------------------------------//
            // get (ibf jbf | M | X(B) )
            if(atomA != atomB){
              double* ibfjbf_M_X_buffer_B = allocate<double>(dfnbfB);
              VectorMap ibfjbf_M_X_B(ibfjbf_M_X_buffer_B, dfnbfB);
              munu_M_X->retrieve_pair_subblock(
                  0, ibf,     // index in unit basis, index in mu
                  ints_type_idx,
                  jbf,        jbf+1,             // nu start, nu fence
                  dfbfoffB,   dfbfoffB + dfnbfB, // X start,  X fence
                  ibfjbf_M_X_buffer_B
              );
              ibfjbf_M_X.tail(dfnbfB) = ibfjbf_M_X_B;
              deallocate(ibfjbf_M_X_buffer_B);
            }
            //----------------------------------------//
            // Now solve A * x = b, with A = ( X_(ab) | M | Y_(ab) ) and b = ( ibf jbf | M | Y_(ab) )
            // The result will be dfpart_size rows and 1 column and needs to be stored in the big C
            Eigen::VectorXd Cpart(dfpart_size);
            Cpart = decomps[atomAB]->solve(ibfjbf_M_X);
            #endif
            //----------------------------------------//
            // and store in the big C
            std::string dfkey = ParsedDensityFittingKey::key(braspace->id(), ketspace->id(), dfspace->id(), metric_key);
            const std::shared_ptr<Eigen::VectorXd> Cpart = df_info->runtime()->get(dfkey, ibf, jbf);
            C.block(
                dfbfoffA, ibf*ketnbf + jbf,
                dfnbfA, 1
            ) = Cpart->head(dfnbfA);
            if(atomA != atomB) {
              C.block(
                  dfbfoffB, ibf*ketnbf + jbf,
                  dfnbfB, 1
              ) = Cpart->tail(dfnbfB);
            }
            //----------------------------------------//
            timer_exit(3);
          } // end loop over basis functions (jbf)
        } // end loop over basis functions (ibf)
        // clean up memory
        for(DMap_iter it=decomps.begin(); it != decomps.end(); ++it){
          delete it->second;
        }
      }
      timer_change("05 - global sum C", 2);
      msg->sum(C.data(), dfnbf * branbf * ketnbf);
      timer_exit(2);
      /*****************************************************************************************/
      /*=======================================================================================*/
      /* Get ( mu nu | g | X )                                 		                        {{{1 */
      //----------------------------------------//
      // < mu X | coulomb | nu >
      // in chemists notation: ( mu nu | coulomb | X )
      timer_change("03 - get (mu nu | g | X)", 1);
      timer_enter("01 - setup", 2);
      EigenMatrix munu_g_X(branbf*ketnbf, dfnbf);
      unsigned int g_type_idx = TwoBodyOperSetDescr::instance(TwoBodyOperSet::ERI)->opertype(metric_oper);
      timer_change("02 - compute", 2);
      Ref<DistArray4> munu_g_X_D4 = munu_M_X;
      // Only need to recompute if we're using a non-coulomb kernel
      if(noncoulomb_kernel){
        if (munu_M_X->data_persistent()) munu_M_X->deactivate();
        munu_M_X = 0;
        const std::string munu_g_X_key = ParsedTwoBodyThreeCenterIntKey::key(
            braspace->id(),
            dfspace->id(),
            ketspace->id(),
            "ERI", ""
        );
        const Ref<TwoBodyThreeCenterMOIntsTransform>& munu_g_X_tform = int3c_rtime->get(munu_g_X_key);
        munu_g_X_tform->compute();
        munu_g_X_D4 = munu_g_X_tform->ints_acc(); munu_g_X_D4->activate();
      }
      timer_change("03 - transfer to eigen", 2);
      double* ibfjbf_g_X_buffer = allocate<double>(dfnbf);
      VectorMap ibfjbf_g_X(ibfjbf_g_X_buffer, dfnbf);
      for(int ibf = 0; ibf < branbf; ++ibf){
        for(int jbf = 0; jbf < ketnbf; ++jbf){
          munu_g_X_D4->retrieve_pair_subblock(
              0, ibf,      // unit basis index, mu_index
              g_type_idx,
              jbf, jbf+1,  // nu_start, nu_fence
              0,   dfnbf,  // X_start, X_fence
              ibfjbf_g_X_buffer
          );
          munu_g_X.row(ibf*ketnbf + jbf) = ibfjbf_g_X;
        }
      }
      deallocate(ibfjbf_g_X_buffer);
      if (munu_g_X_D4->data_persistent()) munu_g_X_D4->deactivate();
      munu_g_X_D4 = 0;
      timer_exit(2);
      /*****************************************************************************************/
      /*=======================================================================================*/
      /* Get ( X | g | Y )                                     		                        {{{1 */
      timer_change("04 - get (X | g | Y)", 1);
      // hopefully this doesn't do a copy?
      timer_enter("01 - setup", 2);
      RefSCMatrix coulomb_2c_ints = metric_2c_ints.pointer();
      timer_change("02 - compute", 2);
      if(noncoulomb_kernel){
        const std::string coulomb2c_key = ParsedTwoBodyTwoCenterIntKey::key(
            dfspace->id(),
            dfspace->id(),
            "ERI", ""
        );
        coulomb_2c_ints = int2c_rtime->get(metric2c_key);
      }
      timer_change("03 - transfer", 2);
      double* coulomb_2c_ints_ptr = allocate<double>(dfnbf*dfnbf);
      coulomb_2c_ints.convert(coulomb_2c_ints_ptr);
      EigenMatrixMap X_g_Y(coulomb_2c_ints_ptr, dfnbf, dfnbf);
      timer_exit(2);
      /*****************************************************************************************/
      /*=======================================================================================*/
      /* Compute J                                             		                        {{{1 */
      timer_change("05 - compute J", 1);
      Eigen::VectorXd J(branbf*ketnbf);
      double *P_ptr = allocate<double>(branbf*ketnbf);
      {
        RefSymmSCMatrix Ptmp = P.copy();
        Ptmp.convert2RefSCMat().convert(P_ptr);
      }
      VectorMap D(P_ptr, branbf*ketnbf);
      Eigen::VectorXd CD = C * D;
      Eigen::VectorXd Dg = munu_g_X.transpose() * D;
      // C is X by (mu nu)
      // CD is X (by 1)
      // Dg is X (by 1)
      Eigen::VectorXd GT(dfnbf);
      GT = X_g_Y * CD;
      //----------------------------------------//
      J = C.transpose() * Dg;
      J += munu_g_X * CD;
      J -= CD.transpose() * X_g_Y * C;
      //----------------------------------------//
      // Include exact semidiagonal integral contributions
      //   if that option is enabled
      Ref<Integral> integral = int3c_rtime->factory()->integral();
      integral->set_basis(brabs, obs, ketbs, obs);
      Ref<TwoBodyInt> eri_eval = integral->electron_repulsion();
      const double* buffer = eri_eval->buffer();
      if(df_info->params()->exact_diag_J()){
        for(int atomA = 0; atomA < brabs->ncenter(); ++atomA){
          const int nshA = brabs->nshell_on_center(atomA);
          const int nbfA = brabs->nbasis_on_center(atomA);
          const int shoffA = brabs->shell_on_center(atomA, 0);
          const int bfoffA = brabs->shell_to_function(shoffA);
          const int dfnshA = dfbs->nshell_on_center(atomA);
          const int dfnbfA = dfbs->nbasis_on_center(atomA);
          const int dfshoffA = dfbs->shell_on_center(atomA, 0);
          const int dfbfoffA = dfbs->shell_to_function(dfshoffA);
          //----------------------------------------//
          for(int atomB = 0; atomB < ketbs->ncenter(); ++atomB){
            const int nshB = ketbs->nshell_on_center(atomB);
            const int nbfB = ketbs->nbasis_on_center(atomB);
            const int shoffB = ketbs->shell_on_center(atomB, 0);
            const int bfoffB = ketbs->shell_to_function(shoffB);
            const int dfnshB = dfbs->nshell_on_center(atomB);
            const int dfnbfB = dfbs->nbasis_on_center(atomB);
            const int dfshoffB = dfbs->shell_on_center(atomB, 0);
            const int dfbfoffB = dfbs->shell_to_function(dfshoffB);
            //----------------------------------------//
            int dfpart_size = dfnbfA;
            if(atomA != atomB) dfpart_size += dfnbfB;
            Eigen::VectorXd dt(dfpart_size);
            dt = Eigen::VectorXd::Zero(dfpart_size);
            Eigen::VectorXd ct(dfpart_size);
            ct = Eigen::VectorXd::Zero(dfpart_size);
            const double perm_fact =  (atomA == atomB ? 1.0 : 2.0);
            //----------------------------------------//
            for(int ishA = shoffA; ishA < shoffA + nshA; ++ishA){
              const int inbfA = brabs->shell(ishA).nfunction();
              const int ibfoffA = brabs->shell_to_function(ishA);
              for(int jshB = shoffB; jshB < shoffB + nshB; ++jshB){
                const int jnbfB = ketbs->shell(jshB).nfunction();
                const int jbfoffB = ketbs->shell_to_function(jshB);
                //----------------------------------------//
                for(int ibfA = 0; ibfA < inbfA; ++ibfA){
                  const int mu = ibfoffA + ibfA;
                  for(int jbfB = 0; jbfB < jnbfB; ++jbfB){
                    const int nu = jbfoffB + jbfB;
                    const int munu = mu*ketnbf + nu;
                    //----------------------------------------//
                    for(int XshA = dfshoffA; XshA < dfshoffA + dfnshA; ++XshA){
                      const int XnbfA = dfbs->shell(XshA).nfunction();
                      const int XbfoffA = dfbs->shell_to_function(XshA);
                      for(int XbfA = 0; XbfA < XnbfA; ++XbfA){
                        const int X = XbfoffA + XbfA;
                        dt(XbfoffA - dfbfoffA + XbfA) += perm_fact * munu_g_X(munu, X) * D(munu);
                        ct(XbfoffA - dfbfoffA + XbfA) += perm_fact * C(X, munu) * D(munu);
                      }
                    }
                    //----------------------------------------//
                    if(atomA != atomB){
                      for(int XshB = dfshoffB; XshB < dfshoffB + dfnshB; ++XshB){
                        const int XnbfB = dfbs->shell(XshB).nfunction();
                        const int XbfoffB = dfbs->shell_to_function(XshB);
                        for(int XbfB = 0; XbfB < XnbfB; ++XbfB){
                          const int X = XbfoffB + XbfB;
                          dt(XbfoffB - dfbfoffB + dfnbfA + XbfB) += perm_fact * munu_g_X(munu, X) * D(munu);
                          ct(XbfoffB - dfbfoffB + dfnbfA + XbfB) += perm_fact * C(X, munu) * D(munu);
                        }
                      }
                    }
                    //----------------------------------------//
                  }
                }
                //----------------------------------------//
              }
            }
            Eigen::VectorXd gt(dfpart_size);
            gt = Eigen::VectorXd::Constant(dfpart_size, 9e99);
            gt.head(dfnbfA) = X_g_Y.block(
                dfbfoffA, dfbfoffA,
                dfnbfA, dfnbfA
            ) * ct.head(dfnbfA);
            if(atomA != atomB){
              gt.head(dfnbfA) += X_g_Y.block(
                  dfbfoffA, dfbfoffB,
                  dfnbfA, dfnbfB
              ) * ct.tail(dfnbfB);
              gt.tail(dfnbfB) = X_g_Y.block(
                  dfbfoffB, dfbfoffB,
                  dfnbfB, dfnbfB
              ) * ct.tail(dfnbfB);
              gt.tail(dfnbfB) += X_g_Y.block(
                  dfbfoffB, dfbfoffA,
                  dfnbfB, dfnbfA
              ) * ct.head(dfnbfA);
            }
            //----------------------------------------//
            for(int ishA = shoffA; ishA < shoffA + nshA; ++ishA){
              const int inbfA = brabs->shell(ishA).nfunction();
              const int ibfoffA = brabs->shell_to_function(ishA);
              for(int jshB = shoffB; jshB < shoffB + nshB; ++jshB){
                const int jnbfB = ketbs->shell(jshB).nfunction();
                const int jbfoffB = ketbs->shell_to_function(jshB);
                //----------------------------------------//
                for(int kshA = shoffA; kshA < shoffA + nshA; ++kshA){
                  const int knbfA = obs->shell(kshA).nfunction();
                  const int kbfoffA = obs->shell_to_function(kshA);
                  for(int lshB = shoffB; lshB < shoffB + nshB; ++lshB){
                    const int lnbfB = obs->shell(lshB).nfunction();
                    const int lbfoffB = obs->shell_to_function(lshB);
                    //----------------------------------------//
                    // Compute shell (ishA jshB | kshA lshB) and add contributions to J
                    eri_eval->compute_shell(ishA, jshB, kshA, lshB);
                    int buff_off = 0;
                    for(int ibfA = 0; ibfA < inbfA; ++ibfA){
                      const int mu = ibfoffA + ibfA;
                      for(int jbfB = 0; jbfB < jnbfB; ++jbfB){
                        const int nu = jbfoffB + jbfB;
                        const int munu = mu*ketnbf + nu;
                        for(int kbfA = 0; kbfA < knbfA; ++kbfA){
                          const int rho = kbfoffA + kbfA;
                          for(int lbfB = 0; lbfB < lnbfB; ++lbfB){
                            const int sigma = lbfoffB + lbfB;
                            const int rhosigma = rho*obsnbf + sigma;
                            //                    (mu nu | rho sigma) * D^(rho sigma)
                            J(munu) += perm_fact * buffer[buff_off++] * D(rhosigma);
                          }
                        }
                      }
                    }
                    //----------------------------------------//
                  } // end loop over lshB
                } // end loop over kshA
                //----------------------------------------//
                for(int ibfA = 0; ibfA < inbfA; ++ibfA){
                  const int mu = ibfoffA + ibfA;
                  for(int jbfB = 0; jbfB < jnbfB; ++jbfB){
                    const int nu = jbfoffB + jbfB;
                    const int munu = mu*ketnbf + nu;
                    J(munu) -= dt.head(dfnbfA).transpose() * C.col(munu).segment(dfbfoffA,  dfnbfA);
                    J(munu) -= munu_g_X.row(munu).segment(dfbfoffA, dfnbfA) * ct.head(dfnbfA);
                    J(munu) += gt.head(dfnbfA).transpose() * C.col(munu).segment(dfbfoffA, dfnbfA);
                    if(atomA != atomB){
                      J(munu) -= dt.tail(dfnbfB).transpose() * C.col(munu).segment(dfbfoffB, dfnbfB);
                      J(munu) -= munu_g_X.row(munu).segment(dfbfoffB, dfnbfB) * ct.tail(dfnbfB);
                      J(munu) += gt.tail(dfnbfB).transpose() * C.col(munu).segment(dfbfoffB, dfnbfB);
                    }
                  }
                }
                //----------------------------------------//
              }
            }
            //----------------------------------------//
          }
        }
      }
      //----------------------------------------//
      /*****************************************************************************************/ #endif //1}}}
      /*=======================================================================================*/
      /* BEGIN NEW APPROACH                                                                    */
      /*=======================================================================================*/
      /* Loop over basis function pairs to form C_tilde and d_tilde                       {{{1 */ #if 1 // begin fold
      timer_change("02 - form C_tilde", 1);
      timer_enter("misc", 2);
      //----------------------------------------//
      Eigen::VectorXd J(branbf*ketnbf);
      J = Eigen::VectorXd::Zero(branbf*ketnbf);
      double *P_ptr = allocate<double>(obsnbf*obsnbf);
      {
        RefSymmSCMatrix Ptmp = P.copy();
        Ptmp.convert2RefSCMat().convert(P_ptr);
      }
      VectorMap D(P_ptr, obsnbf*obsnbf);
      //----------------------------------------//
      // Get the munu_g_X key and the munu_g_X_tform for later (we will need it now to
      //   determine data locality)
      timer_change("01 - compute (mu nu | g | X )", 2);
      const std::string munu_g_X_key = ParsedTwoBodyThreeCenterIntKey::key(
          braspace->id(),
          dfspace->id(),
          ketspace->id(),
          "ERI", ""
      );
      //----------------------------------------//
      const Ref<TwoBodyThreeCenterMOIntsTransform>& munu_g_X_tform = int3c_rtime->get(munu_g_X_key);
      munu_g_X_tform->compute();
      Ref<DistArray4> munu_g_X = munu_g_X_tform->ints_acc(); munu_g_X->activate();
      //----------------------------------------//
      timer_change("02 - compute Ctilde", 2);
      timer_enter("misc", 3);
      Eigen::VectorXd Ctilde(dfnbf);
      Ctilde = Eigen::VectorXd::Zero(dfnbf);
      std::string dfkey = ParsedDensityFittingKey::key(
          braspace->id(),
          ketspace->id(),
          dfspace->id(),
          metric_key
      );
      for(int mu = 0; mu < obsnbf; ++mu){
        if(not munu_g_X->is_local(0, mu))
          continue;
        //----------------------------------------//
        const int ishA = obs->function_to_shell(mu);
        const int atomA = obs->shell_to_center(ishA);
        const int dfnshA = dfbs->nshell_on_center(atomA);
        const int dfnbfA = dfbs->nbasis_on_center(atomA);
        const int dfshoffA = dfbs->shell_on_center(atomA, 0);
        const int dfbfoffA = dfbs->shell_to_function(dfshoffA);
        //----------------------------------------//
        for(int nu = 0; nu < obsnbf; ++nu){
          const int jshB = obs->function_to_shell(nu);
          const int atomB = obs->shell_to_center(jshB);
          const int dfnshB = dfbs->nshell_on_center(atomB);
          const int dfnbfB = dfbs->nbasis_on_center(atomB);
          const int dfshoffB = dfbs->shell_on_center(atomB, 0);
          const int dfbfoffB = dfbs->shell_to_function(dfshoffB);
          //----------------------------------------//
          // Assume locality is the same for munu_g_X and the coefficients
          timer_change("01 - get coefficients", 3);
          std::shared_ptr<Eigen::VectorXd> Cpart = df_rtime->get(dfkey, mu, nu);
          timer_change("02 - contract coefficients with D", 3);
          Ctilde.segment(dfbfoffA, dfnbfA) += D(mu*obsnbf + nu) * Cpart->head(dfnbfA);
          if(atomA != atomB){
            Ctilde.segment(dfbfoffB, dfnbfB) += D(mu*obsnbf + nu) * Cpart->tail(dfnbfB);
          }
          //----------------------------------------//
          timer_change("misc", 3);
        } // end loop over nu
      } // end loop over mu
      timer_change("03 - global sum C_tilde", 3);
      msg->sum(Ctilde.data(), dfnbf);
      timer_exit(3);
      timer_exit(2);
      /*****************************************************************************************/ #endif //1}}}
      /*=======================================================================================*/
      /* Form gtilde                                           		                        {{{1 */ #if 1 // begin fold
      timer_change("03 - Form gtilde", 1);
      timer_enter("01 - compute (X | g | Y)", 2);
      Eigen::VectorXd gtilde(dfnbf);
      {
        const std::string coulomb2c_key = ParsedTwoBodyTwoCenterIntKey::key(
            dfspace->id(),
            dfspace->id(),
            "ERI", ""
        );
        RefSCMatrix coulomb_2c_ints = int2c_rtime->get(coulomb2c_key);
        //----------------------------------------//
        timer_change("02 - transfer (X | g | Y)", 2);
        double* coulomb_2c_ints_ptr = allocate<double>(dfnbf*dfnbf);
        coulomb_2c_ints.convert(coulomb_2c_ints_ptr);
        EigenMatrixMap X_g_Y(coulomb_2c_ints_ptr, dfnbf, dfnbf);
        //----------------------------------------//
        timer_change("03 - contract with Ctilde", 2);
        gtilde = X_g_Y * Ctilde;
        deallocate(coulomb_2c_ints_ptr);
      }
      //----------------------------------------//
      timer_exit(2);
      /*****************************************************************************************/ #endif //1}}}
      /*=======================================================================================*/
      /* Loop over basis function pairs in (mu nu | g | X)     		                        {{{1 */ #if 1 // begin fold
      timer_change("03 - loop over (mu nu | g | X)", 1);
      unsigned int g_type_idx = TwoBodyOperSetDescr::instance(TwoBodyOperSet::ERI)->opertype(metric_oper);
      Eigen::VectorXd gpart(dfnbf);
      Eigen::VectorXd dtilde(dfnbf);
      dtilde = Eigen::VectorXd::Zero(dfnbf);
      for(int mu = 0; mu < obsnbf; ++mu){
        if(not munu_g_X->is_local(0, mu))
          continue;
        //----------------------------------------//
        const int ishA = obs->function_to_shell(mu);
        const int atomA = obs->shell_to_center(ishA);
        const int dfnshA = dfbs->nshell_on_center(atomA);
        const int dfnbfA = dfbs->nbasis_on_center(atomA);
        const int dfshoffA = dfbs->shell_on_center(atomA, 0);
        const int dfbfoffA = dfbs->shell_to_function(dfshoffA);
        //----------------------------------------//
        for(int nu = 0; nu < obsnbf; ++nu){
          const int jshB = obs->function_to_shell(nu);
          const int atomB = obs->shell_to_center(jshB);
          const int dfnshB = dfbs->nshell_on_center(atomB);
          const int dfnbfB = dfbs->nbasis_on_center(atomB);
          const int dfshoffB = dfbs->shell_on_center(atomB, 0);
          const int dfbfoffB = dfbs->shell_to_function(dfshoffB);
          //----------------------------------------//
          munu_g_X->retrieve_pair_subblock(
              0, mu,      // unit basis index, mu_index
              g_type_idx,
              nu, nu+1,  // nu_start, nu_fence
              0,  dfnbf,  // X_start, X_fence
              gpart.data()
          );
          //----------------------------------------//
          // dtilde contribution
          dtilde += D(mu*obsnbf + nu) * gpart;
          //----------------------------------------//
          // J contribution from second term
          J(mu*obsnbf + nu) += Ctilde.transpose() * gpart;
        }
      }
      // Global sum dtilde
      timer_enter("global sum dtilde", 2);
      msg->sum(dtilde.data(), dfnbf);
      timer_exit(2);
      /*****************************************************************************************/ #endif //1}}}
      /*=======================================================================================*/
      /* Loop over basis function pairs for first and third term contributions to J       {{{1 */ #if 1 // begin fold
      timer_change("04 - contributions to J", 1);
      // Combine dtilde and gtilde beforehand
      dtilde -= gtilde;
      // now the first and third terms can be done together
      for(int mu = 0; mu < obsnbf; ++mu){
        if(not munu_g_X->is_local(0, mu))
          continue;
        //----------------------------------------//
        const int ishA = obs->function_to_shell(mu);
        const int atomA = obs->shell_to_center(ishA);
        const int dfnshA = dfbs->nshell_on_center(atomA);
        const int dfnbfA = dfbs->nbasis_on_center(atomA);
        const int dfshoffA = dfbs->shell_on_center(atomA, 0);
        const int dfbfoffA = dfbs->shell_to_function(dfshoffA);
        //----------------------------------------//
        for(int nu = 0; nu < obsnbf; ++nu){
          const int jshB = obs->function_to_shell(nu);
          const int atomB = obs->shell_to_center(jshB);
          const int dfnshB = dfbs->nshell_on_center(atomB);
          const int dfnbfB = dfbs->nbasis_on_center(atomB);
          const int dfshoffB = dfbs->shell_on_center(atomB, 0);
          const int dfbfoffB = dfbs->shell_to_function(dfshoffB);
          //----------------------------------------//
          std::shared_ptr<Eigen::VectorXd> Cpart = df_rtime->get(dfkey, mu, nu);
          //----------------------------------------//
          // First and third term contributions to J added together
          J(mu*obsnbf + nu) += Cpart->head(dfnbfA).transpose() * dtilde.segment(dfbfoffA, dfnbfA);
          if(atomA != atomB){
            J(mu*obsnbf + nu) += Cpart->tail(dfnbfB).transpose() * dtilde.segment(dfbfoffB, dfnbfB);
          }
        } // end loop over nu
      } // end loop over mu
      timer_enter("global sum J", 2);
      msg->sum(J.data(), branbf*ketnbf);
      timer_exit(2);
      /*****************************************************************************************/ #endif //1}}}
      /*=======================================================================================*/
      /* Transfer J to a RefSCMatrix                           		                        {{{1 */ #if 1 // begin fold
      timer_change("05 - transfer J to RefSCMatrix", 1);
      Ref<Integral> localints = int3c_rtime->factory()->integral()->clone();
      localints->set_basis(brabs);
      Ref<PetiteList> brapl = localints->petite_list();
      localints->set_basis(ketbs);
      Ref<PetiteList> ketpl = localints->petite_list();
      RefSCDimension bradim = brapl->AO_basisdim();
      RefSCDimension ketdim = ketpl->AO_basisdim();
      RefSCMatrix result(
          bradim,
          ketdim,
          brabs->so_matrixkit()
      );
      result.assign(J.data());
      /*****************************************************************************************/ #endif //1}}}
      /*=======================================================================================*/
      /* Cleanup stuff                                         		                        {{{1 */ #if 1 // begin fold
      //----------------------------------------//
      timer_change("06 - cleanup", 1);
      deallocate(P_ptr);
      timer_exit(1);
      tim.exit();
      return result;
      /*****************************************************************************************/ #endif //1}}}
      /*=======================================================================================*/
    }

    RefSCMatrix exchange_df(const Ref<DensityFittingInfo>& df_info,
                            const RefSymmSCMatrix& P,
                            SpinCase1 spin,
                            const Ref<GaussianBasisSet>& brabs,
                            const Ref<GaussianBasisSet>& ketbs,
                            const Ref<GaussianBasisSet>& obs,
                            const Ref<FockBuildRuntime::PSqrtRegistry>& psqrtregistry) {

      Ref<MessageGrp> msg = MessageGrp::get_default_messagegrp();

      Timer tim("exchange(DF)");
      const double wall_time_start = tim.wall_time("exchange(DF)");

      int me = msg->me();
      int nproc = msg->n();
      //ExEnv::out0() << endl << indent
      //         << "Entered exchange(DF) matrix evaluator" << endl;
      //ExEnv::out0() << incindent;

      //////////////////////////////////////////////////////////////
      //
      // Evaluation of the exchange matrix proceeds as follows:
      // compute S so that S*S^t = P (S is a "square root" of the density)
      // density-fit (a S| -> C_aS^Mu
      // contract C_aS^Mu with cC_bS^Mu to yield K_ab
      //
      /////////////////////////////////////////////////////////////////////////////////

      const Ref<AOSpaceRegistry> ao_registry = df_info->runtime()->moints_runtime()->factory()->ao_registry();
      const Ref<OrbitalSpace>& braspace = ao_registry->value(brabs);
      const Ref<OrbitalSpace>& ketspace = ao_registry->value(ketbs);
      const Ref<OrbitalSpace>& dfspace = ao_registry->value(df_info->params()->basis());
      const Ref<OrbitalSpace>& obs_space = ao_registry->value(obs);

      const Ref<DensityFittingRuntime>& df_rtime = df_info->runtime();
      const Ref<TwoBodyThreeCenterMOIntsRuntime>& int3c_rtime = df_rtime->moints_runtime()->runtime_3c();

      // get square root of P
      Ref<OrbitalSpace> Sspace;
      if (psqrtregistry->key_exists(P) == false) {

        RefDiagSCMatrix Pevals = P.eigvals();
        RefSCMatrix Pevecs = P.eigvecs();
        const double Peval_max = Pevals->maxabs();
        const double Peval_threshold = Peval_max * 1e-8;
        const int nao = obs_space->rank();
        int Srank = 0;
        for(int ao=0; ao<nao; ++ao) {
          const double value = Pevals(ao);
          assert(value > -1e-8); // Negative eigenvalues? BAD
          if (value > Peval_threshold)
            ++Srank;
        }

        // handle the case of zero electrons
        if (Srank == 0) {
          Ref<Integral> localints = int3c_rtime->factory()->integral()->clone();
          localints->set_basis(brabs);
          Ref<PetiteList> brapl = localints->petite_list();
          localints->set_basis(ketbs);
          Ref<PetiteList> ketpl = localints->petite_list();
          RefSCDimension bradim = brapl->AO_basisdim();
          RefSCDimension ketdim = ketpl->AO_basisdim();
          RefSCMatrix result(bradim,
                             ketdim,
                             brabs->so_matrixkit());
          result.assign(0.0);
          //ExEnv::out0() << decindent;
          //ExEnv::out0() << indent << "Exited exchange(DF) matrix evaluator" << endl;
          tim.exit();
          return result;
        }

        RefSCDimension Sdim = new SCDimension(Srank, 1); Sdim->blocks()->set_subdim(0, new SCDimension(Srank));
        RefSCMatrix S = Pevecs.kit()->matrix(obs_space->coefs().coldim(), Sdim);
        {  // compute S from Pevecs
          // since both are blocked matrices with 1 block, must operate on the blocks directly
          RefSCMatrix Sblk = S.block(0);
          RefSCMatrix Ublk = Pevecs.block(0);
          int s = 0;
          for (int ao = 0; ao < nao; ++ao) {
            const double peval = Pevals(ao);
            if (peval > Peval_threshold) {
              Sblk.assign_subblock(Ublk.get_subblock(0, nao - 1, ao, ao)
                                   * sqrt(peval), 0, nao - 1, s, s);
              ++s;
            }
          }
        }

        // make an orbital space out of sqrt(P)
        const std::string skey = ParsedOrbitalSpaceKey::key(std::string("dd"),spin);
        Sspace = new OrbitalSpace(skey, prepend_spincase(spin,"sqrt(P)"),
                                  obs_space->coefs() * S,
                                  obs_space->basis(),
                                  obs_space->integral());
        psqrtregistry->add(P, Sspace);
      }
      else { // Sspace is in registry
        Sspace = psqrtregistry->value(P);
      }

      /// add sqrt(P) space to the OrbitalSpaceRegistry (it will be removed after the computation)
      Ref<OrbitalSpaceRegistry> oreg = df_info->runtime()->moints_runtime()->factory()->orbital_registry();
      oreg->add(Sspace->id(), Sspace);

      // fit (bra S|
      const std::string C_key = ParsedDensityFittingKey::key(braspace->id(),
                                                             Sspace->id(),
                                                             dfspace->id(),
                                                             df_info->params()->kernel_key());
      Ref<DistArray4> C = df_rtime->get(C_key);  C->activate();
//      {
//        Ref<DistArray4> Ctmp = permute23(C);
//        if (C->data_persistent()) C->deactivate();
//        C = Ctmp;  C->activate();
//      }

      // get (ket S| integrals
      const std::string cC_key = ParsedTwoBodyThreeCenterIntKey::key(ketspace->id(),
                                                                     dfspace->id(),
                                                                     Sspace->id(),
                                                                     "ERI","");
      const Ref<TwoBodyThreeCenterMOIntsTransform>& cC_tform = int3c_rtime->get(cC_key);
      cC_tform->compute();
      Ref<DistArray4> cC = cC_tform->ints_acc();  cC->activate();

      const int ndf = dfspace->rank();
      const blasint nsdf = Sspace->rank() * ndf;
      const unsigned int nbra = braspace->rank();
      const unsigned int nket = ketspace->rank();
      const unsigned int nbraket = nbra * nket;
      double* K = new double[nbraket];
      memset(K,0,sizeof(double)*nbraket);

      {
        // Compute the number of tasks that have full access to the integrals
        // and split the work among them
        // assume that C and cC have the same access
        vector<int> proc_with_ints;
        int nproc_with_ints = C->tasks_with_access(proc_with_ints);

        if (C->has_access(me)) {

          // figure out how big of a block of a/b can be held in memory
          // each block is nsdf doubles big
          const size_t nbytes_per_block = nsdf * sizeof(double);
          const size_t max_num_SR_blocks = ConsumableResources::get_default_instance()->memory() / (2 * nbytes_per_block); // 2 to account for bra and ket storage
          if (max_num_SR_blocks < 1)
            throw MemAllocFailed("not enough memory to hold 1 block, increase memory in ConsumableResources",
                                 __FILE__, __LINE__, 2 * nbytes_per_block);
          const size_t max_blk_size = max_num_SR_blocks; // based on how much memory we have
          const size_t min_blk_size = std::min(nbra,nket)/static_cast<size_t>(floor(sqrt(nproc_with_ints))); // based on the amount of parallelism

          // emphasize distributed concurrency over minimizing local bandwidth
          // i.e. prefer smaller blocks but more units of work
          const size_t blk_size = std::min(max_blk_size, min_blk_size);
          const size_t nblk_a = (nbra + blk_size - 1) / blk_size;
          const size_t nblk_b = (nket + blk_size - 1) / blk_size;
          const size_t blk_size_a = (nbra + nblk_a - 1) / nblk_a;
          const size_t blk_size_b = (nket + nblk_b - 1) / nblk_b;

          std::vector<double> C_a_pR_bufs(blk_size_a * nsdf);
          std::vector<double> cC_b_pR_bufs(blk_size_b * nsdf);
          ConsumableResources::get_default_instance()->consume_memory( (blk_size_a+blk_size_b) * nbytes_per_block );
          // will read a blocks in the innermost loop, but lazily
          bool read_a_bufs = false;
          // will hold one block of the result
          std::vector<double> K_blk(blk_size_a * blk_size_b);

          int a_offset = 0;
          for(int blk_a = 0; blk_a<nblk_a; ++blk_a, a_offset+=blk_size_a) {

            const size_t bsize_a = nbra > (a_offset + blk_size_a) ? blk_size_a : nbra - a_offset;

            int b_offset = 0;
            for(int blk_b = 0; blk_b<nblk_b; ++blk_b, b_offset+=blk_size_b) {

              const size_t task_id = blk_a * nblk_b + blk_b;
              const int proc = task_id % nproc_with_ints;
              if (proc != proc_with_ints[me])
                continue;

              const size_t bsize_b = nket > (b_offset + blk_size_b) ? blk_size_b : nket - b_offset;

              if (not read_a_bufs) {
                for(int a=0; a<bsize_a; ++a)
                  C->retrieve_pair_block(0, a_offset + a, 0,
                                         &C_a_pR_bufs[a*nsdf]);
                read_a_bufs = true;
              }

              for(int b=0; b<bsize_b; ++b)
                cC->retrieve_pair_block(0, b_offset + b,
                                        TwoBodyOper::eri,
                                        &cC_b_pR_bufs[b*nsdf]);

              C_DGEMM('n', 't', bsize_a, bsize_b, nsdf,
                      1.0, &C_a_pR_bufs[0], nsdf,
                      &cC_b_pR_bufs[0], nsdf, 0.0, &K_blk[0], bsize_b);

              for(int b=0; b<bsize_b; ++b)
                cC->release_pair_block(0, b_offset + b, TwoBodyOper::eri);

              for(int a=0, ab=0; a<bsize_a; ++a) {
                const int aa = a_offset + a;
                for(int b=0; b<bsize_b; ++b, ++ab) {
                  const int bb = b_offset + b;
                  K[aa * nket + bb] += K_blk[ab];
                }
              }

            } // end of b block loop

            if (read_a_bufs) {
              for(int a=0; a<bsize_a; ++a)
                C->release_pair_block(0, a_offset + a, 0);
              read_a_bufs = false;
            }

          } // end of a block loop


          ConsumableResources::get_default_instance()->release_memory( (blk_size_a+blk_size_b) * nbytes_per_block );
          //C_a_pR_bufs.~vector();
          //cC_b_pR_bufs.~vector();
        }

        // sum all contributions
        msg->sum(K, nbraket);
      }

      if (cC->data_persistent()) cC->deactivate();
      if (C->data_persistent()) C->deactivate();

      // remove P-depenedent entries from registries
      psqrtregistry->remove(P);
      oreg->remove(Sspace->id());
      df_rtime->remove_if(Sspace->id());
      int3c_rtime->remove_if(Sspace->id());

      Ref<Integral> localints = int3c_rtime->factory()->integral()->clone();
      localints->set_basis(brabs);
      Ref<PetiteList> brapl = localints->petite_list();
      localints->set_basis(ketbs);
      Ref<PetiteList> ketpl = localints->petite_list();
      RefSCDimension bradim = brapl->AO_basisdim();
      RefSCDimension ketdim = ketpl->AO_basisdim();
      RefSCMatrix result(bradim,
                         ketdim,
                         brabs->so_matrixkit());

      result.assign(K);

      tim.exit();
      //ExEnv::out0() << decindent;
      //ExEnv::out0() << indent << "Exited exchange(DF) matrix evaluator ("
      //    << (tim.wall_time("exchange(DF)") - wall_time_start) << " sec)" << endl;

      return result;
    }

    RefSCMatrix exchange_df_local(const Ref<DensityFittingInfo>& df_info,
                            const RefSymmSCMatrix& P,
                            SpinCase1 spin,
                            const Ref<GaussianBasisSet>& brabs,
                            const Ref<GaussianBasisSet>& ketbs,
                            const Ref<GaussianBasisSet>& obs,
                            const Ref<FockBuildRuntime::PSqrtRegistry>& psqrtregistry) {
      /*=======================================================================================*/
      /* Setup and stuff                                       		                        {{{1 */ #if 1 // begin fold
      Timer tim("exchange(DF local)");
      timer_enter("01 - setup", 1);
      //----------------------------------------//
      // Only closed shell is implemented so far...
      if(spin != AnySpinCase1)
        throw FeatureNotImplemented("open shell local DF exchange", __FILE__, __LINE__);
      //----------------------------------------//
      Ref<MessageGrp> msg = MessageGrp::get_default_messagegrp();
      int me = msg->me();
      int nproc = msg->n();
      //ExEnv::out0() << endl << indent
      //         << "Entered Coulomb(DF) matrix evaluator" << endl;
      //ExEnv::out0() << incindent;
      //----------------------------------------//
      const Ref<DensityFittingRuntime>& df_rtime = df_info->runtime();
      //----------------------------------------//
      const Ref<AOSpaceRegistry> ao_registry = df_rtime->moints_runtime()->factory()->ao_registry();
      const Ref<OrbitalSpace>& braspace = ao_registry->value(brabs);
      const Ref<OrbitalSpace>& ketspace = ao_registry->value(ketbs);
      const Ref<GaussianBasisSet>& dfbs = df_info->params()->basis();
      assert(dfbs.nonnull());
      const Ref<OrbitalSpace>& dfspace = ao_registry->value(dfbs);
      const Ref<OrbitalSpace>& obs_space = ao_registry->value(obs);
      //----------------------------------------//
      const blasint dfnbf = dfspace->rank();
      int branbf = brabs->nbasis();
      int ketnbf = ketbs->nbasis();
      int obsnbf = obs->nbasis();
      //----------------------------------------//
      const Ref<TwoBodyThreeCenterMOIntsRuntime>& int3c_rtime = df_rtime->moints_runtime()->runtime_3c();
      const Ref<TwoBodyTwoCenterMOIntsRuntime>& int2c_rtime = df_rtime->moints_runtime()->runtime_2c();
      //----------------------------------------//
      double *P_ptr = allocate<double>(branbf*ketnbf);
      {
        RefSymmSCMatrix Ptmp = P.copy();
        Ptmp.convert2RefSCMat().convert(P_ptr);
      }
      EigenMatrixMap D(P_ptr, branbf, ketnbf);
      //----------------------------------------//
      /*****************************************************************************************/ #endif //1}}}
      /*=======================================================================================*/
      /* Get the coefficients                                  		                        {{{1 */ #if 1 // begin fold
      timer_change("02 - get coefficients", 1);
      //----------------------------------------//
      // Setup
      timer_enter("01 - (mu nu | M | X) setup", 2);
      std::string metric_key = df_info->params()->kernel_key();
      const bool noncoulomb_kernel = (metric_key.find("exp") != std::string::npos);
      std::string params_key = df_info->params()->intparams_key();
      std::string operset_key = noncoulomb_kernel ? "G12'" : "ERI";
      TwoBodyOper::type metric_oper =
          noncoulomb_kernel ? TwoBodyOper::r12_0_g12 : TwoBodyOper::eri;
      unsigned int ints_type_idx = TwoBodyOperSetDescr::instance(
          noncoulomb_kernel ? TwoBodyOperSet::G12NC : TwoBodyOperSet::ERI
      )->opertype(metric_oper);
      //----------------------------------------//
      // Go ahead and compute them all for now.  It will be faster to
      //   compute only the ones we need and save them for later, but
      //   that requires significant infrastructure changes.
      const std::string munu_M_X_key = ParsedTwoBodyThreeCenterIntKey::key(obs_space->id(),
                                                                           dfspace->id(),
                                                                           obs_space->id(),
                                                                           operset_key, params_key);
      timer_change("02 - (mu nu | M | X) compute", 2);
      const Ref<TwoBodyThreeCenterMOIntsTransform>& munu_M_X_tform = int3c_rtime->get(munu_M_X_key);
      munu_M_X_tform->compute();
      Ref<DistArray4> munu_M_X = munu_M_X_tform->ints_acc(); munu_M_X->activate();
      //----------------------------------------//
      // < X Y | metric | >
      // in chemists notation: ( X | metric | Y )
      timer_change("03 - compute (X | M | Y)", 2);
      const std::string metric2c_key = ParsedTwoBodyTwoCenterIntKey::key(
          dfspace->id(),
          dfspace->id(),
          operset_key, params_key
      );
      RefSCMatrix metric_2c_ints = int2c_rtime->get(metric2c_key);
      // hopefully this doesn't do a copy?
      RefSCMatrix coulomb_2c_ints = metric_2c_ints.pointer();
      timer_change("04 - compute (X | g | Y)", 2);
      if(noncoulomb_kernel){
        const std::string coulomb2c_key = ParsedTwoBodyTwoCenterIntKey::key(
            dfspace->id(),
            dfspace->id(),
            "ERI", ""
        );
        coulomb_2c_ints = int2c_rtime->get(metric2c_key);
      }
      timer_change("05 - transfer (X | g | Y)", 2);
      double* coulomb_2c_ints_ptr = allocate<double>(dfnbf*dfnbf);
      coulomb_2c_ints.convert(coulomb_2c_ints_ptr);
      EigenMatrixMap X_g_Y(coulomb_2c_ints_ptr, dfnbf, dfnbf);
      //----------------------------------------//
      // loop over atom pairs...
      timer_change("06 - loop over basis functions", 2);
      timer_enter("01 - setup", 3)
      EigenMatrix d_mu_siX(branbf, obsnbf*dfnbf);
      EigenMatrix gt_nu_siX(branbf, obsnbf*dfnbf);
      d_mu_siX = EigenMatrix::Zero(branbf, obsnbf*dfnbf);
      gt_nu_siX = EigenMatrix::Zero(branbf, obsnbf*dfnbf);
      if(munu_M_X->has_access(me)){
        DecompositionMap decomps;
        for(int ibf = 0; ibf < branbf; ++ibf){
          if(not munu_M_X->is_local(0, ibf))
            continue;
          //----------------------------------------//
          // Convenience variables for atom A
          const int ishA = brabs->function_to_shell(ibf);
          const int atomA = brabs->shell_to_center(ishA);
          const int nshA = brabs->nshell_on_center(atomA);
          const int nbfA = brabs->nbasis_on_center(atomA);
          const int shoffA = brabs->shell_on_center(atomA, 0);
          const int bfoffA = brabs->shell_to_function(shoffA);
          const int dfnshA = dfbs->nshell_on_center(atomA);
          const int dfnbfA = dfbs->nbasis_on_center(atomA);
          const int dfshoffA = dfbs->shell_on_center(atomA, 0);
          const int dfbfoffA = dfbs->shell_to_function(dfshoffA);
          //----------------------------------------//
          // TODO Permutational symmetry
          for(int jbf = 0; jbf < ketnbf; ++jbf){
            //----------------------------------------//
            // Convenience variables for atom B
            const int jshB = ketbs->function_to_shell(jbf);
            const int atomB = ketbs->shell_to_center(jshB);
            const int nshB = ketbs->nshell_on_center(atomB);
            const int nbfB = ketbs->nbasis_on_center(atomB);
            const int shoffB = ketbs->shell_on_center(atomB, 0);
            const int bfoffB = ketbs->shell_to_function(shoffB);
            const int dfnshB = dfbs->nshell_on_center(atomB);
            const int dfnbfB = dfbs->nbasis_on_center(atomB);
            const int dfshoffB = dfbs->shell_on_center(atomB, 0);
            const int dfbfoffB = dfbs->shell_to_function(dfshoffB);
            int dfpart_size = dfnbfA;
            if(atomA != atomB) dfpart_size += dfnbfB;
            IntPair atomAB(atomA, atomB);
            /*---------------------------------------------------------------*/
            /*  Get and Decompose ( X_(ab) | M | Y_(ab) )               {{{2 */ #if 2 // begin fold
            //===============================================================//
            if(decomps.count(atomAB) == 0) {
              timer_change("02 - get ( X_(ab) | M | Y_(ab) ) block", 3);
              double* metric_2c_part_ptr = allocate<double>(dfpart_size*dfpart_size);
              double** blockptr = allocate<double*>(dfpart_size);
              //---------------------------------------------------------------//
              // Get the block corresponding to (A|m|A)
              for(int idfpart = 0; idfpart < dfnbfA; ++idfpart) {
                blockptr[idfpart] = &(metric_2c_part_ptr[idfpart * dfpart_size]);
              }
              {
                RefSCMatrix block = metric_2c_ints.get_subblock(
                    dfbfoffA, dfbfoffA+dfnbfA-1,
                    dfbfoffA, dfbfoffA+dfnbfA-1
                );
                block.convert(blockptr);
              }
              // if atomA == atomB, we're done.  Otherwise,
              if(atomA != atomB){
                // Get the block corresponding to (A|m|B)
                for(int idfpart = 0; idfpart < dfnbfA; ++idfpart) {
                  blockptr[idfpart] = &(metric_2c_part_ptr[idfpart * dfpart_size + dfnbfA]);
                }
                {
                  RefSCMatrix block = metric_2c_ints.get_subblock(
                      dfbfoffA, dfbfoffA+dfnbfA-1,
                      dfbfoffB, dfbfoffB+dfnbfB-1
                  );
                  block.convert(blockptr);
                }
                //---------------------------------------------------------------//
                // Get the block corresponding to (B|m|A)
                for(int idfpart = 0; idfpart < dfnbfB; ++idfpart) {
                  blockptr[idfpart] = &(metric_2c_part_ptr[(idfpart+dfnbfA) * dfpart_size]);
                }
                {
                  RefSCMatrix block = metric_2c_ints.get_subblock(
                      dfbfoffB, dfbfoffB+dfnbfB-1,
                      dfbfoffA, dfbfoffA+dfnbfA-1
                  );
                  block.convert(blockptr);
                }
                //---------------------------------------------------------------//
                // Get the block corresponding to (B|m|B)
                for(int idfpart = 0; idfpart < dfnbfB; ++idfpart) {
                  blockptr[idfpart] = &(metric_2c_part_ptr[(idfpart+dfnbfA) * dfpart_size + dfnbfA]);
                }
                {
                  RefSCMatrix block = metric_2c_ints.get_subblock(
                      dfbfoffB, dfbfoffB+dfnbfB-1,
                      dfbfoffB, dfbfoffB+dfnbfB-1
                  );
                  block.convert(blockptr);
                }
              }
              //-----------------------------------------------------------------//
              EigenMatrixMap X_m_Y(metric_2c_part_ptr, dfpart_size, dfpart_size);
              //for(int X = 0; X < dfpart_size; ++X){
              //  for(int Y = 0; Y < dfpart_size; ++Y){
              //    if(fabs(double(X_m_Y(X,Y))) < 1e-12){
              //      X_m_Y(X,Y) = 0;
              //    }
              //  }
              //}
              timer_change("03 - decompose ( X(ab) | M | Y(ab) )", 3);
              decomps[atomAB] = new Decomposition(X_m_Y);
              /*---------------------------------------------------------------*/
              timer_change("misc", 3);
              deallocate(blockptr);
              deallocate(metric_2c_part_ptr);
              /*---------------------------------------------------------------*/
            }
            /*****************************************************************/ #endif //2}}}
            /*---------------------------------------------------------------*/
            timer_change("04 - solve equations", 3);
            //----------------------------------------//
            // set up the ij_M_X(AB) vector
            timer_enter("01 - get ibfjbf_M_X", 4);
            Eigen::VectorXd ibfjbf_M_X(dfpart_size);
            //----------------------------------------//
            // get ij_M_X(A)
            double* ibfjbf_M_X_buffer_A = allocate<double>(dfnbfA);
            VectorMap ibfjbf_M_X_A(ibfjbf_M_X_buffer_A, dfnbfA);
            munu_M_X->retrieve_pair_subblock(
                0, ibf,     // index in unit basis, index in mu
                ints_type_idx,
                jbf,        jbf+1,             // nu start, nu fence
                dfbfoffA,   dfbfoffA + dfnbfA, // X start,  X fence
                ibfjbf_M_X_buffer_A
            );
            ibfjbf_M_X.head(dfnbfA) = ibfjbf_M_X_A;
            deallocate(ibfjbf_M_X_buffer_A);
            //----------------------------------------//
            // get ij_M_X(B)
            if(atomA != atomB){
              double* ibfjbf_M_X_buffer_B = allocate<double>(dfnbfB);
              VectorMap ibfjbf_M_X_B(ibfjbf_M_X_buffer_B, dfnbfB);
              munu_M_X->retrieve_pair_subblock(
                  0, ibf,     // index in unit basis, index in mu
                  ints_type_idx,
                  jbf,        jbf+1,             // nu start, nu fence
                  dfbfoffB,   dfbfoffB + dfnbfB, // X start,  X fence
                  ibfjbf_M_X_buffer_B
              );
              ibfjbf_M_X.tail(dfnbfB) = ibfjbf_M_X_B;
              deallocate(ibfjbf_M_X_buffer_B);
            }
            //----------------------------------------//
            timer_change("02 - get Cpart", 4);
            // Now solve A * x = b, with A = ( X_(ab) | M | Y_(ab) )
            //   and b = ( ibf jbf | M | Y_(ab) )
            // The result will be dfpart_size rows and 1 column
            Eigen::VectorXd Cpart(dfpart_size);
            Cpart = decomps[atomAB]->solve(ibfjbf_M_X);
            //----------------------------------------//
            // contract with X_g_Y to form gt_mu_siX
            timer_change("03 - get gt_mu_siX", 4);
            gt_nu_siX.block(
                ibf, jbf*dfnbf,
                1, dfnbf
            ) += Cpart.head(dfnbfA).transpose() * X_g_Y.middleRows(dfbfoffA, dfnbfA);
            if(atomA != atomB){
              gt_nu_siX.block(
                  ibf, jbf*dfnbf,
                  1, dfnbf
              ) += Cpart.tail(dfnbfB).transpose() * X_g_Y.middleRows(dfbfoffB, dfnbfB);
            }
            //----------------------------------------//
            // contract with density and put it into d_mu_siX
            timer_change("04 - get d_mu_siX", 4);
            for(int sigma = 0; sigma < obsnbf; ++sigma){
              d_mu_siX.block(
                  ibf, sigma*dfnbf + dfbfoffA,
                  1, dfnbfA
              ) += D(jbf, sigma) * Cpart.head(dfnbfA).transpose();
              if(atomA != atomB){
                d_mu_siX.block(
                    ibf, sigma*dfnbf + dfbfoffB,
                    1, dfnbfB
                ) += D(jbf, sigma) * Cpart.tail(dfnbfB).transpose();
              }
            }
            //----------------------------------------//
            timer_exit(4);
            timer_change("01 - setup", 3);
          } // end loop over basis functions jbf
        } // end loop over basis functions ibf
        for(DMap_iter it=decomps.begin(); it != decomps.end(); ++it){
          delete it->second;
        }
      } // end if has access(me)
      timer_exit(3);
      timer_change("07 - global sum C", 2);
      msg->sum(d_mu_siX.data(), branbf * obsnbf * dfnbf);
      msg->sum(gt_nu_siX.data(), branbf * obsnbf * dfnbf);
      timer_exit(2);
      /*****************************************************************************************/ #endif //1}}}
      /*=======================================================================================*/
      /* Get ( mu nu | g | X )                                		                        {{{1 */ #if 1 // begin fold
      //----------------------------------------//
      // < mu X | coulomb | nu >
      // in chemists notation: ( mu nu | coulomb | X )
      timer_change("03 - get (mu nu | g | X)", 1);
      timer_enter("01 - setup", 2);
      EigenMatrix K_tilde(branbf, ketnbf);
      K_tilde = EigenMatrix::Zero(branbf, ketnbf);
      //EigenMatrix munu_g_X(branbf*ketnbf, dfnbf);
      unsigned int g_type_idx = TwoBodyOperSetDescr::instance(TwoBodyOperSet::ERI)->opertype(metric_oper);
      timer_change("02 - compute", 2);
      Ref<DistArray4> munu_g_X_D4 = munu_M_X;
      // Only need to recompute if we're using a non-coulomb kernel
      if(noncoulomb_kernel){
        if (munu_M_X->data_persistent()) munu_M_X->deactivate();
        munu_M_X = 0;
        const std::string munu_g_X_key = ParsedTwoBodyThreeCenterIntKey::key(
            braspace->id(),
            dfspace->id(),
            ketspace->id(),
            "ERI", ""
        );
        const Ref<TwoBodyThreeCenterMOIntsTransform>& munu_g_X_tform = int3c_rtime->get(munu_g_X_key);
        munu_g_X_tform->compute();
        munu_g_X_D4 = munu_g_X_tform->ints_acc(); munu_g_X_D4->activate();
      }
      timer_change("03 - compute K_tilde", 2);
      timer_enter("misc", 3);
      if(munu_g_X_D4->has_access(me)){
        double* ibfjbf_g_X_buffer = allocate<double>(dfnbf);
        VectorMap ibfjbf_g_X(ibfjbf_g_X_buffer, dfnbf);
        for(int ibf = 0; ibf < branbf; ++ibf){
          if(not munu_g_X_D4->is_local(0, ibf))
            continue;
          for(int jbf = 0; jbf < ketnbf; ++jbf){
            timer_change("01 - get pair subblock", 3);
            munu_g_X_D4->retrieve_pair_subblock(
                0, ibf,      // unit basis index, mu_index
                g_type_idx,
                jbf, jbf+1,  // nu_start, nu_fence
                0,   dfnbf,  // X_start, X_fence
                ibfjbf_g_X_buffer
            );
            timer_change("02 - get K_tilde contribution", 3);
            K_tilde.row(ibf) += d_mu_siX.middleCols(jbf*dfnbf, dfnbf) * ibfjbf_g_X;
            timer_change("misc", 3);
          }
        }
        deallocate(ibfjbf_g_X_buffer);
      }
      timer_change("03 - global sum K_tilde", 3);
      msg->sum(K_tilde.data(), branbf*ketnbf);
      timer_change("misc", 3);
      if(munu_g_X_D4->data_persistent()) munu_g_X_D4->deactivate();
      munu_g_X_D4 = 0;
      munu_M_X = 0;
      timer_exit(3);
      timer_exit(2);
      /*****************************************************************************************/ #endif //1}}}
      /*=======================================================================================*/
      /* Compute K                                             		                        {{{1 */ #if 1 // begin fold
      timer_change("04 - compute K", 1);
      timer_enter("01 - two body part", 2);
      K_tilde -= 0.5 * (d_mu_siX * gt_nu_siX.transpose());
      timer_change("02 - symmetrize", 2);
      EigenMatrix K(branbf, ketnbf);
      K = K_tilde + K_tilde.transpose();
      //----------------------------------------//
      // Transfer K to a RefSCMatrix
      timer_change("03 - transfer to RefSCMatrix", 2);
      Ref<Integral> localints = int3c_rtime->factory()->integral()->clone();
      localints->set_basis(brabs);
      Ref<PetiteList> brapl = localints->petite_list();
      localints->set_basis(ketbs);
      Ref<PetiteList> ketpl = localints->petite_list();
      RefSCDimension bradim = brapl->AO_basisdim();
      RefSCDimension ketdim = ketpl->AO_basisdim();
      RefSCMatrix result(
          bradim,
          ketdim,
          brabs->so_matrixkit()
      );
      result.assign(K.data());
      timer_exit(2);
      /*****************************************************************************************/ #endif //1}}}
      /*=======================================================================================*/
      /* Cleanup stuff                                         		                        {{{1 */ #if 1 // begin fold
      timer_change("05 - cleanup", 1);
      deallocate(P_ptr);
      deallocate(coulomb_2c_ints_ptr);
      timer_exit(1);
      tim.exit();
      return result;
      /*****************************************************************************************/ #endif //1}}}
      /*=======================================================================================*/
    }

}} // namespace sc::detail

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
