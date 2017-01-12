/*
 * curvy_steps.hpp
 *
 *  Created on: Aug 14, 2013
 *      Author: drewlewis
 */

#ifndef CURVY_STEPS_HPP_
#define CURVY_STEPS_HPP_

#include "common.hpp"
#include "tiledarray_fock.hpp"

namespace mpqc {
namespace tests {
    using Array2 = TA::Array<double,2>;
    using Array3 = TA::Array<double,3>;
    using Array4 = TA::Array<double,4>;

namespace curvy_steps {

    void
    Purify(Array2 &R, const Array2 &S){

        // Purify an approximate density matrix so that it meets the requirements
        // of a density matirx mainly that
        // R = R^T
        // Trace(S.R) = Number of electrons * 0.5
        // All eigenvalues of the matrix are 1.
        //lets force symmetry via averaging.  R_{ij} = (R_{ij} + R_{ji}) / 2
        R("i,j") = (R("i,j") + R("j,i")) * 0.5;

        // This routine is essentially R = 3*R^2 - 2 * R^3 it only converges if R is
        // close enough to the correct answer.
        double rdiff;
        do {
            const Array2 RS = R("i,k") * S("k,j");

            /// Peform the actual iteration step.
            const Array2 RSR = RS("i,k") * R("k,j");
            const Array2 Rnew = 3 * RSR("i,j") - 2 * RS("i,m") * RSR("m,j");
            const Array2 Rnorm_temp = Rnew("i,j") - R("i,j");

            rdiff = TA::expressions::norminf(Rnorm_temp("i,j"));
            R = Rnew;
        } while(rdiff > 1e-8);
    }

    Array2
    Scom(const Array2 &R, const Array2 &X, const Array2 &S){

        //Compute the S commutator of 2 matrices [R,X]_s = R.S.X - X.S.R
        Array2 RSX = R("i,k") * S("k,n") * X("n,j");

        // RSX_{ij} = - RSX_{ji}
        return RSX("i,j") + RSX("j,i");
    }

    void
    RotateR(Array2 &R, const Array2 &X, const Array2 &S,
            const size_t order = 4){

        // This routine approximates the rotation of the density matrix which
        // looks like R(X) = e^{-XS} R e^{SX} via the BCH approximation
        // R(X) = R + [R,X]_s + 1/(2!) [[R,X]_s,X]_s + 1/(3!)[[[R,X]_s,X]_s,X]_s . . .
        Array2 LHS = R("i,j");

        //Avoid computing factorial pieces
        double Fac[] = {1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800,
                        39916800};

        //Compute the first S comumator because we always want at least 1st order
        Array2 steps = Scom(R, X, S);
        for(std::size_t i = 0; i < order; ++i){
            //add the next order approximation to R and then compute the one for the
            //next round.  This will always compute one more step than asked for right
            //now, but it only returns the number of steps asked for.
            LHS("i,j") = LHS("i,j") + 1/Fac[i] * steps("i,j");
            steps = Scom(steps, X, S);
        }
        R = LHS;
    }

    void
    conjgradmat(const Array2 &F, Array2 &X, const Array2 &D,
                const Array2 &Esymm, const Array2 &S,
                const Array2 &Eanti, const std::size_t max_iter = 100) {
        //Solve problem Eanti = F.X.D + D.X.F - 0.5( S.X.Esymm + Esymm.X.S)
        //Using a modified conjugate gradient method.

        Array2 FXD = X; // = F("i,n") * X("n,m") * D("m,j") but X = 0
        Array2 SXEs = X; // = S("i,n") * X("n,m") * Esymm("m,j") but X = 0
        Array2 RHS = X; // (FXD("i,j") - FXD("j,i"))

        // Gradient = -Residual
        Array2 R = Eanti; // (Eanti("i,j") - RHS("i,j")); RHS = 0

        // Search direction is the opposite of gradient
        Array2 P_i = R;

        //Scalar product of the residual
        double Rnorm2 = TA::expressions::norm2(R("i,j"));
        double Rsold = Rnorm2 * Rnorm2;


        Array2 d_old = P_i;

        std::size_t count = 0;
        do{
            F.world().gop.fence();
            Array2 FPD = F("i,n") * P_i("n,m") * D("m,j");
            Array2 SPEs = S("i,n") * P_i("n,m") * Esymm("m,j");

            //Precompute what traditioally is A * p_i in regular cg method
            Array2 Ap_i = (FPD("i,j") - FPD("j,i") -
            0.5 * (SPEs("i,j") - SPEs("j,i")) );

            //Forming alpha the optimized length to travel along p_i
            double alpha = Rsold / TA::expressions::dot(P_i("i,j"), Ap_i("i,j"));

            // that the new step
            Array2 d = alpha * P_i("i,j");

            //Updating the unknown
            X("i,j") = X("i,j") +  d("i,j");

            //Calculating the new residual
            R("i,j") = R("i,j") - alpha * Ap_i("i,j");
            double Rsnew = TA::expressions::norm2(R("i,j"));
            Rsnew = Rsnew * Rsnew;

            //Finding a new search direction
            P_i("i,j") = R("i,j") + (Rsnew/Rsold) * P_i("i,j");
            Rsold = Rsnew;

            ++count;
        }while( std::sqrt(Rsold) >= 1e-5 && count < max_iter);
            if(count == max_iter){
                std::cout << "Convergence Factor in Conjgradmat = "
                                << std::sqrt(Rsold) << std::endl;
                //conjgradmat(F, X, D, Esymm, S, Eanti);
                std::cout << "The maximum number of iterations was reached"
                          " in conjugate gradient" << std::endl;
        }
        F.world().gop.fence();
    }

    void Density_Update(Array2 &R, const Array2 &S,
       const Array2 &F, const std::size_t iter,
       const std::size_t rotation_order = 2){
      // Optimize the one electron density matrix to solve for the Hartree Fock
      // Energy using the second order approach to modifying R (The denisty matrix)
      // The equations being solved are
      // F.Rn.S - S.Rn.F = F.X.S.Rn.S + S.Rn.S.X.F - 0.5( S.X.Esymm + Esymm.X.S)
      // Esymm = F.Rn.S + S.Rn.F
      //
      // For more through explination see P.477 of Molecular Electronic-Structure Theory

      ///Precomputing contractions
      /// FRS is anti symmetric so FRS(i,j) = -FRS(j,i) = -SRF(i,j)
      Array2 FRS = F("i,n") * R("n,m") * S("m,j");

      //Forming the RHS Eanti
      Array2 Eanti = FRS("i,j") - FRS("j,i");
      //Forming Esymm from above
      Array2 Esymm = FRS("i,j") + FRS("j,i");

      //Pre computing S.R.S into SRS
      Array2 SRS = S("i,n") * R("n,m") * S("m,j");

      //Making the X matrix which is an anti-symmetric matrix that performs the
      //Rotations on R.
      Array2 XD(S.world(), S.trange());
      XD.set_all_local(0.0);

      //Solving for XD using a modified conjugate gradient method. `
      conjgradmat(F, XD, SRS, Esymm, S, Eanti, iter * 20);

      /// For now this scales XD so that is is not too large.  If XD gets too large
      /// Then the purification step is not able to bring R back to idompotency.
      double scale = 0.1/TA::expressions::norminf(XD("i,j"));
      XD("i,j") = XD("i,j") * std::min(1.0,scale);

      XD.world().gop.fence();
      //Rotating the Rm to Rn with X
      RotateR(R,XD,S,rotation_order);

      XD.world().gop.fence();
      //Purifying Rn so that it meets the criteria for a valid density matrix
      Purify(R, S);

      R.world().gop.fence();
    }
} // namespace curvy_steps
} // namespace tests
} // namespace mpqc


#endif /* CURVY_STEPS_HPP_ */
