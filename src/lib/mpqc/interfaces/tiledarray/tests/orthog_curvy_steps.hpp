//
// orthog_curvy_steps.hpp
//
// Copyright (C) 2013 Drew Lewis
//
// Authors: Drew Lewis
// Maintainer: Drew Lewis and Edward Valeev
//
// This file is part of the MPQC Toolkit.
//
// The MPQC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The MPQC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the MPQC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef MPQC_ORTHOG_CURVY_STEPS_HPP
#define MPQC_ORTHOG_CURVY_STEPS_HPP

#include <mpqc/math/matrix.hpp>
#include "purification.hpp"

namespace mpqc {
    void conjugategradient(const Matrix &F, Matrix &X, const Matrix &D,
                           const Matrix &Esymm, const Matrix &Eanti,
                           performance &data){
        Matrix FXD = X; // X guess is 0 so anything times X is 0.
        Matrix RHS = X;

        // Gradient
        Matrix R = Eanti;
        Matrix P_i = R;

        double Rnorm2 = R.lpNorm<2>();
        double Rsold = Rnorm2 * Rnorm2;

        Matrix d_old = P_i;

        int iter = 0;
        while(std::sqrt(Rsold) > 1e-6 && iter < 4){
            Matrix FPD = F * P_i * D; data.multiplies++;
            Matrix PEs = P_i * Esymm; data.multiplies++;

            Matrix Ap_i = (FPD - FPD.transpose()) - 0.5*(PEs - PEs.transpose());
            data.adds += 3;

            double alpha = Rsold / (Matrix(P_i * Ap_i).lpNorm<2>());
            // Not doing a multiply here since this could be optimized away.

            Matrix d = alpha * P_i;

            X = X + d; data.adds++;

            R = R - alpha * Ap_i; data.adds++;

            double Rsnew = R.lpNorm<2>();
            Rsnew = Rsnew * Rsnew;

            P_i = R + (Rsnew/Rsold) * P_i; data.adds++;

            Rsold = Rsnew;
            ++iter;
        }
        std::cout << "Finished CG" << std::endl;
    }

    Matrix commuter(const Matrix &D, const Matrix &X, performance &data){
        Matrix DX = D * X; data.multiplies++; data.adds++; // Adds is for the return statement
        return DX + DX.transpose();
    }

    void BCH_rotate(Matrix &D, const Matrix &X, std::size_t order,
                    performance &data){
       Matrix LHS = D;
       double Fac[] = {1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800,
                       39916800};

       Matrix steps = commuter(D,X,data);
       for(auto i = 0; i < order; ++i){
           LHS = LHS + 1/Fac[i] * steps; data.adds++;
           steps = commuter(steps, X, data);
       }
       D = LHS;
    }

    void Purify(Matrix &D, performance &data){
        std::cout << "Dpre = \n" << D << std::endl;
        Matrix D2 = D*D; data.multiplies++;
        while((D-D2).lpNorm<Eigen::Infinity>() > 1e-6){
            D = 3 * D2 - 2 * D2 * D; data.adds++; data.multiplies++;
            D2 = D*D;
        }
        std::cout << "Dpure = \n" << D << std::endl;
        std::cout << "Finish pure" << std::endl;
    }

    performance curvy_step(const Matrix &F, Matrix &D){
        performance data={0,0,0};

        Matrix FD = F * D; data.multiplies++;
        Matrix Eanti = FD - FD.transpose(); data.adds++;
        Matrix Esymm = FD + FD.transpose(); data.adds++;

        Matrix X = Eigen::MatrixXd::Zero(F.cols(), F.rows());
        conjugategradient(F,X,D,Esymm,Eanti,data);

        double scale = 0.001/X.lpNorm<Eigen::Infinity>();
        X = scale * X;

        std::size_t rotation_order = 2;
        BCH_rotate(D, X, rotation_order, data);

        Purify(D, data);

        return data;
    }
}




#endif /* MPQC_ORTHOG_CURVY_STEPS_HPP */
