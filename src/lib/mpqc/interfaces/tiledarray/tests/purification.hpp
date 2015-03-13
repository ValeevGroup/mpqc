//
// purification.hpp
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

#ifndef MPQC_INTERFACES_PURIFICATION_HPP
#define MPQC_INTERFACES_PURIFICATION_HPP

#include <mpqc/math/matrix.hpp>
#include <tiledarray.h>
#include <array>

namespace mpqc {
namespace purification {

    using TArray2 = TA::Array<double, 2>;
    // Holds data we want for each type of algorithm
    struct performance {
        std::size_t multiplies;
        std::size_t adds;
        std::size_t traces;
    };

    // Compute Eigenvalue Estimates first element is largest eigenvalue.
    std::array<double ,2> BoundsTA(const TArray2 &F){
        // make madness::Group all process that contain important row.
        // define local reduction reduces tiles into the vector
        // global reduction to finish the vector off.
        // or maybe an all reduce
        // Final

        double min = 10, max = -10;
        return {{max, min}};
    }

    std::array<double, 2> Bounds(const Matrix &F){
        auto n = F.rows();
        double min = 10;
        for(auto i = 0; i < n; ++i){
            double row_sum = 0;
            for(auto j = 0; j < n; ++j){
                row_sum += (j != i) ?  -std::abs(F(i,j)) : F(i,i);
            }
            min = (row_sum < min) ? row_sum : min;
        }

        double max = -10;
        for(auto i = 0; i < n; ++i){
            double row_sum = 0;
            for(auto j = 0; j < n; ++j){
                row_sum += (j != i) ? std::abs(F(i,j)) : F(i,i);
            }
            max = (row_sum > max) ? row_sum : max;
        }
        return {{max, min}};
    }

    // Performs Trace Resetting Purification and returns number of multiplies and additions
    performance tr_pure(const TArray2 &F, TArray2 &D, std::size_t occ){
        performance data;
        data.multiplies = 0;
        data.adds = 1;
        data.traces = 1; // Start at 1 since last while loop won't be counted
                         // and has one of each.
        auto E = Bounds(F);

        Matrix I = Eigen::MatrixXd::Identity(F.rows(), F.cols());

        D = (E[0] * I - F)/(E[0] - E[1]);
        data.adds++;

        Matrix IDD = (I - D) * D;
        data.adds++; data.multiplies++;

        while((IDD).lpNorm<Eigen::Infinity>() > 1e-6){ // While D^2 != D
            data.adds++;  // Update counters

            int gamma = (occ - D.trace() > 0) ? 1 : -1;
            data.traces++;

            D = D + gamma * IDD;
            data.adds++;
            D = 0.5 * Matrix(D + D.transpose());
            IDD = (I - D) * D;
            data.adds++; data.multiplies++;

        }
        return data;
    }

    Matrix Fpure(Matrix &D){
        Matrix I = Eigen::MatrixXd::Identity(D.rows(), D.cols());
        return Matrix(D * D * D * (4 * I - 3 * D));
    }

    Matrix Gpure(Matrix &D){
        Matrix I = Eigen::MatrixXd::Identity(D.rows(), D.cols());
        return Matrix( D * D * (I - D) * (I - D));
    }

    performance fourth_order_pure(const Matrix &F, Matrix &D, std::size_t occ){
        performance data;
        data.multiplies = 0;
        data.adds = 1;
        data.traces = 1; // Start at 1 since last while loop won't be counted
                         // and has one of each.
        auto E = Bounds(F);

        Matrix I = Eigen::MatrixXd::Identity(F.rows(), F.cols());

        D = (E[0] * I - F)/(E[0] - E[1]);
        data.adds++;

        Matrix D2 = D*D;
        data.multiplies++;

        while((D - D2).lpNorm<Eigen::Infinity>() > 1e-9){ // While D^2 != D
            data.adds++; // Update counters

            Matrix D3 = D2 * D;
            Matrix ID = I - D;
            data.adds++; data.multiplies++;

            double gamma = (double(occ) - Matrix(D3 * (4 * I - 3 * D)).trace())/
                            (Matrix(D2 * ID * ID).trace());
            data.adds += 1; data.multiplies += 3; data.traces += 2;
            if(gamma >= 6) {
                D = 2 * D - D2;
                data.adds++;
            }
            else if (gamma < 0){
                D = D2;
            }
            else {
                D = D3 * (4 * I - 3 * D) + gamma * (D2 * ID * ID);
                data.adds += 4; data.multiplies += 3;
            }
            D2 = D * D;
            data.multiplies++;
        }

        return data;
    }

    performance canonical_pure(const Matrix &F, Matrix &D, std::size_t occ){
        performance data;
        data.multiplies = 1;
        data.adds = 1;
        data.traces = 1; // Start at 1 since last while loop won't be counted
                         // and has one of each.
        auto n = F.rows();
        auto E = Bounds(F);

        Matrix I = Eigen::MatrixXd::Identity(n, n);
        double mu = F.trace()/n;
        data.traces++;
        double lambda = std::min(double(occ)/(E[0]-mu),
                                 double(n-occ)/(mu - E[1]));
        D = (lambda/double(n)) * (mu * I - F) + double(occ)/double(n) * I;
        data.adds += 2;

        Matrix D2 = D * D;
        data.multiplies++;

        while((D - D2).lpNorm<Eigen::Infinity>() > 1e-6
                        ){ // While D^2 != D
            data.adds++; // Update counters

            Matrix D3 = D2 * D;
            data.multiplies++;

            double c = (Matrix(D2 - D3).trace())/
                            (Matrix(D - D2).trace());
            data.adds += 2; data.traces +=2;
            if(c >= 0.5){
                D = (1.0/c) * ((1 + c) * D2 - D3);
                data.adds += 1;
            } else {
                D = 1.0/(1.0 - c) * ((1 - 2*c) * D + (1 + c)*D2 - D3);
                data.adds += 2;
            }

            D2 = D * D;
            data.multiplies++;
        }

        return data;
    }
} // namespace purification
} // namespace mpqc






#endif /* MPQC_INTERFACES_PURIFICATION_HPP */
