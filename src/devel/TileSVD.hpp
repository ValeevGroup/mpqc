/*
 *  This file is a part of TiledArray.
 *  Copyright (C) 2014  Virginia Tech
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef TILESVD_HPP
#define TILESVD_HPP

#include<Eigen/Dense>
#include<iostream>

class SVDTile {

 public:
  // Build with an Eigen matrix
  explicit SVDTile(Eigen::MatrixXd& mat, double cut = 1e-6):
    svd_(mat, Eigen::ComputeThinU | Eigen::ComputeThinV), cut_(cut)
  {
    nvals = std::count_if(svd_.singularValues().data(),
                  svd_.singularValues().data()+svd_.singularValues().rows(),
              [=](double x){return x > cut;}
            );
  }

  // Get the V matrix
  Eigen::MatrixXd Vt(){
    //Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(nvals, svd_.matrixV().rows());
    //for(std::size_t i = 0; i < nvals; ++i){
    //  temp.row(i) = svd_.matrixV().transpose().row(i);
    //}
    return svd_.matrixV().transpose().topRows(nvals);
  }

  // Get the singular values
  Eigen::MatrixXd SVals(){
    //Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(nvals, 1);
    //for(std::size_t i = 0; i < nvals; ++i){
    //  temp(i) = svd_.singularValues()(i);
    //}
    return svd_.singularValues().topRows(nvals);
  }

  // Get the U matrix
  Eigen::MatrixXd U(){
    //Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(svd_.matrixU().rows(), nvals);
    //for(std::size_t i = 0; i < nvals; ++i){
    //  temp.col(i) = svd_.matrixU().col(i);
    //}
    return svd_.matrixU().leftCols(nvals);
  }

  Eigen::MatrixXd U_contract(){
    Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(U().rows(), U().cols());
    for(auto i = 0; i < U().cols(); ++i){
      temp.col(i) = SVals()(i) * U().col(i);
    }
    return temp;
  }

  Eigen::MatrixXd Vt_contract(){
    Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(Vt().rows(), Vt().cols());
    for(auto i = 0; i < Vt().rows(); ++i){
      temp.row(i) = SVals()(i) * Vt().row(i);
    }
    return temp;
  }

  Eigen::MatrixXd operator*(SVDTile &rhs){
    return U_contract() * Vt() * rhs.U() * rhs.Vt_contract();
  }

  std::size_t nelements(){
    return nvals * (U().rows() + Vt().cols() + 1);
  }

private:

  Eigen::JacobiSVD<Eigen::MatrixXd> svd_;
  double cut_;
  std::size_t nvals = 0;
};


#endif // TILESVD_HPP
