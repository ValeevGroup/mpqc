//
// eigen.h
//
// Copyright (C) 2013 MPQC Developers
//
// Author: David Hollman <dhollman@vt.edu>
// Maintainer: DSH, EV
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

/*******************************************************************************
 * This file contains various methods and functions for interfacing with the
 * Eigen C++ library.
 *******************************************************************************/

#ifndef _math_mmisc_eigen_h
#define _math_mmisc_eigen_h

#include <Eigen/Dense>
#include <util/misc/sharedptr.h>

namespace sc {

    class StateOut;
    class StateIn;

    void ToStateOut(const Eigen::MatrixXd &m, StateOut &so, int &count);
    void ToStateOut(const Eigen::VectorXd &m, StateOut &so, int &count);
    void ToStateOut(const Eigen::RowVectorXd &m, StateOut &so, int &count);

    void FromStateIn(Eigen::MatrixXd &m, StateIn &si, int &count);
    void FromStateIn(Eigen::VectorXd &m, StateIn &si, int &count);
    void FromStateIn(Eigen::RowVectorXd &m, StateIn &si, int &count);

    template <typename T>
    void ToStateOut(const std::shared_ptr<T>& val, StateOut &so, int &count){
      ToStateOut(*val, so, count);
    }

    template <typename T>
    void FromStateIn(std::shared_ptr<T>& val, StateIn &si, int &count){
      FromStateIn(*val, si, count);
    }

}

#endif
