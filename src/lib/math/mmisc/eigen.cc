//
// eigen.cc
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

#include <math/mmisc/eigen.h>
#include <util/state/stateout.h>
#include <util/state/statein.h>

void sc::ToStateOut(const Eigen::MatrixXd &m, StateOut &so, int &count)
{
  count += so.put(m.rows());
  count += so.put(m.cols());
  count += so.put_array_double(m.data(), m.rows() * m.cols());
}

void sc::ToStateOut(const Eigen::VectorXd &m, StateOut &so, int &count)
{
  count += so.put(m.rows());
  count += so.put_array_double(m.data(), m.rows());
}

void sc::ToStateOut(const Eigen::RowVectorXd &m, StateOut &so, int &count)
{
  count += so.put(m.rows());
  count += so.put_array_double(m.data(), m.rows());
}


void sc::FromStateIn(Eigen::MatrixXd &m, StateIn &si, int &count){
  int rows = 0;
  int cols = 0;
  count += si.get(rows);
  count += si.get(cols);
  m.resize(rows, cols);
  count += si.get_array_double(m.data(), rows*cols);
}

void sc::FromStateIn(Eigen::VectorXd &m, StateIn &si, int &count){
  int rows = 0;
  count += si.get(rows);
  m.resize(rows);
  count += si.get_array_double(m.data(), rows);
}

void sc::FromStateIn(Eigen::RowVectorXd &m, StateIn &si, int &count){
  int rows = 0;
  count += si.get(rows);
  m.resize(rows);
  count += si.get_array_double(m.data(), rows);
}
