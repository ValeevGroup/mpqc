//
// twobodygrid.cc
//
// Copyright (C) 2004 Edward Valeev
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

#include <cmath>
#include <stdexcept>
#include <util/misc/formio.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/mbptr12/twobodygrid.h>

using namespace std;
using namespace sc;

/*---------------
  TwoBodyGrid
 ---------------*/
static ClassDesc TwoBodyGrid_cd(
  typeid(TwoBodyGrid),"TwoBodyGrid",1,"virtual public SavableState",
  0, create<TwoBodyGrid>, create<TwoBodyGrid>);

TwoBodyGrid::TwoBodyGrid(StateIn& si) : SavableState(si), O_(0.0)
{
  si.get(name_);
  int npts; si.get(npts); r1_.resize(npts); r2_.resize(npts);
  for(int pt=0; pt<npts; pt++) {
    si.get(r1_[pt].x());
    si.get(r1_[pt].y());
    si.get(r1_[pt].z());
  }
  for(int pt=0; pt<npts; pt++) {
    si.get(r2_[pt].x());
    si.get(r2_[pt].y());
    si.get(r2_[pt].z());
  }
  si.get(O_.x());
  si.get(O_.y());
  si.get(O_.z());
}

TwoBodyGrid::TwoBodyGrid(const Ref<KeyVal>& keyval)
{
  name_ = keyval->stringvalue("name",KeyValValuestring("two-body grid"));

  // Default is to assume Cartesian coordinates
  bool polar = keyval->booleanvalue("polar",KeyValValueboolean((int)false));
  // 2D grid is default
  unsigned int ndim = keyval->intvalue("ndim",KeyValValueint(2));

  bool O_is_given = keyval->exists("origin");
  if (O_is_given) {
    const int dim = keyval->count("origin");
    if (dim != 3)
      throw std::runtime_error("TwoBodyGrid::TwoBodyGrid() -- keyword origin must be an array of 3 elements");
    for(int xyz=0; xyz<3; xyz++)
      O_.elem(xyz) = keyval->doublevalue("origin",xyz);
  }
  else
    O_ = 0.0;

  bool r1_is_given = keyval->exists("r1");
  bool r2_is_given = keyval->exists("r2");
  if (r1_is_given == false)
    throw std::runtime_error("TwoBodyGrid::TwoBodyGrid() -- keyword r1 must be given");
  if (r2_is_given == false)
    throw std::runtime_error("TwoBodyGrid::TwoBodyGrid() -- keyword r2 must be given");

  const int nelem1 = keyval->count("r1");
  const int nelem2 = keyval->count("r2");
  if (nelem1 == 0)
    throw std::runtime_error("TwoBodyGrid::TwoBodyGrid() -- keyword r1 must be an array of 3-dimensional vectors");
  if (nelem2 == 0)
    throw std::runtime_error("TwoBodyGrid::TwoBodyGrid() -- keyword r2 must be an array of 3-dimensional vectors");
  if (ndim == 1 && nelem1 != nelem2)
    throw InputError("TwoBodyGrid::TwoBodyGrid -- ndim is 1 but number of elements in r1 and r2 do not match",__FILE__,__LINE__);

  const int nelem = (ndim == 2) ? nelem1 * nelem2 : nelem1;
  
  r1_.resize(nelem);
  r2_.resize(nelem);

  std::vector<SCVector3> r1(nelem1);
  std::vector<SCVector3> r2(nelem2);
  for(int i=0; i<nelem1; i++) {
    SCVector3 R1;
    const int dim = keyval->count("r1",i);
    if (dim != 3) {
      std::ostringstream oss;
      oss << "value for keyword r1:" << i << " must be an array of 3 elements";
      throw InputError(oss.str().c_str(),__FILE__,__LINE__);
    }
    for(int xyz=0; xyz<3; xyz++)
      R1.elem(xyz) = keyval->doublevalue("r1",i,xyz);
    if (polar) {
      R1.spherical_to_cartesian(r1[i]);
      r1[i] += O_;
    }
    else
      r1[i] = R1 + O_;
  }
  for(int i=0; i<nelem2; i++) {
    SCVector3 R2;
    const int dim = keyval->count("r2",i);
    if (dim != 3) {
      std::ostringstream oss;
      oss << "value for keyword r2:" << i << " must be an array of 3 elements";
      throw InputError(oss.str().c_str(),__FILE__,__LINE__);
    }
    for(int xyz=0; xyz<3; xyz++)
      R2.elem(xyz) = keyval->doublevalue("r2",i,xyz);
    if (polar) {
      R2.spherical_to_cartesian(r2[i]);
      r2[i] += O_;
    }
    else
      r2[i] = R2 + O_;
  }

  if (ndim == 1) {
    for(int i=0; i<nelem1; i++) {
      r1_[i] = r1[i];
      r2_[i] = r2[i];
    }
  }
  else {
    int ij = 0;
    for(int i=0; i<nelem1; i++)
      for(int j=0; j<nelem2; j++, ++ij) {
	r1_[ij] = r1[i];
	r2_[ij] = r2[j];
      }
  }
}

TwoBodyGrid::~TwoBodyGrid()
{
}

void
TwoBodyGrid::save_data_state(StateOut& so)
{
  so.put(name_);
  const int npts = r1_.size();
  so.put(npts);
  for(int pt=0; pt<npts; pt++) {
    so.put(r1_[pt].x());
    so.put(r1_[pt].y());
    so.put(r1_[pt].z());
  }
  for(int pt=0; pt<npts; pt++) {
    so.put(r2_[pt].x());
    so.put(r2_[pt].y());
    so.put(r2_[pt].z());
  }
  so.put(O_.x());
  so.put(O_.y());
  so.put(O_.z());
}

const std::string&
TwoBodyGrid::name() const { return name_; }

int
TwoBodyGrid::nelem() const { return r1_.size(); }

const SCVector3&
TwoBodyGrid::origin() const { return O_; }


SCVector3
TwoBodyGrid::xyz1(int i, const SCVector3& O) const
{
  return r1_[i] - O;
}

SCVector3
TwoBodyGrid::xyz2(int i, const SCVector3& O) const
{
  return r2_[i] - O;
}

SCVector3
TwoBodyGrid::rtp1(int i, const SCVector3& O) const
{
  SCVector3 RO = r1_[i] - O;
  double r = RO.norm();
  double theta = acos(RO.z()/r);
  double phi = acos(RO.x()/(r*sin(theta)));
  return SCVector3(r,theta,phi);
}

SCVector3
TwoBodyGrid::rtp2(int i, const SCVector3& O) const
{
  SCVector3 RO = r2_[i] - O;
  double r = RO.norm();
  double theta = acos(RO.z()/r);
  double phi = acos(RO.x()/(r*sin(theta)));
  return SCVector3(r,theta,phi);
}

void
TwoBodyGrid::print(std::ostream& os) const
{
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
