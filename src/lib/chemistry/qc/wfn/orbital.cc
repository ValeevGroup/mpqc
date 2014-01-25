//
// orbital.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
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
#include <cassert>

#include <util/misc/formio.h>

#include <math/scmat/local.h>
#include <math/scmat/vector3.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/molshape.h>
#include <chemistry/qc/wfn/orbital.h>
#include <util/misc/scexception.h>
#include <sstream>

using namespace sc;

namespace {
  /** a simple function calculating the number of digits in an integer in base 10
   */
  int IntDigitNum(int II)
  {
    int NumDigit = 0;
    while(II > 0)
    {
      NumDigit += 1;
      II /= 10;
    }
    return NumDigit;
  }
}

static ClassDesc Orbital_cd(
  typeid(Orbital),"Orbital",1,"public Volume",
  0, create<Orbital>, 0);

Orbital::Orbital(const Ref<KeyVal> &keyval):
  Volume(keyval)
{
  wfn_ << keyval->describedclassvalue("wfn");
  orbital_ = keyval->intvalue("orbital");
}

Orbital::Orbital(const Ref<OneBodyWavefunction>& wfn, int orbital):
  Volume(),
  wfn_(wfn),
  orbital_(orbital)
{
}

Orbital::~Orbital() {}

void
Orbital::compute()
{
  SCVector3 r;
  get_x(r);
  if (value_needed()) {
      set_value(fabs(wfn_->orbital(r, orbital_)));
      set_actual_value_accuracy(desired_value_accuracy());
    }
  if (gradient_needed() || hessian_needed()) {
      ExEnv::err0() << indent
           << "Orbital::compute(): gradient & hessian not implemented\n";
      abort();
    }
}

// make a wild guess about the bounding box
void
Orbital::boundingbox(double valuemin,
                     double valuemax,
                     SCVector3& p1, SCVector3& p2)
{
  Molecule& mol = *wfn_->molecule();

  if (mol.natom() == 0) {
      for (int i=0; i<3; i++) p1[i] = p2[i] = 0.0;
    }

  int i;
  for (i=0; i<3; i++) p1[i] = p2[i] = mol.r(0,i);
  for (i=1; i<mol.natom(); i++) {
      for (int j=0; j<3; j++) {
          if (mol.r(i,j) < p1[i]) p1[i] = mol.r(i,j);
          if (mol.r(i,j) > p2[i]) p2[i] = mol.r(i,j);
        }
    }
  for (i=0; i<3; i++) {
      p1[i] = p1[i] - 3.0;
      p2[i] = p2[i] + 3.0;
    }
}





/////////////////////////////////////////////////////////////////////////////
// WriteOrbital

static ClassDesc WriteOrbital_cd(
    typeid(WriteOrbital),"WriteOrbital",1,
    "public WriteGrid", 0, create<WriteOrbital>, 0);

WriteOrbital::WriteOrbital(const Ref<KeyVal> &keyval):
    WriteGrid(WriteOrbital::process_keyval_for_base_class(keyval))
{
  obwfn_ << keyval->describedclassvalue("obwfn");
  if (obwfn_.null()) {
      InputError ex("valid \"obwfn\" missing",
                    __FILE__, __LINE__, "obwfn", "(null)", class_desc());
      try {
          ex.elaborate()
              << "WriteOrbital KeyVal ctor requires"
              << " that \"obwfn\" specifies an object"
              << " of type Wavefunction" << std::endl;
        }
      catch (...) {}
      throw ex;
    }

  orbital_ = keyval->intvalue("orbital", KeyValValueint(-1));
  if (orbital_ < 0 || orbital_ >= obwfn_->oso_dimension().n())
  {
    char buff[IntDigitNum(orbital_)];
    sprintf(buff, "%d", orbital_);
    throw InputError("invalid value", __FILE__, __LINE__, "orbital", buff, class_desc());
  }
}

Ref<KeyVal>
WriteOrbital::process_keyval_for_base_class(Ref<KeyVal> kv) {
  Ref<OneBodyWavefunction> obwfn; obwfn << kv->describedclassvalue("obwfn");
  if (obwfn.null()) { // WriteOrbital constructor will fail anyway, do nothing
    return kv;
  }
  else {
    // if grid is not given, create automatically
    if (kv->exists("grid") == false) {
      const Ref<VDWShape> vdwshape = new VDWShape(obwfn->molecule());
      SCVector3 gmin, gmax;
      vdwshape->boundingbox(-1.0, 1.0, gmin, gmax);
      SCVector3 gorigin = gmin;
      SCVector3 gsize = gmax - gmin;
      const double resolution = 0.2;
      int n[3]; for(int i=0; i<3; ++i) n[i] = int(std::ceil(gsize[i] / 0.2));
      SCVector3 axis0(gsize[0]/n[0], 0.0, 0.0);
      SCVector3 axis1(0.0, gsize[1]/n[1], 0.0);
      SCVector3 axis2(0.0, 0.0, gsize[2]/n[2]);
      Ref<Grid> grid = new Grid(n[0], n[1], n[2], gorigin, axis0, axis1, axis2);

      Ref<AssignedKeyVal> akv = new AssignedKeyVal;
      akv->assign("grid", grid.pointer());

      Ref<AggregateKeyVal> aggkv = new AggregateKeyVal(kv, akv);
      return aggkv;
    }
    else // grid given, do nothing
      return kv;
  }
  return kv; // unreachable
}

WriteOrbital::~WriteOrbital() {}

void
WriteOrbital::label(char* buffer)
{
  sprintf(buffer, "WriteOrbital");
}

Ref<Molecule>
WriteOrbital::get_molecule()
{
  return obwfn_->molecule();
}

double
WriteOrbital::calculate_value(SCVector3 point)
{
	return obwfn_->orbital(point, orbital_);
}

void
WriteOrbital::initialize() {}

/////////////////////////////////////////////////////////////////////////////
// WriteOrbitals

static ClassDesc WriteOrbitals_cd(
    typeid(WriteOrbitals),"WriteOrbitals",1,
    "public WriteVectorGrid", 0, create<WriteOrbitals>, 0);

WriteOrbitals::WriteOrbitals(const Ref<KeyVal> &keyval):
    WriteVectorGrid(WriteOrbital::process_keyval_for_base_class(keyval))
{
  obwfn_ << keyval->describedclassvalue("obwfn");
  if (obwfn_.null()) {
      InputError ex("valid \"obwfn\" missing",
                    __FILE__, __LINE__, "obwfn", "(null)", class_desc());
      try {
          ex.elaborate()
              << "WriteOrbitals KeyVal ctor requires"
              << " that \"obwfn\" specifies an object"
              << " of type Wavefunction" << std::endl;
        }
      catch (...) {}
      throw ex;
    }

  int first = keyval->intvalue("first", KeyValValueint(1));
  const int last = keyval->intvalue("last", KeyValValueint(obwfn_->oso_dimension().n()));
  MPQC_ASSERT(first < obwfn_->oso_dimension().n());
  MPQC_ASSERT(last <= obwfn_->oso_dimension().n());
  const int nmo = last - first + 1;
  for(int o=0; o<nmo; ++o)
    omap_.map.push_back(first++);
}

WriteOrbitals::WriteOrbitals(const Ref<OrbitalSpace> & orbs,
                             const std::vector<int>& labels,
                             const Ref<sc::Grid> & grid,
                             std::string gridformat,
                             std::string gridfile) : WriteVectorGrid(grid,
                                                                     gridformat,
                                                                     gridfile),
                                                     orbs_(orbs)
{
  if (labels.empty() == false) {
    MPQC_ASSERT(orbs_->rank() == labels.size());
    omap_.map = labels;
  }
  else {
    omap_.map.resize(orbs->rank());
    for(int o=0; o<orbs->rank(); ++o)
      omap_.map[o] = o + 1;
  }
}

WriteOrbitals::~WriteOrbitals() {}


void
WriteOrbitals::label(char* buffer)
{
  sprintf(buffer, "WriteOrbitals");
}

Ref<Molecule>
WriteOrbitals::get_molecule()
{
  return obwfn_.nonnull() ? obwfn_->molecule() : orbs_->basis()->molecule();
}

void
WriteOrbitals::calculate_values(const std::vector<SCVector3>& Points, std::vector<double>& Values)
{
  if (obwfn_.nonnull())
    obwfn_->orbitals(Points, Values, omap_.map.front()-1, omap_.map.back()-1, true);
  else
    Wavefunction::orbitals(orbs_, Points, Values);
}

void
WriteOrbitals::initialize() {}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
