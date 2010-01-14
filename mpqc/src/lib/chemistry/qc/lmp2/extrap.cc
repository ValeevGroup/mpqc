
/*
 * Copyright 2009 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This file is a part of the MPQC LMP2 library.
 *
 * The MPQC LMP2 library is free software: you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 */

#include <math.h>

#include <chemistry/qc/lmp2/extrap.h>

namespace sc {

namespace sma2 {

/////////////////////////////////////////////////////////////////////
// Array24SCExtrapData

Array24SCExtrapData::Array24SCExtrapData(Array<2>&array2,bool distrib2,
                                         Array<4>&array4,bool distrib4,
                                         const sc::Ref<sc::MessageGrp> &msg):
  distrib2_(distrib2),
  distrib4_(distrib4),
  msg_(msg)
{
  array2_.init_blocks(array2);
  array2_("p","q") = array2("p","q");
  array4_.init_blocks(array4);
  array4_("p","q","r","s") = array4("p","q","r","s");
}


Array24SCExtrapData::~Array24SCExtrapData()
{
}

void
Array24SCExtrapData::save_data_state(sc::StateOut&)
{
  throw std::runtime_error("Array24SCExtrapData::save_data_state: not implemented");
}

sc::SCExtrapData*
Array24SCExtrapData::copy()
{
  return new Array24SCExtrapData(array2_,distrib2_,array4_,distrib4_,msg_);
}

void
Array24SCExtrapData::zero()
{
  array2_.zero();
  array4_.zero();
}

void
Array24SCExtrapData::accumulate_scaled(double f,
                                       const sc::Ref<sc::SCExtrapData>& ardat)
{
  sc::Ref<Array24SCExtrapData> array
      = sc::require_dynamic_cast<Array24SCExtrapData*>
      (ardat.pointer(), "Array24SCExtrapError::scalar_product");
  array2_("p","q") += f * array->array2_("p","q");
  array4_("p","q","r","s") += f * array->array4_("p","q","r","s");
}

void
Array24SCExtrapData::update(Array<2> &array2,Array<4> &array4)
{
  array2("p","q") = array2_("p","q");
  array4("p","q","r","s") = array4_("p","q","r","s");
}

/////////////////////////////////////////////////////////////////////
// Array24SCExtrapError

Array24SCExtrapError::Array24SCExtrapError(Array<2>& array2,bool distrib2,
                                           Array<4>& array4,bool distrib4,
                                           const sc::Ref<sc::MessageGrp> &msg):
  distrib2_(distrib2),
  distrib4_(distrib4),
  msg_(msg)
{
  array2_.init_blocks(array2);
  array2_("p","q") = array2("p","q");
  array4_.init_blocks(array4);
  array4_("p","q","r","s") = array4("p","q","r","s");
}

Array24SCExtrapError::~Array24SCExtrapError()
{
}

void
Array24SCExtrapError::save_data_state(sc::StateOut&)
{
  throw std::runtime_error("Array24SCExtrapError::save_data_state: not implemented");
}

double
Array24SCExtrapError::error()
{
  double r2 = array2_("p","q") * array2_("p","q");
  double r4 = array4_("p","q","r","s") * array4_("p","q","r","s");
  if (distrib2_) msg_->sum(r2);
  if (distrib4_) msg_->sum(r4);
  return sqrt(r2+r4);
}

double
Array24SCExtrapError::scalar_product(const sc::Ref<sc::SCExtrapError>&eerr)
{
  sc::Ref<Array24SCExtrapError> aeerr
      = sc::require_dynamic_cast<Array24SCExtrapError*>
      (eerr.pointer(), "Array24SCExtrapError::scalar_product");
  double r2 = array2_("p","q") * aeerr->array2_("p","q");
  double r4 = array4_("p","q","r","s") * aeerr->array4_("p","q","r","s");
  if (distrib2_) msg_->sum(r2);
  if (distrib4_) msg_->sum(r4);
  return r2 + r4;
}

/////////////////////////////////////////////////////////////////////
// Array2SCExtrapData

Array2SCExtrapData::Array2SCExtrapData(Array<2>&array,bool distrib,
                                       const sc::Ref<sc::MessageGrp> &msg):
  distrib_(distrib),
  msg_(msg)
{
  array_.init_blocks(array);
  array_("p","q") = array("p","q");
}


Array2SCExtrapData::~Array2SCExtrapData()
{
}

void
Array2SCExtrapData::save_data_state(sc::StateOut&)
{
  throw std::runtime_error("Array2SCExtrapData::save_data_state: not implemented");
}

sc::SCExtrapData*
Array2SCExtrapData::copy()
{
  return new Array2SCExtrapData(array_,distrib_,msg_);
}

void
Array2SCExtrapData::zero()
{
  array_.zero();
}

void
Array2SCExtrapData::accumulate_scaled(double f,
                                      const sc::Ref<sc::SCExtrapData>& ardat)
{
  sc::Ref<Array2SCExtrapData> array
      = sc::require_dynamic_cast<Array2SCExtrapData*>
      (ardat.pointer(), "Array2SCExtrapError::scalar_product");
  array_("p","q") += f * array->array_("p","q");
}

void
Array2SCExtrapData::update(Array<2> &array)
{
  array("p","q") = array_("p","q");
}

/////////////////////////////////////////////////////////////////////
// Array2SCExtrapError

Array2SCExtrapError::Array2SCExtrapError(Array<2>& array,bool distrib,
                                         const sc::Ref<sc::MessageGrp> &msg):
  distrib_(distrib),
  msg_(msg)
{
  array_.init_blocks(array);
  array_("p","q") = array("p","q");
}

Array2SCExtrapError::~Array2SCExtrapError()
{
}

void
Array2SCExtrapError::save_data_state(sc::StateOut&)
{
  throw std::runtime_error("Array2SCExtrapError::save_data_state: not implemented");
}

double
Array2SCExtrapError::error()
{
  double r = array_("p","q") * array_("p","q");
  if (distrib_) msg_->sum(r);
  return sqrt(r);
}

double
Array2SCExtrapError::scalar_product(const sc::Ref<sc::SCExtrapError>&eerr)
{
  sc::Ref<Array2SCExtrapError> aeerr
      = sc::require_dynamic_cast<Array2SCExtrapError*>
      (eerr.pointer(), "Array2SCExtrapError::scalar_product");
  double r = array_("p","q") * aeerr->array_("p","q");
  if (distrib_) msg_->sum(r);
  return r;
}

/////////////////////////////////////////////////////////////////////
// Array4SCExtrapData

Array4SCExtrapData::Array4SCExtrapData(Array<4>&array,bool distrib,
                                       const sc::Ref<sc::MessageGrp> &msg):
  distrib_(distrib),
  msg_(msg)
{
  array_.init_blocks(array);
  array_("p","q","r","s") = array("p","q","r","s");
}


Array4SCExtrapData::~Array4SCExtrapData()
{
}

void
Array4SCExtrapData::save_data_state(sc::StateOut&)
{
  throw std::runtime_error("Array4SCExtrapData::save_data_state: not implemented");
}

sc::SCExtrapData*
Array4SCExtrapData::copy()
{
  return new Array4SCExtrapData(array_,distrib_,msg_);
}

void
Array4SCExtrapData::zero()
{
  array_.zero();
}

void
Array4SCExtrapData::accumulate_scaled(double f,
                                      const sc::Ref<sc::SCExtrapData>& ardat)
{
  sc::Ref<Array4SCExtrapData> array
      = sc::require_dynamic_cast<Array4SCExtrapData*>
      (ardat.pointer(), "Array4SCExtrapError::scalar_product");
  array_("p","q","r","s") += f * array->array_("p","q","r","s");
}

void
Array4SCExtrapData::update(Array<4> &array)
{
  array("p","q","r","s") = array_("p","q","r","s");
}

/////////////////////////////////////////////////////////////////////
// Array4SCExtrapError

Array4SCExtrapError::Array4SCExtrapError(Array<4>& array,bool distrib,
                                         const sc::Ref<sc::MessageGrp> &msg):
  distrib_(distrib),
  msg_(msg)
{
  array_.init_blocks(array);
  array_("p","q","r","s") = array("p","q","r","s");
}

Array4SCExtrapError::~Array4SCExtrapError()
{
}

void
Array4SCExtrapError::save_data_state(sc::StateOut&)
{
  throw std::runtime_error("Array4SCExtrapError::save_data_state: not implemented");
}

double
Array4SCExtrapError::error()
{
  double r = array_("p","q","r","s") * array_("p","q","r","s");
  if (distrib_) msg_->sum(r);
  return sqrt(r);
}

double
Array4SCExtrapError::scalar_product(const sc::Ref<sc::SCExtrapError>&eerr)
{
  sc::Ref<Array4SCExtrapError> aeerr
      = sc::require_dynamic_cast<Array4SCExtrapError*>
      (eerr.pointer(), "Array4SCExtrapError::scalar_product");
  double r = array_("p","q","r","s") * aeerr->array_("p","q","r","s");
  if (distrib_) msg_->sum(r);
  return r;
}

/////////////////////////////////////////////////////////////////////
// Array6SCExtrapData

Array6SCExtrapData::Array6SCExtrapData(Array<6>&array,bool distrib,
                                       const sc::Ref<sc::MessageGrp> &msg):
  distrib_(distrib),
  msg_(msg)
{
  array_.init_blocks(array);
  array_("p","q","r","s","t","u") = array("p","q","r","s","t","u");
}


Array6SCExtrapData::~Array6SCExtrapData()
{
}

void
Array6SCExtrapData::save_data_state(sc::StateOut&)
{
  throw std::runtime_error("Array6SCExtrapData::save_data_state: not implemented");
}

sc::SCExtrapData*
Array6SCExtrapData::copy()
{
  return new Array6SCExtrapData(array_,distrib_,msg_);
}

void
Array6SCExtrapData::zero()
{
  array_.zero();
}

void
Array6SCExtrapData::accumulate_scaled(double f,
                                      const sc::Ref<sc::SCExtrapData>& ardat)
{
  sc::Ref<Array6SCExtrapData> array
      = sc::require_dynamic_cast<Array6SCExtrapData*>
      (ardat.pointer(), "Array6SCExtrapError::scalar_product");
  array_("p","q","r","s","t","u") += f * array->array_("p","q","r","s","t","u");
}

void
Array6SCExtrapData::update(Array<6> &array)
{
  array("p","q","r","s","t","u") = array_("p","q","r","s","t","u");
}

/////////////////////////////////////////////////////////////////////
// Array6SCExtrapError

Array6SCExtrapError::Array6SCExtrapError(Array<6>& array,bool distrib,
                                         const sc::Ref<sc::MessageGrp> &msg):
  distrib_(distrib),
  msg_(msg)
{
  array_.init_blocks(array);
  array_("p","q","r","s","t","u") = array("p","q","r","s","t","u");
}

Array6SCExtrapError::~Array6SCExtrapError()
{
}

void
Array6SCExtrapError::save_data_state(sc::StateOut&)
{
  throw std::runtime_error("Array6SCExtrapError::save_data_state: not implemented");
}

double
Array6SCExtrapError::error()
{
  double r = array_("p","q","r","s","t","u") * array_("p","q","r","s","t","u");
  if (distrib_) msg_->sum(r);
  return sqrt(r);
}

double
Array6SCExtrapError::scalar_product(const sc::Ref<sc::SCExtrapError>&eerr)
{
  sc::Ref<Array6SCExtrapError> aeerr
      = sc::require_dynamic_cast<Array6SCExtrapError*>
      (eerr.pointer(), "Array6SCExtrapError::scalar_product");
  double r = array_("p","q","r","s","t","u") * aeerr->array_("p","q","r","s","t","u");
  if (distrib_) msg_->sum(r);
  return r;
}

}

}
