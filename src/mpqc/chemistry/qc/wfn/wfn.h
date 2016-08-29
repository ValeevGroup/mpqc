/*
 * wfn.h
 *
 *  Created on: Apr 27, 2016
 *      Author: Drew Lewis
 */

#ifndef MPQC_CHEMISTRY_QC_WFN_WFN_H_
#define MPQC_CHEMISTRY_QC_WFN_WFN_H_

#include <mpqc/util/keyval/keyval.hpp>
#include <mpqc/chemistry/qc/wfn/wfn_world.h>
#include <mpqc/chemistry/qc/properties/propertybase.h>

#include <memory>
#include <functional>

namespace mpqc {
namespace qc {

class PropertyBase;

class Wfn : public DescribedClass {
 public:
  using ArrayType = WfnWorld::ArrayType;

 private:
  /*! Pointer to to the WfnWorld
   *
   * \note No need to make this shared Wfn is just a member of the world it
   *lives in so no ownership here.
   *
   * \warning Wfn should never delete or allocate this pointer.
   */
  WfnWorld* wfn_world_ = nullptr;

 public:
  Wfn(KeyVal const& kv) : wfn_world_(kv.value<WfnWorld*>("wfn_world")){};

  virtual void compute(PropertyBase* pb) = 0;

  WfnWorld* wfn_world() { return wfn_world_; }
};

}  // namespace qc
}  // namespace mpqc

#endif  // MPQC_CHEMISTRY_QC_WFN_WFN_H_
