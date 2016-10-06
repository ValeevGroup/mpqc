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

/**
 * Wfn is a wave function that lives in a WfnWorld
 */
class Wfn : public DescribedClass {
 public:
  using ArrayType = WfnWorld::ArrayType;

 private:
  /*! Pointer to the WfnWorld
   *
   * \note No need to make this shared Wfn is just a member of the world it
   *lives in so no ownership here.
   *
   * \warning Wfn should never delete or allocate this pointer.
   *
   * \note by chong I changed this to shared pointer, for example, MP2 and HF
   *          will share the same wfn_world
   */
  std::shared_ptr<WfnWorld> wfn_world_;

 public:
  /**
   * KeyVal constructor
   *
   * keywords
   *
   * @param  wfn_world,  WfnWorld object, if not provided, will check "$:wfn_world" and "..:wfn_world"
   *
   */
  Wfn(const KeyVal& kv);
  virtual ~Wfn();

  virtual void compute(PropertyBase* pb) = 0;
  virtual double value() = 0;

  std::shared_ptr<WfnWorld> wfn_world() { return wfn_world_; }
};

}  // namespace qc
}  // namespace mpqc

#endif  // MPQC_CHEMISTRY_QC_WFN_WFN_H_
