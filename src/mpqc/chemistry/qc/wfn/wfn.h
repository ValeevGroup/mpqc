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
#include <memory>

namespace mpqc {
namespace qc {

class PropertyBase;

class Wfn : public DescribedClass {
 private:
  std::shared_ptr<WfnWorld> wfn_world_;
  std::vector<std::shared_ptr<PropertyBase>> properties_;

 public:
  Wfn(KeyVal const &kv) 
      : wfn_world_(std::make_shared<WfnWorld>(kv.keyval("wfn_world"))){};

  virtual void visit(PropertyBase *) = 0;

  std::shared_ptr<WfnWorld> wfn_world() { return wfn_world_; }
};

}  // namespace qc
}  // namespace mpqc

#endif  // MPQC_CHEMISTRY_QC_WFN_WFN_H_
