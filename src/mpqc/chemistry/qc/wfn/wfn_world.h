/*
 * wfn_world.h
 *
 *  Created on: Aug 18, 2016
 *      Author: Drew Lewis
 */
#ifndef MPQC_CHEMISTRY_QC_WFN_WFN_WORLD_H_
#define MPQC_CHEMISTRY_QC_WFN_WFN_WORLD_H_

#include <mpqc/util/keyval/keyval.hpp>

namespace mpqc {
namespace qc {

class WfnWorld : public DescribedClass {
 public:
  WfnWorld(KeyVal const &kv);
};

}  // namespace qc
}  // namespace mpqc

#endif  // MPQC_CHEMISTRY_QC_WFN_WFN_WORLD_H_
