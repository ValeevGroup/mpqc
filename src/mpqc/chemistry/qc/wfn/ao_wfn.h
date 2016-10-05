/*
 * ao_wfn.h
 *
 *  Created on: Aug 17, 2016
 *      Author: Drew Lewis
 */

#ifndef MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_
#define MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_

#include <mpqc/chemistry/qc/wfn/wfn.h>
#include <mpqc/util/keyval/keyval.hpp>

namespace mpqc {
namespace qc {

class AOWfn : public Wfn {
 public:
  AOWfn(KeyVal const &kv);
  ~AOWfn();

  void compute(PropertyBase *pb) override;
};

}  // namespace qc
}  // namespace mpqc

#endif  // MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_
