/*
 * ao_wfn.h
 *
 *  Created on: Apr 27, 2016
 *      Author: Drew Lewis
 */

#ifndef MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_
#define MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_

#include <mpqc/util/keyval/keyval.hpp>
#include <mpqc/chemistry/qc/wfn/wfn.h>

namespace mpqc {
namespace qc {

class AOWfn : public Wfn {

 public:
  AOWfn(KeyVal const &kv);

  void visit(PropertyBase *pb) override; 
};

}  // namespace qc
}  // namespace mpqc

#endif  // MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_
