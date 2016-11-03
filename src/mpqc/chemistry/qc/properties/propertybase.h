/*
 * property_base.h
 *
 *  Created on: Aug 18, 2016
 *      Author: Drew Lewis
 */

#ifndef MPQC_CHEMISTRY_QC_WFN_PROPERTY_BASE_H_
#define MPQC_CHEMISTRY_QC_WFN_PROPERTY_BASE_H_

#include <mpqc/util/keyval/keyval.hpp>
#include <mpqc/chemistry/qc/wfn/wfn_forward.h>

namespace mpqc {
namespace qc {

class PropertyBase : public DescribedClass {
 public:
  virtual void apply(Wavefunction *) = 0;
};

}  // namespace qc
}  // namespace mpqc

#endif  // MPQC_CHEMISTRY_QC_WFN_PROPERTY_BASE_H_
