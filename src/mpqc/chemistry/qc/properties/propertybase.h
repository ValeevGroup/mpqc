/*
 * property_base.h
 *
 *  Created on: Aug 18, 2016
 *      Author: Drew Lewis
 */

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_PROPERTYBASE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_PROPERTYBASE_H_

#include "mpqc/util/keyval/keyval.h"
#include "mpqc/chemistry/qc/wfn/wfn_forward.h"

namespace mpqc {
namespace lcao {

class PropertyBase : public DescribedClass {
 public:
  virtual void apply(Wavefunction *) = 0;
};

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_PROPERTYBASE_H_
