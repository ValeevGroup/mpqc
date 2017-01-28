#ifndef SRC_MPQC_CHEMISTRY_QC_BASIS_LINKAGE_H_
#define SRC_MPQC_CHEMISTRY_QC_BASIS_LINKAGE_H_

#include "mpqc/util/keyval/forcelink.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

class AtomicBasis;
mpqc::detail::ForceLink<AtomicBasis> fl1;

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc

#endif /* SRC_MPQC_CHEMISTRY_QC_BASIS_LINKAGE_H_ */
