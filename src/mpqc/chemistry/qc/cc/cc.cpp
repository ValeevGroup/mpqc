//
// Created by Chong Peng on 6/6/16.
//

#include <mpqc/chemistry/qc/cc/ccsd.h>
//#include <mpqc/chemistry/qc/cc/ccsd_t.h>
//#include <mpqc/chemistry/qc/cc/dbccsd.h>

MPQC_CLASS_EXPORT_KEY2("CCSD", mpqc::cc::CCSD<TA::TensorD, TA::SparsePolicy>);

namespace mpqc {
namespace cc {

template class CCSD<TA::TensorD, TA::SparsePolicy>;

// template class CCSD_T<TA::TensorD, TA::SparsePolicy>;

// template class DBCCSD<TA::TensorD, TA::SparsePolicy>;
}
}
