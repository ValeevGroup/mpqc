/*
 * ccsd_pno.cpp
 *
 *  Created on: Jan 4, 2017
 *      Author: jinmei
 */
#include "ccsd_pno.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 0
template class mpqc::lcao::CCSD_PNO<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("CCSD-PNO", mpqc::lcao::CCSD_PNO<TA::TensorD, TA::DensePolicy>);
#elif TA_DEFAULT_POLICY == 1
template class mpqc::lcao::CCSD_PNO<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("CCSD-PNO", mpqc::lcao::CCSD_PNO<TA::TensorD, TA::SparsePolicy>);
#endif
