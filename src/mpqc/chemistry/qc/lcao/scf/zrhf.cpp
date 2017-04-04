
#include "mpqc/chemistry/qc/lcao/scf/zrhf.h"
#include "mpqc/util/keyval/forcelink.h"

#include "mpqc/chemistry/qc/lcao/scf/decomposed_rij.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_cond_ortho.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_soad.h"

#include <clocale>
#include <sstream>

#include <tiledarray.h>

#include "mpqc/chemistry/qc/lcao/wfn/wfn.h"
#include "mpqc/util/external/madworld/parallel_file.h"
#include "mpqc/util/external/madworld/parallel_print.h"
#include "mpqc/util/keyval/keyval.h"
#include "mpqc/util/misc/formio.h"

#if TA_DEFAULT_POLICY == 0

#elif TA_DEFAULT_POLICY == 1

template class mpqc::lcao::zRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("zRHF", mpqc::lcao::zRHF<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::DFzRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("DF-zRHF", mpqc::lcao::DFzRHF<TA::TensorD, TA::SparsePolicy>);

#endif


