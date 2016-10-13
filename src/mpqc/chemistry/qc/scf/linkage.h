//
// Created by Chong Peng on 10/5/16.
//

#ifndef MPQC_CHEMISTRY_QC_SCF_LINKAGE_H
#define MPQC_CHEMISTRY_QC_SCF_LINKAGE_H

#include <mpqc/chemistry/qc/scf/rhf.h>
#include <mpqc/util/keyval/forcelink.h>

namespace mpqc{
namespace qc{

mpqc::detail::ForceLink<mpqc::scf::RHF> fl1;
mpqc::detail::ForceLink<mpqc::scf::RIRHF> fl2;


}
}

#endif //MPQC_CHEMISTRY_QC_SCF_LINKAGE_H
