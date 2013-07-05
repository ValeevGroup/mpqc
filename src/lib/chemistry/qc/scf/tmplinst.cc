
#ifdef HAVE_CONFIG_H
#include <mpqc_config.h>
#endif

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION

#include <util/misc/formio.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/ltbgrad.h>

#include <chemistry/qc/scf/clhftmpl.h>
#include <chemistry/qc/scf/hsoshftmpl.h>
#include <chemistry/qc/scf/osshftmpl.h>
#include <chemistry/qc/scf/tchftmpl.h>
#include <chemistry/qc/scf/uhftmpl.h>

using namespace sc;

template class GBuild<LocalCLHFContribution>;
template class GBuild<LocalCLHFEnergyContribution>;

template class LocalGBuild<LocalCLHFContribution>;
template class LocalGBuild<LocalCLHFEnergyContribution>;

template class TBGrad<LocalCLHFGradContribution>;
template class LocalTBGrad<LocalCLHFGradContribution>;

///////////////////////////////////////////////////////////////////////////

template class GBuild<LocalHSOSContribution>;
template class GBuild<LocalHSOSEnergyContribution>;
template class LocalGBuild<LocalHSOSContribution>;
template class LocalGBuild<LocalHSOSEnergyContribution>;

template class TBGrad<LocalHSOSGradContribution>;
template class LocalTBGrad<LocalHSOSGradContribution>;

///////////////////////////////////////////////////////////////////////////

template class GBuild<LocalOSSContribution>;
template class GBuild<LocalOSSEnergyContribution>;
template class LocalGBuild<LocalOSSContribution>;
template class LocalGBuild<LocalOSSEnergyContribution>;

template class TBGrad<LocalOSSGradContribution>;
template class LocalTBGrad<LocalOSSGradContribution>;

///////////////////////////////////////////////////////////////////////////

template class GBuild<LocalTCContribution>;
template class GBuild<LocalTCEnergyContribution>;
template class LocalGBuild<LocalTCContribution>;
template class LocalGBuild<LocalTCEnergyContribution>;

template class TBGrad<LocalTCGradContribution>;
template class LocalTBGrad<LocalTCGradContribution>;

///////////////////////////////////////////////////////////////////////////

template class GBuild<LocalUHFContribution>;
template class GBuild<LocalUHFEnergyContribution>;
template class LocalGBuild<LocalUHFContribution>;
template class LocalGBuild<LocalUHFEnergyContribution>;

template class TBGrad<LocalUHFGradContribution>;
template class LocalTBGrad<LocalUHFGradContribution>;

#endif
