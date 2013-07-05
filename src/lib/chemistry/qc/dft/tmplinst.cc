
#ifdef HAVE_CONFIG_H
#include <mpqc_config.h>
#endif

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION

#include <util/misc/formio.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/ltbgrad.h>

#include <chemistry/qc/dft/clkstmpl.h>
#include <chemistry/qc/dft/ukstmpl.h>
#include <chemistry/qc/dft/hsoskstmpl.h>

using namespace sc;

template class GBuild<LocalUKSContribution>;
template class GBuild<LocalUKSEnergyContribution>;
template class LocalGBuild<LocalUKSContribution>;
template class LocalGBuild<LocalUKSEnergyContribution>;

///////////////////////////////////////////////////////////////////////////

template class GBuild<LocalCLKSContribution>;
template class GBuild<LocalCLKSEnergyContribution>;

template class LocalGBuild<LocalCLKSContribution>;
template class LocalGBuild<LocalCLKSEnergyContribution>;

///////////////////////////////////////////////////////////////////////////

template class GBuild<LocalHSOSKSContribution>;
template class GBuild<LocalHSOSKSEnergyContribution>;

template class LocalGBuild<LocalHSOSKSContribution>;
template class LocalGBuild<LocalHSOSKSEnergyContribution>;

#endif
