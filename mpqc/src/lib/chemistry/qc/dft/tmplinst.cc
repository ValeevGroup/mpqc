
#ifdef __GNUC__

#include <util/misc/formio.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/ltbgrad.h>

#include <chemistry/qc/dft/clkstmpl.h>
#include <chemistry/qc/dft/ukstmpl.h>


template class GBuild<LocalUKSContribution>;
template class GBuild<LocalUKSEnergyContribution>;
template class LocalGBuild<LocalUKSContribution>;
template class LocalGBuild<LocalUKSEnergyContribution>;

template class TBGrad<LocalUKSGradContribution>;
template class LocalTBGrad<LocalUKSGradContribution>;

///////////////////////////////////////////////////////////////////////////

template class GBuild<LocalCLKSContribution>;
template class GBuild<LocalCLKSEnergyContribution>;

template class LocalGBuild<LocalCLKSContribution>;
template class LocalGBuild<LocalCLKSEnergyContribution>;

template class TBGrad<LocalCLKSGradContribution>;
template class LocalTBGrad<LocalCLKSGradContribution>;

#endif
