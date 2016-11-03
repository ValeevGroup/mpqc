#include <mpqc/chemistry/qc/properties/propertybase.h>
#include <mpqc/chemistry/qc/wfn/ao_wfn.h>

namespace mpqc{
namespace qc{


template class AOWavefunction<TA::TensorD, TA::SparsePolicy>;
template class AOWavefunction<TA::TensorD, TA::DensePolicy>;

}//namespace qc
}//namespace mpqc

