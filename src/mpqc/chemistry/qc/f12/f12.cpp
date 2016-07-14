//
// Created by Chong Peng on 4/29/16.
//


#include <mpqc/chemistry/qc/f12/ccsdf12.h>
#include <mpqc/chemistry/qc/f12/mp2f12.h>
#include <mpqc/chemistry/qc/f12/dbmp2f12.h>
#include <mpqc/chemistry/qc/f12/dbccsdf12.h>

namespace mpqc{
namespace f12{


    template class MP2F12<TA::TensorD>;

    template class CCSDF12<TA::TensorD>;

    template class DBMP2F12<TA::TensorD>;

    template class DBCCSDF12<TA::TensorD>;

}
}
