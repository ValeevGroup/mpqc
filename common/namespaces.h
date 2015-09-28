#pragma once
#ifndef TCC_COMMON_NAMESPACES_H
#define TCC_COMMON_NAMESPACES_H

namespace TiledArray {}
namespace Eigen {}
namespace madness {}

namespace tcc {
namespace integrals {}
namespace tensor {}
}


namespace TA = TiledArray;
namespace Eig = Eigen;
namespace mad = madness;

namespace tints = tcc::integrals;
namespace ttensor = tcc::tensor;

// MPQC FWD DECS
namespace mpqc {
    namespace integrals {}
}

namespace mpqc_ints = mpqc::integrals;

#endif // TCC_COMMON_NAMESPACES_H
