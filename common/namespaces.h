#pragma once
#ifndef TCC_COMMON_NAMESPACES_H
#define TCC_COMMON_NAMESPACES_H

namespace TiledArray {}
namespace Eigen {}
namespace madness {}


namespace TA = TiledArray;
namespace Eig = Eigen;
namespace mad = madness;

// MPQC FWD DECS
namespace mpqc {
    namespace integrals {}
    namespace molecule {}
    namespace mol = molecule;
}

namespace mpqc_ints = mpqc::integrals;

#endif // TCC_COMMON_NAMESPACES_H
