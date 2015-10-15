#pragma once
#ifndef MPQC_SCF_SCFBASE_H
#define MPQC_SCF_SCFBASE_H

#include "../common/typedefs.h"
#include "../molecule/molecule_fwd.h"
#include "../basis/basis_fwd.h"

namespace mpqc {
namespace scf {

    /*! \brief Scf base class that defines the scf interface.
     */
    class ScfBase {
        public:
            virtual ~ScfBase() = default;

            virtual double energy() const = 0;

            virtual MatrixD fock_matrix() const = 0;

            virtual MatrixD density_matrix() const = 0;

            virtual MatrixD coeff_matrix() const = 0;

            virtual basis::Basis basis() const = 0;

            virtual molecule::Molecule const & molecule() const = 0;
    };

} // namespace scf
} // namespace mpqc

#endif // MPQC_SCF_SCFBASE_H
