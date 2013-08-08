/*
 * dmm_scf.h
 *
 *  Created on: Jul 22, 2013
 *      Author: drewlewis
 */

#ifndef DMM_SCF_H
#define DMM_SCF_H

#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include <util/misc/exenv.h>
#include <tiled_array.h>

namespace TA = TiledArray;
namespace sc {
namespace dmm {
    class DMMHF {
    public:

        DMMHF(const Ref<KeyVal>&);
        ~DMMHF();

        void print(std::ostream&o=ExEnv::out0()) const;

    private:

        Ref<Molecule> mol_;
        Ref<GaussianBasisSet> ao_basis_;
        madness::world &world;

    };
} // namespace dmm
} // namespace sc

#endif


#endif /* DMM_SCF_H_ */
