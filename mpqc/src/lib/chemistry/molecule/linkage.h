
#ifndef _chemistry_molecule_linkage_h
#define _chemistry_molecule_linkage_h

#ifndef __PIC__

#include <chemistry/molecule/coor.h>

#include <math/scmat/linkage.h>
#include <math/optimize/linkage.h>

const ClassDesc &molecule_force_link_a_ = RedundMolecularCoor::class_desc_;
const ClassDesc &molecule_force_link_b_ = CartMolecularCoor::class_desc_;
const ClassDesc &molecule_force_link_c_ = SymmMolecularCoor::class_desc_;

#endif /* __PIC__ */

#endif
