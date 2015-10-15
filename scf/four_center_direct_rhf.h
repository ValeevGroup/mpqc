#pragma once
#ifndef MPQC_SCF_FOURCENTERDIRECTRHF_H
#define MPQC_SCF_FOURCENTERDIRECTRHF_H

#include "scf_base.h"
#include "../include/tiledarray.h"
#include "../integrals/direct_tile.h"
#include "../include/libint.h"
#include "../integrals/integral_engine_pool.h"

#include "../integrals/task_integrals_common.h"




namespace mpqc {
namespace scf {

class FourCenterDirectRHF : public ScfBase {
  private:
    using array_type = TAArray<2, SpPolicy>;
    array_type F_;
    array_type D_;
    array_type S_;
    array_type H_;

    using TwoEngine = libint2::TwoBodyEngine<libint2::Coulomb>;
    static std::shared_ptr<integrals::ShrPool<TwoEngine>> eri_engine_pool_;

  public:

    
};

} // namespace scf
} // namespace mpqc
#endif // MPQC_SCF_FOURCENTERDIRECTRHF_H
