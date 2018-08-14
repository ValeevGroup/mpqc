//
// Created by Eduard Valeyev on 8/13/18.
//

#include <clocale>
#include <madness/world/worldgop.h>
#include <tiledarray.h>

#include "mpqc/mpqc_init.h"
#include "mpqc/util/keyval/keyval.h"
#include "mpqc/util/core/exenv.h"
#include "mpqc/chemistry/qc/lcao/factory/ao_factory.h"

#include "mpqc/chemistry/qc/lcao/factory/linkage.h"

madness::World* world_ptr = nullptr;

int main( int argc, char* argv[] )
{
  world_ptr = &madness::initialize(argc, argv);
  assert(world_ptr != nullptr);
  mpqc::initialize(argc, argv, *world_ptr);

  std::shared_ptr<mpqc::KeyVal> kv =
      mpqc::MPQCInit::instance().make_keyval(*world_ptr, "aoints.json");
  kv->assign("world", world_ptr);   // set "$:world" keyword to &world to define the default execution context for this input
  TA::set_default_world(*world_ptr);  // must specify default world to avoid madness::World::get_default() getting called

  // madness::World::get_default() getting called
  using mpqc::lcao::gaussian::AOFactory;
  auto aofactory = std::make_shared<AOFactory<TA::TensorD,TA::SparsePolicy>>(*kv);

  auto eri4 = aofactory->compute(L"(μ ν| G |κ λ)");
  auto eri3 = aofactory->compute(L"( Κ | G |κ λ)");
  auto eri2 = aofactory->compute(L"( Κ | G | Λ )");
  auto eri2_sqrt_inv = aofactory->compute(L"( Κ | G | Λ )[inv_sqr]");

  // test eri2_sqrt_inv
  {
    auto identity = mpqc::math::create_diagonal_matrix(eri2, 1.0);
    TA::TSpArrayD zero;
    zero("m,l") = identity("m,l") - eri2("m,n") * eri2_sqrt_inv("n,k") * eri2_sqrt_inv("k,l");
    mpqc::ExEnv::out0() << "testing V^{-1/2}: this should be zero = " << TA::norm2(zero) << std::endl;
  }

  // test DF reconstruction of eri4
  {
    TA::TSpArrayD B;
    B("K,p,q") = eri2_sqrt_inv("K,X") * eri3("X,p,q");
    TA::TSpArrayD dferror4;
    dferror4("p,q,r,s") = eri4("p,q,r,s") - B("K,p,q") * B("K,r,s");
    mpqc::ExEnv::out0() << "testing DF: this should be small (and decrease with increasing DF basis) = " << TA::norm2(dferror4) << std::endl;
  }
  
  mpqc::finalize();
  madness::finalize();

  return 0;
}
