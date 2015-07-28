//
// Created by Chong Peng on 7/27/15.
//

#ifndef TILECLUSTERCHEM_INTEGRAL_GENERATOR_H
#define TILECLUSTERCHEM_INTEGRAL_GENERATOR_H

#include <libint2/engine.h>
#include "../include/tiledarray.h"
#include <TiledArray/tensor/tensor_map.h>
#include "../common/namespaces.h"
#include "../basis/cluster_shells.h"
#include "../integrals/integral_engine_pool.h"

namespace tcc{
  namespace  cc {


    //IntegralGenerator for two body AO integrals, use to work with LazyIntegral
    template<libint2::MultiplicativeSphericalTwoBodyKernel Kernel>
    class TwoBodyIntGenerator {

    public:
      typedef libint2::TwoBodyEngine<Kernel> Engine;
      typedef tcc::integrals::EnginePool<Engine> EnginePool;


      TwoBodyIntGenerator() = delete;
      TwoBodyIntGenerator(TwoBodyIntGenerator const &) = default;
      TwoBodyIntGenerator& operator=(TwoBodyIntGenerator const &) = default;


      TwoBodyIntGenerator(const std::shared_ptr<EnginePool> &pool,
                          const std::shared_ptr<std::vector<basis::ClusterShells>> &cluster_shells)
              : pool_(pool), cluster_shells_(cluster_shells) { }

      TA::Tensor<double> compute(const TA::Range& r, const std::array<std::size_t, 4>& index){

      // create the tile
      TA::Tensor<double> tile(r);

      // get the shells in each dimension
      auto shells0 = (*cluster_shells_)[index[0]].flattened_shells();
      auto shells1 = (*cluster_shells_)[index[1]].flattened_shells();
      auto shells2 = (*cluster_shells_)[index[2]].flattened_shells();
      auto shells3 = (*cluster_shells_)[index[3]].flattened_shells();

      // total number of shells in each dimension
      std::size_t nshells0 = shells0.size();
      std::size_t nshells1 = shells1.size();
      std::size_t nshells2 = shells2.size();
      std::size_t nshells3 = shells3.size();

      // total number of basis function in each dimension
      std::size_t nfunctions0 = (*cluster_shells_)[index[0]].flattened_nfunctions();
      std::size_t nfunctions1 = (*cluster_shells_)[index[1]].flattened_nfunctions();
      std::size_t nfunctions2 = (*cluster_shells_)[index[2]].flattened_nfunctions();
      std::size_t nfunctions3 = (*cluster_shells_)[index[3]].flattened_nfunctions();
      std::size_t nfunctions23 = nfunctions2*nfunctions3;
      std::size_t nfunctions123 = nfunctions1*nfunctions23;


      // get the engine
      auto engine = pool_->local();


      // bf, to track the position in total basis function
      std::size_t bf0, bf1, bf2, bf3;
      bf0 = bf1 = bf2 = bf3 = 0l;

      // compute
      for (auto s0=0l; s0!=nshells0; ++s0){

        std::size_t ns0 = shells0[s0].size();
        bf1 = 0l;

        for (auto s1=0l; s1!=nshells1; ++s1){

          std::size_t ns1 = shells1[s1].size();
          bf2 = 0l;

          for (auto s2=0l; s2!=nshells2; ++s2){

            std::size_t ns2 = shells2[s2].size();
            bf3 = 0l;

            for (auto s3=0l; s3!=nshells3; ++s3){

              std::size_t ns3 = shells3[s3].size();

              // compute shell pair
              // (s0 s1|s2, s3)
              const auto* buf = engine.compute(shells0[s0],shells1[s1],shells2[s2],shells3[s3]);

              //store it into tile
              auto lowbound = {bf0,bf1,bf2,bf3};
              auto upbound = {bf0+ns0,bf1+ns1,bf2+ns2,bf3+ns3};
              auto view = tile.block(lowbound, upbound);
              auto map = TA::make_map(buf, lowbound, upbound);
              view = map;

              bf3 += ns3;
            }
            bf2 += ns2;
          }
          bf1 += ns1;
        }
        bf0 += ns0;
      }

      return tile;
    }

    private:
      std::shared_ptr<EnginePool> pool_;
      std::shared_ptr<std::vector<tcc::basis::ClusterShells>> cluster_shells_;
    };

  }
}

#endif //TILECLUSTERCHEM_INTEGRAL_GENERATOR_H
