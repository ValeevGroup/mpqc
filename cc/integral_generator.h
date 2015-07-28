//
// Created by Chong Peng on 7/27/15.
//

#ifndef TILECLUSTERCHEM_INTEGRAL_GENERATOR_H
#define TILECLUSTERCHEM_INTEGRAL_GENERATOR_H

#include <libint2/engine.h>
#include "../include/tiledarray.h"
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
      auto shells0 = (*cluster_shells_)[0].shells(index[0]);
      auto shells1 = (*cluster_shells_)[1].shells(index[1]);
      auto shells2 = (*cluster_shells_)[2].shells(index[2]);
      auto shells3 = (*cluster_shells_)[3].shells(index[3]);

      // total number of shells in each dimension
      std::size_t nshells0 = shells0.size();
      std::size_t nshells1 = shells1.size();
      std::size_t nshells2 = shells2.size();
      std::size_t nshells3 = shells3.size();

      // total number of basis function in each dimension
      std::size_t nfunctions0 = (*cluster_shells_)[0].nfunctions(index[0]);
      std::size_t nfunctions1 = (*cluster_shells_)[1].nfunctions(index[1]);
      std::size_t nfunctions2 = (*cluster_shells_)[2].nfunctions(index[2]);
      std::size_t nfunctions3 = (*cluster_shells_)[3].nfunctions(index[3]);
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

        for (auto s1=0l; s1!=nshells1; ++s1){

          std::size_t ns1 = shells1[s1].size();

          for (auto s2=0l; s2!=nshells2; ++s2){

            std::size_t ns2 = shells2[s2].size();

            for (auto s3=0l; s3!=nshells3; ++s3){

              std::size_t ns3 = shells3[s3].size();

              // compute shell pair
              // (s0 s1|s2, s3)
              const auto* buf = engine.compute(shells0[s0],shells1[s1],shells2[s2],shells3[s3]);

              //store it into tile
              std::size_t ns23 = ns2*ns3;
              std::size_t ns123 = ns1*ns23;
              for (auto ins0 = 0l; ins0 < ns0; ins0++){
                for (auto ins1 = 0l; ins1 < ns1; ins1++){
                  for (auto ins2 = 0l; ins2 < ns2; ins2++){

                    // position in the buffer
                    std::size_t buf_dest = ins0*ns123 + ins1*ns23 + ins2*ns3;

                    // position in the whole tile
                    std::size_t tile_dest = (bf0+ins0)*nfunctions123 +(bf1+ins1)*nfunctions23
                                            + (bf2+ins2)*nfunctions3 + bf3;

                    std::copy(buf+buf_dest, buf+buf_dest+ns3, tile.data()+tile_dest);
                  }
                }
              }

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
