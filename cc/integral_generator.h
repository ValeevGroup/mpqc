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

namespace tcc {
    namespace  cc {


        //IntegralGenerator for two body integrals, use to work with LazyIntegral
        // it used libint2 for generation of integrals
        template<libint2::MultiplicativeSphericalTwoBodyKernel Kernel>
        class TwoBodyIntGenerator {

        public:
            typedef libint2::TwoBodyEngine<Kernel> Engine;
            typedef tcc::integrals::EnginePool<Engine> EnginePool;


            TwoBodyIntGenerator() = default;

            TwoBodyIntGenerator(TwoBodyIntGenerator const &) = default;

            TwoBodyIntGenerator &operator=(
                    TwoBodyIntGenerator const &) = default;


            TwoBodyIntGenerator(const std::shared_ptr<EnginePool> &pool,
                                const std::shared_ptr<std::vector<basis::ClusterShells>> &cluster_shells)
                    : pool_(pool), cluster_shells_(cluster_shells) { }


            void set_pool(const std::shared_ptr<EnginePool> &pool){
                pool_ = pool;
            }

            void set_shell(const std::shared_ptr<std::vector<basis::ClusterShells>> &cluster_shells){
                cluster_shells_ = cluster_shells;
            }

            // return a tile of integral in chemical notation
            // (block1 block2| block3 block4)
            TA::Tensor <double> compute(const TA::Range &r,
                                        const std::vector<std::size_t> &index) {

                // assert size of index is 4
                assert(index.size() == 4);
                // create the tile
                TA::Tensor <double> tile(r);

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

                // get the engine
                auto engine = pool_->local();

                // bf, to track the position in total basis function
                std::size_t bf0, bf1, bf2, bf3 = 0;
                auto lo0 = r.lobound_data()[0];
                auto lo1 = r.lobound_data()[1];
                auto lo2 = r.lobound_data()[2];
                auto lo3 = r.lobound_data()[3];

                // compute
                bf0 = lo0;
                for (auto s0 = 0l; s0 != nshells0; ++s0) {

                    std::size_t ns0 = shells0[s0].size();
                    bf1 = lo1;

                    for (auto s1 = 0l; s1 != nshells1; ++s1) {

                        std::size_t ns1 = shells1[s1].size();
                        bf2 = lo2;

                        for (auto s2 = 0l; s2 != nshells2; ++s2) {

                            std::size_t ns2 = shells2[s2].size();
                            bf3 = lo3;

                            for (auto s3 = 0l; s3 != nshells3; ++s3) {

                                std::size_t ns3 = shells3[s3].size();

                                // compute shell pair
                                // (s0 s1|s2, s3)
                                const auto *buf = engine.compute(shells0[s0],
                                                                 shells1[s1],
                                                                 shells2[s2],
                                                                 shells3[s3]);

                                auto lowbound = {bf0, bf1, bf2, bf3};
                                auto upbound = {bf0 + ns0, bf1 + ns1, bf2 + ns2,
                                                    bf3 + ns3};
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
