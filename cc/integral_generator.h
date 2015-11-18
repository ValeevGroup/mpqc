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
#include "../integrals/task_integral_kernels.h"
#include "../integrals/integral_engine_pool.h"

namespace mpqc {
    namespace  cc {


        //IntegralGenerator for two body integrals, use to work with LazyIntegral
        // it used libint2 for generation of integrals
        template<libint2::MultiplicativeSphericalTwoBodyKernel Kernel>
        class TwoBodyIntGenerator {

        public:
            typedef libint2::TwoBodyEngine<Kernel> Engine;
            typedef mpqc::integrals::EnginePool<Engine> EnginePool;


            TwoBodyIntGenerator() = default;

            TwoBodyIntGenerator(TwoBodyIntGenerator const &) = default;

            TwoBodyIntGenerator &operator=(
                    TwoBodyIntGenerator const &) = default;


            TwoBodyIntGenerator(const std::shared_ptr<EnginePool> &pool,
                                const std::shared_ptr<std::vector<ShellVec>> &cluster_shells,
                                const std::shared_ptr<integrals::Screener> screen)
                    : pool_(pool), shells_(cluster_shells), screener_(screen) { }

            void set_pool(const std::shared_ptr<EnginePool> &pool){
                pool_ = pool;
            }

            void set_shell(const std::shared_ptr<std::vector<ShellVec>> &cluster_shells){
                shells_ = cluster_shells;
            }

            void set_screener(const std::shared_ptr<integrals::Screener> & screen) {
                screener_ = screen;
            }

// return a tile of integral in chemical notation
            // (block1 block2| block3 block4)
            TA::Tensor <double> compute(const TA::Range &r,
                                        const std::vector<std::size_t> &index) {

                // assert size of index is 4
                assert(index.size() == 4);
//                // create the tile
//                TA::Tensor <double> tile(r);
                auto range(r);

//                // get the engine
                auto engine = pool_->local();
                // get the shells in each dimension
                integrals::detail::VecArray<4> array_shells;
                array_shells[0] = &(*shells_)[index[0]];
                array_shells[1] = &(*shells_)[index[1]];
                array_shells[2] = &(*shells_)[index[2]];
                array_shells[3] = &(*shells_)[index[3]];

                auto tile = integrals::detail::integral_kernel(engine,std::move(range),array_shells, *screener_);

                return tile;
            }

        private:
            std::shared_ptr<EnginePool> pool_;
            std::shared_ptr<std::vector<ShellVec>> shells_;
            std::shared_ptr<integrals::Screener> screener_;
        };

    }


    // this part of code is not used and tested!
    // 2015 Septemper Chong Peng
    // IntegralGenerator for two electron AO integrals, use to work with LazyIntegral
    // it used three center integral for generation of integrals
    //
    template <typename Policy>
    class TwoElectronIntDFGenerator{

    public:

        typedef TA::Array <double, 3, TA::Tensor<double>, Policy> TArray3;
        typedef TA::Array <double, 4, TA::Tensor<double>, Policy> TArray4;

        TwoElectronIntDFGenerator() = default;

        TwoElectronIntDFGenerator(TwoElectronIntDFGenerator const &) = default;

        TwoElectronIntDFGenerator(const TArray3& Xpq): Xpq_(Xpq)
        {
            Xlobound_ = Xpq.trange().data().front().tiles().first;
            Xupbound_ = Xpq.trange().data().front().tiles().second;
        }


        TwoElectronIntDFGenerator &operator=(
                TwoElectronIntDFGenerator const &) = default;


        // use chemical notation
        // (pq|rs)
        TA::Tensor<double> compute(const TA::Range &r,
                                   const std::vector<std::size_t> &index){
            // assert size of index is 4
            assert(index.size() == 4);
            // create the tile
            TA::Tensor <double> tile(r);

            //construct lowbound and upbound
            std::vector<std::size_t> bralow;
            std::vector<std::size_t> braup;
            std::vector<std::size_t> ketlow;
            std::vector<std::size_t> ketup;

            bralow.reserve(3);
            ketlow.reserve(3);
            braup.reserve(3);
            ketup.reserve(3);

            // X dimension is all tiles
            bralow[0] = Xlobound_;
            ketlow[0] = Xlobound_;
            braup[0] = Xupbound_;
            ketup[0] = Xupbound_;

            for(auto i = 1; i < 3; ++i){
                bralow[i] = index[i-1];
                ketlow[i] = index[i+1];
                braup[i] = bralow[i] + 1;
                ketlow[i] = ketlow[i] + 1;
            }

            // do the block contraction
            TArray4 ao_block;
            ao_block("p,q,r,s") = Xpq_("X,p,q").block(bralow,braup)*Xpq_("X,r,s").block(ketlow,ketup);

            // get the tile
            auto future_tile = ao_block.find({0,0,0,0});
            tile = future_tile.get();
            return tile;

        }

    private:
        // three center integral
        TArray3 Xpq_;

        // lowbound for X, should be 0
        std::size_t Xlobound_;
        // upbound for X, should be the max
        std::size_t Xupbound_;
    };

}

#endif //TILECLUSTERCHEM_INTEGRAL_GENERATOR_H
