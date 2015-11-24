//
// Created by Chong Peng on 7/27/15.
//

#ifndef TILECLUSTERCHEM_LAZY_TILE_H
#define TILECLUSTERCHEM_LAZY_TILE_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include "integral_generator.h"

static auto DIRECTAOTWOELECTONINTEGRAL = std::make_shared<mpqc::cc::TwoBodyIntGenerator<libint2::Coulomb>>();

namespace mpqc {
    namespace cc {

        // lazy integral interface
        // Integral tile of a DIM-order TA::Array that's "evaluated" when needed
        // by calling IntegralGenerator.compute(TA::Range range_, std::vector<std::size_t> index_)
        // check integral_generator.h TwoBodyIntGenerator class for an example of IntegralGenerator template

        template<unsigned int DIM, typename IntegralGenerator>
        class LazyTile {

        public:
            typedef double value_type;
            typedef TA::Tensor<double> eval_type;
            typedef TA::Range range_type;

            /// Default constructor
            LazyTile() = default;

            /// Copy constructor
            LazyTile(const LazyTile &other) :
                    range_(other.range_), index_(other.index_), integral_generator_(other.integral_generator_) { }

            /// Assignment operator
            LazyTile &operator=(const LazyTile &other) {
                if (this == &other){
                    return *this;
                }
                range_ = other.range_;
                index_ = other.index_;
                integral_generator_ = other.integral_generator_;
                return *this;
            }

            /// Constructor
            LazyTile(range_type range,
                     const std::vector<std::size_t> &index,
                     std::shared_ptr<IntegralGenerator> integral_generator) :
                    range_(range), index_(index), integral_generator_(integral_generator)
            {
                assert(index.size() == DIM);
            }

            // Convert lazy tile to data tile
            operator TA::Tensor<double>() const {
                return integral_generator_->compute(range_, index_);
            }

            template<typename Archive>
            typename std::enable_if<madness::archive::is_output_archive<Archive>::value>::type
            serialize(Archive &ar) {
                ar & range_;
                ar & index_;
            }

            template<typename Archive>
            typename std::enable_if<madness::archive::is_input_archive<Archive>::value>::type
            serialize(Archive &ar) {
                ar & range_;
                ar & index_;
                integral_generator_ = DIRECTAOTWOELECTONINTEGRAL;
            }

        private:
            range_type range_;
            std::vector<std::size_t> index_;
            std::shared_ptr<IntegralGenerator> integral_generator_;

        };

        // direct two electron integral
        typedef mpqc::cc::LazyTile<4, TwoBodyIntGenerator < libint2::Coulomb>> LazyTwoElectronTile;
        typedef TA::Array<double, 4, LazyTwoElectronTile, TA::DensePolicy> DirectTwoElectronDenseArray;
        typedef TA::Array<double, 4, LazyTwoElectronTile, TA::SparsePolicy> DirectTwoElectronSparseArray;

        // function to make direct two electron dense TArray
        DirectTwoElectronDenseArray make_lazy_two_electron_dense_array(
                madness::World &world, mpqc::basis::Basis &basis,
                const TA::TiledRange &trange, const int screen_option) {

            // compute cluster shell
            auto cluster_shells = basis.cluster_shells();
            auto p_cluster_shells = std::make_shared<std::vector<ShellVec>>(cluster_shells);

            // compute engine pull
            auto two_body_coulomb_engine = mpqc::integrals::make_2body(basis);
            auto p_engine_pool = std::make_shared<mpqc::integrals::EnginePool<libint2::TwoBodyEngine<libint2::Coulomb>>>(
                    two_body_coulomb_engine);

            // compute screener
            std::shared_ptr<integrals::Screener> p_screen;
            if(screen_option == 1){
                auto screen_builder = integrals::init_schwarz_screen(1e-10);
                p_screen = std::make_shared<integrals::Screener>(screen_builder(world, p_engine_pool, basis));

                if(world.rank() == 0){
                    std::cout << "schwarz screen" << std::endl;
                }
            }
            else if(screen_option == 2){
                auto screen_builder = integrals::init_qqr_screen{};
                p_screen = std::make_shared<integrals::Screener>(screen_builder(world, p_engine_pool, basis));
                if(world.rank() == 0){
                    std::cout << "qqr screen" << std::endl;
                }
            }
            else{
                p_screen = std::make_shared<integrals::Screener>(integrals::Screener{});
                if(world.rank() == 0){
                    std::cout << "no screen" << std::endl;
                }
            }

            DIRECTAOTWOELECTONINTEGRAL->set_pool(p_engine_pool);
            DIRECTAOTWOELECTONINTEGRAL->set_shell(p_cluster_shells);
            DIRECTAOTWOELECTONINTEGRAL->set_screener(p_screen);

            DirectTwoElectronDenseArray lazy_two_electron(world, trange);

            // set the functor of tile
            DirectTwoElectronDenseArray::iterator it = lazy_two_electron.begin();
            DirectTwoElectronDenseArray::iterator end = lazy_two_electron.end();
            for (; it != end; ++it) {
                TA::Range range = lazy_two_electron.trange().make_tile_range(
                        it.ordinal());
                auto index = it.index();
                lazy_two_electron.set(index, LazyTwoElectronTile(range, index,
                                                                 DIRECTAOTWOELECTONINTEGRAL));

            }
            return lazy_two_electron;
        }


        // TODO update and test sparse direct array
        // function to make direct two electron Sparse TArray
        DirectTwoElectronSparseArray make_lazy_two_electron_sparse_array(
                madness::World &world, const mpqc::basis::Basis &basis,
                const TA::TiledRange &trange, const int screen_option) {

            // compute cluster shell
            auto cluster_shells = basis.cluster_shells();
            auto p_cluster_shells = std::make_shared<std::vector<ShellVec>>(cluster_shells);

            // compute engine pull
            auto two_body_coulomb_engine = mpqc::integrals::make_2body(basis);
            auto p_engine_pool = std::make_shared<mpqc::integrals::EnginePool<libint2::TwoBodyEngine<libint2::Coulomb>>>(
                    two_body_coulomb_engine);

            // compute screener
            std::shared_ptr<integrals::Screener> p_screen;
            if(screen_option == 1){
                auto screen_builder = integrals::init_schwarz_screen(1e-10);
                p_screen = std::make_shared<integrals::Screener>(screen_builder(world, p_engine_pool, basis));

                if(world.rank() == 0){
                    std::cout << "schwarz screen" << std::endl;
                }
            }
            else if(screen_option == 2){
                auto screen_builder = integrals::init_qqr_screen{};
                p_screen = std::make_shared<integrals::Screener>(screen_builder(world, p_engine_pool, basis));
                if(world.rank() == 0){
                    std::cout << "qqr screen" << std::endl;
                }
            }
            else{
                p_screen = std::make_shared<integrals::Screener>(integrals::Screener{});
                if(world.rank() == 0){
                    std::cout << "no screen" << std::endl;
                }
            }
            DIRECTAOTWOELECTONINTEGRAL->set_pool(p_engine_pool);
            DIRECTAOTWOELECTONINTEGRAL->set_shell(p_cluster_shells);
            DIRECTAOTWOELECTONINTEGRAL->set_screener(p_screen);

            // make shape
            TA::TensorF tile_norms(trange.tiles(),0.0);
            auto t_volume = trange.tiles().volume();
            auto pmap = TA::SparsePolicy::default_pmap(world,t_volume);

            auto compute_tile = [=](int64_t ord, TA::Range range, TA::TensorF* ptr_tile_norm){
                auto index = trange.tiles().idx(ord);
                auto ta_tile = DIRECTAOTWOELECTONINTEGRAL->compute(range,index);

                const auto tile_volume = ta_tile.range().volume();
                const auto tile_norm = ta_tile.norm();

                if(tile_norm >= tile_volume*TA::SparseShape<float>::threshold()){
                    (*ptr_tile_norm)[ord] = tile_norm;
                }

            };
            // need to work with this later
            for(auto const &ord: *pmap){
                auto range = trange.make_tile_range(ord);
                world.taskq.add(compute_tile,ord,range,&tile_norms);
            }
            world.gop.fence();

            TA::SparseShape<float> shape(world, tile_norms,trange);

            DirectTwoElectronSparseArray lazy_two_electron(world, trange, shape, pmap);

            // set the functor of tile
            DirectTwoElectronSparseArray::iterator it = lazy_two_electron.begin();
            DirectTwoElectronSparseArray::iterator end = lazy_two_electron.end();
            for (; it != end; ++it) {
                TA::Range range = lazy_two_electron.trange().make_tile_range(
                        it.ordinal());
                auto index = it.index();
                lazy_two_electron.set(index, LazyTwoElectronTile(range, index,
                                                                 DIRECTAOTWOELECTONINTEGRAL));

            }
            return lazy_two_electron;
        }

        // direct ao integral using DF three center integral
        typedef mpqc::cc::LazyTile<4, TwoElectronIntDFGenerator < TA::DensePolicy>> LazyTwoElectronDFDenseTile;
        typedef TA::Array<double, 4, LazyTwoElectronDFDenseTile, TA::DensePolicy> DirectTwoElectronDFDenseArray;
    } // namespace cc
} // namespace mpqc



#endif //TILECLUSTERCHEM_LAZY_INTEGRAL_H
