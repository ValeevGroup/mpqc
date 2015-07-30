//
// Created by Chong Peng on 7/27/15.
//

#ifndef TILECLUSTERCHEM_LAZY_INTEGRAL_H
#define TILECLUSTERCHEM_LAZY_INTEGRAL_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include "lazy_integral.h"

static auto DIRECTAOTWOELECTONINTEGRAL = std::make_shared<tcc::cc::TwoBodyIntGenerator<libint2::Coulomb>>();

namespace tcc{
  namespace cc{

    // lazy integral interface
    // Integral tile of a DIM-order TA::Array that's "evaluated" when needed
    // by calling IntegralGenerator.compute(TA::Range range_, std::vector<std::size_t> index_)
    // check integral_generator.h TwoBodyIntGenerator class for an example of IntegralGenerator template

    template <unsigned int DIM, typename IntegralGenerator>
    class LazyIntegral {

    public:
      typedef double value_type;
      typedef TA::Tensor<double> eval_type;
      typedef TA::Range range_type;

      /// Default constructor
      LazyIntegral() = default;

      /// Copy constructor
      LazyIntegral(const LazyIntegral& other) :
              range_(other.range_), index_(other.index_), integral_generator_(other.integral_generator_)

      { }

      /// Assignment operator
      LazyIntegral& operator= (const LazyIntegral& other)
      {
        range_ = other.range_;
        index_ = other.index_;
        integral_generator_ = other.integral_generator_;
        return *this;
      }

      /// Constructor
      LazyIntegral(range_type range,
                  const std::vector<std::size_t>& index,
                  std::shared_ptr<IntegralGenerator>  integral_generator):
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
      serialize(Archive& ar){
          ar &range_;
          ar &index_;
      }

        template<typename Archive>
        typename std::enable_if<madness::archive::is_input_archive<Archive>::value>::type
        serialize(Archive& ar){
            ar &range_;
            ar &index_;
            integral_generator_ = DIRECTAOTWOELECTONINTEGRAL;
        }
    private:
      range_type range_;
      std::vector<std::size_t> index_;
      std::shared_ptr<IntegralGenerator> integral_generator_;

    };

  // direct two electron integral
  typedef tcc::cc::LazyIntegral<4, TwoBodyIntGenerator<libint2::Coulomb>> LazyTwoElectronTile;
  typedef TA::Array<double, 4, LazyTwoElectronTile, TA::DensePolicy> DirectTwoElectronDenseArray;

    // function to make direct two electron dense TArray
    DirectTwoElectronDenseArray make_lazy_two_electron_array(
            madness::World &world, const tcc::basis::Basis &basis,
            const TA::TiledRange &trange) {

      auto p_cluster_shells = std::make_shared<std::vector<tcc::basis::ClusterShells>>(
              basis.cluster_shells());

      auto two_body_coulomb_engine = tcc::integrals::make_2body(basis);

      auto p_engine_pool = std::make_shared<tcc::integrals::EnginePool<libint2::TwoBodyEngine<libint2::Coulomb>>>(
              two_body_coulomb_engine);

        DIRECTAOTWOELECTONINTEGRAL->set_pool(p_engine_pool);
        DIRECTAOTWOELECTONINTEGRAL->set_shell(p_cluster_shells);

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

      // direct ao integral using DF three center integral
      typedef tcc::cc::LazyIntegral<4, TwoElectronIntDFGenerator<TA::DensePolicy>> LazyTwoElectronDFDenseTile;
      typedef TA::Array<double, 4, LazyTwoElectronDFDenseTile, TA::DensePolicy> DirectTwoElectronDFDenseArray;
  } // namespace cc
} // namespace tcc



#endif //TILECLUSTERCHEM_LAZY_INTEGRAL_H
