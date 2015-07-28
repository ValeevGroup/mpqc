//
// Created by Chong Peng on 7/27/15.
//

#ifndef TILECLUSTERCHEM_LAZY_INTEGRAL_H
#define TILECLUSTERCHEM_LAZY_INTEGRAL_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"


namespace tcc{
  namespace cc{

    /// Integral tile of a DIM-order TA::Array that's "evaluated" when needed
    // by calling IntegralGenerator.compute(range_, index_)

    template <unsigned int DIM, typename IntegralGenerator>
    class LazyIntegral {

    public:
      typedef double value_type;
      typedef TA::Tensor<double> eval_type;
      typedef TA::Range range_type;

      /// Default constructor
      LazyIntegral() {}

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
                  const std::array<std::size_t, DIM>& index,
                  std::shared_ptr<IntegralGenerator>  integral_generator):
              range_(range), index_(index), integral_generator_(integral_generator)
              { }

      // Convert lazy tile to data tile
      operator TA::Tensor<double>() const {
        return integral_generator_->compute(range_, index_);
      }

      template<typename Archive>
      void serialize(Archive& ar){
        assert(false);
      }

    private:
      range_type range_;
      std::array<std::size_t, DIM> index_;
      std::shared_ptr<IntegralGenerator> integral_generator_;

    };
  }
}



#endif //TILECLUSTERCHEM_LAZY_INTEGRAL_H
