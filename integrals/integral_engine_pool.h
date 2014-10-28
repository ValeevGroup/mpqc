//
// integral_engine_pool.h
//
// Copyright (C) 2014 Drew Lewis
// Maintainer Drew Lewis
//
// Based on file integral_engine_pool.hpp from mpqc
//

#include "../include/tbb.h"
#include "../include/libint.h"

namespace tcc {
namespace integrals {

template <typename E>
class IntegralEnginePool {
  public:
    using EngType = E;
    /// Don't allow copies or default initialization.
    IntegralEnginePool() = delete;
    IntegralEnginePool(IntegralEnginePool const &) = delete;
    IntegralEnginePool& operator=(IntegralEnginePool const&) = delete;

    IntegralEnginePool& operator=(IntegralEnginePool &&) = default;
    IntegralEnginePool(IntegralEnginePool &&a) = default;

    /// Initialize class with engine, engines_ needs lambda due to way 
    /// tbb::enumerable_thread_specific constructor works. 
    explicit IntegralEnginePool(E e)
        : engine_(std::move(e)),
          engines_(engine_) {}

    /// Get reference to thread local engine. 
    E& local(){ return engines_.local(); }

  private:
    E engine_;
    tbb::enumerable_thread_specific<E> engines_;
};

// Type traits for the engines in the pool
template<typename Pool> struct IntegralPoolTypeTraits {
    constexpr unsigned long order() const;
};
template<>
struct IntegralPoolTypeTraits<libint2::OneBodyEngine>{
    constexpr unsigned long order() const {return 2ul;}
};
template<>
struct IntegralPoolTypeTraits<libint2::TwoBodyEngine<libint2::Coulomb>>{
    constexpr unsigned long order() const {return 4ul;}
};

} // namespace integrals
} // namespac tcc
