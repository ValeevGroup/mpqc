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
class EnginePool {
  public:
    using EngType = E;
    /// Don't allow copies or default initialization.
    EnginePool() = delete;
    EnginePool(EnginePool const &) = delete;
    EnginePool &operator=(EnginePool const &) = delete;

    EnginePool &operator=(EnginePool &&) = default;
    EnginePool(EnginePool &&a) = default;

    /// Initialize class with engine, engines_ needs lambda due to way
    /// tbb::enumerable_thread_specific constructor works.
    explicit EnginePool(E e)
        : engine_(std::move(e)), engines_(engine_) {}

    /// Get reference to thread local engine.
    E &local() { return engines_.local(); }

  private:
    E engine_;
    tbb::enumerable_thread_specific<E> engines_;
};

template <typename E>
EnginePool<E> make_pool(E e) {
    return EnginePool<E>(std::move(e));
}

// Constexpr function to return the order of the integral pools.
template <typename Pool>
constexpr unsigned long pool_order();

template <>
constexpr unsigned long pool_order<EnginePool<libint2::OneBodyEngine>>() {
    return 2ul;
}

template <>
constexpr unsigned long
pool_order<EnginePool<libint2::TwoBodyEngine<libint2::Coulomb>>>() {
    return 4ul;
}

} // namespace integrals
} // namespac tcc
