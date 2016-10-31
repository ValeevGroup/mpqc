#pragma once
#ifndef MPQC_INTEGRALS_INTEGRALENGINEPOOL_H
#define MPQC_INTEGRALS_INTEGRALENGINEPOOL_H
//
// integral_engine_pool.h
//
// Copyright (C) 2014 Drew Lewis
// Maintainer Drew Lewis
//
// Based on file integral_engine_pool.hpp from mpqc
//

#include <tbb/enumerable_thread_specific.h>

#include <libint2/engine.h>

namespace mpqc {
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
std::shared_ptr<EnginePool<E>> make_pool(E e) {
    return std::make_shared<EnginePool<E>>(std::move(e));
}

} // namespace integrals
} // namespac mpqc

#endif /* end of include guard: MPQC_INTEGRALS_INTEGRALENGINEPOOL_H */
