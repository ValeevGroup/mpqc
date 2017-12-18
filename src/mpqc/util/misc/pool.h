
#ifndef MPQC4_SRC_MPQC_UTIL_MISC_POOL_H_
#define MPQC4_SRC_MPQC_UTIL_MISC_POOL_H_

#include <memory>

#include <tbb/enumerable_thread_specific.h>

namespace mpqc {
namespace utility {

/// A pool of thread-specific objects
template <typename Item>
class TSPool {
  public:
    /// Don't allow copies or default initialization.
    TSPool() = delete;
    TSPool(TSPool const &) = delete;
    TSPool &operator=(TSPool const &) = delete;

    TSPool &operator=(TSPool &&) = default;
    TSPool(TSPool &&a) = default;

    /// Initializes the pool with a single @c Item
    explicit TSPool(Item e)
        : item_(std::move(e)), items_(item_) {}
    /// @return reference to the thread-local @c Item instance.
    Item &local() { return items_.local(); }

  private:
    Item item_;
    tbb::enumerable_thread_specific<Item> items_;
};

template <typename Item>
std::shared_ptr<TSPool<Item>> make_pool(Item e) {
    return std::make_shared<TSPool<Item>>(std::move(e));
}

}  // namespace utility
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_INTEGRAL_ENGINE_POOL_H_
