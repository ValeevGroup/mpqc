#include "mpqc/chemistry/qc/lcao/scf/pbc/util.h"

namespace mpqc {
namespace pbc {
namespace detail {

shellpair_list_t parallel_compute_shellpair_list(madness::World &world,
                                                 const Basis &basis1,
                                                 const Basis &basis2,
                                                 double threshold,
                                                 double target_precision) {
  using ::mpqc::lcao::gaussian::make_engine_pool;
  using ::mpqc::lcao::gaussian::detail::to_libint2_operator;
  // initialize engine
  auto engine_pool = make_engine_pool(
      libint2::Operator::overlap, utility::make_array_of_refs(basis1, basis2),
      libint2::BraKet::x_x);

  std::mutex mx;
  shellpair_list_t result;

  const auto &shv1 = basis1.flattened_shells();
  const auto &shv2 = basis2.flattened_shells();
  const auto nsh1 = shv1.size();
  const auto nsh2 = shv2.size();

  auto compute = [&](int64_t input_s1) {

    auto n1 = shv1[input_s1].size();
    const auto engine_precision = target_precision;
    auto engine = engine_pool->local();
    engine.set_precision(engine_precision);
    const auto &buf = engine.results();

    for (auto s2 = 0l; s2 != nsh2; ++s2) {
      auto on_same_center = (shv1[input_s1].O == shv2[s2].O);
      bool significant = on_same_center;
      if (!on_same_center) {
        auto n2 = shv2[s2].size();
        engine.compute1(shv1[input_s1], shv2[s2]);
        Eigen::Map<const RowMatrixXd> buf_mat(buf[0], n1, n2);
        auto norm = buf_mat.norm();
        significant = (norm >= threshold);
      }

      if (significant) {
        mx.lock();
        result[input_s1].emplace_back(s2);
        mx.unlock();
      }
    }
  };

  for (auto s1 = 0l; s1 != nsh1; ++s1) {
    result.insert(std::make_pair(s1, std::vector<size_t>()));
    world.taskq.add(compute, s1);
  }
  world.gop.fence();

  engine_pool.reset();

  // resort shell list in increasing order
  for (auto s1 = 0l; s1 != nsh1; ++s1) {
    auto &list = result[s1];
    std::sort(list.begin(), list.end());
  }

  return result;
}

}  // namespace detail
}  // namespace pbc
}  // namespace mpqc
