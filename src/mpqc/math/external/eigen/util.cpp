#include "mpqc/math/external/eigen/util.h"

namespace mpqc {
namespace detail {

void sort_eigen(VectorZ &eigVal, MatrixZ &eigVec) {
  auto val = eigVal.real();

  // Sort by ascending eigenvalues
  std::vector<std::pair<double, int>> sortedVal;
  sortedVal.reserve(val.size());
  for (auto i = 0; i != val.size(); ++i) {
    auto pair = std::make_pair(val(i), i);
    sortedVal.push_back(pair);
  };
  std::sort(sortedVal.begin(), sortedVal.end());

  // Build sorted eigenvalues and eigenvectors
  VectorZ sortedEigVal(eigVal);
  MatrixZ sortedEigVec(eigVec);
  for (auto i = 0; i != val.size(); ++i) {
    sortedEigVal(i) = eigVal(sortedVal[i].second);
    sortedEigVec.col(i) = eigVec.col(sortedVal[i].second);
  }

  eigVal = sortedEigVal;
  eigVec = sortedEigVec;
}

}  // namespace detail
}  // namespace mpqc

