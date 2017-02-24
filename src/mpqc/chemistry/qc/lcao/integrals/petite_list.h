#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_INTEGRALS_PETITE_LIST_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_INTEGRALS_PETITE_LIST_H_

#include <stdexcept>

namespace mpqc {
namespace lcao {
namespace utility {

class PetiteList {
 public:
  PetiteList() = default;
  virtual ~PetiteList() = default;

  enum class Symmetry {
    i,
    ab = i,
    aa,
    abc = i,
    aab = aa,
    abcd = i,
    aabb,
    aaaa
  };

  virtual bool is_canonical(long idx0) const = 0;
  virtual bool is_canonical(long idx0, long idx1) const = 0;
  virtual bool is_canonical(long idx0, long idx1, long idx2) const = 0;
  virtual bool is_canonical(long idx0, long idx1, long idx2,
                            long idx3) const = 0;

  virtual int64_t multiplicity(long idx0) const = 0;
  virtual int64_t multiplicity(long idx0, long idx1) const = 0;
  virtual int64_t multiplicity(long idx0, long idx1, long idx2) const = 0;
  virtual int64_t multiplicity(long idx0, long idx1, long idx2,
                               long idx3) const = 0;
};

namespace detail {

template <PetiteList::Symmetry TC>
struct TupleCanonicalizer;

template <>
struct TupleCanonicalizer<PetiteList::Symmetry::i> {
  template <typename... Is>
  static bool is_canonical(Is... idxs) {
    return true;
  }
  template <typename... Is>
  static int64_t multiplicity(Is... idxs) {
    return 1;
  }
};

template <>
struct TupleCanonicalizer<PetiteList::Symmetry::aaaa> {
  template <typename I>
  static bool is_canonical(I idx0) {
    return true;
  }
  template <typename I>
  static bool is_canonical(I idx0, I idx1) {
    return idx0 >= idx1;
  }
  template <typename I>
  static bool is_canonical(I idx0, I idx1, I idx2) {
    return is_canonical(idx0, idx1);
  }
  template <typename I>
  static bool is_canonical(I idx0, I idx1, I idx2, I idx3) {
    if (idx0 < idx1) return false;
    if (idx2 < idx3) return false;
    if (idx0 == idx2 && idx1 < idx3) return false;
    return true;
  }
  template <typename I>
  static int64_t multiplicity(I idx0) {
    throw std::logic_error("underfined TupleCanonicalizer::multiplicity");
  }
  template <typename I>
  static int64_t multiplicity(I idx0, I idx1) {
    throw std::logic_error("underfined TupleCanonicalizer::multiplicity");
  }
  template <typename I>
  static int64_t multiplicity(I idx0, I idx1, I idx2) {
    throw std::logic_error("underfined TupleCanonicalizer::multiplicity");
  }
  template <typename I>
  static int64_t multiplicity(I idx0, I idx1, I idx2, I idx3) {
    if (is_canonical(idx0, idx1, idx2, idx3)) {
      return ((idx0 == idx1) ? 1 : 2) * ((idx0 == idx1) ? 1 : 2) *
             ((idx0 == idx2 && idx1 == idx3) ? 1 : 2);
    } else
      return 0;
  }
};

}  // namespace detail

template <PetiteList::Symmetry TC>
class SymmPetiteList : public PetiteList {
 public:
  bool is_canonical(long idx0) const override {
    return detail::TupleCanonicalizer<TC>::is_canonical(idx0);
  }
  bool is_canonical(long idx0, long idx1) const override {
    return detail::TupleCanonicalizer<TC>::is_canonical(idx0, idx1);
  }
  bool is_canonical(long idx0, long idx1, long idx2) const override {
    return detail::TupleCanonicalizer<TC>::is_canonical(idx0, idx1, idx2);
  }
  bool is_canonical(long idx0, long idx1, long idx2, long idx3) const override {
    return detail::TupleCanonicalizer<TC>::is_canonical(idx0, idx1, idx2, idx3);
  }
  int64_t multiplicity(long idx0) const override {
    return detail::TupleCanonicalizer<TC>::multiplicity(idx0);
  }
  int64_t multiplicity(long idx0, long idx1) const override {
    return detail::TupleCanonicalizer<TC>::multiplicity(idx0, idx1);
  }
  int64_t multiplicity(long idx0, long idx1, long idx2) const override {
    return detail::TupleCanonicalizer<TC>::multiplicity(idx0, idx1, idx2);
  }
  int64_t multiplicity(long idx0, long idx1, long idx2,
                       long idx3) const override {
    return detail::TupleCanonicalizer<TC>::multiplicity(idx0, idx1, idx2, idx3);
  }
};

}  // namespace utility
}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_INTEGRALS_PETITE_LIST_H_
