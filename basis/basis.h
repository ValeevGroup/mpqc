#pragma once
#ifndef MPQC_BASIS_BASIS_H
#define MPQC_BASIS_BASIS_H

#include "basis_fwd.h"
#include "../molecule/molecule_fwd.h"

#include <vector>
#include <iosfwd>
#include <memory>

#include "../common/typedefs.h"

#include "../include/tiledarray.h"
#include "../include/libint.h"


namespace mpqc {
namespace basis {

class Basis {
  public:
    Basis();
    ~Basis();
    Basis(Basis const &);
    Basis(Basis &&);
    Basis &operator=(Basis const &);
    Basis &operator=(Basis &&);

    Basis(std::vector<ShellVec> cs);

    Basis join(const Basis &basis);

    std::vector<ShellVec> const &cluster_shells() const;

    TiledArray::TiledRange1 create_trange1() const;

    int64_t max_nprim() const;
    int64_t max_am() const;
    int64_t nfunctions() const;
    int64_t nshells() const;
    int64_t nclusters() const { return shells_.size(); };
    std::vector<Shell> flattened_shells() const;

    template <typename Archive>
    typename std::
          enable_if<madness::archive::is_output_archive<Archive>::value>::type
          serialize(Archive &ar) {
        auto nvecs = shells_.size();
        ar &nvecs;
        for(auto const &v : shells_){
            ar &v;
        }
    }

    template <typename Archive>
    typename std::
          enable_if<madness::archive::is_input_archive<Archive>::value>::type
          serialize(Archive &ar) {
        auto nvecs = 0;
        ar & nvecs;

        for(auto i = 0; i < nvecs; ++i){
            ShellVec tmp;
            ar & tmp;
            shells_.emplace_back(std::move(tmp));
        }
    }

  private:
    std::vector<ShellVec> shells_;
};

std::ostream &operator<<(std::ostream &, Basis const &);

template <typename Archive>
void serialize(Archive &ar, libint2::Shell const &s){
    ar & s.alpha;
    ar & s.contr;
    ar & s.O;
    ar & s.max_ln_coeff;
}

template <typename Archive>
void serialize(Archive &ar, libint2::Shell::Contraction const &c){
    ar & c.l;
    ar & c.pure;
    ar & c.coeff;
}

/*! \brief reblock allows for reblocking a basis
 *
 * \warning If reblocking a basis with the intent to use it with tensors
 * computed with the old basis you must be careful not to reorder the shells.
 *
 * \param op should be a function that takes a std::vector<Shell> and returns
 * a std::vector<std::vector<Shell>> for use in initializing a Basis.
 */
template <typename Op, typename... Args>
Basis reblock(Basis const &basis, Op op, Args... args) {
    return Basis(op(basis.flattened_shells(), args...));
}

} // namespace basis
} // namespace mpqc
#endif /* end of include guard: MPQC_BASIS_BASIS_H */
