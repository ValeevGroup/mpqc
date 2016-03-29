//
// Created by Chong Peng on 2/9/16.
//

#ifndef MPQC_SCF_SCF_H
#define MPQC_SCF_SCF_H


#include "../common/typedefs.h"
#include "../include/tiledarray.h"
#include "builder.h"
#include "density_builder.h"
#include "../utility/json_handling.h"

#include <memory>

namespace mpqc {
namespace scf {

class ClosedShellSCF {
  public:
    using array_type = TA::TSpArrayD;

  protected:
    array_type H_;
    array_type S_;
    array_type F_;
    array_type F_diis_;
    array_type D_;
    array_type C_;
    TiledArray::DIIS<array_type> diis_;

    std::unique_ptr<FockBuilder> f_builder_;
    std::unique_ptr<DensityBuilder> d_builder_;

    std::vector<double> scf_times_;
    std::vector<double> d_times_;
    std::vector<double> build_times_;

    double repulsion_;

  public:
    ClosedShellSCF() = default;

    template <typename FBuilder, typename DBuilder>
    ClosedShellSCF(array_type const &H, array_type const &S, double rep,
                   FBuilder f_builder, DBuilder d_builder,
                   array_type const &F_guess = array_type{})
            : H_(H),
              S_(S),
              f_builder_(make_unique<FBuilder>(std::move(f_builder))),
              d_builder_(make_unique<DBuilder>(std::move(d_builder))),
              repulsion_(rep) {

        if (F_guess.is_initialized()) {
            F_ = F_guess;
        } else {
            F_ = H_;
        }

        F_diis_ = F_;
        compute_density();
    }

    ClosedShellSCF(array_type const &H, array_type const &S, double rep,
                   std::unique_ptr<FockBuilder> &&f_builder,
                   std::unique_ptr<DensityBuilder> &&d_builder,
                   array_type const &F_guess = array_type{})
            : H_(H),
              S_(S),
              f_builder_(std::move(f_builder)),
              d_builder_(std::move(d_builder)),
              repulsion_(rep) {

        if (F_guess.is_initialized()) {
            F_ = F_guess;
        } else {
            F_ = H_;
        }

        F_diis_ = F_;
        compute_density();
    }

    inline array_type const &overlap() const { return S_; }
    inline array_type const &fock() const { return F_; }
    inline array_type const &density() const { return D_; }
    inline array_type const &coefficents() const { return C_; }

    double energy() const;

    /*! Function to compute the density to the desired accuracy.
     *
     * Takes some form of integral and does the scf iterations.  The place to
     *specialized is in build_fock.
     *
     * returns true if the calculation converged to the desired threshold in
     *fewer than max_iters
     */
    bool solve(int64_t max_iters, double thresh);

    virtual rapidjson::Value results(rapidjson::Document &) const;

  private:
    void compute_density();
    void build_F();
};

} // end of namespace scf
} // end of namespace mpqc


#endif // MPQC_SCF_SCF_H
