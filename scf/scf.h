//
// Created by Chong Peng on 2/9/16.
//

#ifndef MPQC_SCF_SCF_H
#define MPQC_SCF_SCF_H


#include "../common/typedefs.h"
#include "../include/tiledarray.h"
#include "builder.h"

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
    array_type D_;
    array_type C_;
    TiledArray::DIIS<array_type> diis_;

    std::unique_ptr<FockBuilder> builder_;

    std::vector<double> scf_times_;

    int64_t occ_;
    double repulsion_;

    double TcutC_ = 0;

  public:
    ClosedShellSCF() = default;

    template <typename Builder>
    ClosedShellSCF(array_type const &H, array_type const &S, int64_t occ,
                   double rep, Builder builder,
                   array_type const &F_guess = array_type{}, double TcutC = 0)
            : H_(H),
              S_(S),
              builder_(make_unique<Builder>(std::move(builder))),
              occ_(occ),
              repulsion_(rep),
              TcutC_(TcutC) {

        if (F_guess.is_initialized()) {
            F_ = F_guess;
        } else {
            F_ = H_;
        }

        compute_density(occ_);
    }

    inline array_type const &overlap() const { return S_; }
    inline array_type const &fock() const { return F_; }
    inline array_type const &density() const { return D_; }
    inline array_type const &coefficents() const { return C_; }

    double energy();

    /*! Function to compute the density to the desired accuracy.
     *
     * Takes some form of integral and does the scf iterations.  The place to
     *specialized is in build_fock.
     *
     * returns true if the calculation converged to the desired threshold in
     *fewer than max_iters
     */
    bool solve(int64_t max_iters, double thresh);

  private:
    void compute_density(int64_t occ);
    void build_F();
};

} // end of namespace scf
} // end of namespace mpqc


#endif // MPQC_SCF_SCF_H
