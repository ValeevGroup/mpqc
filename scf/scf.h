//
// Created by Chong Peng on 2/9/16.
//

#ifndef TILECLUSTERCHEM_SCF_H
#define TILECLUSTERCHEM_SCF_H


#include "../integrals/integrals.h"
#include "../include/tiledarray.h"
#include "../ta_routines/array_to_eigen.h"

namespace mpqc{

class RHF {
public:
    using array_type = TA::TSpArrayD;

    RHF(array_type const &H, array_type const &S, int64_t occ, double rep)
    : H_(H), S_(S), occ_(occ), repulsion_(rep) {
        F_ = H_;
        compute_density(occ_);
    }

    const array_type get_overlap() const {
        return S_;
    }

    const array_type get_fock() const {
        return F_;
    }

    virtual void solve(int64_t max_iters, double thresh, const array_type &eri4);

    double energy() {
        return repulsion_
               + D_("i,j").dot(F_("i,j") + H_("i,j"), D_.get_world()).get();
    }
private:

    template<typename Integral>
    array_type compute_j(Integral const &eri4);

    template<typename Integral>
    array_type compute_k(Integral const &eri4);

    void compute_density(int64_t occ);

    template <typename Integral>
    void form_fock(Integral const &eri4);

protected:
    array_type H_;
    array_type S_;
    array_type F_;
    array_type D_;
    TiledArray::DIIS<array_type> diis_;

    std::vector<double> k_times_;
    std::vector<double> j_times_;
    std::vector<double> scf_times_;

    int64_t occ_;
    double repulsion_;
};



class DFRHF : public RHF {

public:
    using array_type = RHF::array_type ;

    DFRHF(array_type const &H, array_type const &F_guess, array_type const &S,
    array_type const &L_invV, int64_t occ, double rep)
    : RHF(H,S,occ,rep), L_invV_(L_invV){
        F_ = F_guess;
        compute_density(occ_);
    }

    void solve(int64_t max_iters, double thresh, const array_type &eri3) override;

private:
    void compute_density(int64_t occ);

    template <typename Integral>
    void form_fock(Integral const &eri3);

private:
    array_type C_;
    array_type L_invV_;
    std::vector<double> w_times_;
};

}//end of namespace mpqc


#endif //TILECLUSTERCHEM_SCF_H
