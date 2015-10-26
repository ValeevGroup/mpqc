//
// Created by Chong Peng on 10/14/15.
//

#ifndef TILECLUSTERCHEM_ATOMIC_INTEGRAL_H
#define TILECLUSTERCHEM_ATOMIC_INTEGRAL_H

#include <string>



#include"../common/namespaces.h"
#include "../include/tiledarray.h"
#include "../basis/basis.h"
#include "../expression/formula.h"

namespace mpqc{
namespace integrals{

    //TODO Atomic Integral class
    template<typename Tile, typename Policy>
    class AtomicIntegral {
    public:
        using TArray2 = TA::Array <double, 2, Tile, Policy>;
        using TArray3 = TA::Array <double, 3, Tile, Policy>;
        using TArray4 = TA::Array <double, 4, Tile, Policy>;


        AtomicIntegral(std::shared_ptr<basis::Basis> obs,
                       std::shared_ptr<basis::Basis> dfbs = nullptr,
                       std::shared_ptr<basis::Basis> auxbs = nullptr) :
                obs_(obs), dfbs_(dfbs), auxbs_(auxbs)  { }

    private:

        TArray2 compute_one_electron(const std::wstring& formula);
        TArray4 compute_two_electron(const std::wstring& formula);

        TArray2 compute_one_electron(const Formula& formula);
        TArray4 compute_two_electron(const Formula& formula);

    private:

        std::shared_ptr<basis::Basis> obs_;
        std::shared_ptr<basis::Basis> dfbs_;
        std::shared_ptr<basis::Basis> auxbs_;

    };

    template <typename Tile, typename Policy>
    typename AtomicIntegral<Tile,Policy>::TArray2 AtomicIntegral<Tile,Policy>::compute_one_electron(const std::wstring &formula) {
        Formula tmp(formula);
        return compute_one_electron(std::move(tmp));
    }

    template <typename Tile, typename Policy>
    typename AtomicIntegral<Tile,Policy>::TArray2 AtomicIntegral<Tile,Policy>::compute_one_electron(const Formula& formula) {
    }

    template <typename Tile, typename Policy>
    typename AtomicIntegral<Tile,Policy>::TArray4 AtomicIntegral<Tile,Policy>::compute_two_electron(const std::wstring &formula) {
        Formula tmp(formula);
        return compute_one_electron(std::move(tmp));
    }

    template <typename Tile, typename Policy>
    typename AtomicIntegral<Tile,Policy>::TArray4 AtomicIntegral<Tile,Policy>::compute_two_electron(const Formula& formula) {
    }


    //TODO R12 Integral
    //TODO geminal parametors

    // R12 Atomic Integral

//    template<typename Tile, typename Policy>
//    class R12AtomicIntegral : public AtomicIntegral<Tile,Policy>{
//    public:
//        typedef TA::Array <double, 2, Tile, Policy> TArray2;
//        typedef TA::Array <double, 3, Tile, Policy> TArray3;
//        typedef TA::Array <double, 4, Tile, Policy> TArray4;
//
//    private:
//         geminal parametors
//    };

}
}


#endif //TILECLUSTERCHEM_ATOMIC_INTEGRAL_H
