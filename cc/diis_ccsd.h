//
// Created by Chong Peng on 7/7/15.
//

#ifndef TILECLUSTERCHEM_DIIS_CCSD_H
#define TILECLUSTERCHEM_DIIS_CCSD_H

#include <cmath>

#include "../include/tiledarray.h"
#include "../common/namespaces.h"

namespace tcc {

  namespace cc{

    template<typename T, typename Tile, typename  Policy>
    using TArray2 = TA::Array<T, 2, Tile, Policy>;

    template<typename T, typename Tile, typename  Policy>
    using TArray4 = TA::Array<T, 4, Tile, Policy>;


    template <typename T, typename Tile, typename Policy>
    struct T1T2 {

      typedef typename TArray2<T,Tile,Policy>::element_type element_type;

      T1T2(TArray2<T,Tile,Policy> &t1 , TArray4<T,Tile,Policy> &t2)
              : first(t1), second(t2) { }

      TArray2 <T, Tile, Policy>& first;
      TArray4 <T, Tile, Policy>& second;

      double norm(){
        const std::string var2 = TA::detail::dummy_annotation(2ul);
        const std::string var4 = TA::detail::dummy_annotation(4ul);

        auto n1 = first(var2).norm().get();
        auto n2 = second(var4).norm().get();

        return double(std::sqrt(n1*n1+n2*n2));
      }
    };


  template<typename T, typename Tile, typename Policy>
  inline size_t size(const T1T2<T,Tile,Policy> &a) {
    // this is the number of tiles
    if (a.first.size() && a.second.size() > 0) // assuming dense shape
      return a.first.trange().elements().volume() +
             a.second.trange().elements().volume();
    else
      return 0;
  };


  template<typename T, typename Tile, typename Policy>
  inline void zero(T1T2<T,Tile,Policy> &a) {
    const std::string var2 = TA::detail::dummy_annotation(2ul);
    const std::string var4 = TA::detail::dummy_annotation(4ul);

    typedef typename T1T2<T,Tile,Policy>::element_type element_type;
    a.first(var2) = element_type(0) * (a.first(var2));
    a.second(var4) = element_type(0) * (a.second(var4));

  }

  template<typename T, typename Tile, typename Policy>
  inline typename T1T2<T,Tile,Policy>::element_type
  dot_product(const T1T2<T,Tile,Policy> &a,
              const T1T2<T,Tile,Policy> &b) {

    const std::string var2 = TA::detail::dummy_annotation(2ul);
    const std::string var4 = TA::detail::dummy_annotation(4ul);

    return a.first(var2).dot(b.first(var2)).get() +
           a.second(var4).dot(b.second(var4)).get();

  };

  template<typename T, typename Tile, typename Policy>
  inline void axpy(T1T2<T,Tile,Policy>  &y,
                   typename T1T2<T,Tile,Policy>::element_type a,
                   const T1T2<T,Tile,Policy>  &x) {

    const std::string var2 = TA::detail::dummy_annotation(2ul);
    const std::string var4 = TA::detail::dummy_annotation(4ul);

    y.first(var2) = y.first(var2) + a * x.first(var2);
    y.second(var4) = y.second(var4) + a * x.second(var4);
  };

  } // namespace cc
} //namespace tcc

#endif //TILECLUSTERCHEM_DIIS_CCSD_H
