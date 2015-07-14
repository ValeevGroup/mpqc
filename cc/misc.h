//
// Created by Chong Peng on 7/7/15.
//

#ifndef TILECLUSTERCHEM_MISC_H
#define TILECLUSTERCHEM_MISC_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"

namespace tcc {


  template<typename T, typename Tile, typename  Policy>
  using TArray2 = TA::Array<T, 2, Tile, Policy>;

  template<typename T, typename Tile, typename  Policy>
  using TArray4 = TA::Array<T, 4, Tile, Policy>;

  template<typename T, typename Tile, typename Policy>
  inline size_t size(const std::pair<TArray2<T,Tile,Policy>, TArray4<T,Tile,Policy>> &a) {
    // this is the number of tiles
    if (a.first.size() && a.second.size() > 0) // assuming dense shape
      return a.first.trange().elements().volume() +
             a.second.trange().elements().volume();
    else
      return 0;
  };


  template<typename T, typename Tile, typename Policy>
  inline void zero(const std::pair<TArray2<T,Tile,Policy>, TArray4<T,Tile,Policy>> &a) {
    const std::string var2 = TA::detail::dummy_annotation(2ul);
    const std::string var4 = TA::detail::dummy_annotation(4ul);

    typedef typename TArray2<T,Tile,Policy>::element_type element_type;
    a.first(var2) = element_type(0) * a.first(var2);
    a.second(var4) = element_type(0) * a.second(var4);

  }

  template<typename T, typename Tile, typename Policy>
  inline typename TArray2<T,Tile,Policy>::element_type
  dot_product(const std::pair<TArray2<T,Tile,Policy>, TArray4<T,Tile,Policy>> &a,
              const std::pair<TArray2<T,Tile,Policy>, TArray4<T,Tile,Policy>> &b) {

    const std::string var2 = TA::detail::dummy_annotation(2ul);
    const std::string var4 = TA::detail::dummy_annotation(4ul);

    return a.first(var2).dot(b.first(var2)).get() +
           a.second(var4).dot(b.second(var4)).get();

  };

  template<typename T, typename Tile, typename Policy>
  inline void axpy(const std::pair<TArray2<T,Tile,Policy>, TArray4<T,Tile,Policy>> &y,
                   typename TArray2<T,Tile,Policy>::element_type a,
                   const std::pair<TArray2<T,Tile,Policy>, TArray4<T,Tile,Policy>> &x) {

    const std::string var2 = TA::detail::dummy_annotation(2ul);
    const std::string var4 = TA::detail::dummy_annotation(4ul);

    y.first(var2) = y.first(var2) + a * x.first(var2);
    y.second(var4) = y.second(var4) + a * x.second(var4);
  };

}

#endif //TILECLUSTERCHEM_MISC_H
