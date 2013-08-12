/*
 * integrals.hpp
 *
 *  Created on: Aug 5, 2013
 *      Author: drewlewis
 */

#ifndef mpqc_interfaces_tiledarray_integrals_hpp
#define mpqc_interfaces_tiledarray_integrals_hpp

#include "tile_ints.hpp"


namespace TA = TiledArray;
namespace mpqc {

    template<std::size_t N, typename IntEngPool>
    std::array<tiling::TRange1, N>
    get_blocking(const IntEngPool &pool, const TRange1Gen &trange1gen){

        std::array<tiling::TRange1, N> blocking;

        for(auto i = 0; i < N; ++i){
            blocking[i] = trange1gen(pool.instance()->basis(i));
        }

        return blocking;
    }

    template <typename IntEngPool>
    TA::Array<double, EngineTypeTraits<typename IntEngPool::engine_type>::ncenters >
    Integrals(madness::World &world, const IntEngPool &pool,
              const TRange1Gen &trange1gen){

        typedef typename IntEngPool::engine_type engine_type;
        const size_t rank = EngineTypeTraits<engine_type>::ncenters;

        std::array<tiling::TRange1, rank> blocking =
                        get_blocking<rank>(pool, trange1gen);

        TA::TiledRange trange(blocking.begin(), blocking.end());

        TA::Array<double, rank> array(world, trange);

        fill_tiles(array, pool);

        return array;
    }


} // namespace mpqc


#endif /* mpqc_interfaces_tiledarray_integrals_hpp */
