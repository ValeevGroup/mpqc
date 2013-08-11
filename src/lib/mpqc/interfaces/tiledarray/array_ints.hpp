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

    template<std::size_t N>
    TA::TiledRange
    make_TiledRange(std::array<tiling::TRange1, N> &blocking){
        return TA::TiledRange(blocking.begin(), blocking.end());
    }

    template <typename IntEngPool>
    TA::Array<double, EngineTypeTraits<typename IntEngPool::engine_type>::ncenters >
    Integrals(madness::World &world, const IntEngPool &pool,
              const TRange1Gen &trange1gen){

        typedef typename IntEngPool::engine_type engine_type;
        const size_t engine_ncenters = EngineTypeTraits<engine_type>::ncenters;
        /*
        static_assert(N == engine_ncenters,
                      "Rank of Tensor does not corresond "
                      "to Rank of this type of integral Engine");
                      */

        std::array<tiling::TRange1, engine_ncenters> blocking =
                        get_blocking<engine_ncenters>(pool, trange1gen);

        TA::Array<double,engine_ncenters> array(world, make_TiledRange(blocking));

        fill_tiles(array, pool);

        return array;
    }


} // namespace mpqc


#endif /* mpqc_interfaces_tiledarray_integrals_hpp */
