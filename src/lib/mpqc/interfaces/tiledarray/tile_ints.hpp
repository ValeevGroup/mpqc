/*
 * tile_ints.hpp
 *
 *  Created on: Aug 7, 2013
 *      Author: drewlewis
 */

#ifndef TILE_INTS_HPP_
#define TILE_INTS_HPP_

#include <tiled_array.h>
#include <mpqc/interfaces/tiledarray/trange1.hpp>
#include <chemistry/qc/basis/integral.h>
#include <boost/ref.hpp>
#include <mpqc/utility/foreach.hpp>
#include <mpqc/integrals/integrals.hpp>

namespace mpqc{

    template <typename IntegralEngine>
    struct EngineTypeTraits;
    template <> struct EngineTypeTraits<sc::Ref<sc::TwoBodyInt> > {
        static const size_t ncenters = 4u;
    };
    template <> struct EngineTypeTraits<sc::Ref<sc::TwoBodyThreeCenterInt> > {
        static const size_t ncenters = 3u;
    };
    template <> struct EngineTypeTraits<sc::Ref<sc::TwoBodyTwoCenterInt> > {
        static const size_t ncenters = 2u;
    };
    template <> struct EngineTypeTraits<sc::Ref<sc::OneBodyInt> > {
        static const size_t ncenters = 2u;
    };
    template <> struct EngineTypeTraits<sc::Ref<sc::OneBodyOneCenterInt> > {
        static const size_t ncenters = 1u;
    };

    namespace int_details {
        typedef std::vector<int> shell_range;

        template<typename RefEngine>
        void
        get_integrals(const std::vector<shell_range> &s,
                      RefEngine &engine,
                      TensorRef<double, 2> &tile_map){
            mpqc::integrals::evaluate(engine, s[0], s[1], tile_map);
        }

        template<typename RefEngine>
        void
        get_integrals(const std::vector<shell_range> &s,
                      RefEngine &engine,
                      TensorRef<double, 3> &tile_map){
            mpqc::integrals::evaluate(engine, s[0], s[1], s[2], tile_map);
        }

        template<typename RefEngine>
        void
        get_integrals(const std::vector<shell_range> &s,
                      RefEngine &engine,
                      TensorRef<double, 4> &tile_map){
            mpqc::integrals::evaluate(engine, s[0], s[1], s[2], s[3], tile_map);
        }


    } // namespace int_details

    template<typename Tile, typename RefEngine>
    void
    get_integrals(Tile &tile, RefEngine &engine){
        constexpr std::size_t rank = EngineTypeTraits<RefEngine>::ncenters;
        typedef int_details::shell_range shell_range;

        std::size_t t_end = tile.size() - 1;

        std::vector<shell_range> ranges(rank);

        for(auto i = 0; i < rank; ++i){

            std::size_t first_func = tile.range().idx(0)[i];
            std::size_t last_func = tile.range().idx(t_end)[i];
            std::size_t first_shell =
                           engine->basis(i)->function_to_shell(first_func);
            std::size_t last_shell =
                           engine->basis(i)->function_to_shell(last_func);

            for(auto j = first_shell; j <  last_shell + 1; ++j){
                ranges[i].push_back(j);
            }
        }

        std::size_t dim[rank];
        for(auto i = 0; i < rank; ++i){
            dim[i] = tile.range().size()[i];
        }

        TensorRef<double, rank> tile_map(tile.data(), dim);

        int_details::get_integrals(ranges, engine, tile_map);

    }

    template<typename RefPool, typename It, class A>
    void make_integral_task(It first, It last, const A &array,
                       RefPool pool);

    template<typename RefPool, typename It, class A>
    void
    integral_task(It first, It last, A &array,
                  RefPool pool){

        while(last - first > 5){
            It middle = first;
            std::advance(middle, std::distance(first,last)/2);
            make_integral_task(middle, last, array, pool);
            last = middle;
        }
        typedef typename boost::unwrap_reference<RefPool>::type PoolType;

        typename PoolType::engine_type engine = pool.get().instance();

        for(; first != last; ++first){
            typename A::value_type tile(
                            array.trange().make_tile_range(*first));
            get_integrals(tile, engine);

            array.set(*first, tile);
        }
    }

    template<typename RefPool, typename It, class A>
    void
    make_integral_task(It first, It last, const A &array,
                       RefPool pool){
        array.get_world().taskq.add(&integral_task<RefPool, It, A>, first,
                                    last, array, pool);
    }

    template<typename IntEngPool, class A>
    void fill_tiles(A &array, const IntEngPool &pool){
        auto begin = array.get_pmap()->begin();
        auto end = array.get_pmap()->end();

        make_integral_task(begin, end, array, boost::cref(pool));
    }

} // namespace mpqc



#endif /* TILE_INTS_HPP_ */
