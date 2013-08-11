/*
 * tilings.hpp
 *
 *  Created on: Aug 5, 2013
 *      Author: drewlewis
 */

#ifndef mpqc_interface_tiledarray_trange1_hpp
#define mpqc_interface_tiledarray_trange1_hpp

#include <tiled_array.h>
#include <chemistry/qc/basis/basis.h>
#include <vector>

namespace mpqc {


    namespace tiling {

        typedef sc::Ref<sc::GaussianBasisSet> RefBasis;
        typedef TiledArray::TiledRange1 TRange1;

        std::vector<std::size_t> vec_init(std::size_t size_guess){
            std::vector<std::size_t> vec_guess;
            vec_guess.reserve(size_guess); //  Reserve space
            vec_guess.push_back(0); // Tile 0 always starts at 0
            return vec_guess;
        }
        TRange1 tile_by_shell(const RefBasis &basis){

            std::size_t nshell = basis->nshell();
            std::size_t nbasis = basis->nbasis();

            std::vector<std::size_t> tilesizes = vec_init(nshell);


            if(nshell != 0){
                for(std::size_t i = 0; i < nshell; ++i){
                    // Compute one past the last funciton in the shell.
                    // If i is the last shell then one past the last function
                    // on the shell is one past the last function in the basis.
                    std::size_t shell_end = (i == nshell - 1) ? nbasis :
                                    basis->shell_to_function(i + 1);

                    // Push back the last function on the shell.
                    // This way if shell one has 1 function and shell two has 3
                    // functions . The vector will look like vector[0] = |0,1) and
                    // vector[1] = |1,4)
                    tilesizes.push_back(shell_end);
                }
            }
            else {
                std::cout << "ADD SOME EXCEPTION STUFF HERE PRONTO" << std::endl;
            }

            return TRange1(tilesizes.begin(), tilesizes.end());
        }

        TRange1 tile_by_atom(const RefBasis & basis){

            std::size_t ncenters = basis->ncenter();
            std::vector<std::size_t> tilesizes = vec_init(ncenters);

            if(ncenters != 0){
                for(std::size_t i = 0; i < ncenters; ++i){
                    tilesizes.push_back(tilesizes[i] +
                                        basis->nbasis_on_center(i)
                                        );
                }
            }
            else {
                std::cout << "ADD SOME EXCEPTION STUFF HERE PRONTO" << std::endl;
            }

            return TRange1(tilesizes.begin(), tilesizes.end());
        }


    } // namespace tiling

    typedef tiling::TRange1 (*TRange1Gen)(const tiling::RefBasis &);

} // namespace mpqc



#endif /* mpqc_interface_tiledarray_trange1_hpp */
