//
// Created by Chong Peng on 10/1/15.
//

#ifndef TILECLUSTERCHEM_UTILITY_H_H
#define TILECLUSTERCHEM_UTILITY_H_H

#include <cmath>

namespace mpqc{
    namespace  cc{

        // reblock based on blocksize
        std::vector<std::vector<libint2::Shell>>
                reblock_basis(std::vector<libint2::Shell> shells, std::size_t blocksize){

            std::vector<std::vector<libint2::Shell>> result;

            std::vector<libint2::Shell> tmp;
            std::size_t tmp_size = 0;
            for (auto& shell : shells){

                std::size_t shell_size = shell.size();
                tmp_size += shell_size;

                // if current size is not greater than 2/3 of blocksize
                if(3*shell_size < 5*blocksize){
                    tmp.push_back(shell);
                }else{
                    result.push_back(tmp);
                    tmp = std::vector<libint2::Shell>();
                    tmp.push_back(shell);
                    tmp_size = shell_size;
                }
            }
            return result;
        }

        // compute the min and max block size in TiledRange1
        std::pair<std::size_t, std::size_t>
                minmax_blocksize(TiledArray::TiledRange1 tr1){

            std::vector<std::size_t> block_sizes;
            for (auto block = tr1.begin(); block != tr1.end(); ++block){
                auto block_size = block->second - block->first;
                block_sizes.push_back(block_size);
            }
            auto minmax_block_size = std::minmax_element(block_sizes.begin(),block_sizes.end());

            std::size_t min_block_size = *(minmax_block_size.first);
            std::size_t max_block_size = *(minmax_block_size.second);
            auto result = std::make_pair(min_block_size,max_block_size);
            return result;
        };

        // create matrix d("a,b,i,j) = 1/(ei + ej - ea - eb)
        template<typename Tile, typename Policy>
        void create_d_abij(TA::Array<double, 4, Tile, Policy> &abij,
                           const Eigen::VectorXd &ens, std::size_t n_occ) {
            typedef typename TA::Array<double, 4, Tile, Policy>::range_type range_type;
            typedef typename TA::Array<double, 4, Tile, Policy>::iterator iterator;

            auto make_tile = [&ens, n_occ](range_type &range) {

                auto result_tile = Tile(range);

                // compute index
                const auto a0 = result_tile.range().lobound()[0];
                const auto an = result_tile.range().upbound()[0];
                const auto b0 = result_tile.range().lobound()[1];
                const auto bn = result_tile.range().upbound()[1];
                const auto i0 = result_tile.range().lobound()[2];
                const auto in = result_tile.range().upbound()[2];
                const auto j0 = result_tile.range().lobound()[3];
                const auto jn = result_tile.range().upbound()[3];

                auto tile_idx = 0;
                typename Tile::value_type tmp = 1.0;
                for (auto a = a0; a < an; ++a) {
                    const auto e_a = ens[a + n_occ];
                    for (auto b = b0; b < bn; ++b) {
                        const auto e_b = ens[b + n_occ];
                        for (auto i = i0; i < in; ++i) {
                            const auto e_i = ens[i];
                            for (auto j = j0; j < jn; ++j, ++tile_idx) {
                                const auto e_j = ens[j];
                                const auto e_iajb = e_i + e_j - e_a - e_b;
                                const auto result_abij = tmp / (e_iajb);
                                result_tile[tile_idx] = result_abij;
                            }
                        }
                    }
                }
                return result_tile;
            };

            for (iterator it = abij.begin(); it != abij.end(); ++it) {

                madness::Future<Tile> tile = abij.get_world().taskq.add(
                        make_tile,
                        abij.trange().make_tile_range(it.ordinal()));

                *it = tile;
            }

        }


        // create matrix d("a,i") = 1/(ei - ea)
        template<typename Tile, typename Policy>
        void create_d_ai(TA::Array<double, 2, Tile, Policy> &f_ai, const Eigen::VectorXd &ens, int n_occ) {
            typedef typename TA::Array <double, 2, Tile, Policy>::range_type range_type;
            typedef typename TA::Array <double, 2, Tile, Policy>::iterator iterator;

            auto make_tile = [&ens, n_occ](range_type &range) {

                auto result_tile = Tile(range);
                const auto a0 = result_tile.range().lobound()[0];
                const auto an = result_tile.range().upbound()[0];
                const auto i0 = result_tile.range().lobound()[1];
                const auto in = result_tile.range().upbound()[1];

                auto ai = 0;
                typename Tile::value_type tmp = 1.0;
                for (auto a = a0; a < an; ++a) {
                    const auto e_a = ens[a + n_occ];
                    for (auto i = i0; i < in; ++i, ++ai) {
                        const auto e_i = ens[i];
                        const auto e_ia = e_i - e_a;
                        const auto result_ai = tmp / (e_ia);
                        result_tile[ai] = result_ai;
                    }
                }
                return result_tile;
            };

            for (iterator it = f_ai.begin(); it != f_ai.end(); ++it) {

                madness::Future<Tile> tile = f_ai.get_world().taskq.add(
                        make_tile,
                        f_ai.trange().make_tile_range(it.ordinal()));

                *it = tile;
            }

        }

        // print progress
        void print_progress(int lowprogress, int upprogress, int total){
            int divide = 10;
            if(total < 10){
                divide = total;
                int percent = 100*double(upprogress)/divide;
                std::cout << percent << "% done." << std::endl;
            }else{
                int increase = std::round(double(total)/divide);
                std::vector<int> progress_points;
                for(int i = 0; i < total; i += increase){
                   progress_points.push_back(i);
                }

                for (int i = lowprogress; i < upprogress; i++){
                    if (std::find(progress_points.begin(),progress_points.end(),i) != progress_points.end()){
                        int percent = 100*double(i)/total;
                        std::cout << percent << "% done." << std::endl;
                    }
                }
            }

        }
    }
}


#endif //TILECLUSTERCHEM_UTILITY_H_H
