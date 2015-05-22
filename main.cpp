#include "common/namespaces.h"
#include "common/typedefs.h"
#include "include/tiledarray.h"
#include "array_ops/deep_filter.h"

#include <utility>

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);

    TA::TiledRange1 tr1(0, 1, 2, 4, 8);
    std::vector<decltype(tr1)> tr1s = {tr1, tr1};
    TA::TiledRange tr(tr1s.begin(), tr1s.end());

    TA::Tensor<float> tile_norms(tr.tiles(), 4.0);

    TA::SparseShape<float> shape(world, tile_norms, tr);

    TA::Array<double, 2, TA::Tensor<double>, TA::SparsePolicy> T(world, tr,
                                                                 shape);
    auto pmap = T.get_pmap();
    const auto end = pmap->end();
    for(auto it = pmap->begin(); it != end; ++it){
        auto range = tr.make_tile_range(*it);
        auto idx = tr.tiles().idx(*it);
        if(idx[0] < 2){
            if(idx[1] < 2){
                T.set(*it, TA::Tensor<double>(range, 1));
            } else {
                T.set(*it, TA::Tensor<double>(range, 2));
            }
        } else {
            if(idx[1] < 2){
                T.set(*it, TA::Tensor<double>(range, 2));
            } else {
                T.set(*it, TA::Tensor<double>(range, 4));
            }
        }
    }
    T.truncate();
    std::cout << "Tensor = \n" << T << "\n" << std::endl;

    auto o_range = std::make_pair(0ul, 2ul);
    std::array<decltype(o_range), 2> oo_range = {{o_range, o_range}};
    auto T_oo = tcc::array_ops::deep_filter(T, oo_range);
    std::cout << "T_oo = \n" << T_oo << "\n" << std::endl;

    auto v_range = std::make_pair(2ul, 8ul);
    std::array<decltype(o_range), 2> ov_range = {{o_range, v_range}};
    auto T_ov = tcc::array_ops::deep_filter(T, ov_range);
    std::cout << "T_ov = \n" << T_ov << "\n" << std::endl;

    std::array<decltype(o_range), 2> vo_range = {{v_range, o_range}};
    auto T_vo = tcc::array_ops::deep_filter(T, vo_range);
    std::cout << "T_vo = \n" << T_vo << "\n" << std::endl;

    std::array<decltype(o_range), 2> vv_range = {{v_range, v_range}};
    auto T_vv = tcc::array_ops::deep_filter(T, vv_range);
    std::cout << "T_vv = \n" << T_vv << "\n" << std::endl;

    madness::finalize();
    return 0;
}
