#include "../common/typedefs.h"

#include "../include/tiledarray.h"
#include "../include/libint.h"

#include "../utility/time.h"

#include "../integrals/task_integrals_helper.h"

using namespace mpqc;
namespace intd = integrals::detail;

namespace tcc_time = tcc::utility::time;

using coulomb_engine = libint2::TwoBodyEngine<libint2::Coulomb>;
std::vector<libint2::Atom> read_geometry(const std::string &filename);

using ShellVec = std::vector<libint2::Shell>;

TA::TensorD
integral_kernel_meh(coulomb_engine &eng, TA::Range &&rng, ShellVec const &sh0,
                    ShellVec const &sh1, ShellVec const &sh2);

TA::TensorD integral_kernel_better(coulomb_engine &eng, TA::Range &&rng,
                                   ShellVec const &sh0, ShellVec const &sh1,
                                   ShellVec const &sh2);
TA::TensorD
integral_kernel_best(coulomb_engine &eng, TA::Range &&rng, ShellVec const &sh0,
                     ShellVec const &sh1, ShellVec const &sh2);

TA::Range
make_range(libint2::BasisSet const &dfbs, libint2::BasisSet const &obs) {
    return TA::Range(dfbs.nbf(), obs.nbf(), obs.nbf());
}


int main(int argc, char **argv) {
    const auto filename
          = (argc > 1)
                  ? argv[1]
                  : throw std::invalid_argument("There was not a file name");
    const std::string basisname
          = (argc > 2)
                  ? argv[2]
                  : throw std::invalid_argument("There was not a basis name");
    const std::string dfbs_basisname
          = (argc > 3)
                  ? argv[3]
                  : throw std::invalid_argument("There was not a basis name");

    std::vector<libint2::Atom> atoms = read_geometry(filename);
    libint2::BasisSet obs(basisname, atoms);
    libint2::BasisSet dfbs(dfbs_basisname, atoms);

    libint2::init();
    auto max_nprim = std::max(obs.max_nprim(), dfbs.max_nprim());
    auto max_l = std::max(obs.max_l(), dfbs.max_l());
    coulomb_engine engine(max_nprim, max_l, 0);

    auto meh0 = tcc_time::now();
    double norm = 0.0;
    for (auto i = 0; i < 10; ++i) {
        auto tensor = integral_kernel_meh(engine, make_range(dfbs, obs),
                                              dfbs, obs, obs);
        norm += tensor.norm();
    }
    auto meh1 = tcc_time::now();
    auto meh_time = tcc_time::duration_in_s(meh0, meh1) / 10;
    std::cout << "Meh time: " << meh_time << ", norm: " << norm / 10
              << std::endl;

    auto better0 = tcc_time::now();
    norm = 0.0;
    for (auto i = 0; i < 10; ++i) {
        auto tensor = integral_kernel_better(engine, make_range(dfbs, obs),
                                              dfbs, obs, obs);
        norm += tensor.norm();
    }
    auto better1 = tcc_time::now();
    auto better_time = tcc_time::duration_in_s(better0, better1) / 10;
    std::cout << "Better time: " << better_time << ", norm: " << norm / 10
              << std::endl;

    auto best0 = tcc_time::now();
    norm = 0.0;
    for (auto i = 0; i < 10; ++i) {
        auto tensor = integral_kernel_best(engine, make_range(dfbs, obs),
                                              dfbs, obs, obs);
        norm += tensor.norm();
    }
    auto best1 = tcc_time::now();
    auto best_time = tcc_time::duration_in_s(best0, best1) / 10;
    std::cout << "Best time: " << best_time << ", norm: " << norm / 10
              << std::endl;

    libint2::cleanup();

    return 0;
}

std::vector<libint2::Atom> read_geometry(const std::string &filename) {

    std::ifstream is(filename);
    assert(is.good());

    std::ostringstream oss;
    oss << is.rdbuf();
    std::istringstream iss(oss.str());

    if (filename.rfind(".xyz") != std::string::npos)
        return libint2::read_dotxyz(iss);
    else
        throw "only .xyz files are accepted";
}

TA::TensorD
integral_kernel_meh(coulomb_engine &eng, TA::Range &&rng, ShellVec const &sh0,
                    ShellVec const &sh1, ShellVec const &sh2) {

    std::size_t bf0, bf1, bf2 = 0;
    auto lo0 = rng.lobound_data()[0];
    auto lo1 = rng.lobound_data()[1];
    auto lo2 = rng.lobound_data()[2];

    auto tile = TA::TensorD(std::move(rng));

    // compute
    bf0 = lo0;
    for (auto s0 : sh0) {

        std::size_t ns0 = s0.size();
        bf1 = lo1;

        for (auto s1 : sh1) {

            std::size_t ns1 = s1.size();
            bf2 = lo2;

            for (auto s2 : sh2) {
                std::size_t ns2 = s2.size();

                const auto *buf
                      = eng.compute(s0, libint2::Shell::unit(), s1, s2);

                auto lowbound = {bf0, bf1, bf2};
                auto upbound = {bf0 + ns0, bf1 + ns1, bf2 + ns2};
                auto view = tile.block(lowbound, upbound);
                auto map = TA::make_map(buf, lowbound, upbound);
                view = map;

                bf2 += ns2;
            }
            bf1 += ns1;
        }
        bf0 += ns0;
    }

    return tile;
}

TA::TensorD integral_kernel_better(coulomb_engine &eng, TA::Range &&rng,
                                   ShellVec const &sh0, ShellVec const &sh1,
                                   ShellVec const &sh2) {

    auto const &lobound = rng.lobound();
    unsigned long lb [] = {lobound[0], lobound[1], lobound[2]};
    unsigned long ub [3];

    for(auto i = 0; i < 33; ++i){
        ub[i] = lb[i];
    }

    auto tile = TA::TensorD(std::move(rng));

    for (auto const &s0 : sh0) {
        const auto ns0 = s0.size();
        ub[0] += ns0;

        lb[1] = ub[1] = lobound[1];
        for (auto const &s1 : sh1) {
            const auto ns1 = s1.size();
            ub[1] += ns1;

            lb[2] = ub[2] = lobound[2];
            for (auto const &s2 : sh2) {
                const auto ns2 = s2.size();
                ub[2] += ns2;

                auto map
                      = TA::make_map(intd::shell_set(eng, s0, s1, s2), 
                              {lb[0], lb[1], lb[2]},
                              {ub[0], ub[1], ub[2]});
                tile.block({lb[0], lb[1], lb[2]}, {ub[0], ub[1], ub[2]}) = map;

                lb[2] = ub[2];
            }
            lb[1] = ub[1];
        }
        lb[0] = ub[0];
    }

    return tile;
}

TA::TensorD
integral_kernel_best(coulomb_engine &eng, TA::Range &&rng, ShellVec const &sh0,
                     ShellVec const &sh1, ShellVec const &sh2) {

    auto const &lobound = rng.lobound();
    std::array<unsigned long, 3> lb = {{lobound[0], lobound[1], lobound[2]}};
    std::array<unsigned long, 3> ub = lb;

    auto tile = TA::TensorD(std::move(rng));

    // init map
    const double dummy = 0.0;
    auto map = TA::make_map(&dummy, {0, 0, 0}, {1, 1, 1});

    for (auto const &s0 : sh0) {
        const auto ns0 = s0.size();
        ub[0] += ns0;

        lb[1] = ub[1] = lobound[1];
        for (auto const &s1 : sh1) {
            const auto ns1 = s1.size();
            ub[1] += ns1;

            lb[2] = ub[2] = lobound[2];
            for (auto const &s2 : sh2) {
                const auto ns2 = s2.size();
                ub[2] += ns2;

                TA::remap(map, intd::shell_set(eng, s0, s1, s2), lb, ub);
                tile.block(lb, ub) = map;

                lb[2] = ub[2];
            }
            lb[1] = ub[1];
        }
        lb[0] = ub[0];
    }

    return tile;
}
