#define CATCH_CONFIG_RUNNER

#include <clocale>
#include <madness/world/world.h>
#include <madness/world/worldgop.h>
#include "../../../tests/unit/catch.hpp"

#include "mpqc/util/keyval/keyval.h"
#include "mpqc_init.h"
#include "mpqc_task.h"
#include "mpqc/util/misc/exenv.h"

#include "mpqc/chemistry/qc/scf/linkage.h"

madness::World* world_ptr = nullptr;

int main( int argc, char* argv[] )
{
  Catch::Session session;
  // global setup...
  std::setlocale(LC_ALL,"en_US.UTF-8");

  world_ptr = &madness::initialize(argc, argv);
  assert(world_ptr != nullptr);
  mpqc::initialize(argc, argv, nullptr, *world_ptr);
  mpqc::FormIO::set_printnode(-1); // disable all output

  int result = session.run( argc, argv );

  return result;
}

TEST_CASE("multiple tasks", "[multiworld]") {
  auto& world = *world_ptr;
  {
    using namespace madness;
    using namespace mpqc;

    // make subworlds
    const auto rank = world.rank();
    const bool odd = rank & 0x1;
    std::vector<int> ranks;
    for(auto r=odd?1:0; r<world.size(); r+=2)
      ranks.push_back(r);
    SafeMPI::Group group = world.mpi.comm().Get_group().Incl(ranks.size(), &ranks[0]);
    SafeMPI::Intracomm comm = world.mpi.comm().Create(group);
    World my_world(comm);
    world.gop.fence();

    double energy[] = {0.0, 0.0};
    {
      const std::string srcdir = std::string(SRCDIR);
      const std::string ifile =
          srcdir +
          std::string(odd ? "/mpqc_test_odd.json" : "/mpqc_test_even.json");
      std::shared_ptr<mpqc::KeyVal> kv =
          mpqc::MPQCInit::instance().make_keyval(my_world, ifile);
      kv->assign("file_prefix", srcdir);
      mpqc::MPQCTask task(my_world, kv);
      task.run();
      if (rank == 0 || rank == 1) energy[rank] = kv->value<double>("wfn:energy");
    }
    world.gop.sum(energy, sizeof(energy)/sizeof(double));

    // check energies on rank 0 only
    if (rank == 0) {
      auto eref0 = -75.8234375775808;
      auto eref1 = -75.2522350404768;
      REQUIRE(std::abs(energy[0] - eref0) < 1e-11);
      if (world.size() > 1)
        REQUIRE(std::abs(energy[1] - eref1) < 1e-11);
    }
  }
}
