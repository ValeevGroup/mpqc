#define CATCH_CONFIG_RUNNER

#include <clocale>
#include <madness/world/world.h>
#include <madness/world/worldgop.h>
#include "../../../tests/unit/catch.hpp"

#include "mpqc/util/keyval/keyval.h"
#include "mpqc_init.h"
#include "mpqc_task.h"
#include "mpqc/chemistry/qc/scf/linkage.h"

int main( int argc, char* argv[] )
{
  Catch::Session session;
  // global setup...
  std::setlocale(LC_ALL,"en_US.UTF-8");

  auto &world = madness::initialize(argc, argv);
  mpqc::initialize(argc, argv);

  int result = session.run( argc, argv );

  return result;
}

TEST_CASE("multiple tasks", "[multiworld]") {
  auto& world = madness::World::get_default();
  if (world.size() == 2) {
    using namespace madness;

    // make subworlds
    const auto rank = world.rank();
    const bool odd = rank & 0x1;
    int even_ranks[] = {0};
    int odd_ranks[] =  {1};
    SafeMPI::Group group = world.mpi.comm().Get_group().Incl(1, odd ? odd_ranks : even_ranks);
    SafeMPI::Intracomm comm = world.mpi.comm().Create(group);
    World my_world(comm);
    world.gop.fence();

    double energy[2];
    {
      const std::string srcdir = std::string(SRCDIR);
      const std::string ifile =
          srcdir +
          std::string(odd ? "mpqc_test_odd.json" : "mpqc_test_even.json");
      std::shared_ptr<mpqc::KeyVal> kv =
          mpqc::MPQCInit::instance().make_keyval(my_world, ifile);
      kv->assign("file_prefix", srcdir);
      mpqc::MPQCTask task(my_world, kv);
      task.run();
      energy[rank] = kv->value<double>("wfn:energy");
    }

    if (rank == 0) {
      madness::World::get_default().gop.broadcast(energy);

      // check energies
      REQUIRE(energy[0] == 0.0);
      REQUIRE(energy[1] == 0.0);
    }
  }
}
