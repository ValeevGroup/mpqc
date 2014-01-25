
#include <mpqc_config.h>
#include <util/group/pregtime.h>
#include <util/misc/exenv.h>
#if defined(HAVE_MPI)
#  include <util/group/messmpi.h>
#endif


void
dotest(const sc::Ref<sc::MessageGrp> &grp,
       bool node_independent)
{
  sc::Ref<sc::RegionTimer> rt
      = new sc::ParallelRegionTimer(grp->clone(),"prttest",1,1);

  sc::Timer t1(rt,"t1_a");

  if (node_independent || grp->me() == 0) {
      sc::Timer t2(rt,"t2_a");
      sc::Timer t3(rt,"t3_a");
      t3.reset("t3_b");
      t3.reset("t3_c");
      t3.reset("t3_d");
      t3.reset();
      t2.reset("t2_b");
      if (node_independent || grp->me() == 0) {
          sc::Timer t4(rt,"t4_a");
          t4.reset("t4_b");
          t4.reset();
        }
      t2.reset("t2_c");
      t2.reset();
    }

  t1.reset();

  rt->print(sc::ExEnv::outn());
}

int
main(int argc, char**argv)
{
  sc::Ref<sc::MessageGrp> grp;
#if defined(HAVE_MPI) && defined(ALWAYS_USE_MPI)
  grp = new sc::MPIMessageGrp(&argc, &argv);
#endif
  if (grp == 0) grp = sc::MessageGrp::initial_messagegrp(argc, argv);
  if (grp)
      sc::MessageGrp::set_default_messagegrp(grp);
  else
      grp = sc::MessageGrp::get_default_messagegrp();

  std::cout << grp->me() << ": using " << grp->class_name() << std::endl;

  dotest(grp, true);
  dotest(grp, false);

  grp = 0;
  sc::MessageGrp::set_default_messagegrp(grp);

  return 0;
}
