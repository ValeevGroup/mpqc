
#include "actmsg.h"

#include <util/group/pregtime.h>

using namespace sc;

static Ref<ThreadGrp>
init_thread(const Ref<KeyVal>& keyval, int &argc, char **&argv)
{
  ///////////////////////////////////////////////////////////////
  // Initialize the thread group

  // get the thread group.  first try the commandline and environment
  Ref<ThreadGrp> thread = ThreadGrp::initial_threadgrp(argc, argv);
  
  // if we still don't have a group, try reading the thread group
  // from the input
  if (thread == 0) {
    thread << keyval->describedclassvalue("thread");
  }

  if (thread)
    ThreadGrp::set_default_threadgrp(thread);
  else
    thread = ThreadGrp::get_default_threadgrp();

  std::cout << "nthreads = " << thread->nthread() << std::endl;

  return thread;
}

static Ref<MessageGrp>
init_message(const Ref<KeyVal>& keyval,int &argc, char **&argv)
{
  Ref<MessageGrp> grp;

  grp << keyval->describedclassvalue("message");

  if (grp == 0) grp = MessageGrp::initial_messagegrp(argc, argv);

  if (grp == 0) grp = MessageGrp::get_default_messagegrp();

  MessageGrp::set_default_messagegrp(grp);
  
  RegionTimer::set_default_regiontimer(
      new ParallelRegionTimer(grp,"actmsgtest",1,0));

  SCFormIO::set_printnode(0);

  SCFormIO::setindent(std::cout, 2);
  SCFormIO::setindent(std::cerr, 2);
  
  return grp;
}

static void
clean_up(void)
{
  MemoryGrp::set_default_memorygrp(0);
  MessageGrp::set_default_messagegrp(0);
  ThreadGrp::set_default_threadgrp(0);
//   SCMatrixKit::set_default_matrixkit(0);
//   Integral::set_default_integral(0);
  RegionTimer::set_default_regiontimer(0);
}


int
main(int argc, char** argv)
{
  Ref<KeyVal> keyval = new ParsedKeyVal("actmsgtest.in");

  /////////////////////////////////////////////////////////
  // Parallel initialization:

  Ref<ThreadGrp> thr = init_thread(keyval,argc,argv);
  Ref<MessageGrp> msg = init_message(keyval,argc,argv);
  Ref<ActiveMessageGrp> amsggrp = new ActiveMessageGrp(msg,thr);

  Ref<ActiveMessageEcho> amsg = new ActiveMessageEcho(msg->me());

  Ref<StateSend> out = amsggrp->get_statesend();

  amsggrp->activate();
  for (int i=0; i<msg->n(); i++) {
      int remote = (1+i+msg->me())%msg->n();
      std::cout << msg->me()
                << " sending active message to " << remote << std::endl;
      amsggrp->send(remote,out,amsg);
    }
  amsggrp->deactivate();

  clean_up();

  return 0;
}
