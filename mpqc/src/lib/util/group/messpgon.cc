
#include <stdio.h>
extern "C" {
#include <nx.h>
void gsync(void);
void crecv(long typesel, char *buf, long count);
void csend(long type, char *buf, long count, long node, long ptype);
}
#include <util/group/messpgon.h>

#define CLASSNAME ParagonMessageGrp
#define PARENTS public intMessageGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
ParagonMessageGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  intMessageGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

ParagonMessageGrp::ParagonMessageGrp(const RefKeyVal&)
{
  initialize();
}

ParagonMessageGrp::ParagonMessageGrp()
{
  initialize();
}

void
ParagonMessageGrp::initialize()
{
  int nprocs = numnodes();
  int mynodeid = mynode();
  intMessageGrp::initialize(mynodeid, nprocs, 30);
}

int
ParagonMessageGrp::basic_probe(int type)
{
  return iprobe(type);
}

ParagonMessageGrp::~ParagonMessageGrp()
{
}

void
ParagonMessageGrp::sync()
{
  gsync();
}

void
ParagonMessageGrp::basic_recv(int type, void* buf, int bytes)
{
  crecv(type, (char*) buf, bytes);
}

void
ParagonMessageGrp::basic_send(int dest, int type, void* buf, int bytes)
{
  csend(type, (char*) buf, bytes, dest, 0);
}
 
int
ParagonMessageGrp::last_source()
{
  return infonode();
}

int
ParagonMessageGrp::last_size()
{
  return infocount();
}

int
ParagonMessageGrp::last_type()
{
  return msgtype_typ(infotype());
}
