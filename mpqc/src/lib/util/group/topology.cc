#ifdef __GNUG__
#pragma implementation
#endif

#include <util/group/message.h>
#include <util/group/topology.h>

#define CLASSNAME GlobalMsgIter
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
GlobalMsgIter::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

GlobalMsgIter::GlobalMsgIter(int nproc, int me, int root)
{
  nproc_ = nproc;
  me_ = me;
  root_ = root;
  forwards();
}

#define CLASSNAME MachineTopology
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
MachineTopology::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

MachineTopology::MachineTopology()
{
}

MachineTopology::MachineTopology(const RefKeyVal&)
{
}

MachineTopology::~MachineTopology()
{
}
