
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/scmat/block.h>

/////////////////////////////////////////////////////////////////////////////
// SCMatrixBlock member functions

SavableState_REF_def(SCMatrixBlock);

#define CLASSNAME SCMatrixBlock
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SCMatrixBlock::SCMatrixBlock()
{
}

void *
SCMatrixBlock::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCMatrixBlock::~SCMatrixBlock()
{
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixRectBlock member functions

SavableState_REF_def(SCMatrixRectBlock);

#define CLASSNAME SCMatrixRectBlock
#define PARENTS public SCMatrixBlock
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

SCMatrixRectBlock::SCMatrixRectBlock(int is, int ie, int js, int je):
  istart(is),
  jstart(js),
  iend(ie),
  jend(je)
{
  data = new double[(ie-is)*(je-js)];
}

SCMatrixRectBlock::SCMatrixRectBlock(StateIn&s):
  SavableState(s,SCMatrixRectBlock::class_desc_)
{
  s.get(istart);
  s.get(jstart);
  s.get(iend);
  s.get(jend);
  s.get(data);
}

void
SCMatrixRectBlock::save_data_state(StateOut&s)
{
  s.put(istart);
  s.put(jstart);
  s.put(iend);
  s.put(jend);
  s.put(data,(iend-istart)*(jend-jstart));
}

void *
SCMatrixRectBlock::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCMatrixBlock::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCMatrixRectBlock::~SCMatrixRectBlock()
{
  delete[] data;
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixLTriBlock member functions

SavableState_REF_def(SCMatrixLTriBlock);

#define CLASSNAME SCMatrixLTriBlock
#define PARENTS public SCMatrixBlock
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

SCMatrixLTriBlock::SCMatrixLTriBlock(int s,int e):
  start(s),
  end(e)
{
  data = new double[((e-s)*(e-s+1))/2];
}

SCMatrixLTriBlock::SCMatrixLTriBlock(StateIn&s):
  SavableState(s,SCMatrixLTriBlock::class_desc_)
{
  s.get(start);
  s.get(end);
  s.get(data);
}

void
SCMatrixLTriBlock::save_data_state(StateOut&s)
{
  s.put(start);
  s.put(end);
  s.put(data,((end-start)*(end-start+1))/2);
}

void *
SCMatrixLTriBlock::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCMatrixBlock::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCMatrixLTriBlock::~SCMatrixLTriBlock()
{
  delete[] data;
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixDiagBlock member functions

SavableState_REF_def(SCMatrixDiagBlock);

#define CLASSNAME SCMatrixDiagBlock
#define PARENTS public SCMatrixBlock
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

SCMatrixDiagBlock::SCMatrixDiagBlock(int s, int e):
  istart(s),
  iend(e),
  jstart(s)
{
  data = new double[e-s];
}

SCMatrixDiagBlock::SCMatrixDiagBlock(int is, int ie,int js):
  istart(is),
  iend(ie),
  jstart(js)
{
  data = new double[ie-is];
}

SCMatrixDiagBlock::SCMatrixDiagBlock(StateIn&s):
  SavableState(s,SCMatrixDiagBlock::class_desc_)
{
  s.get(istart);
  s.get(jstart);
  s.get(iend);
  s.get(data);
}

void
SCMatrixDiagBlock::save_data_state(StateOut&s)
{
  s.put(istart);
  s.put(jstart);
  s.put(iend);
  s.put(data,iend-istart);
}

void *
SCMatrixDiagBlock::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCMatrixBlock::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCMatrixDiagBlock::~SCMatrixDiagBlock()
{
  delete[] data;
}

/////////////////////////////////////////////////////////////////////////////
// SCVectorSimpleBlock member functions

SavableState_REF_def(SCVectorSimpleBlock);

#define CLASSNAME SCVectorSimpleBlock
#define PARENTS public SCMatrixBlock
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

SCVectorSimpleBlock::SCVectorSimpleBlock(int s, int e):
  istart(s),
  iend(e)
{
  data = new double[e-s];
}

SCVectorSimpleBlock::SCVectorSimpleBlock(StateIn&s):
  SavableState(s,SCVectorSimpleBlock::class_desc_)
{
  s.get(istart);
  s.get(iend);
  s.get(data);
}

void
SCVectorSimpleBlock::save_data_state(StateOut&s)
{
  s.put(istart);
  s.put(iend);
  s.put(data,iend-istart);
}

void *
SCVectorSimpleBlock::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCMatrixBlock::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCVectorSimpleBlock::~SCVectorSimpleBlock()
{
  delete[] data;
}
