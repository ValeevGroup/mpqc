
#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/scf/effh.h>

///////////////////////////////////////////////////////////////////////////
// GSGeneralEffH

#define CLASSNAME GSGeneralEffH
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public AccumEffectiveH
//#include <util/state/statei.h>
#include <util/class/classi.h>
void *
GSGeneralEffH::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = AccumEffectiveH::_castdown(cd);
  return do_castdowns(casts,cd);
}

void
GSGeneralEffH::init()
{
  coef(0,0,0) =  1.0;
  coef(0,1,0) =  2.0; coef(0,1,1) = 1.0;
  coef(0,2,0) =  1.0; coef(0,2,1) = 0.0; coef(0,2,2) = 1.0;

  coef(1,0,0) =  0.0;
  coef(1,1,0) = -1.0; coef(1,1,1) = 0.0;
  coef(1,2,0) =  0.0; coef(1,2,1) = 1.0; coef(1,2,2) = 0.0;
}

GSGeneralEffH::GSGeneralEffH(StateIn& s) :
  AccumEffectiveH(s)
{
  init();
}

GSGeneralEffH::GSGeneralEffH(const RefKeyVal& keyval) :
  AccumEffectiveH(keyval)
{
  init();
}

GSGeneralEffH::~GSGeneralEffH()
{
}

///////////////////////////////////////////////////////////////////////////
// GSHighSpinEffH

#define CLASSNAME GSHighSpinEffH
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public AccumEffectiveH
#include <util/class/classi.h>
void *
GSHighSpinEffH::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = AccumEffectiveH::_castdown(cd);
  return do_castdowns(casts,cd);
}

void
GSHighSpinEffH::init()
{
  coef(0,0,0) =  2.0;
  coef(0,1,0) =  2.0; coef(0,1,1) = 2.0;
  coef(0,2,0) =  1.0; coef(0,2,1) = 0.0; coef(0,2,2) = 2.0;

  coef(1,0,0) = -1.0;
  coef(1,1,0) = -1.0; coef(1,1,1) = -1.0;
  coef(1,2,0) =  0.0; coef(1,2,1) =  1.0; coef(1,2,2) = -1.0;
}

GSHighSpinEffH::GSHighSpinEffH(StateIn& s) :
  AccumEffectiveH(s)
{
  init();
}

GSHighSpinEffH::GSHighSpinEffH(const RefKeyVal& keyval) :
  AccumEffectiveH(keyval)
{
  init();
}

GSHighSpinEffH::~GSHighSpinEffH()
{
}

///////////////////////////////////////////////////////////////////////////
// TestEffH

#define CLASSNAME TestEffH
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public AccumEffectiveH
#include <util/class/classi.h>
void *
TestEffH::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = AccumEffectiveH::_castdown(cd);
  return do_castdowns(casts,cd);
}

void
TestEffH::init()
{
  coef(0,0,0) =  0.0;
  coef(0,1,0) =  2.0; coef(0,1,1) = 0.0;
  coef(0,2,0) =  1.0; coef(0,2,1) = 0.0; coef(0,2,2) = 0.0;

  coef(1,0,0) =  1.0;
  coef(1,1,0) = -1.0; coef(1,1,1) = 1.0;
  coef(1,2,0) =  0.0; coef(1,2,1) = 1.0; coef(1,2,2) = 1.0;
}

TestEffH::TestEffH(StateIn& s) :
  AccumEffectiveH(s)
{
  init();
}

TestEffH::TestEffH(const RefKeyVal& keyval) :
  AccumEffectiveH(keyval)
{
  init();
}

TestEffH::~TestEffH()
{
}

///////////////////////////////////////////////////////////////////////////
// PsiEffH

#define CLASSNAME PsiEffH
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS public AccumEffectiveH
#include <util/class/classi.h>
void *
PsiEffH::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = AccumEffectiveH::_castdown(cd);
  return do_castdowns(casts,cd);
}

void
PsiEffH::init()
{
  coef(0,0,0) =  1.0;
  coef(0,1,0) =  2.0; coef(0,1,1) = 0.0;
  coef(0,2,0) =  1.0; coef(0,2,1) = 0.0; coef(0,2,2) = 0.0;

  coef(1,0,0) =  0.0;
  coef(1,1,0) = -1.0; coef(1,1,1) = 1.0;
  coef(1,2,0) =  0.0; coef(1,2,1) = 1.0; coef(1,2,2) = 1.0;
}

PsiEffH::PsiEffH(StateIn& s) :
  AccumEffectiveH(s)
{
  init();
}

PsiEffH::PsiEffH(const RefKeyVal& keyval) :
  AccumEffectiveH(keyval)
{
  init();
}

PsiEffH::~PsiEffH()
{
}
