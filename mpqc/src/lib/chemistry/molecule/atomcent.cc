
#ifdef __GNUC__
#pragma implementation
#endif

#include <string.h>

#include "atomcent.h"

DescribedClass_REF_def(AtomicCenter);

#define CLASSNAME AtomicCenter
#define PARENTS public SavableState
#define HAVE_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
AtomicCenter::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

AtomicCenter::AtomicCenter() :
  label_(0)
{
}

AtomicCenter::AtomicCenter(const AtomicCenter&ac) :
  p(ac.p),element_(ac.element_),label_(0)
{
  if (ac.label_) {
    label_ = new char[strlen(ac.label_)+1];
    strcpy(label_,ac.label_);
  }
}

AtomicCenter::AtomicCenter(const char*symbol,double x,double y,double z, 
                           const char *lab) :
  p(x,y,z),
  element_(symbol),
  label_(0)
{
  if (lab) {
    label_ = new char[strlen(lab)+1];
    strcpy(label_,lab);
  }
}

AtomicCenter::~AtomicCenter()
{
  if (label_) delete[] label_; label_=0;
}

AtomicCenter& AtomicCenter::operator=(const AtomicCenter&ac)
{
  p = ac.p;
  element_ = ac.element_;
  if (ac.label_) {
    label_ = new char[strlen(ac.label_)+1];
    strcpy(label_,ac.label_);
  }

  return *this;
}

void AtomicCenter::save_data_state(StateOut& so)
{
  p.save_object_state(so);
  element_.save_object_state(so);
  so.putstring(label_);
}

AtomicCenter::AtomicCenter(StateIn& si):
  SavableState(si),
  p(si),
  element_(si)
{
  si.getstring(label_);
}

void AtomicCenter::print(FILE*fp)
{
  fprintf(fp,"%2s",element().symbol());
  point().print(fp);
}
