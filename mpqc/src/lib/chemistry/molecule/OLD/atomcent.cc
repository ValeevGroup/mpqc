
#include "atomcent.h"

DescribedClass_IMPL(AtomicCenter,1,"","")
SavableState_IMPL(AtomicCenter)
void * AtomicCenter::_castdown(const ClassDesc *cd)
{
  if(&class_desc_ == cd) return this;
  return 0;
  }

AtomicCenter::AtomicCenter()
{
}

AtomicCenter::AtomicCenter(const AtomicCenter&ac) :
p(ac.p),element_(ac.element_)
{
}

AtomicCenter::AtomicCenter(const char*symbol,double x,double y,double z):
element_(symbol),
p(x,y,z)
{
}

AtomicCenter::~AtomicCenter()
{
}

AtomicCenter& AtomicCenter::operator=(const AtomicCenter&ac)
{
  p = ac.p;
  element_ = ac.element_;
  return *this;
}

void AtomicCenter::save_data_state(StateOut& so)
{
  so.put(p); so.put(element_);
}

void AtomicCenter::restore_data_state(int v, StateIn& si)
{
  si.get(p); si.get(element_);
}

void AtomicCenter::print(FILE*fp)
{
  fprintf(fp,"%2s",element().symbol());
  point().print(fp);
}
