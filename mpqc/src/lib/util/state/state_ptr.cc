
#ifdef __GNUC__
#pragma implementation
#endif

#include <util/class/class.h>
#include "state.h"
#include "stateptrImplSet.h"
#include "statenumImplSet.h"

// Returns the number of the new object, if this object is new.
// If the object is found then zero is returned.
int StateIn::getpointer(void**p)
{
  int refnum;
  get(refnum);
  if (refnum == 0) {
    *p = 0;
    //printf("StateOut::getpointer: pointer is 0\n");
    return 0;
    }
  StateDataNum num(refnum);
  Pix ind = (ps_?ps_->seek(num):0);
  //printf("StateOut::getpointer: looking for %d and got %d\n",refnum,(int)ind);
  if (ind == 0) {
    *p = 0;
    return refnum;
    }
  else {
    *p = ((*ps_)(ind)).ptr();
    //printf("StateOut::getpointer: pointer is made 0x%x\n",*p);
    return 0;
    }
  }

void StateIn::nextobject(int objnum)
{
  _nextobject = objnum;
}

void StateIn::havepointer(void*p)
{
  if (_nextobject) {
      havepointer(_nextobject,p);
      _nextobject = 0;
    }
}

void StateIn::havepointer(int objnum,void*p)
{
  StateDataNum num(objnum,p);
  if (ps_) ps_->add(num);
}

// Returns 0 if the object has already been written.
// Returns 1 if the object must yet be written.
int StateOut::putpointer(void*p)
{
  if (p == 0) {
    put(0);
    return 0;
    }
  StateDataPtr dp(p);
  Pix ind = (ps_?ps_->seek(dp):0);
  //printf("StateOut::putpointer: ind = %d for 0x%x\n",(int)ind,p);
  if (ind == 0) {
      if (ps_) {
          dp.assign_num(next_pointer_number++);
          ps_->add(dp);
        }
      put(dp.num());
      return 1;
    }
  else {
      put((*ps_)(ind).num());
      return 0;
    }
}
