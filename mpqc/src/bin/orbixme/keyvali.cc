
#include <strstream.h>
#include <iostream.h>
#include "keyvali.h"

#include <util/class/class.h>
#include <util/keyval/keyval.h>

C_KeyValCreatableImpl::C_KeyValCreatableImpl()
{
  dc_ = 0;
}

void
C_KeyValCreatableImpl::keyval_create(const char *s, CORBA_Environment &IT_env)
{
  cout << "create: " << s << endl;

  istrstream in(s);

  RefKeyVal keyval = new ParsedKeyVal(in);

  RefDescribedClass dc = keyval->describedclassvalue("object");
  if (dc.null()) {
      cout << "create failed" << endl;
      return;
    }
  else {
      cout << "created a \"" << dc->class_name() << "\"" << endl;
    }

  dc_ = dc.pointer();
  dc_->reference();
}

C_KeyValCreatableImpl::~C_KeyValCreatableImpl()
{
  dc_->dereference();
  if (dc_->nreference() == 0) delete dc_;
}

