
#include "mkclasses.h"

MkClassesOptions::MkClassesOptions()
{
  gensrc_ = false;
  geninc_ = false;
}

MkClassesOptions::MkClassesOptions(const MkClassesOptions &o)
{
  *this = o;
}

MkClassesOptions::~MkClassesOptions()
{
}

