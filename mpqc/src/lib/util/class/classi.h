//
// classi.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef CLASSNAME
#define CLASSNAME you_forgot_to_define_CLASSNAME
#endif
#ifndef VERSION
#define VERSION 1
#endif
#ifndef PARENTS
#define PARENTS
#endif

#define stringize_(arg) # arg
// this has () so commas are treated correctly
#define stringize2(arg) stringize_((arg))
#define stringize(arg) stringize_(arg)

ClassDesc CLASSNAME::class_desc_(stringize(CLASSNAME),
                                 VERSION,
                                 stringize2(PARENTS),
#ifdef HAVE_CTOR
                                 &CLASSNAME::create,
#else
                                 0,
#endif
#ifdef HAVE_KEYVAL_CTOR
                                 &CLASSNAME::create,
#else
                                 0,
#endif
#ifdef HAVE_STATEIN_CTOR
                                 &CLASSNAME::create
#else
                                 0
#endif
                                 );
const ClassDesc* CLASSNAME::static_class_desc()
{
  return &CLASSNAME::class_desc_;
}
const ClassDesc* CLASSNAME::class_desc() const
{
    return &CLASSNAME::class_desc_;
}
CLASSNAME*
CLASSNAME::castdown(DescribedClass*p)
{
  if (!p) return 0;
  return (CLASSNAME*) p->_castdown(CLASSNAME::static_class_desc());
}
CLASSNAME*
CLASSNAME::require_castdown(DescribedClass*p,const char * errmsg,...)
{
  if (!p) return 0;
  CLASSNAME* t = (CLASSNAME*) p->_castdown(CLASSNAME::static_class_desc());
  if (!t) {
      va_list args;
      va_start(args,errmsg);
      fprintf(stderr,"A required castdown failed in: ");
      vfprintf(stderr,errmsg,args);
      fprintf(stderr,"\nwanted type \"%s\" but got \"%s\"\n",
              stringize(CLASSNAME),p?p->class_name():"(null)");
      fflush(stderr);
      va_end(args);
      ::abort();
  }
  return t;
}
CLASSNAME*
CLASSNAME::castdown(const RefDescribedClass&p)
{
  if (p.null()) return 0;
  return (CLASSNAME*) p->_castdown(CLASSNAME::static_class_desc());
}
void *
CLASSNAME::do_castdowns(void**casts,const ClassDesc*cd)
{
  if (cd == &CLASSNAME::class_desc_) {
      return this;
    }
  void* p = 0;
  const ParentClasses& parents = CLASSNAME::class_desc_.parents();
  int n = parents.n();
  for (int i=0; i<n; i++) {
      if (!parents[i].access() == ParentClass::Private) {
          void * tmp = casts[i];
          if (!tmp) continue;
          if (p && tmp != p) {
              cerr << CLASSNAME::class_desc_.name()
                   << ": castdown to " << cd->name()
                   << " ambiguous (from "
                   << CLASSNAME::class_desc_.name() << ")" << endl
                   << " tmp = 0x" << setbase(16) << (long)tmp
                   << " p = 0x" << setbase(16) << (long)p << endl;
            }
          p = tmp;
        }
    }
  return p;
}
#ifdef HAVE_CTOR
DescribedClass*
CLASSNAME::create()
{
  return (DescribedClass*) new CLASSNAME();
}
#endif
#ifdef HAVE_KEYVAL_CTOR
DescribedClass*
CLASSNAME::create(const RefKeyVal& keyval)
{
  return (DescribedClass*) new CLASSNAME(keyval);
}
#endif
#ifdef HAVE_STATEIN_CTOR
DescribedClass*
CLASSNAME::create(StateIn& statein)
{
  return (DescribedClass*) new CLASSNAME(statein);
}
#endif

#undef CLASSNAME
#undef PARENTS
#undef VERSION
#undef HAVE_CTOR
#undef HAVE_KEYVAL_CTOR
#undef HAVE_STATEIN_CTOR

#undef stringize_
#undef stringize2
#undef stringize

