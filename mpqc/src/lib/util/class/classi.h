
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

//#define POINTER_TO_CASTDOWN__(arg1,arg2) arg1 ## arg2
//#define POINTER_TO_CASTDOWN_(arg1,arg2) POINTER_TO_CASTDOWN__(arg1,arg2)
//#define POINTER_TO_CASTDOWN POINTER_TO_CASTDOWN_(CLASSNAME,pointer_to_castdown)

ClassDesc CLASSNAME::class_desc_(stringize(CLASSNAME),
                                 VERSION,
                                 stringize2(PARENTS),
#ifdef HAVE_CTOR
                                 CLASSNAME::create,
#else
                                 0,
#endif
#ifdef HAVE_KEYVAL_CTOR
                                 CLASSNAME::create,
#else
                                 0,
#endif
#ifdef HAVE_STATEIN_CTOR
                                 CLASSNAME::create
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
      va_end(args);
      abort();
  }
  return t;
}
CLASSNAME*
CLASSNAME::castdown(RefDescribedClass&p)
{
  if (p.null()) return 0;
  return (CLASSNAME*) p->_castdown(CLASSNAME::static_class_desc());
}
void *
//CLASSNAME::do_castdowns(POINTER_TO_CASTDOWN*casts,const ClassDesc*cd)
CLASSNAME::do_castdowns(void**casts,const ClassDesc*cd)
{
  if (cd == &CLASSNAME::class_desc_) {
      return this;
    }
  void* p = 0;
  for (int i=0; i<CLASSNAME::class_desc_.parents().n(); i++) {
      if (!CLASSNAME::class_desc_.parents()[i].access()
          == ParentClass::Private) {
          //void * tmp = casts[i](cd);
          void * tmp = casts[i];
          if (!tmp) continue;
          if (p && tmp != p) {
              fprintf(stderr,"%s: castdown to %s ambiguous (from %s)\n",
                      CLASSNAME::class_desc_.name(),
                      cd->name(),
                      CLASSNAME::class_desc_.name());
              fprintf(stderr," tmp = 0x%x p = 0x%x\n",tmp,p);
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
#undef HAVE_CTOR
#endif
#ifdef HAVE_KEYVAL_CTOR
DescribedClass*
CLASSNAME::create(KeyVal& keyval)
{
  return (DescribedClass*) new CLASSNAME(keyval);
}
#undef HAVE_KEYVAL_CTOR
#endif
#ifdef HAVE_STATEIN_CTOR
DescribedClass*
CLASSNAME::create(StateIn& statein)
{
  return (DescribedClass*) new CLASSNAME(statein);
}
#undef HAVE_STATEIN_CTOR
#endif

#undef CLASSNAME
#undef PARENTS
#undef VERSION

#undef stringize_
#undef stringize

//#undef POINTER_TO_CASTDOWN__
//#undef POINTER_TO_CASTDOWN_
//#undef POINTER_TO_CASTDOWN
