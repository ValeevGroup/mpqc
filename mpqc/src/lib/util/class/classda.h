
#ifndef CLASSNAME
#define CLASSNAME you_forgot_to_define_CLASSNAME
#endif

  private:
    static ClassDesc class_desc_;
    //typedef void* CLASSNAME::(* POINTER_TO_CASTDOWN ) (const ClassDesc*);
    void* CLASSNAME::do_castdowns(void**,
                                  const ClassDesc*cd);
  public:
    void* _castdown(const ClassDesc*);
    static CLASSNAME* require_castdown(DescribedClass*p,const char*,...);
    static CLASSNAME* castdown(DescribedClass*p);
    static CLASSNAME* castdown(RefDescribedClass&p);
    static const ClassDesc* static_class_desc();
    const ClassDesc* class_desc() const;
#ifdef HAVE_CTOR
#undef HAVE_CTOR
#endif
#ifdef HAVE_KEYVAL_CTOR
#undef HAVE_KEYVAL_CTOR
#endif
#ifdef HAVE_STATEIN_CTOR
#undef HAVE_STATEIN_CTOR
#endif
  private:

#undef HAVE_KEYVAL_CTOR
#undef HAVE_STATEIN_CTOR
#undef HAVE_CTOR
#undef CLASSNAME
