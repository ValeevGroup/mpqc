
#ifndef CLASSNAME
#define CLASSNAME you_forgot_to_define_CLASSNAME
#endif

//#define POINTER_TO_CASTDOWN__(arg1,arg2) arg1 ## arg2
//#define POINTER_TO_CASTDOWN_(arg1,arg2) POINTER_TO_CASTDOWN__(arg1,arg2)
//#define POINTER_TO_CASTDOWN POINTER_TO_CASTDOWN_(CLASSNAME,pointer_to_castdown)

  public:
     // this is public to make in easier to force linkage
    static ClassDesc class_desc_;
  private:
    //typedef void* CLASSNAME::(* POINTER_TO_CASTDOWN ) (const ClassDesc*);
    void* CLASSNAME::do_castdowns(//void* CLASSNAME::(*casts[])(ClassDesc*),
                                  //POINTER_TO_CASTDOWN *,
                                  void**,
                                  const ClassDesc*cd);
  public:
    void* _castdown(const ClassDesc*);
    static CLASSNAME* castdown(DescribedClass*p);
    static CLASSNAME* castdown(RefDescribedClass&p);
    inline static const ClassDesc* static_class_desc()
      { return &CLASSNAME::class_desc_; }
    inline const ClassDesc* class_desc() const { return &class_desc_; }
#ifdef HAVE_CTOR
    static DescribedClass* create();
#undef HAVE_CTOR
#endif
#ifdef HAVE_KEYVAL_CTOR
    static DescribedClass* create(KeyVal&);
#undef HAVE_KEYVAL_CTOR
#endif
#ifdef HAVE_STATEIN_CTOR
    static DescribedClass* create(StateIn&);
#undef HAVE_STATEIN_CTOR
#endif
  private:

#undef HAVE_KEYVAL_CTOR
#undef HAVE_STATEIN_CTOR
#undef HAVE_CTOR
#undef CLASSNAME

//#undef POINTER_TO_CASTDOWN__
//#undef POINTER_TO_CASTDOWN_
//#undef POINTER_TO_CASTDOWN
