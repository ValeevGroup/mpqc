
#ifdef __GNUG__
#pragma interface
#endif

template <class T>
class  Ref  : private RefBase {
  private:
    T* p;
  public:
    T* operator->() const { return p; }
    T* pointer() const { return p; }
    REF_TYPE_CAST_DEC(T);
    T& operator *() const { return *p; };
    int null() const { return p == 0; }
    int nonnull() const { return p != 0; }
    int operator!=(const Ref<T> &a) const { return p != a.p; }
    int compare(const Ref<T> &a) const
    { return ((p==a.p)?0:((p<a.p)?-1:1)); }
    int operator==(const Ref<T> &a) const { return p == a.p; }
    int operator!=(const  T * a) const { return p != a; }
    int operator==(const  T * a) const { return p == a; }
    int operator>=(const  Ref<T> &a) const { return p >= a.p; }
    int operator<=(const  Ref<T> &a) const { return p <= a.p; }
    int operator>(const  Ref<T> &a) const { return p > a.p; }
    int operator<(const  Ref<T> &a) const { return p < a.p; }
    Ref(): p(0) {}
    Ref(T*a): p(a)
    {
      if (p) {
          if (DO_REF_CHECK_STACK(p)) {
              DO_REF_UNMANAGE(p);
              warn_ref_to_stack();
            }
          p->reference();
        }
    }
    Ref(const Ref<T> &a): p(a.p)
    {
      if (p) p->reference();
    }
    ~ Ref ()
    {
      clear();
    }
    void clear()
    {
      if (p && p->dereference()<=0) {
          delete p;
        }
      p = 0;
    }
    Ref<T>& operator=(const Ref<T> & c)
    {
      if (c.p) c.p->reference();
      clear();
      p=c.p;
      return *this;
    }
    Ref<T>& operator=(T* cr)
    {
      assign_pointer(cr);
      return *this;
    }
    void assign_pointer(T* cr)
    {
      if (cr) {
          if (DO_REF_CHECK_STACK(cr)) {
              DO_REF_UNMANAGE(cr);
              warn_ref_to_stack();
            }
          cr->reference();
        }
      clear();
      p = cr;
    }
    void check_pointer() const
    {
      if (p && p->nreference() <= 0) {
          warn_bad_ref_count();
        }
    }
    void ref_info(FILE*fp) const
    {
      RefBase::ref_info(p,fp);
    }
    void warn(const char*s) const { RefBase::warn(s); }
};
