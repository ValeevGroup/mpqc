
#ifdef __GNUG__
#pragma interface
#endif

template <class T>
class DCRef : public RefDescribedClassBase {
  private:
    T* p;
  public:
    DescribedClass* parentpointer() const { return p; }
    DCRef(const RefDescribedClassBase&a) {
        p = T::castdown(a.parentpointer());
        if (p) p->reference();
      }
    DCRef<T>& operator=(const RefDescribedClassBase&a) {
        T* cr = T::castdown(a.parentpointer());
        if (cr) cr->reference();
        clear();
        p = cr;
        return *this;
      }
    T* operator->() const { return p; }
    T* pointer() const { return p; }
    REF_TYPE_CAST_DEC(T);
    T& operator *() const { return *p; };
    int null() const { return p == 0; }
    int nonnull() const { return p != 0; }
    int compare(const DCRef<T> &a) const {
        return ((p==a.p)?0:((p<a.p)?-1:1));
      }
    int operator==(const DCRef<T> &a) const { return p == a.p; }
    int operator!=(const DCRef<T> &a) const { return p != a.p; }
    int operator>=(const DCRef<T> &a) const { return p >= a.p; }
    int operator<=(const DCRef<T> &a) const { return p <= a.p; }
    int operator> (const DCRef<T> &a) const { return p >  a.p; }
    int operator< (const DCRef<T> &a) const { return p <  a.p; }
    int operator==(const DescribedClass*a) const {return parentpointer() == a;}
    int operator!=(const DescribedClass*a) const {return parentpointer() != a;}
    int operator>=(const DescribedClass*a) const {return parentpointer() >= a;}
    int operator<=(const DescribedClass*a) const {return parentpointer() <= a;}
    int operator> (const DescribedClass*a) const {return parentpointer() >  a;}
    int operator< (const DescribedClass*a) const {return parentpointer() <  a;}
    int operator==(const RefDescribedClassBase &a) const {
        return parentpointer() == a.parentpointer(); }
    int operator!=(const RefDescribedClassBase &a) const {
        return parentpointer() != a.parentpointer(); }
    int operator>=(const RefDescribedClassBase &a) const {
        return parentpointer() >= a.parentpointer(); }
    int operator<=(const RefDescribedClassBase &a) const {
        return parentpointer() <= a.parentpointer(); }
    int operator> (const RefDescribedClassBase &a) const {
        return parentpointer() >  a.parentpointer(); }
    int operator< (const RefDescribedClassBase &a) const {
        return parentpointer() <  a.parentpointer(); }
    DCRef(): p(0) {}
    DCRef(T*a): p(a)
    {
      if (p) {
          if (DO_REF_CHECK_STACK(p)) {
              DO_REF_UNMANAGE(p);
              warn_ref_to_stack();
            }
          p->reference();
        }
    }
    DCRef(const DCRef<T> &a): p(a.p)
    {
      if (p) p->reference();
    }
    ~ DCRef ()
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
    DCRef<T>& operator=(const DCRef<T> & c)
    {
      if (c.p) c.p->reference();
      clear();
      p=c.p;
      return *this;
    }
    DCRef<T>& operator=(T* cr)
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
      RefDescribedClassBase::ref_info(p,fp);
    }
    void warn(const char*s) const { RefDescribedClassBase::warn(s); }
};
