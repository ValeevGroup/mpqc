
#ifdef __GNUG__
#pragma interface
#endif

template <class T>
class DCRef : public DCRefBase {
  protected:
    T* p;
  public:
    DescribedClass* parentpointer() const { return p; }
    DCRef(const DCRefBase&a) {
        p = T::castdown(a.parentpointer());
        reference(p);
      }
    DCRef<T>& operator=(const DCRefBase&a) {
        T* cr = T::castdown(a.parentpointer());
        reference(cr);
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
        return (eq(p,a.p)?0:(lt(p,a.p)?-1:1));
      }
    int operator==(const DCRef<T> &a) const { return eq(p,a.p); }
    int operator!=(const DCRef<T> &a) const { return ne(p,a.p); }
    int operator>=(const DCRef<T> &a) const { return ge(p,a.p); }
    int operator<=(const DCRef<T> &a) const { return le(p,a.p); }
    int operator> (const DCRef<T> &a) const { return gt(p,a.p); }
    int operator< (const DCRef<T> &a) const { return lt(p,a.p); }
    DCRef(): p(0) {}
    DCRef(T*a): p(a)
    {
      reference(p);
    }
    DCRef(const DCRef<T> &a): p(a.p)
    {
      reference(p);
    }
    ~ DCRef ()
    {
      clear();
    }
    void clear()
    {
      dereference(p);
      p = 0;
    }
    DCRef<T>& operator=(const DCRef<T> & c)
    {
      reference(c.p);
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
      reference(cr);
      clear();
      p = cr;
    }
};
