
#ifndef _util_container_ref_h
#define _util_container_ref_h

#include <stdio.h>
#include <stdlib.h>

// The macro version of the reference counting class
#define REF_dec_custom(T,custom)					      \
class  Ref ## T  {							      \
  private:								      \
    T* p;								      \
  public:								      \
    inline T* operator->() { return p; }				      \
    inline const T* operator->() const { return p; }			      \
    inline T* pointer() { return p; }					      \
    inline const T* pointer() const { return p; }			      \
    inline operator T*() { return p; };					      \
    inline const operator T*() const { return p; };			      \
    inline T& operator *() { return *p; };				      \
    inline const T& operator *() const { return *p; };			      \
    Ref ## T ();							      \
    Ref ## T (T*a);							      \
    Ref ## T ( Ref ## T &a);						      \
    ~Ref ## T ();							      \
    inline int null() { return p == 0; }				      \
    inline int nonnull() { return p != 0; }				      \
    Ref ## T  operator=( Ref ## T & c);					      \
    Ref ## T  operator=(T* cr);						      \
    void assign_pointer(T* cr);						      \
    inline int operator==( Ref ## T &a) { return p == a.p; }		      \
    inline int operator==( T * a) { return p == a; }			      \
    inline int operator>=( Ref ## T &a) { return p >= a.p; }		      \
    inline int operator<=( Ref ## T &a) { return p <= a.p; }		      \
    inline int operator>( Ref ## T &a) { return p > a.p; }		      \
    inline int operator<( Ref ## T &a) { return p < a.p; }		      \
    void  ref_info(FILE*fp=stdout);					      \
    custom								      \
}

#define REF_dec(T) REF_dec_custom(T,)

#define REF_def(T)							      \
Ref ## T :: Ref ## T (): p(0) {}					      \
Ref ## T :: Ref ## T (T*a): p(a) { if (p) p->count++; }			      \
Ref ## T :: Ref ## T ( Ref ## T &a): p(a.p) { if (p) p->count++; }	      \
Ref ## T :: ~Ref ## T () { if (p && --p->count<=0) delete p; }		      \
Ref ## T  Ref ## T :: operator=( Ref ## T & c)				      \
{									      \
  if (   p								      \
         && --p->count <= 0						      \
         && p != c.p)							      \
    delete p;								      \
  p=c.p;								      \
  if (p) p->count++;							      \
  return *this;								      \
}									      \
Ref ## T  Ref ## T :: operator=(T* cr)					      \
{									      \
  if (p && --p->count <= 0 && p != cr) delete p;			      \
  p=cr;									      \
  if (p) p->count++;							      \
  return *this;								      \
}									      \
void Ref ## T :: assign_pointer(T* cr)					      \
{									      \
  if (p && --p->count <= 0 && p != cr) delete p;			      \
  p=cr;									      \
  if (p) p->count++;							      \
}									      \
void Ref ## T :: ref_info(FILE*fp=stdout)				      \
{									      \
  if (nonnull()) fprintf(fp,"count = %d\n",p->count);			      \
  else fprintf(fp,"reference is null\n");				      \
}

//template <class T> class Ref;
class RefCount {
    //template <class T> friend class Ref;
  public: // should be private, but template friends don't work
    int count;
  protected:
    inline RefCount():count(0) {};
  public:
    inline ~RefCount() {};
};

//template <class T> class Ref;
class VRefCount {
    //template <class T> friend class Ref;
  public: // should be private, but template friends don't work
    int count;
  protected:
    inline VRefCount():count(0) {};
  public:
    virtual ~VRefCount();
};

// template <class T>
// class Ref {
//   private:
//     T* p;
//   public:
//     inline T* operator->() const { return p; }
//     inline T* pointer() const { return p; }
// 
//     Ref(): p(0) {}
//     Ref(T*a): p(a) { if (p) p->count++; }
//     Ref(Ref&a): p(a.p) { if (p) p->count++; }
//     ~Ref() { if (p && --p->count<=0) delete p; }
// 
//     inline int null() { return p == 0; }
//     inline int nonnull() { return p != 0; }
// 
//     Ref operator=(Ref& c)
//       {
//       if (   p
//           && --p->count <= 0
//           && p != c.p)
//         delete p;
//       p=c.p;
//       if (p) p->count++;
//       return *this;
//       }
//     Ref operator=(T* cr)
//       {
//       if (p && --p->count <= 0 && p != cr) delete p;
//       p=cr;
//       if (p) p->count++;
//       return *this;
//       }
//     void assign_pointer(T* cr)
//       {
//       if (p && --p->count <= 0 && p != cr) delete p;
//       p=cr;
//       if (p) p->count++;
//       }
// 
//     inline int operator==(Ref&a) { return p == a.p; }
//     inline int operator>=(Ref&a) { return p >= a.p; }
//     inline int operator<=(Ref&a) { return p <= a.p; }
//     inline int operator>(Ref&a) { return p > a.p; }
//     inline int operator<(Ref&a) { return p < a.p; }
// 
//     void ref_info(FILE*fp=stdout)
//       {
//       if (nonnull()) fprintf(fp,"count = %d\n",p->count);
//       else fprintf(fp,"reference is null\n");
//       }
// 
//     // Doesn't work under gcc 2.3.1
//     //template <class U> operator Ref<U>() { Ref<U> r(p); return r; }
//     // A hack to compensate for the above:
//     //class GuiObject;
//     //operator Ref<GuiObject>() { Ref<GuiObject> r(p); return r; }
//     // but, for now I won't allow it in case the template member function
//     // really isn't allowed by the language.  Type conversions must now
//     // be done manually with the assistance of:
//     inline operator T*() { return p; }
//     // example of use of the above:
//     // Ref<GuiCascadeButton> cb;
//     // Ref<GuiObject> o = Ref<GuiObject>((GuiCascadeButton*)cb);
//   };

#endif

