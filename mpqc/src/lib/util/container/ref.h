
#ifndef _util_container_ref_h
#define _util_container_ref_h

#include <stdio.h>
#include <stdlib.h>

extern "C" void * sbrk(int);

// The macro version of the reference counting class
#define REF_dec_custom(T,custom)					      \
class  Ref ## T  {							      \
  private:								      \
    T* p;								      \
  public:								      \
    T* operator->();							      \
    const T* operator->() const;					      \
    T* pointer();							      \
    const T* pointer() const;						      \
    operator T*();							      \
    const operator T*() const;						      \
    T& operator *();							      \
    const T& operator *() const;					      \
    Ref ## T ();							      \
    Ref ## T (T*a);							      \
    Ref ## T ( Ref ## T &a);						      \
    ~Ref ## T ();							      \
    int null();								      \
    int nonnull();							      \
    Ref ## T  operator=( Ref ## T & c);					      \
    Ref ## T  operator=(T* cr);						      \
    void assign_pointer(T* cr);						      \
    int operator==( Ref ## T &a);					      \
    int operator==( T * a);						      \
    int operator>=( Ref ## T &a);					      \
    int operator<=( Ref ## T &a);					      \
    int operator>( Ref ## T &a);					      \
    int operator<( Ref ## T &a);					      \
    void  ref_info(FILE*fp=stdout);					      \
    void warn(const char *);						      \
    void clear();							      \
    void check_pointer();						      \
    custom								      \
}

#define REF_dec(T) REF_dec_custom(T,)

#ifndef REF_CHECK_POINTER
#define REF_CHECK_POINTER 1
#endif

#ifndef REF_CHECK_STACK
#define REF_CHECK_STACK 1
#endif

#if REF_CHECK_STACK
#include <unistd.h>
#endif

#define REF_def(T)							      \
T* Ref ## T :: operator->() { return p; }				      \
const T* Ref ## T :: operator->() const { return p; }			      \
T* Ref ## T :: pointer() { return p; }					      \
const T* Ref ## T :: pointer() const { return p; }			      \
Ref ## T :: operator T*() { return p; };				      \
Ref ## T :: const operator T*() const { return p; };			      \
T& Ref ## T :: operator *() { return *p; };				      \
const T& Ref ## T :: operator *() const { return *p; };			      \
int Ref ## T :: null() { return p == 0; }				      \
int Ref ## T :: nonnull() { return p != 0; }				      \
int Ref ## T :: operator==( Ref ## T &a) { return p == a.p; }		      \
int Ref ## T :: operator==( T * a) { return p == a; }			      \
int Ref ## T :: operator>=( Ref ## T &a) { return p >= a.p; }		      \
int Ref ## T :: operator<=( Ref ## T &a) { return p <= a.p; }		      \
int Ref ## T :: operator>( Ref ## T &a) { return p > a.p; }		      \
int Ref ## T :: operator<( Ref ## T &a) { return p < a.p; }		      \
Ref ## T :: Ref ## T (): p(0) {}					      \
Ref ## T :: Ref ## T (T*a): p(a)					      \
{									      \
  if (REF_CHECK_STACK && (void*) p > sbrk(0)) {				      \
      warn("Ref" # T ": creating a reference to stack data");		      \
    }									      \
  if (p) p->reference();						      \
  if (REF_CHECK_POINTER) check_pointer();				      \
}									      \
Ref ## T :: Ref ## T ( Ref ## T &a): p(a.p)				      \
{									      \
  if (p) p->reference();						      \
  if (REF_CHECK_POINTER) check_pointer();				      \
}									      \
Ref ## T :: ~Ref ## T ()						      \
{									      \
  clear();								      \
}									      \
void									      \
Ref ## T :: clear()							      \
{									      \
  if (REF_CHECK_POINTER) check_pointer();				      \
  if (p && p->dereference()<=0) {					      \
      if (REF_CHECK_STACK && (void*) p > sbrk(0)) {			      \
          warn("Ref" # T ": skipping delete of object on the stack");	      \
        }								      \
      else {								      \
           delete p;							      \
         }								      \
    }									      \
  p = 0;								      \
}									      \
void									      \
Ref ## T :: warn ( const char * msg)					      \
{									      \
  fprintf(stderr,"WARNING: %s\n",msg);					      \
}									      \
Ref ## T  Ref ## T :: operator=( Ref ## T & c)				      \
{									      \
  if (c.p) c.p->reference();						      \
  clear();								      \
  p=c.p;								      \
  if (REF_CHECK_POINTER) check_pointer();				      \
  return *this;								      \
}									      \
Ref ## T  Ref ## T :: operator=(T* cr)					      \
{									      \
  if (cr) cr->reference();						      \
  clear();								      \
  p = cr;								      \
  if (REF_CHECK_POINTER) check_pointer();				      \
  return *this;								      \
}									      \
void Ref ## T :: assign_pointer(T* cr)					      \
{									      \
  if (cr) cr->reference();						      \
  clear();								      \
  p = cr;								      \
  if (REF_CHECK_POINTER) check_pointer();				      \
}									      \
void Ref ## T :: check_pointer()					      \
{									      \
  if (p && p->nreference() <= 0) {					      \
      warn("Ref" # T ": bad reference count in referenced object\n");	      \
    }									      \
}									      \
void Ref ## T :: ref_info(FILE*fp=stdout)				      \
{									      \
  if (nonnull()) fprintf(fp,"nreference() = %d\n",p->nreference());	      \
  else fprintf(fp,"reference is null\n");				      \
}

//template <class T> class Ref;
class RefCount {
    //template <class T> friend class Ref;
  public: // should be private, but template friends don't work
    int _reference_count_;
  protected:
    RefCount();
    RefCount(const RefCount&);
    RefCount& operator=(const RefCount&);
  public:
    ~RefCount();
    int reference();
    int nreference();
    int dereference();
};

//template <class T> class Ref;
class VRefCount {
    //template <class T> friend class Ref;
  public: // should be private, but template friends don't work
    int _reference_count_;
  protected:
    VRefCount();
    VRefCount(const VRefCount&);
    VRefCount& operator=(const VRefCount&);
  public:
    virtual ~VRefCount();
    int reference();
    int nreference();
    int dereference();
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

#ifdef INLINE_FUNCTIONS
#include <util/container/ref_i.h>
#endif

#endif

