
#ifndef _util_container_ref_h
#define _util_container_ref_h

#include <stdio.h>
#include <stdlib.h>

#ifndef SUNMOS
#ifndef REF_CHECK_POINTER
#define REF_CHECK_POINTER 1
#endif
#endif

#ifndef REF_CHECK_STACK
#define REF_CHECK_STACK 1
#endif

#if REF_CHECK_STACK
#include <unistd.h>
extern "C" void * sbrk(int);
#define DO_REF_CHECK_STACK(p) (((void*) (p) > sbrk(0)) && (p)->managed())
#else
#define DO_REF_CHECK_STACK(p) (0)
#endif

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
    Ref ## T& operator=( Ref ## T & c);					      \
    Ref ## T& operator=(T* cr);						      \
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
  if (DO_REF_CHECK_STACK(p)) {						      \
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
      if (DO_REF_CHECK_STACK(p)) {					      \
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
Ref ## T& Ref ## T :: operator=( Ref ## T & c)				      \
{									      \
  if (c.p) c.p->reference();						      \
  clear();								      \
  p=c.p;								      \
  if (REF_CHECK_POINTER) check_pointer();				      \
  return *this;								      \
}									      \
Ref ## T& Ref ## T :: operator=(T* cr)					      \
{									      \
  assign_pointer(cr);				\
  return *this;								      \
}									      \
void Ref ## T :: assign_pointer(T* cr)					      \
{						\
  if (cr) {					\
      if (DO_REF_CHECK_STACK(cr)) {		\
          warn("Ref" # T ": creating a reference to stack data"); \
        }					\
      cr->reference();				\
    }						\
  clear();					\
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

class RefCount {
  private:
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
    int managed();
    void unmanage();
};

class VRefCount {
  private:
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

    // unmanaged objects always return 1 for reference counts
    int managed();
    void unmanage();
};

#ifdef INLINE_FUNCTIONS
#include <util/container/ref_i.h>
#endif

#endif

