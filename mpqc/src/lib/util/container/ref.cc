
#ifdef __GNUC__
#pragma implementation
#endif

#include "ref.h"

#define BAD_REF 1

static void unmanage(int&);
static int managed(int);

// The reference count is coded in such way that a bad references
// can be detected.
static int nreference(int rc)
{
  if (!managed(rc)) return 1;
  int nref1 = rc&0x7fff;
  int nref2 = (rc>>16)&0x7fff;
  if (nref1 != nref2) {
      fprintf(stderr,"WARNING: ::nreference: bad reference count\n");
      return -1;
    }
  return nref1;
}
static int refcount(int nref)
{
  if (nref > 0x7fff) {
      fprintf(stderr,"WARNING: ::refcount: too many references to object\n");
      nref = 0x7fff;
    }
  int rc = nref | (nref<<16);
  return rc;
}
static int reference(const int& constrc)
{
  // cast away the constness:
  int& rc = (int) constrc;
  if (!managed(rc)) return 1;
  int nref = nreference(rc);
  nref++;
  rc = refcount(nref);
  return nref;
}
static int dereference(int& rc)
{
  if (!managed(rc)) return 1;
  int nref = nreference(rc);
  nref--;
  if (nref < 0) {
      fprintf(stderr,"ref.cc:dereference: attempted to delete a reference"
              " to an unreferenced object\n");
    }
  rc = refcount(nref);
  return nref;
}
static void
unmanage(int&rc)
{
  if (!managed(rc)) return;
  if (nreference(rc)) {
      fprintf(stderr,"[V]RefCount::unmanage: cannot unmanage a managed"
              " object with a nonzero nreference\n");
      abort();
    }
  rc = 0x7fffffff;
}
int
managed(int rc)
{
  return rc != 0x7fffffff;
}

int
VRefCount::nreference() const
{
  return ::nreference(_reference_count_);
}
int
VRefCount::reference() const
{
  return ::reference(_reference_count_);
}
int
VRefCount::dereference()
{
  return ::dereference(_reference_count_);
}
int
VRefCount::managed() const
{
  return ::managed(_reference_count_);
}
void
VRefCount::unmanage()
{
  ::unmanage(_reference_count_);
}

int
RefCount::nreference() const
{
  return ::nreference(_reference_count_);
}
int
RefCount::reference() const
{
  return ::reference(_reference_count_);
}
int
RefCount::dereference()
{
  return ::dereference(_reference_count_);
}
int
RefCount::managed() const
{
  return ::managed(_reference_count_);
}
void
RefCount::unmanage()
{
  ::unmanage(_reference_count_);
}

// The DTORs try to catch bugs by setting the ref count to BAD_REF.  In
// principle, then, the Ref class and the RefCount classes can try
// to detect multiple deletes and warn the user.  Unfortunately,
// with at least the L486 arch, the BAD_REF get blasted away right after
// the DTOR is called.

VRefCount::~VRefCount()
{
  if (_reference_count_ == BAD_REF) {
      fprintf(stderr,"WARNING: VRefCount: deleting a deleted object\n");
    }
  else if (managed() && nreference()) {
      fprintf(stderr,"WARNING: VRefCount: deleting a referenced object\n");
    }
  _reference_count_ = BAD_REF;
}

RefCount::~RefCount()
{
  if (_reference_count_ == BAD_REF) {
      fprintf(stderr,"WARNING: RefCount: deleting a deleted object\n");
    }
  else if (managed() && nreference()) {
      fprintf(stderr,"WARNING: RefCount: deleting a referenced object\n");
    }
  _reference_count_ = BAD_REF;
}
