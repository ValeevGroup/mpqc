
#include "ref.h"

#define BAD_REF 1

// The reference count is coded in such way that a bad references
// can be detected.
static int nreference(int rc)
{
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
static int reference(int& rc)
{
  int nref = nreference(rc);
  nref++;
  rc = refcount(nref);
  return nref;
}
static int dereference(int& rc)
{
  int nref = nreference(rc);
  nref--;
  if (nref < 0) {
      fprintf(stderr,"ref.cc:dereference: attempted to delete a reference"
              " to an unreferenced object\n");
    }
  rc = refcount(nref);
  return nref;
}

int
VRefCount::nreference()
{
  return ::nreference(_reference_count_);
}
int
VRefCount::reference()
{
  return ::reference(_reference_count_);
}
int
VRefCount::dereference()
{
  return ::dereference(_reference_count_);
}

int
RefCount::nreference()
{
  return ::nreference(_reference_count_);
}
int
RefCount::reference()
{
  return ::reference(_reference_count_);
}
int
RefCount::dereference()
{
  return ::dereference(_reference_count_);
}

// The DTORs try to catch bugs by setting the ref count to BAD_REF.  In
// principle, then, the Ref class and the RefCount classes can try
// to detect multiple deletes and warn the user.  Unfortunately,
// with at least the L486 arch, the BAD_REF get blasted away right after
// the DTOR is called.

VRefCount::~VRefCount()
{
  if (_reference_count_ > 0) {
      fprintf(stderr,"WARNING: VRefCount: deleting a referenced object\n");
    }
  else if (_reference_count_ < 0) {
      fprintf(stderr,"WARNING: VRefCount: deleting a deleted object\n");
    }
  _reference_count_ = BAD_REF;
}

RefCount::~RefCount()
{
  if (_reference_count_ > 0) {
      fprintf(stderr,"WARNING: RefCount: deleting a referenced object\n");
    }
  else if (_reference_count_ < 0) {
      fprintf(stderr,"WARNING: RefCount: deleting a deleted object\n");
    }
  _reference_count_ = BAD_REF;
}
