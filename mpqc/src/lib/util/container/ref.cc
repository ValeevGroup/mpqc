
#ifdef __GNUC__
#pragma implementation
#endif

#include <util/container/ref.h>

void
VRefCount::error(const char * w) const
{
  fprintf(stderr,"VRefCount: ERROR: %s\n",w);
  abort();
}

void
VRefCount::too_many_refs() const
{
  error("Too many refs.");
}
void
VRefCount::not_enough_refs() const
{
  error("Ref count dropped below zero.");
}

#if REF_CHECKSUM
void
VRefCount::bad_checksum() const
{
  error("Bad checksum.");
}
#endif

VRefCount::~VRefCount()
{
#if REF_CHECKSUM
  if (_reference_count_ == 0 && _checksum_ == 0) {
      error("Deleting a deleted or overwritten object.");
    }
#endif
#if REF_MANAGE
  if (managed() && nreference()) {
      error("Deleting a referenced object.");
    }
#endif
#if REF_CHECKSUM
  _reference_count_ = 0;
  _checksum_ = 0;
#endif
}

///////////////////////////////////////////////////////////////////////

void
RefBase::warn ( const char * msg) const
{
  fprintf(stderr,"WARNING: %s\n",msg);
}
void
RefBase::warn_ref_to_stack() const
{
  warn("Ref: creating a reference to stack data");
}
void
RefBase::warn_skip_stack_delete() const
{
  warn("Ref: skipping delete of object on the stack");
}
void
RefBase::warn_bad_ref_count() const
{
  warn("Ref: bad reference count in referenced object\n");
}
void
RefBase::ref_info(VRefCount*p,FILE*fp) const
{
  if (p) fprintf(fp,"nreference() = %d\n",p->nreference());
  else fprintf(fp,"reference is null\n");
}
