
#ifdef __GNUC__
#pragma interface
#endif

#ifdef INLINE_FUNCTIONS
#define INLINE inline
#else
#define INLINE
#endif

INLINE
Result::Result(Compute*c):
  _compute(0),_computed(0),_c(c)
{
  c->add(this);
}

INLINE int&
Result::compute()
{
  return _compute;
}

INLINE int
Result::compute(int c)
{
  int r = _compute; _compute = c; return r;
}

INLINE int&
Result::computed()
{
  return _computed;
}

INLINE int
Result::needed()
{
  return _compute && (!_computed);
}

#undef INLINE
