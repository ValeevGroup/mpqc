
#ifdef __GNUC__
#pragma interface
#endif

#ifdef INLINE_FUNCTIONS
#define INLINE inline
#else
#define INLINE
#endif

///////////////////////////////////////////////////////////////////////////

INLINE int
CartesianIter::n()
{
  return ((l_>=0)?((((l_)+2)*((l_)+1))>>1):0);
}

INLINE int
CartesianIter::a()
{
  return a_;
}

INLINE int
CartesianIter::b()
{
  return b_;
}

INLINE int
CartesianIter::c()
{
  return c_;
}

INLINE int
CartesianIter::l()
{
  return l_;
}

INLINE int
CartesianIter::l(int i)
{
  return i ? (i==1 ? b_ : c_) : a_;
}

INLINE int
CartesianIter::bfn()
{
  return bfn_;
}

///////////////////////////////////////////////////////////////////////////

INLINE void
RedundantCartesianIter::start()
{
  if (l_==0)
    done_ = 1;
  else
    done_ = 0;

  for (int i=0; i<l_; i++)
    axis_[i] = 0;
}

INLINE void
RedundantCartesianIter::next()
{
  for (int i=0; i<l_; i++) {
    if (axis_[i] == 2)
      axis_[i] = 0;
    else {
      axis_[i]++;
      return;
    }
  }
  done_ = 1;
}

INLINE 
RedundantCartesianIter::operator int()
{
  return !done_;
}

INLINE int
RedundantCartesianIter::l()
{
  return l_;
}

INLINE int
RedundantCartesianIter::a()
{
  return l(0);
}

INLINE int
RedundantCartesianIter::b()
{
  return l(1);
}

INLINE int
RedundantCartesianIter::c()
{
  return l(2);
}

INLINE int
RedundantCartesianIter::l(int axis)
{
  int i;
  int r = 0;
  for (i=0; i<l_; i++) if (axis_[i]==axis) r++;
  return r;
}

INLINE int
RedundantCartesianIter::axis(int i)
{
  return axis_[i];
}

///////////////////////////////////////////////////////////////////////////

INLINE
RedundantCartesianSubIter::operator int()
{
  return !done_;
}

INLINE
RedundantCartesianSubIter::l()
{
  return l_;
}

INLINE
RedundantCartesianSubIter::a()
{
  return e_[0];
}

INLINE
RedundantCartesianSubIter::b()
{
  return e_[1];
}

INLINE
RedundantCartesianSubIter::c()
{
  return e_[2];
}

INLINE
RedundantCartesianSubIter::l(int i)
{
  return e_[i];
}

INLINE
RedundantCartesianSubIter::axis(int i)
{
  return axis_[i];
}
