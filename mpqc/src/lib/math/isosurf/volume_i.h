
#ifdef __GNUC__
#pragma interface
#endif

#ifdef INLINE_FUNCTIONS
#define INLINE inline
#else
#define INLINE
#endif

INLINE double&
Volume::interpolation_accuracy()
{
  return _interp_acc;
};
