
#ifdef ALLOC_BOUND_GLOBALS
#define EXTERN
#else
#define EXTERN extern
#endif

typedef signed char int_bound_t;

EXTERN int_bound_t int_Q;
EXTERN int_bound_t int_R;
EXTERN int_bound_t *int_Qvec;
EXTERN int_bound_t *int_Rvec;

#undef EXTERN

