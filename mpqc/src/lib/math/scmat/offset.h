
#ifndef _math_scmat_offset_h
#define _math_scmat_offset_h

#ifdef __c_plus_plus

static inline int
i_offset(int i)
{
  return ((i*(i+1)) >> 1);
}

static inline int
ij_offset(int i, int j)
{
  return (i>j) ? (((i*(i+1)) >> 1) + j) : (((j*(j+1)) >> 1) + i);
}

static inline int
igtj_offset(int i, int j)
{
  return ((i*(i+1)) >> 1) + j;
}

#else

#define i_offset(i) (((i)*((i)+1))>>1)
#define ij_offset(i,j) (((i)>(j))?(i_offset(i)+(j)):(i_offset(j)+(i)))
#define igtj_offset(i,j) (i_offset(i)+(j))

#endif

#endif
