
/* bitarray.h -- definition of the BitArray Class
 *
 *      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
 *      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
 *      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
 *      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
 *      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
 *      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
 *
 *  Author:
 *      E. T. Seidl
 *      Bldg. 12A, Rm. 2033
 *      Computer Systems Laboratory
 *      Division of Computer Research and Technology
 *      National Institutes of Health
 *      Bethesda, Maryland 20892
 *      Internet: seidl@alw.nih.gov
 *      July, 1993
 */

#ifndef _math_topology_bitarray_h
#define _math_topology_bitarray_h

#include <string.h>
#include <stdlib.h>
#include <math.h>

//
// class BitArray is used as the lower triangle of a boolean matrix.  rather
// than storing an int or a char, just use one bit for each, so instead
// of n(n+1)/2 bytes of storage you have n(n+1)/16 bytes.  A further savings
// of n bits could be obtained by setting the diagonal to always true or always
// false depending on the application, but this would probably be more
// expensive computationally than it's worth.
//

class BitArray {
  private:
    unsigned char *a;
    int n;
    int nm;
    int na;

    inline unsigned int ioff(unsigned int i,unsigned int j) const {
      return (i>j) ? (((i*(i+1)) >> 1) + j) : (((j*(j+1)) >> 1) + i);
      }

  public:
    BitArray(int =0);
    BitArray(int =0, int =0);
    ~BitArray();

    inline void set(unsigned int i) { a[(i>>3)] |= (1 << (i&7)); }
    inline void set(unsigned int i, unsigned int j) { set(ioff(i,j)); }

    inline int is_set(unsigned int i, unsigned int j) const
      { int ij = ioff(i,j); return (a[(ij>>3)] & (1 << (ij&7))); }
    inline int is_set(unsigned int i) const
      { return (a[(i>>3)] & (1 << (i&7))); }

    inline int operator()(unsigned int i, unsigned int j) const
      { int ij = ioff(i,j); return (a[(ij>>3)] & (1 << (ij&7))); }
    inline int operator()(unsigned int i) const
      { return (a[(i>>3)] & (1 << (i&7))); }
    inline int operator[](unsigned int i) const
      { return (a[(i>>3)] & (1 << (i&7))); }

    inline int dim() const { return na; }
    inline int nrow() const { return nm; }
    inline int ncol() const { return nm; }

    inline int degree(unsigned int i) const {
      int nedge=0;
      for (int j=0; j < nm; j++) if ((*this)(i,j)) nedge++;
      return nedge;
      }
};


inline BitArray::BitArray(int sz)
  : a(0), n(0), nm(0), na(sz)
{
  if(sz) {
    n = sz/8 + ((sz%8)?1:0);
    a = new unsigned char[n];
    memset(a,'\0',n);
    nm = (int)((sqrt((double)1+8*sz) - 1)/2);
  }
}

inline BitArray::BitArray(int i, int j)
  : a(0), n(0), nm(0), na(0)
{
  if (i!=j) {
    fprintf(stderr,"BitArray(int,int): i != j\n");
    exit(1);
  }
  int sz = i*(i+1)/2;

  if(sz) {
    n = sz/8 + ((sz%8)?1:0);
    a = new unsigned char[n];
    memset(a,'\0',n);
    na=sz; nm=i;
  }
}

inline BitArray::~BitArray()
{
  if (a) delete[] a; a=0; n=0;
}

#endif
