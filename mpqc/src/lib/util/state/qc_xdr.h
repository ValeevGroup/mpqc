
/* xdr.h -- definition of the external data representation classes
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
 *      Internet: seidl@janed.com
 *      February, 1993
 */

/*
 * as I understand sun's xdr format, things are translated to big-endian
 * byte ordering, which is the byte ordering on sparcs.  For ints we've
 * got it made because we can just use the standard inet ntoh and hton
 * macros/functions.  We'll have to write our own float and double
 * conversions
 */

#ifndef _libQC_xdr_h
#define _libQC_xdr_h

#ifdef HAVE_CONFIG_H
#include <scconfig.h>
#endif

extern "C" {
#if defined(__GNUC__) || defined(sgi)
#include <limits.h>  // use the gnu limits.h
#else
#include <sys/limits.h>
#endif
}

#if defined(WORDS_BIGENDIAN)
#define BIGENDIAN 1
#else
#define BIGENDIAN 0
#endif

//. \clsnm{QCXDR} is used to convert data to and from QCXDR format.
//. Real XDR is not used---standard floating point and integer
//. data formats must be used.  This only gets around endianness
//. differences.
class QCXDR 
{
  private:
    // do not allow copy constructor or assignment
    QCXDR(const QCXDR&);
    operator=(const QCXDR&);
  public:
    QCXDR();
    //. These perform the translation to/from xdr format
    //. single element translations---the value of the passed variable is not
    //. changed.
    char   translate(char);
    short  translate(unsigned short);
    int    translate(unsigned int);
    long   translate(unsigned long);
    void   translate(float*);
    void   translate(double*);
    void*  translate(void*); // for pointers

    //. Array transformations---the elements of the array are changed.
    void   translate(char*,int);
    void   translate(unsigned short*,int);
    void   translate(unsigned int*,int);
    void   translate(unsigned long*,int);
    void   translate(float*,int);
    void   translate(double*,int);

    //. Functions for actually performing byte swapping.
    char   byte_swap(char);
    short  byte_swap(unsigned short);
    int    byte_swap(unsigned int);
    long   byte_swap(unsigned long);
    void   byte_swap(float*);
    void   byte_swap(double*);
    void*  byte_swap(void*); // for pointers
};

#ifdef USE_INLINE
#include <util/state/qc_xdr_i.h>
#endif

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
// End:
