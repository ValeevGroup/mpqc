
/* state_xdr.cc -- implementation of byte-swaping StateOut classes
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
 *      May, 1993
 */

#include <util/class/class.h>
#include "state.h"

///////////////////////////////////////////////////////////////////

StateOutBinXDR::StateOutBinXDR() :
  StateOutBin()
{
}

StateOutBinXDR::StateOutBinXDR(FILE* fp) :
  StateOutBin(fp)
{
}

StateOutBinXDR::StateOutBinXDR(const char *path, const char * mode) :
  StateOutBin(path,mode)
{
}

StateOutBinXDR::~StateOutBinXDR()
{
}

int StateOutBinXDR::put_array_void(const void*v,int i)
{
  return StateOutBin::put_array_void(v,i);
}

int StateOutBinXDR::put_array_char(const char *p, int size)
{
  translate((char*)p,size);
  int ret = put_array_void((void*)p,size*sizeof(char));
  translate((char*)p,size);
  return ret;
}

int StateOutBinXDR::put_array_int(const int *p, int size)
{
  translate((unsigned int*)p,size);
  int ret = put_array_void((void*)p,size*sizeof(int));
  translate((unsigned int*)p,size);
  return ret;
}

int StateOutBinXDR::put_array_float(const float *p, int size)
{
  translate((float*)p,size);
  int ret = put_array_void((void*)p,size*sizeof(float));
  translate((float*)p,size);
  return ret;
}

int StateOutBinXDR::put_array_double(const double *p, int size)
{
  translate((double*)p,size);
  int ret = put_array_void((void*)p,size*sizeof(double));
  translate((double*)p,size);
  return ret;
}

//////////////////////////////////////////////////////////////

StateInBinXDR::StateInBinXDR() :
  StateInBin()
{
}

StateInBinXDR::StateInBinXDR(FILE* fp) :
  StateInBin(fp)
{
}

StateInBinXDR::StateInBinXDR(const char *path, const char * mode) :
  StateInBin(path,mode)
{
}

StateInBinXDR::~StateInBinXDR()
{
}

int StateInBinXDR::get_array_void(void*v,int i)
{
  return StateInBin::get_array_void(v,i);
}

int StateInBinXDR::get_array_char(char *p, int size)
{
  int ret = get_array_void((void*)p,size*sizeof(char));
  translate(p,size);
  return ret;
}

int StateInBinXDR::get_array_int(int *p, int size)
{
  int ret = get_array_void((void*)p,size*sizeof(int));
  translate((unsigned int*)p,size);
  return ret;
}

int StateInBinXDR::get_array_float(float *p, int size)
{
  int ret = get_array_void((void*)p,size*sizeof(float));
  translate(p,size);
  return ret;
}

int StateInBinXDR::get_array_double(double *p, int size)
{
  int ret = get_array_void((void*)p,size*sizeof(double));
  translate(p,size);
  return ret;
}
