
#include <stdio.h>
#include <util/sgen/sgen.h>
#include "atoms.h"
#include "atomsfwr.h"

/* Generated fwrite function for shell_type: */
int
fwrite_shell_type(_fp,_shell_type,_offset)
FILE *_fp;
shell_type_t *_shell_type;
int *_offset;
{
  typedef int boolean;
  typedef char * string;
  int _offset_init= *_offset;
  fwrite_int(_fp,&(_shell_type->am),_offset,
    sizeof(int)*1);
  fwrite_int(_fp,&(_shell_type->puream),_offset,
    sizeof(int)*1);
  return(*_offset-_offset_init);
  }

/* Generated fwrite function for basis: */
int
fwrite_basis(_fp,_basis,_offset)
FILE *_fp;
basis_t *_basis;
int *_offset;
{
  typedef int boolean;
  typedef char * string;
  int i;
  int _offset_init= *_offset;
  fwrite_int(_fp,&(_basis->n),_offset,
    sizeof(int)*1);
  fwrite_string(_fp,&(_basis->name),_offset,
    sizeof(string)*1);
  if(_basis->n!=0) {
    fwrite_pointer(_fp,&(_basis->shell),_offset,sizeof(shell_t *));
    if(_basis->shell!=NULL) {
      for (i=0; i<_basis->n; i++)  {
        fwrite_shell(_fp,&(_basis->shell[i]),_offset);
        }
      }
    }
  return(*_offset-_offset_init);
  }

/* Generated fwrite function for center: */
int
fwrite_center(_fp,_center,_offset)
FILE *_fp;
center_t *_center;
int *_offset;
{
  typedef int boolean;
  typedef char * string;
  int _offset_init= *_offset;
  fwrite_string(_fp,&(_center->atom),_offset,
    sizeof(string)*1);
  fwrite_double(_fp,&(_center->charge),_offset,
    sizeof(double)*1);
  if(3!=0) {
    fwrite_pointer(_fp,&(_center->r),_offset,sizeof(double *));
    if(_center->r!=NULL) {
      fwrite_double(_fp,(_center->r),_offset,
        sizeof(double)*3);
      }
    }
  fwrite_basis(_fp,&(_center->basis),_offset);
  return(*_offset-_offset_init);
  }

/* Generated fwrite function for centers: */
int
fwrite_centers(_fp,_centers,_offset)
FILE *_fp;
centers_t *_centers;
int *_offset;
{
  typedef int boolean;
  typedef char * string;
  int i;
  int _offset_init= *_offset;
  fwrite_int(_fp,&(_centers->n),_offset,
    sizeof(int)*1);
  if(_centers->n!=0) {
    fwrite_pointer(_fp,&(_centers->center),_offset,sizeof(center_t *));
    if(_centers->center!=NULL) {
      for (i=0; i<_centers->n; i++)  {
        fwrite_center(_fp,&(_centers->center[i]),_offset);
        }
      }
    }
  fwrite_int(_fp,&(_centers->shell_offset),_offset,
    sizeof(int)*1);
  fwrite_int(_fp,&(_centers->prim_offset),_offset,
    sizeof(int)*1);
  fwrite_int(_fp,&(_centers->func_offset),_offset,
    sizeof(int)*1);
  fwrite_int(_fp,&(_centers->nshell),_offset,
    sizeof(int)*1);
  fwrite_int(_fp,&(_centers->nprim),_offset,
    sizeof(int)*1);
  fwrite_int(_fp,&(_centers->nfunc),_offset,
    sizeof(int)*1);
  if(_centers->nshell!=0) {
    fwrite_pointer(_fp,&(_centers->center_num),_offset,sizeof(int *));
    if(_centers->center_num!=NULL) {
      fwrite_int(_fp,(_centers->center_num),_offset,
        sizeof(int)*_centers->nshell);
      }
    }
  if(_centers->nshell!=0) {
    fwrite_pointer(_fp,&(_centers->shell_num),_offset,sizeof(int *));
    if(_centers->shell_num!=NULL) {
      fwrite_int(_fp,(_centers->shell_num),_offset,
        sizeof(int)*_centers->nshell);
      }
    }
  if(_centers->nshell!=0) {
    fwrite_pointer(_fp,&(_centers->func_num),_offset,sizeof(int *));
    if(_centers->func_num!=NULL) {
      fwrite_int(_fp,(_centers->func_num),_offset,
        sizeof(int)*_centers->nshell);
      }
    }
  return(*_offset-_offset_init);
  }
