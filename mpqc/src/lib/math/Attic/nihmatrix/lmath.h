
#ifndef lmath_h
#define lmath_h

#include "nihmatrix.h"

void mathqc_diag(DMatrix&, DVector&, DMatrix&, int =1, double =1.0e-15);
void mathqc_diag(DMatrix*, DVector*, DMatrix*, int =1, double =1.0e-15);

double mathqc_invert(DMatrix&);
double mathqc_invert(DMatrix*);

double mathqc_solve_lin(DMatrix*,DVector*);
double mathqc_solve_lin(DMatrix&,DVector&);

void mathqc_mxm(DMatrix&,int,DMatrix&,int,DMatrix&,int,int,int,int,int);

#endif
