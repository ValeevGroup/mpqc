#ifndef tstfcn_h
#define tstfcn_h

#include <math.h>
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))


// Test functions

void rosen(int, int, ColumnVector, double&, ColumnVector&, SymmetricMatrix&);
void quad(int, int, ColumnVector, double&, ColumnVector&, SymmetricMatrix&);
void quad2(int, int, ColumnVector, double&, ColumnVector&, SymmetricMatrix&);
void erosen(int, int, ColumnVector, double&, ColumnVector&, SymmetricMatrix&);
void ds_ex71(int, int, ColumnVector, double&, ColumnVector&, SymmetricMatrix&);
void ds_ex9(int, int, ColumnVector, double&, ColumnVector&, SymmetricMatrix&);
void wood(int, int, ColumnVector, double&, ColumnVector&, SymmetricMatrix&);


#endif
