
#ifndef _math_nihmatrix_newmat_h
#define _math_nihmatrix_newmat_h

// newmat compat. routines
void Convert(DVector&in,ColumnVector&out);
void Convert(DMatrix&in,SymmetricMatrix&out);
void Convert(ColumnVector&in,DVector&out);
void Convert(SymmetricMatrix&in,DMatrix&out);

#endif
