//$$ newmatio.h           definition file for matrix package input/output

// Copyright (C) 1991,2,3: R B Davies

#ifndef MATRIXIO_LIB
#define MATRIXIO_LIB 0

#include <math/newmat7/newmat.h>

/**************************** input/output *****************************/

ostream& operator<<(ostream&, const BaseMatrix&);

ostream& operator<<(ostream&, const GeneralMatrix&);


#endif
