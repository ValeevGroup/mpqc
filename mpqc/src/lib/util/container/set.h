
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_container_set_h
#define _util_container_set_h

#include <Pix.h>
#include <stdlib.h>
#include <util/container/array.h>

#include <util/container/settmpl.h> // The template set declaration.
#include <util/container/setmacr.h> // The macro set declaration.
#define SET_dec(Type) Set_declare(Type)
#define SET_def(Type)

// This class combines array capabilities with set capabilities.
// It is basically a set with the iseek and operator[](int i) members
// added.  At the moment, it requires that Set use an array internally.
// When Set is improved, Arrayset must be updated.
#include <util/container/asettmpl.h> // The template set declaration.
#include <util/container/asetmacr.h> // The macro set declaration.
#define ARRAYSET_dec(Type) Arrayset_declare(Type)
#define ARRAYSET_def(Type)

// declare sets for the basic types
SET_dec(int);
ARRAYSET_dec(int);
ARRAY_dec(Arraysetint);
SET_dec(double);
ARRAYSET_dec(double);

#endif
