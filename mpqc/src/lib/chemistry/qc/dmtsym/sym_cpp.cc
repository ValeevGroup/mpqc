
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <math/symmetry/pointgrp.h>
#include <util/misc/libmisc.h>

extern "C" {
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>

#include "symm_mac.h"
#include "symm.h"

#include "symmzero.h"
#include "symmallc.h"

#include "syminit.gbl"
}

int
sym_struct_from_pg(const PointGroup& pg, centers_t& centers,
                   sym_struct_t& sym_info)
{
  double_array3_t trans;

  if (sym_make_sym_struct(&centers,&sym_info,pg.symbol(),&trans) < 0) {
    fprintf(stderr,"sym_struct_from_pg: sym_make_sym_struct failed\n");
    return -1;
  }

  return 0;
}
