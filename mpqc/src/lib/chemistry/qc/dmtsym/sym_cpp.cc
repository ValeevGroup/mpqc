
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <math/symmetry/pointgrp.h>
#include <util/misc/libmisc.h>
#include <util/keyval/keyval.h>
#include <util/keyval/ipv2c.h>

extern "C" {
#include <tmpl.h>
#include <math/array/math_lib.h>
}

#include <chemistry/qc/intv2/int_libv2.h>

extern "C" {
#include "symm_mac.h"
#include "symm.h"

#include "symmzero.h"
#include "symmallc.h"

#include "syminit.gbl"
}

////////////////////////////////////////////////////////////////////////////
//
// given a reference to a point group and an initialized centers struct
// containing all atoms (not just the unique ones), fill in all the info
// in sym_info
//
int
sym_struct_from_pg(const PointGroup& pg, centers_t& centers,
                   sym_struct_t& sym_info)
{
  double_array3_t trans;

  CharacterTable ct = pg.char_table();

  if (allocbn_double_array3(&trans,"n1 n2 n3",ct.order(),3,3) != 0) {
    fprintf(stderr,"sym_struct_from_pg: could not alloc trans\n");
    return -1;
  }

  for (int g=0; g < ct.order(); g++) {
    SymmetryOperation so = ct.symm_operation(g);

    for (int i=0; i < 3; i++)
      for (int j=0; j < 3; j++)
        trans.d[g][i][j] = so(i,j);
  }

  char *pgrp = strdup(pg.symbol());

  if (sym_make_sym_struct(&centers,&sym_info,pgrp,&trans) < 0) {
    fprintf(stderr,"sym_struct_from_pg: sym_make_sym_struct failed\n");
    return -1;
  }

  free(pgrp);

  return 0;
}

///////////////////////////////////////////////////////////////////////////
//
// given a keyval, read the input to get the unique centers, and then
// calculate all atoms and place in centers, and also initialize sym_info
//

int
sym_init_centers(KeyVal& keyval, centers_t& centers, sym_struct_t& sym_info)
{
  centers_t unique_centers;
  int errcod;
  int nat;
  int i,j;

  char *point_group = keyval.pcharvalue("symmetry");

  int_read_centers(keyval,unique_centers);

  errcod =
    sym_init_given_centers(&unique_centers,&centers,&sym_info,point_group);

  if (errcod < 0) fprintf(stderr,"sym_init_centers: could not init centers\n");

  free_centers(&unique_centers);

  return errcod;
}
