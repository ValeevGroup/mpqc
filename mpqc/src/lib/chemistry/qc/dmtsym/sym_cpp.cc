
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <math/symmetry/pointgrp.h>
#include <util/misc/libmisc.h>
#include <util/keyval/keyval.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/rot.h>

extern "C" {
#include <tmpl.h>
#include <math/array/math_lib.h>
}

#include <chemistry/qc/intv2/int_libv2.h>

extern "C" {
#include <chemistry/qc/dmtsym/symm_mac.h>
#include <chemistry/qc/dmtsym/symm.h>

#include <chemistry/qc/dmtsym/symmzero.h>
#include <chemistry/qc/dmtsym/symmallc.h>

#include <chemistry/qc/dmtsym/syminit.gbl>
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
  free_double_array3(&trans);

  return 0;
}

///////////////////////////////////////////////////////////////////////////
//
// given a keyval, read the input to get the unique centers, and then
// calculate all atoms and place in centers, and also initialize sym_info
//

int
sym_init_centers(const RefKeyVal& keyval, centers_t& centers,
                 sym_struct_t& sym_info)
{
  centers_t unique_centers;
  int errcod;

  char *point_group = keyval->pcharvalue("symmetry");

  int_read_centers(keyval,unique_centers);

  errcod =
    sym_init_given_centers(&unique_centers,&centers,&sym_info,point_group);

  if (errcod < 0) fprintf(stderr,"sym_init_centers: could not init centers\n");

  free_centers(&unique_centers);

  return errcod;
}

////////////////////////////////////////////////////////////////////////////
//
// given a reference to a basis set, fill in all the info in sym_info
//
int
sym_struct_from_gbs(const RefGaussianBasisSet& gbs, sym_struct_t& sym_info)
{
  int i, g;
  int n, nirr;
  int f_exist=0, g_exist=0, h_exist=0;
  enum pgroups ptgrp_sym;
  
  // first let's see if there are h functions or higher am
  int nsh = gbs->nshell();
  int nat = gbs->ncenter();
  int nshtr = nsh*(nsh+1)/2;
  
  for (i=0; i < nat; i++) {
    for (int s=0; s < gbs->nshell_on_center(i); s++) {
      int maxam = gbs->operator()(i,s).max_am();
      switch(maxam) {
      case 5:
        fprintf(stderr,"sym_struct_from_gbs: cannot handle am >= 5\n");
        return -1;

      case 4:
        g_exist=1;
        break;

      case 3:
        f_exist=1;
        // fall through
      default:
        break;
      }
    }
  }

  PointGroup& pg = gbs->molecule()->point_group();
  CharacterTable ct = pg.char_table();
  PetiteList pl(gbs);
  
  char *point_group = strdup(pg.symbol());
  
  int errcod = sym_parse_symbol(point_group, &g, &ptgrp_sym, &n, &nirr);
  if (errcod || g != ct.order() || nirr != ct.nirrep()) {
    fprintf(stderr,"sym_struct_from_gbs: trouble parsing point group\n");
    return -1;
  }

  free(point_group);

  // now allocate memory for sym_info
  // not much to do for C1
  if (ct.order()==1) {
    errcod =
      allocbn_sym_struct(&sym_info, "point_group pg", pg.symbol(), ptgrp_sym);

  } else if (g_exist) {
    errcod = allocbn_sym_struct(&sym_info,
         "natom g nshell nshtri point_group pg nirrep n1 n2 n2p n3 n3p n4 n4p",
         nat, g, nsh, nshtr, pg.symbol(), ptgrp_sym, nirr,
         3, 6, 5, 10, 7, 15, 9);

  } else if (f_exist) {
    errcod = allocbn_sym_struct(&sym_info,
         "natom g nshell nshtri point_group pg nirrep n1 n2 n2p n3 n3p",
         nat , g, nsh, nshtr, pg.symbol(), ptgrp_sym, nirr,
         3, 6, 5, 10, 7);

  } else {
    errcod = allocbn_sym_struct(&sym_info,
                       "natom g nshell nshtri point_group pg nirrep n1 n2 n2p",
                       nat, g, nsh, nshtr, pg.symbol(),
                       ptgrp_sym, nirr, 3, 6, 5);
  }
    
  if (errcod != 0) {
    fprintf(stderr,"sym_struct_from_gbs:  cannot alloc sym_struct\n");
    return -1;
  }

  zero_sym_struct(&sym_info);

  if (g > 1) {
    if ((errcod = allocbn_char_tab(&sym_info.ct,"nirrep",nirr)) != 0) {
      fprintf(stderr,"sym_struct_from_gbs: "
                     "could not allocate memory for char_tab");
      return -1;
    }

    for (i=0; i < sym_info.natom; i++)
      for (g=0; g < sym_info.g; g++)
        sym_info.atom_map[i][g] = pl.atom_map(i,g);
        
    for (i=0; i < sym_info.nshell; i++) {
      sym_info.p1[i] = pl.in_p1(i);
      
      for (g=0; g < sym_info.g; g++)
        sym_info.shell_map[i][g] = pl.shell_map(i,g);
    }
        
    for (i=0; i < nshtr; i++)
      sym_info.lamij[i] = pl.lambda(i);
    
    int nb=0;
    for (i=0; i < nat; i++) {
      sym_info.first[i] = nb;
      for (int s=0; s < gbs->nshell_on_center(i); s++)
        nb += gbs->operator()(i,s).nfunction();
    }

    // don't mess with sym_info.firstp...it's not used anymore

    // now for the rotation matrices
    for (g=0; g < ct.order(); g++) {
      Rotation rp(1, ct.symm_operation(g));
      Rotation rd(2, ct.symm_operation(g));
      Rotation rf(3, ct.symm_operation(g));
      Rotation rg(4, ct.symm_operation(g));

      for (i=0; i < 3; i++)
        for (int j=0; j < 3; j++)
          sym_info.Rp[g][i][j] = rp(i,j);
      
      for (i=0; i < 6; i++)
        for (int j=0; j < 6; j++)
          sym_info.Rd[g][i][j] = rd(i,j);

      double xynorm = sqrt(3.0);
      double xynormi = 1.0/xynorm;
      sym_info.Rd[g][3][5] *= xynorm;
      sym_info.Rd[g][3][0] *= xynorm;
      sym_info.Rd[g][3][2] *= xynorm;
      sym_info.Rd[g][4][5] *= xynorm;
      sym_info.Rd[g][4][0] *= xynorm;
      sym_info.Rd[g][4][2] *= xynorm;
      sym_info.Rd[g][1][5] *= xynorm;
      sym_info.Rd[g][1][0] *= xynorm;
      sym_info.Rd[g][1][2] *= xynorm;

      sym_info.Rd[g][5][3] *= xynormi;
      sym_info.Rd[g][0][3] *= xynormi;
      sym_info.Rd[g][2][3] *= xynormi;
      sym_info.Rd[g][5][4] *= xynormi;
      sym_info.Rd[g][0][4] *= xynormi;
      sym_info.Rd[g][2][4] *= xynormi;
      sym_info.Rd[g][5][1] *= xynormi;
      sym_info.Rd[g][0][1] *= xynormi;
      sym_info.Rd[g][2][1] *= xynormi;
      
      if (f_exist)
        for (i=0; i < 10; i++)
          for (int j=0; j < 10; j++)
            sym_info.Rf[g][i][j] = rf(i,j);
      
      if (g_exist)
        for (i=0; i < 15; i++)
          for (int j=0; j < 15; j++)
            sym_info.Rg[g][i][j] = rg(i,j);
    }

    // now let's finish up the character table
    for (int ir=0; ir < nirr; ir++) {
      if (errcod=allocbn_character(&sym_info.ct.gamma[ir],"g",ct.order())) {
        fprintf(stderr,"sym_struct_from_gbs: ");
        fprintf(stderr,"could not alloc gamma(%d)\n",ir);
        return -1;
      }

      sym_info.ct.gamma[ir].degen = ct.gamma(ir).degeneracy();
      sym_info.ct.gamma[ir].nrot = ct.gamma(ir).nrot();
      sym_info.ct.gamma[ir].ntrans = ct.gamma(ir).ntrans();
      sym_info.ct.gamma[ir].label = strdup(ct.gamma(ir).symbol());

      for (g=0; g < ct.order(); g++)
        sym_info.ct.gamma[ir].rep[g] = ct.gamma(ir).character(g);
    }
  }
      
  return 0;
}
