/*
 * opdens.c
 *
 * Copyright (C) 1996 Limit Point Systems, Inc.
 *
 * Author: Curtis Janssen <cljanss@ca.sandia.gov>
 * Maintainer: LPS
 *
 * This file is part of the SC Toolkit.
 *
 * The SC Toolkit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * The SC Toolkit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
 * the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * The U.S. Government is granted a limited license as per AL 91-7.
 */

/* $Log$
 * Revision 1.5  1996/10/25 22:47:54  etseidl
 * add copyleft notice
 *
 * Revision 1.4  1996/03/23 02:37:05  cljanss
 * Everything can now be configured with autoconf.
 *
 * Revision 1.3  1995/03/17 01:48:29  cljanss
 * Removed -I. and -I$(SRCDIR) from the default include path in
 * GlobalMakefile to avoid name conflicts with system include files.
 * Modified files under src.lib to include all files relative to src.lib.
 * Makefiles under src.bin need to add the -I. and -I$(SRCDIR) back onto
 * INCLUDE and CXXINCLUDE or make other arrangements.
 *
 * Revision 1.2  1994/08/26  22:48:00  etseidl
 * get rid of rcsids
 *
 * Revision 1.1.1.1  1993/12/29  12:53:00  etseidl
 * SC source tree 0.1
 *
 * Revision 1.2  1992/06/17  21:37:21  jannsen
 * cleaned up for saber-c and converted to ngl loops
 *
 * Revision 1.1.1.1  1992/03/17  14:56:52  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  14:56:51  seidl
 * Initial revision
 *
 * Revision 1.1  1992/02/13  00:46:55  seidl
 * Initial revision
 * */

#include <stdio.h>
#include <tmpl.h>
#include <math/dmt/libdmt.h>

#include <chemistry/qc/dmtqc/opdens.gbl>
#include <chemistry/qc/dmtqc/opdens.lcl>

/* This computes the open-shell density */

GLOBAL_FUNCTION void
dmt_open_density (C, ndoc, nsoc, PO)
dmt_matrix C;
int ndoc;
int nsoc;
dmt_matrix PO;
{
  double *col;
  int i;
  loop_t *loop;
  int iind,isize,jsize;

  dmt_fill (PO, 0.0);
  loop = dmt_ngl_create("%mr",C);
  while(dmt_ngl_next(loop)) {
    dmt_ngl_create_inner(loop,0);
    while(dmt_ngl_next_inner_m(loop,&iind,&isize,&i,&jsize,&col)) {
      if (i >= ndoc && i < ndoc+nsoc) {
        /*Sum contribution from |col| into local |PO| blocks.*/
        int localp = dmt_nlocal (PO);
        int il;
        int mub, nub, mu, nu, mustart, nustart, musize, nusize;
        double *block;
        for (il = 0; il<localp; il++) {
          dmt_get_block (PO, il, &mub, &nub, &block);
          dmt_describe_block (PO, mub, &mustart, &musize);
          dmt_describe_block (PO, nub, &nustart, &nusize);
          for (mu = 0; mu<musize; mu++) {
            for (nu = 0; nu<nusize; nu++) {
              block[mu*nusize+nu] += col[nustart+nu]*col[mustart+mu];
              }
            }
          }
        }
      }
    }
  dmt_ngl_kill(loop);

  }
