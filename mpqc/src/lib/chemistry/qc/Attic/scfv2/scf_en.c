
/* $Log$
 * Revision 1.1  1993/12/29 12:53:04  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:08:31  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:08:29  seidl
 * Initial revision
 *
 * Revision 1.4  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.3  1992/01/13  19:11:28  seidl
 * move test for open-shell out of loop, add comments
 *
 * Revision 1.2  1992/01/09  11:38:02  seidl
 * add parallel code
 *
 * Revision 1.1  1991/12/20  16:23:05  seidl
 * Initial revision
 *
 * Revision 1.2  1991/06/22  08:58:55  seidl
 * op is no longer a pointer
 *
 * Revision 1.1  1991/06/19  13:01:36  seidl
 * Initial revision
 * */

static char rcsid[]="$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <math.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/sym/sym_lib.h>
#include "scf.h"

#include "scf_en.gbl"
#include "scf_en.lcl"

static double twocut=1.0;
static double plimit;

GLOBAL_FUNCTION int
scf_calculate_energy(_irreps,_scf_info,_mfp,_mast,
  _fock,_pmat,_dpmat,_gmato,_pmato,_hcore,iter,_outfile)
scf_irreps_t *_irreps;
scf_struct_t *_scf_info;
scf_mf_ptrs_t *_mfp;
int _mast;
sym_d_vector_t *_fock;
sym_d_vector_t *_pmat;
sym_d_vector_t *_dpmat;
sym_d_vector_t *_gmato;
sym_d_vector_t *_pmato;
sym_d_vector_t *_hcore;
int iter;
FILE *_outfile;
{
  int i,j,k,ij,nn;
  double edif,etot;
  double neelec = 0.0;
  double ir_energy, dtemp;
  double cinext;
  double delta;
  scf_irrep_t *s;

  if(!iter) plimit = pow(10.0,(double) -(_scf_info->convergence));

/* calculate the total electronic energy
 * for closed shell, 
 *    Eelec = sum(i,j) D(i,j)*{H(i,j)+F(i,j)}
 *
 * for open shell, 
 *    Eelec = sum(i,j) D(i,j)*{H(i,j)+F(i,j)} - 1/2 Do(i,j)*Go(i,j)
 */

  delta=0.0;
  for(k=0; k < _irreps->nirrep ; k++) {
    s = &_irreps->ir[k];
    if(nn=s->num_so) {
      ir_energy = 0.0;
      if(!_scf_info->iopen) {
        for(i=ij=0; i < nn ; i++)
          for(j = 0 ; j <= i ; j++,ij++)
            ir_energy += 
              0.5*_pmat->ir[k].d[ij]*(_hcore->ir[k].d[ij]+_fock->ir[k].d[ij]);
         }
       else {
        for(i=ij=0; i < nn ; i++)
          for(j = 0 ; j <= i ; j++,ij++)
            ir_energy +=
              0.5*_pmat->ir[k].d[ij]*(_hcore->ir[k].d[ij]+_fock->ir[k].d[ij])
                                - 0.5*_pmato->ir[k].d[ij]*_gmato->ir[k].d[ij];
        }
      neelec += ir_energy;
      for (i = 0; i < nn*(nn+1)/2 ; i++) {
        dtemp = _dpmat->ir[k].d[i];
        delta += dtemp*dtemp;
        }
      }
    }

  delta = sqrt(delta)/_scf_info->mxcoef2;
  etot = _scf_info->nuc_rep + neelec;
  edif = _scf_info->e_elec - neelec;

#if !defined(I860) && !defined(NCUBE)
  if (!iter) {
    fprintf(_outfile,"\n  iter       total energy       ");
    fprintf(_outfile," delta E         delta P          diiser\n");
    }

  fprintf(_outfile, "%5d %20.10f %15.6e %15.6e %15.6e\n", 
                   iter+1, etot, edif, delta, _scf_info->diis_er);
  fflush(_outfile);
#else
  if(mynode0()==0) {
    if (!iter) {
      fprintf(_outfile,"\n  iter       total energy       ");
      fprintf(_outfile," delta E         delta P          diiser\n");
      }

    fprintf(_outfile, "%5d %20.10f %15.6e %15.6e %15.6e\n", 
                   iter+1, etot, edif, delta, _scf_info->diis_er);
    fflush(_outfile);
    }
#endif

  if ( delta < plimit && delta != 0.0) _scf_info->converged=1;

  _scf_info->e_elec = neelec;

#if 0
 /* this is for TCSCF, it determines when new ci coefficients should
  * be calculated */

  cinext = pow(10.0,-twocut);
  if (delta < cinext && delta && !converged) {
     twocut += incr;
     return(1);
     }
  else return(0);
#endif
  return(0);
  }
