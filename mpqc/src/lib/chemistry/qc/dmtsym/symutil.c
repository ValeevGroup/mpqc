
/* $Log$
 * Revision 1.1  1993/12/29 12:52:58  etseidl
 * Initial revision
 *
 * Revision 1.2  1992/06/17  21:56:07  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.1.1.1  1992/03/17  16:27:56  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:27:54  seidl
 * Initial revision
 *
 * Revision 1.1  1992/01/27  12:52:35  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/21  11:33:57  seidl
 * use libintv2
 *
 * Revision 1.1  1991/12/20  15:45:23  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <util/ipv2/ip_libv2.h>


#include "symm.h"
#include "symm_mac.h"
#include "symerr.gbl"

#include "symutil.gbl"
#include "symutil.lcl"

GLOBAL_FUNCTION int
sym_make_center(_center,_coord,_atom)
center_t *_center;
double *_coord;
int _atom;
{
  int i;
  int errcod;
  char *at_lab,*bset;
  double charge;

  init_center(_center);

  errcod = ip_string("atoms",&at_lab,1,_atom);
  if (errcod!=IPE_OK) return errcod;
  for(i=0; i < strlen(at_lab); i++) at_lab[i] = (char) toupper(at_lab[i]);
  _center->atom = at_lab;

  errcod = ip_data("charges","%lf",&charge,1,_atom);
  if(errcod != IPE_OK) _center->charge = atom_to_an(_center->atom);

  _center->r = (double *) malloc(sizeof(double)*3);
  if (!_center->r) return IPE_MALLOC;
  _center->r[0] = _coord[0];
  _center->r[1] = _coord[1];
  _center->r[2] = _coord[2];

  errcod = ip_count("basis",&i,0);
  if(errcod == IPE_NOT_AN_ARRAY) {
    errcod = ip_string("basis",&bset,0);
    if (errcod!=IPE_OK) return errcod;
    }
  else {
    errcod = ip_string("basis",&bset,1,_atom);
    if (errcod!=IPE_OK) return errcod;
    }

  errcod = read_basis(sym_to_atom(_center->atom),bset,&(_center->basis));
  if (errcod!=IPE_OK) return errcod;
  return IPE_OK;

  }


/* These utility routines assist in reading in the basis functions. */
/* they are ripped out of libint */


#define COUNT 1
#define READ 2

static char *prb = "encountered a problem while reading the basis set:\n ";

LOCAL_FUNCTION int
read_basis(atom,basisname,_basis)
char *atom;
char *basisname;
basis_t *_basis;
{
  int errcod;

  /* Copy the basisname to the basis structure. */
  _basis->name = (char *)malloc(strlen(basisname)+1);
  strcpy(_basis->name,basisname);

  /* If the basis name is "nothing" then return with an empty basis set. */
  if (!strcmp(basisname,"nothing")) {
    _basis->n = 0;
    _basis->shell = NULL;
    return IPE_OK;
    }

  if ((errcod = read_basis_(COUNT,atom,basisname,_basis))!=IPE_OK)
    return errcod;

  _basis->shell = (shell_t *) malloc(sizeof(shell_t)*_basis->n);
  _basis->n = 0;
  return read_basis_(READ,atom,basisname,_basis);
  }

LOCAL_FUNCTION int
read_basis_(what,atom,basisname,basis)
int what;
char *atom;
char *basisname;
basis_t *basis;
{
  int i;
  int errcod;
  char key[KEYWORD_LENGTH];
  char *val;

  i = 0;

  while(1) {
    sprintf(key,":basis:%s:%s:%d",atom,basisname,i);
    errcod = ip_value_v(key,&val,0,NULL);
    if (errcod == IPE_KEY_NOT_FOUND) break;
    sprintf(key,":basis:%s:%s:%d:type",atom,basisname,i);
    errcod = ip_value_v(key,&val,0,NULL);
    if (errcod == IPE_OK) {
      ip_warn("%k has been illegally assigned a value");
      return IPE_VAL_NOT_EXPD;
      }
    else if (errcod == IPE_HAS_NO_VALUE) {
      /* We have found a shell read it in (or count). */
      if (what == READ) {
        sprintf(key,":basis:%s:%s:%d",atom,basisname,i);
        errcod = ip_read_shell_v(key,&basis->shell[basis->n],0,NULL);
        if (errcod != IPE_OK) {
          ip_warn("%scouldn't read a shell: looking for \"%k\"",prb);
          return errcod;
          }
        }
      basis->n++;
      }
    else if (errcod == IPE_KEY_NOT_FOUND) {
      /* Not an atom, see if we have a get. */
      sprintf(key,":basis:%s:%s:%d:get",atom,basisname,i);
      errcod = ip_value_v(key,&val,0,NULL);
      if (errcod != IPE_OK) {
        ip_warn("%scouldn't find a shell or a \"get\" while looking for \"%k\"",
                prb);
        return errcod;
        }
      errcod = read_basis_(what,atom,val,basis);
      if (errcod != IPE_OK) return errcod;
      }
    else {
      return errcod;
      }
    i++;
    }
  /* If i is still zero, then something was missing in the input. */
  if (!i) {
    ip_warn("%scould not find \"%k\"",prb);
    exit(1);
    }
  return IPE_OK;
  }
