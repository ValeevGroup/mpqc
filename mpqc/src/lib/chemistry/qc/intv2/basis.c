
/* These utility routines assist in reading in the basis functions. */

/* $Log$
 * Revision 1.1  1993/12/29 12:53:01  etseidl
 * Initial revision
 *
 * Revision 1.4  1993/04/28  00:30:13  jannsen
 * added int_read_basis global function
 *
 * Revision 1.3  1992/06/17  22:04:23  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/31  01:21:20  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/*
 * Revision 1.6  91/10/31  14:32:24  cljanss
 * The basis set name is kept in the basis structure.
 * 
 * Revision 1.5  91/09/28  19:26:42  cljanss
 * Switch to new file naming scheme
 * 
 * Revision 1.4  91/09/10  19:32:34  cljanss
 * If the basis name is "nothing", an empty basis set
 * will be returned.
 * 
 * Revision 1.3  1991/08/09  16:41:23  cljanss
 * fixed basis set reader, it returned incorrect return codes
 *
 * Revision 1.2  1991/07/16  17:55:27  cljanss
 * slight change to make a character initialization compile on the SGI
 *
 * Revision 1.1  1991/06/16  16:40:07  janssen
 * Initial revision
 * */
static char *rcsid = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <util/sgen/sgen.h>
#include <util/ipv2/ip_libv2.h>

#include "atoms.h"
#include "atomsip.h"
#include "atomsinit.h"

#include "basis.gbl"
#include "basis.lcl"

#include "atominfo.gbl"

#define COUNT 1
#define READ 2

static char *prb = "encountered a problem while reading the basis set:\n ";

GLOBAL_FUNCTION int
int_read_basis(atom,basisname,_basis)
char *atom;
char *basisname;
basis_t *_basis;
{
  read_basis(atom,basisname,_basis);
  }

LOCAL_FUNCTION int
read_basis(atom,basisname,_basis)
char *atom;
char *basisname;
basis_t *_basis;
{
  int errcod;

  /* Copy the basisname to the basis structure. */
  _basis->name = (char *)malloc(strlen(basisname)+1);
  sgen_chkmalloc(_basis->name);
  strcpy(_basis->name,basisname);

  /* If the basis name is "nothing" then return with an empty basis set. */
  if (!strcmp(basisname,"nothing")) {
    _basis->n = 0;
    _basis->shell = NULL;
    return IPE_OK;
    }

  if ((errcod = read_basis_(COUNT,atom,basisname,_basis))!=IPE_OK) return errcod;

  _basis->shell = (shell_t *) malloc(sizeof(shell_t)*_basis->n);
  sgen_chkmalloc(_basis->shell);
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

/* This is declared in ip_atoms.h. */
int
ip_read_center_v(keyword, _center, n, v)
char *keyword;
center_t *_center;
int n;
int *v;
{
  int i;
  int errcod;
  char newkey[KEYWORD_LENGTH],basekey[KEYWORD_LENGTH];
  char *val;

  init_center(_center);
  ip_construct_key_v(keyword,basekey,n,v);

  sprintf(newkey,"%s:%s",basekey,"atom");
  errcod = ip_value_v(newkey,&val,0,NULL);
  if (errcod!=IPE_OK && errcod!=IPE_KEY_NOT_FOUND) return errcod;
  if (!val) return IPE_HAS_NO_VALUE;
  _center->atom = (char *) malloc(strlen(val) + 1);
  sgen_chkmalloc(_center->atom);
  strcpy(_center->atom,val);
  if (errcod!=IPE_OK && errcod!=IPE_KEY_NOT_FOUND) return errcod;

  sprintf(newkey,"%s:%s",basekey,"charge");
  errcod = ip_read_double_v(newkey,&(_center->charge),0,NULL);
  if (errcod!=IPE_OK && errcod!=IPE_KEY_NOT_FOUND) return errcod;
  if (errcod == IPE_KEY_NOT_FOUND)
     _center->charge = atom_to_an(_center->atom);

  _center->r = (double *)malloc(sizeof(double)*3);
  sgen_chkmalloc(_center->r);
  if (!_center->r) return IPE_MALLOC;
  for (i=0; i<3; i++)  {
    sprintf(newkey,"%s:%s:%d",basekey,"r",i);
    errcod = ip_read_double_v(newkey,&(_center->r[i]),0,NULL);
    if (errcod!=IPE_OK && errcod!=IPE_KEY_NOT_FOUND) return errcod;
    }
  sprintf(newkey,"%s:%s",basekey,"basis");

  /* val is a pointer to the name of the basis set. */
  errcod = ip_value_v(newkey,&val,0,NULL);
  if (errcod!=IPE_OK && errcod!=IPE_KEY_NOT_FOUND) return errcod;
  if (!val) return IPE_HAS_NO_VALUE;
  errcod = read_basis(sym_to_atom(_center->atom),val,&(_center->basis));
  if (errcod!=IPE_OK) return errcod;
  return IPE_OK;
  }


/* This will print out a shell type.  The integer am is converted to
 * a character and if puream is nonnull, then the number of functions
 * in the shell is appended to the name. */
void print_shell_type(fp,_shell_type)
FILE *fp;
shell_type_t *_shell_type;
{
  char *amnames = "spdfghijklmnoqrtuvwxyz";

  fprintf(fp," am = %c",amnames[_shell_type->am]);
  if (_shell_type->puream) fprintf(fp,"%d",2*_shell_type->am + 1);
  sgen_print_suppress_indent();
  }

/* This fills in the am and puream members of the shell_type_t variable
 * by looking for the "am" keyword and treating it as a character
 * representation of a shell type. */
int ip_read_shell_type_v(keyword,_shell_type,n,v)
char *keyword;
shell_type_t *_shell_type;
int n;
int *v;
{
  char *val;
  char newkey[KEYWORD_LENGTH], basekey[KEYWORD_LENGTH];

  ip_construct_key_v(keyword,basekey,n,v);

  sprintf(newkey,"%s:am",basekey);
  ip_value_v(newkey,&val,0,NULL);

  if (!val) {
    ip_warn("%scould not find \"%k\"",prb);
    return IPE_KEY_NOT_FOUND;
    }
  if (val[0] == '\0') {
    ip_warn("%sgot an empty \"%k\"",prb);
    return IPE_TYPE;
    }

  if (val[0] == 's' || val[0] == 'S') {
    _shell_type->am = 0;
    }
  else if (val[0] == 'p' || val[0] == 'P') {
    _shell_type->am = 1;
    }
  else if (val[0] == 'd' || val[0] == 'D') {
    _shell_type->am = 2;
    }
  else if (val[0] == 'f' || val[0] == 'F') {
    _shell_type->am = 3;
    }
  else if (val[0] == 'g' || val[0] == 'G') {
    _shell_type->am = 4;
    }
  else if (val[0] == 'h' || val[0] == 'H') {
    _shell_type->am = 5;
    }
  else if (val[0] == 'i' || val[0] == 'I') {
    _shell_type->am = 6;
    }
  else if (val[0] == 'j' || val[0] == 'J') {
    _shell_type->am = 7;
    }
  else if (val[0] == 'k' || val[0] == 'K') {
    _shell_type->am = 8;
    }
  else if (val[0] == 'l' || val[0] == 'L') {
    _shell_type->am = 9;
    }
  else {
    ip_warn("%sillegal \"%k\": \"%s\"",prb,_shell_type->am);
    return IPE_TYPE;
    }

  /* Set the puream flag to 1 if pure angular momentum functions are 
   * desired. */
  if ((_shell_type->am < 2) || (val[1] == '\0')) {
    /* A simple shell was given, no pure am is used. */
    _shell_type->puream = 0;
    }
  else {
    int npure = atoi(&val[1]);

    /* If pure am is desired, d5, f7, g9, etc. must be specified,
     * otherwise return a IPE_TYPE error. */
    if (npure != 2*_shell_type->am + 1) {
      ip_warn("%sillegal \"%k\": \"%s\"",prb,_shell_type->am);
      return IPE_TYPE;
      }
    _shell_type->puream = 1;
    /* Cannot specify pure am for s and p shells. */
    if (_shell_type->am < 2) {
      ip_warn("%sillegal \"%k\": \"%s\"",prb,_shell_type->am);
      return IPE_TYPE;
      }
    }

  return IPE_OK;
  }
