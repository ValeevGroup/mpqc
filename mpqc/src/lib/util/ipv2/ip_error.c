
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <tmpl.h>
#include "ip_types.h"
#include "ip_global.h"

/* Cannot include ip_error.global due to xlc's handling of varargs. */
#if !defined(AIXV3) && !defined(SGI) && !defined(NCUBE)
#include "ip_error.gbl"
#endif
#include "ip_error.lcl"
#include "ip_error.h"

#include "scan.gbl"

char lastkeyword[KEYWORD_LENGTH] = { '\0' };

/* Returns some text for an errcod. */
GLOBAL_FUNCTION char *
ip_error_message(errcod)
int errcod;
{
  static char *ipe_ok = "No problem has been detected.";
  static char *ipe_key_not_found = "No match was found for the given keyword.";
  static char *ipe_out_of_bounds = "An array index is out of bounds.";
  static char *ipe_malloc = "Memory allocation failed.";
  static char *ipe_not_an_array = "An index was given for a scalar quantity.";
  static char *ipe_not_a_scalar = "Expected a scalar, but found an array.";
  static char *ipe_type = "The datum is not of the appropiate type.";
  static char *ipe_has_no_value = "The keyword has no value.";
  static char *ipe_val_not_expd = "A value was not expected for the keyword.";
  static char *huh = "The nature of the problem is unknown.";

  if (errcod == IPE_OK) return ipe_ok;
  if (errcod == IPE_KEY_NOT_FOUND) return ipe_key_not_found;
  if (errcod == IPE_OUT_OF_BOUNDS) return ipe_out_of_bounds;
  if (errcod == IPE_MALLOC) return ipe_malloc;
  if (errcod == IPE_NOT_AN_ARRAY) return ipe_not_an_array;
  if (errcod == IPE_NOT_A_SCALAR) return ipe_not_a_scalar;
  if (errcod == IPE_TYPE) return ipe_type;
  if (errcod == IPE_HAS_NO_VALUE) return ipe_has_no_value;
  if (errcod == IPE_VAL_NOT_EXPD) return ipe_val_not_expd;
  return huh;
  }

GLOBAL_VA_FUNCTION VOID
#ifndef __STDC__
ip_error(msg)
char *msg;
#else
ip_error(char *msg,...)
#endif
{
  va_list args;
  va_start(args,msg);
  fprintf(ip_out,"IP_ERROR: ");
  vfprintf(ip_out,msg,args);
  fprintf(ip_out,"\n");
  va_end(args);
  showpos();
  exit(1);
  }

GLOBAL_VA_FUNCTION VOID
#ifndef __STDC__
ip_warn(msg)
char *msg;
#else
ip_warn(char *msg,...)
#endif
{
  va_list args;
  char *newmsg,*poskey;

  /* If msg has a %k in it, then substitute in the last keyword. */
#if defined(NCUBE)||defined(DEC)||defined(I860)
  /* The NCUBE and DEC are missing the strstr function. */
  for (poskey=msg; *poskey!='\0'; poskey++) {
    if (poskey[0] == '%' && poskey[1] == 'k') break;
    }
  if (poskey[0] != '%') poskey = NULL;
#else
  poskey = strstr(msg,"%k");
#endif
  if (poskey) {
    newmsg = (char *) malloc(strlen(msg)-1 + strlen(lastkeyword));
    strcpy(newmsg,msg);
    newmsg[poskey-msg] = '\0';
    strcat(newmsg,lastkeyword);
    strcat(newmsg,&poskey[2]);
    }
  else {
    newmsg = msg;
    }

  va_start(args,msg);
  fprintf(ip_out,"IP_WARN: ");
  vfprintf(ip_out,newmsg,args);
  fprintf(ip_out,"\n");
  va_end(args);

  if (poskey) free(newmsg);
  }

GLOBAL_FUNCTION VOID
ip_lastkeyword(keyword)
char *keyword;
{
  strcpy(lastkeyword,keyword);
  }

GLOBAL_FUNCTION VOID
ip_lastkeywordtree(kt)
ip_keyword_tree_t *kt;
{
  lastkeyword[0] = '\0';
  ip_lastkeyword_(kt);
  }

LOCAL_FUNCTION VOID
ip_lastkeyword_(kt)
ip_keyword_tree_t *kt;
{
  if (kt->up) ip_lastkeyword_(kt->up);
  strcat(lastkeyword,":");
  strcat(lastkeyword,kt->keyword);
  }
