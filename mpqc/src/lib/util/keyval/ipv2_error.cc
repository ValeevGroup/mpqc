
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <util/keyval/ipv2.h>

/* Returns some text for an errcod. */
const char*
IPV2::error_message(IPV2::Status errcod)
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
  
  if (errcod == OK) return ipe_ok;
  if (errcod == KeyNotFound) return ipe_key_not_found;
  if (errcod == OutOfBounds) return ipe_out_of_bounds;
  if (errcod == Malloc) return ipe_malloc;
  if (errcod == NotAnArray) return ipe_not_an_array;
  if (errcod == NotAScalar) return ipe_not_a_scalar;
  if (errcod == Type) return ipe_type;
  if (errcod == HasNoValue) return ipe_has_no_value;
  if (errcod == ValNotExpd) return ipe_val_not_expd;
  return huh;
}

void
IPV2::error(const char *msg,...)
{
  va_list args;
  va_start(args,msg);
  fprintf(stderr,"IPV2::error: ");
  vfprintf(stderr,msg,args);
  fprintf(stderr,"\n");
  va_end(args);
  showpos();
  exit(1);
}

void
IPV2::warn(const char *msg,...)
{
  va_list args;
  char *newmsg;
  const char *poskey;
  
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
      va_start(args,msg);
      fprintf(stderr,"IPV2::warn: ");
      vfprintf(stderr,newmsg,args);
      fprintf(stderr,"\n");
      va_end(args);
      if (poskey) free(newmsg);
    }
  else {
      va_start(args,msg);
      fprintf(stderr,"IPV2::warn: ");
      vfprintf(stderr,msg,args);
      fprintf(stderr,"\n");
      va_end(args);
    }
}

void
IPV2::ip_lastkeyword(const char*keyword)
{
  strncpy(lastkeyword,keyword, KEYWORD_LENGTH-1);
}

void
IPV2::ip_lastkeywordtree(ip_keyword_tree_t*kt)
{
  lastkeyword[0] = '\0';
  ip_lastkeyword_(kt);
}

void
IPV2::ip_lastkeyword_(ip_keyword_tree_t*kt)
{
  if (kt->up) ip_lastkeyword_(kt->up);
  if (strlen(lastkeyword) + strlen(kt->keyword) + 2 > KEYWORD_LENGTH) {
      cerr << "IPV2: keyword too big" << endl;
      abort();
    }
  strcat(lastkeyword,":");
  strcat(lastkeyword,kt->keyword);
}
