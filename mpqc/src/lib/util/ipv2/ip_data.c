
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <tmpl.h>
#include "ip_types.h"
#include "ip_global.h"
#include "ip_error.h"

#include "ip_error.gbl"

#include "ip_data.gbl"
#include "ip_data.lcl"

#include "ip_cwk.gbl"
#include "ip_karray.gbl"

GLOBAL_VA_FUNCTION int
#ifndef __STDC__
ip_boolean(keyword,boolean,n)
char *keyword;
int *boolean;
int n;
#else
ip_boolean(char *keyword,int *boolean,int n,...)
#endif
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
    return ip_boolean_v(keyword,boolean,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) return IPE_MALLOC;
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = ip_boolean_v(keyword,boolean,n,v);
    free(v);
    return r;
    }
  }

GLOBAL_FUNCTION int
ip_boolean_v(keyword,boolean,n,v)
char *keyword;
int *boolean;
int n;
int *v;
{
  char *val;
  int errcod;
  char copy[10],*s;

  if ((errcod = ip_value_v(keyword,&val,n,v)) !=0) return errcod;

  strncpy(copy,val,10);
  copy[9] = '\0';

  /* Convert the string to uppercase. */
  for (s=copy; *s!='\0'; s++) {
    if (*s>='a' && *s <='z') *s = *s + 'A' - 'a';
    }
  
  if (!strcmp(copy,"YES")) *boolean = 1;
  else if (!strcmp(copy,"NO")) *boolean = 0;
  else if (!strcmp(copy,"1")) *boolean = 1;
  else if (!strcmp(copy,"0")) *boolean = 0;
  else if (!strcmp(copy,"TRUE")) *boolean = 1;
  else if (!strcmp(copy,"FALSE")) *boolean = 0;
  else if (!strcmp(copy,"ON")) *boolean = 1;
  else if (!strcmp(copy,"OFF")) *boolean = 0;
  else return IPE_TYPE;

  return IPE_OK;
  }

/* n should always be zero in this version of libip. */
GLOBAL_VA_FUNCTION int
#ifndef __STDC__
ip_exist(keyword,n)
char *keyword;
int n;
#else
ip_exist(char *keyword,int n,...)
#endif
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
    return ip_exist_v(keyword,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) {
      ip_warn("ip_exist: problem mallocing %d integers",n);
      return 0;
      }
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = ip_exist_v(keyword,n,v);
    free(v);
    return r;
    }
  }

/* n should always be zero in this version of libip. */
GLOBAL_FUNCTION int
ip_exist_v(keyword,n,v)
char *keyword;
int n;
int *v;
{
  if (n!=0) {
    fprintf(ip_out,"ip_exist_v: n must be zero n = %d, v = 0x%x\n", n, v);
    exit(1);
    }
  if (ip_cwk_descend_tree(keyword)) return 1;

  return 0;
  }

GLOBAL_VA_FUNCTION int
#ifndef __STDC__
ip_data(keyword,conv,value,n)
char *keyword;
char *conv;
VOID_PTR value;
int n;
#else
ip_data(char *keyword,char *conv,void *value,int n,...)
#endif
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
    return ip_data_v(keyword,conv,value,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) return IPE_MALLOC;
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = ip_data_v(keyword,conv,value,n,v);
    free(v);
    return r;
    }
  }

GLOBAL_FUNCTION int
ip_data_v(keyword,conv,value,n,v)
char *keyword;
char *conv;
VOID_PTR value;
int n;
int *v;
{
  char *val;
  int errcod;

  if ((errcod = ip_value_v(keyword,&val,n,v))!=0) return errcod;

  if (sscanf(val,conv,value) != 1) return IPE_TYPE;

  return IPE_OK;
  }

GLOBAL_VA_FUNCTION int
#ifndef __STDC__
ip_class(keyword,classname,n)
char *keyword;
char** classname;
int n;
#else
ip_class(char *keyword,char **classname,int n,...)
#endif
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
    return ip_class_v(keyword,classname,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) return IPE_MALLOC;
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = ip_class_v(keyword,classname,n,v);
    free(v);
    return r;
    }
  }

GLOBAL_FUNCTION int
ip_class_v(keyword,classname,n,v)
char *keyword;
char** classname;
int n;
int *v;
{
  ip_keyword_tree_t *kt;
  static char newkey[KEYWORD_LENGTH];
  int errcod;

  if ((errcod = ip_construct_key_v(keyword,newkey,n,v))!=IPE_OK) return errcod;

  kt = ip_cwk_descend_tree(newkey);
  if (!kt) {
    ip_lastkeyword(keyword);
    return IPE_KEY_NOT_FOUND;
    }

  ip_lastkeywordtree(kt);
  *classname = kt->classname;

  return IPE_OK;
  }

GLOBAL_VA_FUNCTION int
#ifndef __STDC__
ip_string(keyword,value,n)
char *keyword;
char **value;
int n;
#else
ip_string(char *keyword,char **value,int n,...)
#endif
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
    return ip_string_v(keyword,value,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) return IPE_MALLOC;
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = ip_string_v(keyword,value,n,v);
    free(v);
    return r;
    }
  }

GLOBAL_FUNCTION int
ip_string_v(keyword,value,n,v)
char *keyword;
char **value;
int n;
int *v;
{
  char *val;
  int errcod;

  if ((errcod = ip_value_v(keyword,&val,n,v))!=0) return errcod;

  *value = (char *) malloc(sizeof(char)*(strlen(val)+1));
  if (! *value) return IPE_MALLOC;
  strcpy(*value,val);
  return IPE_OK;
  }

GLOBAL_VA_FUNCTION int
#ifndef __STDC__
ip_value(keyword,value,n)
char *keyword;
char **value;
int n;
#else
ip_value(char *keyword,char **value,int n,...)
#endif
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
    return ip_value_v(keyword,value,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) return IPE_MALLOC;
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = ip_value_v(keyword,value,n,v);
    free(v);
    return r;
    }
  }

GLOBAL_FUNCTION int
ip_value_v(keyword,value,n,v)
char *keyword;
char **value;
int n;
int *v;
{
  int errcod;
  ip_keyword_tree_t *kt;
  static char newkey[KEYWORD_LENGTH];

  if ((errcod = ip_construct_key_v(keyword,newkey,n,v))!=IPE_OK) return errcod;

  /* Get the kt corresponding to the keyword using the cwk. */
  kt = ip_cwk_descend_tree(newkey);

  if (!kt) {
    ip_lastkeyword(keyword);
    return IPE_KEY_NOT_FOUND;
    }

  ip_lastkeywordtree(kt);

  *value = kt->value;
  if (*value == NULL) return IPE_HAS_NO_VALUE;

  return IPE_OK;
  }


GLOBAL_FUNCTION int
ip_construct_key_v(keyword,newkey,n,v)
char *keyword;
char *newkey;
int n;
int *v;
{
  int i;
  char index[11];
  int count,errcod;

  /* Construct the new keyword. */
  strcpy(newkey,keyword);
  for (i=0; i<n; i++) {
    if ((errcod = ip_count_v(newkey,&count,0,NULL)) != IPE_OK) return errcod;
    if (v[i]<0 || v[i]>=count) return IPE_OUT_OF_BOUNDS;
    sprintf(index,":%d",v[i]);
    strcat(newkey,index);
    }

  return IPE_OK;
  }
