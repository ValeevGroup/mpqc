
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "ipv2.h"

IPV2::Status
IPV2::boolean(const char *keyword,int *boolean,int n,...)
{
  va_list args;
  int i;
  int *v;
  Status r;

  if (n==0) {
    return boolean_v(keyword,boolean,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) return Malloc;
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = boolean_v(keyword,boolean,n,v);
    free(v);
    return r;
    }
  }

IPV2::Status
IPV2::boolean_v(const char* keyword,int* boolean,int n,int* v)
{
  const char *val;
  Status errcod;
  char copy[10],*s;

  if ((errcod = value_v(keyword,&val,n,v)) !=0) return errcod;

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
  else return Type;

  return OK;
  }

/* n should always be zero in this version of libip. */
int
IPV2::exist(const char *keyword,int n,...)
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
    return exist_v(keyword,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) {
      warn("IPV2::exist: problem mallocing %d integers",n);
      return 0;
      }
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = exist_v(keyword,n,v);
    free(v);
    return r;
    }
  }

/* n should always be zero in this version of libip. */
int
IPV2::exist_v(const char* keyword,int n,int* v)
{
  if (n!=0) {
    fprintf(ip_out,"IPV2::exist_v: n must be zero n = %d, v = 0x%x\n", n, v);
    exit(1);
    }
  if (ip_cwk_descend_tree(keyword)) return 1;

  return 0;
  }

IPV2::Status
IPV2::data(const char *keyword,const char *conv,void *value,int n,...)
{
  va_list args;
  int i;
  int *v;
  Status r;

  if (n==0) {
    return data_v(keyword,conv,value,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) return Malloc;
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = data_v(keyword,conv,value,n,v);
    free(v);
    return r;
    }
  }

IPV2::Status
IPV2::data_v(const char* keyword,const char* conv,void* value,int n,int*v)
{
  const char *val;
  Status errcod;

  if ((errcod = value_v(keyword,&val,n,v))!=0) return errcod;

  if (sscanf(val,conv,value) != 1) return Type;

  return OK;
  }

IPV2::Status
IPV2::classname(const char *keyword,const char **name,int n,...)
{
  va_list args;
  int i;
  int *v;
  Status r;

  if (n==0) {
    return classname_v(keyword,name,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) return Malloc;
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = classname_v(keyword,name,n,v);
    free(v);
    return r;
    }
  }

IPV2::Status
IPV2::classname_v(const char* keyword,const char** name,int n,int *v)
{
  ip_keyword_tree_t *kt;
  static char newkey[KEYWORD_LENGTH];
  Status errcod;

  if ((errcod = construct_key_v(keyword,newkey,n,v))!=OK) return errcod;

  kt = ip_cwk_descend_tree(newkey);
  if (!kt) {
    ip_lastkeyword(keyword);
    *name = 0;
    return KeyNotFound;
    }

  ip_lastkeywordtree(kt);
  *name = kt->classname;

  return OK;
  }

IPV2::Status
IPV2::truekeyword(const char *keyword,const char **name,int n,...)
{
  va_list args;
  int i;
  int *v;
  Status r;

  if (n==0) {
    return truekeyword_v(keyword,name,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) return Malloc;
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = truekeyword_v(keyword,name,n,v);
    free(v);
    return r;
    }
  }

static char*
get_name(ip_keyword_tree_t*kt,char*currentname)
{
  if (!kt) return currentname;
  char*newname;
  if (currentname) {
      newname = (char*)malloc(strlen(currentname)+strlen(kt->keyword)+2);
      sprintf(newname,"%s:%s",kt->keyword,currentname);
      free(currentname);
    }
  else newname = strcpy((char*)malloc(strlen(kt->keyword)+1),kt->keyword);
  return get_name(kt->up,newname);
}

char*
IPV2::get_truename(ip_keyword_tree_t*kt)
{
  return get_name(kt,0);
}

IPV2::Status
IPV2::truekeyword_v(const char* keyword,const char** name,int n,int *v)
{
  ip_keyword_tree_t *kt;
  static char newkey[KEYWORD_LENGTH];
  Status errcod;

  if ((errcod = construct_key_v(keyword,newkey,n,v))!=OK) return errcod;

  kt = ip_cwk_descend_tree(newkey);
  if (!kt) {
    ip_lastkeyword(keyword);
    *name = 0;
    return KeyNotFound;
    }

  ip_lastkeywordtree(kt);

  if (!kt->truename) {
      kt->truename = get_truename(kt);
    }
  *name = kt->truename;

  return OK;
  }

IPV2::Status
IPV2::string(const char *keyword,char **value,int n,...)
{
  va_list args;
  int i;
  int *v;
  Status r;

  if (n==0) {
    return string_v(keyword,value,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) return Malloc;
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = string_v(keyword,value,n,v);
    free(v);
    return r;
    }
  }

// if an error is encountered, *value is not modified
IPV2::Status
IPV2::string_v(const char* keyword,char** value,int n,int *v)
{
  const char *val;
  Status errcod;

  if ((errcod = value_v(keyword,&val,n,v))!=0) return errcod;

  char *tmp = new char[strlen(val)+1];
  if (! tmp) return Malloc;
  strcpy(tmp,val);
  *value = tmp;
  return OK;
  }

IPV2::Status
IPV2::value(const char *keyword,const char **value,int n,...)
{
  va_list args;
  int i;
  int *v;
  Status r;

  if (n==0) {
    return value_v(keyword,value,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) return Malloc;
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = value_v(keyword,value,n,v);
    free(v);
    return r;
    }
  }

IPV2::Status
IPV2::value_v(const char* keyword,const char**value,int n,int*v)
{
  Status errcod;
  ip_keyword_tree_t *kt;
  static char newkey[KEYWORD_LENGTH];

  if ((errcod = construct_key_v(keyword,newkey,n,v))!=OK) return errcod;

  /* Get the kt corresponding to the keyword using the cwk. */
  kt = ip_cwk_descend_tree(newkey);

  if (!kt) {
    ip_lastkeyword(keyword);
    return KeyNotFound;
    }

  ip_lastkeywordtree(kt);

  *value = kt->value;
  if (*value == NULL) return HasNoValue;

  return OK;
  }


IPV2::Status
IPV2::construct_key_v(const char* keyword,char *newkey,int n,int*v)
{
  int i;
  char index[11];
  int count;
  Status errcod;

  /* Construct the new keyword. */
  strcpy(newkey,keyword);
  for (i=0; i<n; i++) {
    if ((errcod = count_v(newkey,&count,0,NULL)) != OK) return errcod;
    if (v[i]<0 || v[i]>=count) return OutOfBounds;
    sprintf(index,":%d",v[i]);
    strcat(newkey,index);
    }

  return OK;
  }
