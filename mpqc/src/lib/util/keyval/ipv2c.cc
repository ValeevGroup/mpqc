
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include "ipv2.h"
#include "ipv2c.h"

void
ip_initialize(FILE*in,FILE*out)
{
  if (IPV2::have_global()) ip_done();
  IPV2::set_global(new IPV2);
  ip_read(in,out);
}

void
ip_done()
{
  if (IPV2::have_global()) {
      delete IPV2::global();
      IPV2::set_global(0);
    }
}

void
ip_set_uppercase(int u)
{
  IPV2::global()->set_uppercase(u);
}
void
ip_read(FILE*in,FILE*out)
{
  IPV2::global()->read(in,out);
}
void
ip_append_from_input(const char*prefix,FILE*fp)
{
  IPV2::global()->append_from_input(prefix,fp);
}
const char*
ip_error_message(int err)
{
  return IPV2::global()->error_message((IPV2::Status) err);
}
void
ip_error(const char*fmt,...)
{
  va_list args;
  va_start(args,fmt);
  IPV2::global()->error(fmt,args);
  va_end(args);
}
void
ip_warn(const char*fmt,...)
{
  va_list args;
  va_start(args,fmt);
  IPV2::global()->warn(fmt,args);
  va_end(args);
}
void
ip_cwk_root()
{
  IPV2::global()->cwk_root();
}
void
ip_cwk_clear()
{
  IPV2::global()->cwk_clear();
}
void
ip_cwk_add(const char*cwk)
{
  IPV2::global()->cwk_add(cwk);
}
void
ip_cwk_push()
{
  IPV2::global()->cwk_push();
}
void
ip_cwk_pop()
{
  IPV2::global()->cwk_pop();
}
int
ip_boolean(const char*key,int*val,int n,...)
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
      return ip_boolean_v(key,val,n,0);
    }
  else {
      v = (int *) malloc(sizeof(int)*n);
      if (!v) return IPE_MALLOC;
      va_start(args, n);
      for (i=0; i<n; i++) {
          v[i] = va_arg(args,int);
        }
      va_end(args);
      r = ip_boolean_v(key,val,n,v);
      free(v);
      return r;
    }
}
int
ip_boolean_v(const char*key,int*val,int n,int* v)
{
  return IPV2::global()->boolean_v(key,val,n,v);
}
int
ip_exist(const char*key,int n,...)
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
      return ip_exist_v(key,n,0);
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
      r = ip_exist_v(key,n,v);
      free(v);
      return r;
    }
}
int
ip_exist_v(const char*key,int n,int* v)
{
  return IPV2::global()->exist_v(key,n,v);
}
int
ip_construct_key_v(const char*key,char*s,int n,int*v)
{
  return (int) IPV2::global()->construct_key_v(key,s,n,v);
}
int
ip_count(const char*key,int*count,int n,...)
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
      return ip_count_v(key,count,n,0);
    }
  else {
      v = (int *) malloc(sizeof(int)*n);
      if (!v) return IPE_MALLOC;
      va_start(args, n);
      for (i=0; i<n; i++) {
          v[i] = va_arg(args,int);
        }
      va_end(args);
      r = ip_count_v(key,count,n,v);
      free(v);
      return r;
    }
}
int
ip_count_v(const char*key,int*count,int n,int*v)
{
  return (int) IPV2::global()->count_v(key,count,n,v);
}
int
ip_data(const char*key,const char*fmt,void*dat,int n, ...)
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
      return ip_data_v(key,fmt,dat,n,0);
    }
  else {
      v = (int *) malloc(sizeof(int)*n);
      if (!v) return IPE_MALLOC;
      va_start(args, n);
      for (i=0; i<n; i++) {
          v[i] = va_arg(args,int);
        }
      va_end(args);
      r = ip_data_v(key,fmt,dat,n,v);
      free(v);
      return r;
    }
}
int
ip_data_v(const char*key,const char*fmt,void*dat,int n,int*v)
{
  return (int) IPV2::global()->data_v(key,fmt,dat,n,v);
}
int
ip_class(const char*key,const char**clas,int n, ...)
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
      return ip_class_v(key,clas,n,0);
    }
  else {
      v = (int *) malloc(sizeof(int)*n);
      if (!v) return IPE_MALLOC;
      va_start(args, n);
      for (i=0; i<n; i++) {
          v[i] = va_arg(args,int);
        }
      va_end(args);
      r = ip_class_v(key,clas,n,v);
      free(v);
      return r;
    }
}
int
ip_class_v(const char*key,const char**clas,int n,int*v)
{
  return (int) IPV2::global()->classname_v(key,clas,n,v);
}
int
ip_string(const char*key,char**dat,int n,...)
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
      return ip_string_v(key,dat,n,0);
    }
  else {
      v = (int *) malloc(sizeof(int)*n);
      if (!v) return IPE_MALLOC;
      va_start(args, n);
      for (i=0; i<n; i++) {
          v[i] = va_arg(args,int);
        }
      va_end(args);
      r = ip_string_v(key,dat,n,v);
      free(v);
      return r;
    }
}
int
ip_string_v(const char*key,char**dat,int n,int*v)
{
  int ret = (int) IPV2::global()->string_v(key,dat,n,v);
  // convert new'ed memory into malloc'ed memory
  if ((ret == IPE_OK) && *dat) {
      char* tmp = (char*) malloc(strlen(*dat)+1);
      strcpy(tmp,*dat);
      delete[] *dat;
      *dat = tmp;
    }
  return ret;
}
int
ip_value(const char*key,const char**dat,int n, ...)
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
      return ip_value_v(key,dat,n,0);
    }
  else {
      v = (int *) malloc(sizeof(int)*n);
      if (!v) return IPE_MALLOC;
      va_start(args, n);
      for (i=0; i<n; i++) {
          v[i] = va_arg(args,int);
        }
      va_end(args);
      r = ip_value_v(key,dat,n,v);
      free(v);
      return r;
    }
}
int
ip_value_v(const char*key,const char**dat,int n,int*v)
{
  return (int) IPV2::global()->value_v(key,dat,n,v);
}
