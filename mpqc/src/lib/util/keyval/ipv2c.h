
/*  These give the definitions for global c functions that */
/*  read from IPV2::default().  These are intended to completely */
/*  replace the old util/ipv2 routines, where they are still used. */
/*  In places that the IPV2 routines use new, these routines copy */
/*  the data to malloc'ed data and do a copy.  ip_done() is a no-op, */
/*  the use should now delete the default IPV2 and do a set_default(0). */

#ifndef _util_keyval_ipv2_ipv2c_h
#define _util_keyval_ipv2_ipv2c_h

#include <stdio.h>

#define IPE_OK            0  /* No problem. */
#define IPE_KEY_NOT_FOUND 1  /* The keyword was not found. */
#define IPE_OUT_OF_BOUNDS 2  /* An array subscript was out of bounds. */
#define IPE_MALLOC        3  /* Memory allocation failed. */
#define IPE_NOT_AN_ARRAY  4  /* Gave index for data which isn't an array */
#define IPE_NOT_A_SCALAR  5  /* Didn't give index for data which is an array */
#define IPE_TYPE          6  /* The datum is not of the appropiate type. */
#define IPE_HAS_NO_VALUE  7  /* The keyword has no value. */
#define IPE_VAL_NOT_EXPD  8  /* A value was not expected for the keyword. */

#define KEYWORD_LENGTH 256

#ifdef __cplusplus
extern "C" {
#endif
  void ip_set_uppercase(int);
  void ip_initialize(FILE* in,FILE* out);
  void ip_done();
  void ip_read(FILE* in,FILE* out);
  void ip_append_from_input(const char*,FILE*);
  void ip_done();
  const char* ip_error_message(int);
  void ip_error(const char*,...);
  void ip_warn(const char*,...);
  void ip_cwk_root();
  void ip_cwk_clear();
  void ip_cwk_add(const char*);
  void ip_cwk_push();
  void ip_cwk_pop();
  int ip_count(const char*,int*,int,...);
  int ip_count_v(const char*,int*,int,int*);
  int ip_construct_key_v(const char*,char*,int,int*);
  int ip_boolean(const char*,int*,int,...);
  int ip_boolean_v(const char*,int*,int,int*);
  int ip_exist(const char*,int,...);
  int ip_exist_v(const char*,int,int*);
  int ip_data(const char*,const char*,void*,int,...);
  int ip_data_v(const char*,const char*,void*,int,int*);
    /* the character string produced by class must not be free'ed */
  int ip_class(const char*,const char**,int,...);
  int ip_class_v(const char*,const char**,int,int*);
  int ip_string(const char*,char**,int,...);
  int ip_string_v(const char*,char**,int,int*);
    /* the character string produced by value must not be free'ed */
  int ip_value(const char*,const char**,int,...);
  int ip_value_v(const char*,const char**,int,int*);

#ifdef __cplusplus
};
#endif

#endif
