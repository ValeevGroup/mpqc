
#ifndef _SGEN_H
#define _SGEN_H

#if defined(__cplusplus) || defined(__STDC__)
#define ARGS0 void
#define ARGS1(a) a
#define ARGS2(a,b) a,b
#define ARGS3(a,b,c) a,b,c
#define ARGS4(a,b,c,d) a,b,c,d
#define ARGS5(a,b,c,d,e) a,b,c,d,e
#define ARGS6(a,b,c,d,e,f) a,b,c,d,e,f
#else
#define ARGS0
#define ARGS1(a)
#define ARGS2(a,b)
#define ARGS3(a,b,c)
#define ARGS4(a,b,c,d)
#define ARGS5(a,b,c,d,e)
#define ARGS6(a,b,c,d,e,f)
#endif

#ifdef __cplusplus
extern "C" {
#endif

void sgen_print_suppress_indent( ARGS0 );
void sgen_print_indent( ARGS1(FILE*) );

void sgen_bcast0_set_print( ARGS3(int,int,int) );
void sgen_bcast_print( ARGS4(int,int,int,int) );
void bcast0_sync( ARGS0 );
void sgen_reset_bcast0( ARGS0 );

void sgen_chkmalloc( ARGS1(void*) );

void print_boolean( ARGS2(FILE*,int*) );
void print_char( ARGS2(FILE*,char*) );
void print_double( ARGS2(FILE*,double*) );
void print_float( ARGS2(FILE*,float*) );
void print_int( ARGS2(FILE*,int*) );
void print_long( ARGS2(FILE*,long*) );
void print_string( ARGS2(FILE*,char**) );

int iseq_boolean( ARGS2(int,int) );
int iseq_char( ARGS2(int,int) );
int iseq_double( ARGS2(double,double) );
int iseq_float( ARGS2(double,double) );
int iseq_int( ARGS2(int,int) );
int iseq_long( ARGS2(long,long) );
int iseq_string( ARGS2(char*,char*) );

void rbcast0_boolean( ARGS4(int*,int,int,int) );
void rbcast0_char( ARGS4(char*,int,int,int) );
void rbcast0_double( ARGS4(double*,int,int,int) );
void rbcast0_float( ARGS4(float*,int,int,int) );
void rbcast0_int( ARGS4(int*,int,int,int) );
void rbcast0_long( ARGS4(long*,int,int,int) );
void rbcast0_pointer( ARGS4(void*,int,int,int) );
int * rbcast0_test_pointer( ARGS3(int,int,int) );
void rbcast0_string( ARGS4(char**,int,int,int) );

void sbcast0_boolean( ARGS4(int*,int,int,int) );
void sbcast0_char( ARGS4(char*,int,int,int) );
void sbcast0_double( ARGS4(double*,int,int,int) );
void sbcast0_float( ARGS4(float*,int,int,int) );
void sbcast0_int( ARGS4(int*,int,int,int) );
void sbcast0_long( ARGS4(long*,int,int,int) );
void sbcast0_pointer( ARGS4(void*,int,int,int) );
void sbcast0_string( ARGS4(char**,int,int,int) );

void recv0_boolean( ARGS4(int*,int,int,int) );
void recv0_char( ARGS4(char*,int,int,int) );
void recv0_double( ARGS4(double*,int,int,int) );
void recv0_float( ARGS4(float*,int,int,int) );
void recv0_int( ARGS4(int*,int,int,int) );
void recv0_long( ARGS4(long*,int,int,int) );
void recv0_pointer( ARGS4(void*,int,int,int) );
int * recv0_test_pointer( ARGS3(int,int,int) );
void recv0_string( ARGS4(char**,int,int,int) );

void send0_boolean( ARGS4(int*,int,int,int) );
void send0_char( ARGS4(char*,int,int,int) );
void send0_double( ARGS4(double*,int,int,int) );
void send0_float( ARGS4(float*,int,int,int) );
void send0_int( ARGS4(int*,int,int,int) );
void send0_long( ARGS4(long*,int,int,int) );
void send0_pointer( ARGS4(void*,int,int,int) );
void send0_string( ARGS4(char**,int,int,int) );

void sgen_sndrcv0_set_print( ARGS3(int,int,int) );
void sgen_sndrcv_print( ARGS4(int,int,int,int) );
void sndrcv0_sync( ARGS0 );
void sgen_reset_sndrcv0( ARGS0 );

#ifdef __cplusplus
}
#endif

#undef ARGS0
#undef ARGS1
#undef ARGS2
#undef ARGS3
#undef ARGS4
#undef ARGS5
#undef ARGS6

#endif /* _SGEN_H */

