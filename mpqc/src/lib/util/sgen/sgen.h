
#ifndef _SGEN_H
#define _SGEN_H

#ifdef __STDC__
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

void sgen_print_suppress_indent( ARGS0 );
void sgen_print_indent( ARGS1(FILE*) );

void sgen_bcast0_set_print( ARGS3(int,int,int) );
void sgen_bcast_print( ARGS4(int,int,int,int) );
void bcast0_sync( ARGS0 );
void sgen_reset_bcast0( ARGS0 );

void bread_boolean( ARGS4(int,int*,int*,int) );
void bread_char( ARGS4(int,char*,int*,int) );
void bread_double( ARGS4(int,double*,int*,int) );
void bread_float( ARGS4(int,float*,int*,int) );
void bread_int( ARGS4(int,int*,int*,int) );
void bread_long( ARGS4(int,long*,int*,int) );
void bread_pointer( ARGS4(int,void*,int*,int) );
int * bread_test_pointer( ARGS3(int,int*,int) );
void bread_string( ARGS4(int,char**,int*,int) );

void bwrite_boolean( ARGS4(int,int*,int*,int) );
void bwrite_char( ARGS4(int,char*,int*,int) );
void bwrite_double( ARGS4(int,double*,int*,int) );
void bwrite_float( ARGS4(int,float*,int*,int) );
void bwrite_int( ARGS4(int,int*,int*,int) );
void bwrite_long( ARGS4(int,long*,int*,int) );
void bwrite_pointer( ARGS4(int,void*,int*,int) );
void bwrite_string( ARGS4(int,char**,int*,int) );

void sgen_chkmalloc( ARGS1(void*) );

void fread_boolean( ARGS4(FILE*,int*,int*,int) );
void fread_char( ARGS4(FILE*,char*,int*,int) );
void fread_double( ARGS4(FILE*,double*,int*,int) );
void fread_float( ARGS4(FILE*,float*,int*,int) );
void fread_int( ARGS4(FILE*,int*,int*,int) );
void fread_long( ARGS4(FILE*,long*,int*,int) );
void fread_pointer( ARGS4(FILE*,void*,int*,int) );
int * fread_test_pointer( ARGS3(FILE*,int*,int) );
void fread_string( ARGS4(FILE*,char**,int*,int) );

void fwrite_boolean( ARGS4(FILE*,int*,int*,int) );
void fwrite_char( ARGS4(FILE*,char*,int*,int) );
void fwrite_double( ARGS4(FILE*,double*,int*,int) );
void fwrite_float( ARGS4(FILE*,float*,int*,int) );
void fwrite_int( ARGS4(FILE*,int*,int*,int) );
void fwrite_long( ARGS4(FILE*,long*,int*,int) );
void fwrite_pointer( ARGS4(FILE*,void*,int*,int) );
void fwrite_string( ARGS4(FILE*,char**,int*,int) );

void ip_read_boolean( ARGS4(char*,int*,int,int*) );
void ip_read_char( ARGS4(char*,char*,int,int*) );
void ip_read_double( ARGS4(char*,double*,int,int*) );
void ip_read_float( ARGS4(char*,float*,int,int*) );
void ip_read_int( ARGS4(char*,int*,int,int*) );
void ip_read_long( ARGS4(char*,long*,int,int*) );
void ip_read_string( ARGS4(char*,char**,int,int*) );

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


#undef ARGS0
#undef ARGS1
#undef ARGS2
#undef ARGS3
#undef ARGS4
#undef ARGS5
#undef ARGS6

#endif /* _SGEN_H */

