
/* All of the official picl routines are declared here.  These
 * declarations apply to all of the picl packages. */

#ifndef _PICL_H
#define _PICL_H

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

/* Low level routines. */
int host0( ARGS0 );
void check0( ARGS1(int) );
double clock0( ARGS0 );
void close0( ARGS1(int) );
void load0( ARGS2(char*,int) );
void message0( ARGS1(char*) );
void open0( ARGS3(int*,int*,int*) );
int probe0( ARGS1(int) );
void recv0( ARGS3(void*,int,int) );
void recvinfo0( ARGS3(int*,int*,int*) );
void send0( ARGS4(void*,int,int,int) );
void sync0( ARGS0 );
void who0( ARGS3(int*,int*,int*) );

/* High level routines. */
void setarc0( ARGS4(int*,int*,int*,int*) );
void getarc0( ARGS4(int*,int*,int*,int*) );

void bcast0( ARGS4(void*,int,int,int ) );
void bcast1( ARGS4(void*,int,int,int) );

void gcomb0( ARGS6(void*,int,int,int,int,void(*comb)()) );
void gand0( ARGS5(void*,int,int,int,int) );
void gmax0( ARGS5(void*,int,int,int,int) );
void gmin0( ARGS5(void*,int,int,int,int) );
void gor0( ARGS5(void*,int,int,int,int) );
void gprod0( ARGS5(void*,int,int,int,int) );
void gsum0( ARGS5(void*,int,int,int,int) );
void gxor0( ARGS5(void*,int,int,int,int) );

void barrier0( ARGS0 );

int gray0( ARGS1(int) );
int ginv0( ARGS1(int) );

#undef ARGS0
#undef ARGS1
#undef ARGS2
#undef ARGS3
#undef ARGS4
#undef ARGS5
#undef ARGS6

#endif /* _PICL_H */
