
#ifndef _PICLEXT_H
#define _PICLEXT_H

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

int cubedim( ARGS0 );
int mynode ( ARGS0 );
int numnodes( ARGS0 );
int infocount( ARGS0 );
int mynode0( ARGS0 );
int cubedim0( ARGS0 );
int numnodes0( ARGS0 );

int loop_out_neigh( ARGS0 );
int my_loop_index( ARGS0 );
int loop_in_neigh( ARGS0 );
int loop_neigh( ARGS1(int) );
void gop0( ARGS4(double*,int,int,int) );
void gop1( ARGS5(double*,int,double*,int,int) );
void gop0_sc( ARGS4(signed char*,int,int,int) );

void picl_prober( ARGS0 );

#if 0
void who1( ARGS3(int*,int*,int*) );
#endif

void gdcomb(ARGS5(int,double*,double*,int,int));

#undef ARGS0
#undef ARGS1
#undef ARGS2
#undef ARGS3
#undef ARGS4
#undef ARGS5
#undef ARGS6

#endif /* _PICLEXT_H */
