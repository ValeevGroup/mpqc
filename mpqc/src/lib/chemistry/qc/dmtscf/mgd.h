
#ifndef _dmtscf_mgd_h
#define _dmtscf_mgd_h

#define MIN0(a,b) ((a)<(b)) ? (a) : (b)
#define MAX0(a,b) ((a)>(b)) ? (a) : (b)

#ifndef IOFF
#define IOFF(a,b) ((a)>(b))?(a)*((a)+1)/2+(b):(b)*((b)+1)/2+(a)
#endif

struct mgd {
  int si;
  int sj;
  int sk;
  int sl;
  int isz;
  int jsz;
  int ksz;
  int lsz;
  double *glp;
  double *glpo;
  double *plp;
  double *plpo;
  double *gloc;
  double *gloco;
  double *ploc;
  double *ploco;
} ;

#endif
