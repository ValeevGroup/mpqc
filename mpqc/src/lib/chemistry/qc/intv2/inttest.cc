
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <math/array/math_lib.h>
}

#include <util/keyval/keyval.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/basis/basis.h>

void test_int_shell_1e(const RefKeyVal&, centers_t&,
            void (*int_shell_1e)(centers_t*,centers_t*,double*,int,int),
            int permute);
void test_3_center(const RefKeyVal&, centers_t&);
void test_4_center(const RefKeyVal&, centers_t&);
void test_4der_center(const RefKeyVal&, centers_t&);

#define maxint 9

main()
{
  int ii, i,j,k,l,m,n;

  char *infile = new char[strlen(SRCDIR)+strlen("/inttest.in")+1];
  sprintf(infile,SRCDIR "/inttest.in");

  RefKeyVal pkv(new ParsedKeyVal(infile));
  RefKeyVal keyval(new PrefixKeyVal(":centers :basis",pkv));
  RefKeyVal tkeyval(new PrefixKeyVal(":test", pkv));

  RefGaussianBasisSet basis = keyval->describedclassvalue("basis");
  RefMolecule mol = basis->molecule();
  centers_t *pcenters = basis->convert_to_centers_t();
  centers_t &centers = *pcenters;

  int_normalize_centers(&centers);
  int_initialize_offsets2(&centers,&centers,&centers,&centers);
  
  if (tkeyval->booleanvalue("print_centers")) print_centers(stdout,&centers);
  if (tkeyval->booleanvalue("overlap")) {
      printf("testing overlap:\n");
      int permute = keyval->booleanvalue("permute");
      test_int_shell_1e(tkeyval, centers, int_shell_overlap, permute);
    }
  if (tkeyval->booleanvalue("kinetic")) {
      printf("testing kinetic:\n");
      int permute = keyval->booleanvalue("permute");
      test_int_shell_1e(tkeyval, centers, int_shell_kinetic, permute);
    }
  if (tkeyval->booleanvalue("hcore")) {
      printf("testing hcore:\n");
      int permute = keyval->booleanvalue("permute");
      test_int_shell_1e(tkeyval, centers, int_shell_hcore, permute);
    }
  if (tkeyval->booleanvalue("nuclear")) {
      printf("testing nuclear:\n");
      int permute = keyval->booleanvalue("permute");
      test_int_shell_1e(tkeyval, centers, int_shell_nuclear, permute);
    }
  if (tkeyval->booleanvalue("3")) test_3_center(tkeyval, centers);
  if (tkeyval->booleanvalue("4")) test_4_center(tkeyval, centers);
  if (tkeyval->booleanvalue("4der")) test_4der_center(tkeyval, centers);

  int_done_offsets2(&centers,&centers,&centers,&centers);

  pkv = keyval = 0;

  return 0;
}

void
do_shell_test_1e(centers_t& centers,
                 void (*int_shell_1e)(centers_t*,centers_t*,double*,int,int),
                 int permute, int i, int j, int na, int nb,
                 double *buf, double *pbuf)
{
  int ii = 0;
  pbuf[na*nb] = 123.456;
  (*int_shell_1e)(&centers, &centers, buf, i, j);
  (*int_shell_1e)(&centers, &centers, pbuf, j, i);
  if (pbuf[na*nb] != 123.456) {
      printf("------- buffer overwritten for shells %d %d\n", j, i);
    }
  for (int a=0; a<na; a++) {
      for (int b=0; b<nb; b++) {
          if (fabs(buf[ii] - pbuf[a + na*b]) > 1.0e-13) {
              printf("----- 1e perm failed:"
                     "<%d %d|%d %d>:"
                     " %18.14f != %18.14f "
                     "<%d %d|%d %d>\n",
                     i, a, j, b,
                     buf[ii],
                     pbuf[a + na*b],
                     j, b, i, a);
            }
          if (fabs(buf[ii]) > 1.0e-15) {
              printf(" <(%d %d)|(%d %d)> = %15.11f\n",
                     i,a,j,b, buf[ii]);
            }
          ii++;
        }
    }
}

void
test_int_shell_1e(const RefKeyVal& keyval, centers_t& centers,
                  void (*int_shell_1e)(centers_t*,centers_t*,double*,int,int),
                  int permute)
{
  int flags = 0;
  double *buf = int_initialize_1e(flags, 1, &centers, &centers);
  int maxfunc = int_find_nfuncmax(&centers);
  int size = maxfunc * maxfunc;
  double *pbuf = new double[size+1];

  for (int i=0; i<centers.nshell; i++) {
      int na = INT_SH(&centers,i).nfunc;
      for (int j=0; j<centers.nshell; j++) {
          int nb = INT_SH(&centers,j).nfunc;
          do_shell_test_1e(centers, int_shell_1e, permute,
                           i, j, na, nb, buf, pbuf);

        }
    }
}

void
test_3_center(const RefKeyVal& keyval, centers_t& centers)
{
  int ii, i,j,k,l,m,n;

  int flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2|INT_NOPERM|INT_REDUND;
  double *buffer =
    int_initialize_erep(flags,0,&centers,&centers,&centers,&centers);

  for (i=0; i<centers.nshell; i++) {
      for (j=0; j<centers.nshell; j++) {
          int sh[2], sizes[2];
          sh[0] = i;
          sh[1] = j;
          int_erep2_v(flags,sh,sizes);
          ii = 0;
          for (k=0; k<sizes[0]; k++) {
              for (l=0; l<sizes[1]; l++) {
                  if (INT_NONZERO(buffer[ii]))
                      printf(" ((%d %d)|(%d %d)) = %15.11f\n",
                             sh[0],k,sh[1],l, buffer[ii]);
                  ii++;
                }
            }
        }
    }

  for (i=0; i<centers.nshell; i++) {
      for (j=0; j<centers.nshell; j++) {
          for (m=0; m<centers.nshell; m++) {
              int sh[3], sizes[3];
              sh[0] = i;
              sh[1] = j;
              sh[2] = m;
              int_erep3_v(flags,sh,sizes);
              ii = 0;
              for (k=0; k<sizes[0]; k++) {
                  for (l=0; l<sizes[1]; l++) {
                      for (n=0; n<sizes[2]; n++) {
                          if (INT_NONZERO(buffer[ii]))
                              printf(" ((%d %d)|(%d %d)(%d %d)) = %15.11f\n",
                                     sh[0],k,sh[1],l,sh[2],n, buffer[ii]);
                          ii++;
                        }
                    }
                }
            }
        }
    }

  int_done_erep();
}

void
init_shell_perm(int flags, double *integrals,
                double buff[maxint][maxint][maxint][maxint],
                int sh[4], int sizes[4])
{
  int i, j, k, l;
  int_erep_v(flags|INT_NOPERM, sh, sizes);
  for (i=0; i<sizes[0]; i++) {
      for (j=0; j<sizes[1]; j++) {
          for (k=0; k<sizes[2]; k++) {
              for (l=0; l<sizes[3]; l++) {
                  buff[i][j][k][l] = *integrals++;
                }
            }
        }
    }
}

void
check_shell_perm(int flags, double *integrals,
                 double buff[maxint][maxint][maxint][maxint],
                 int sh[4], int sizes[4], int p0, int p1, int p2, int p3)
{
  int ip[4], p[4];
  int psizes[4];
  int psh[4];
  int index = 0;
  int i[4];
  p[0] = p0;
  p[1] = p1;
  p[2] = p2;
  p[3] = p3;
  ip[p0] = 0;
  ip[p1] = 1;
  ip[p2] = 2;
  ip[p3] = 3;
  psh[0] = sh[p0];
  psh[1] = sh[p1];
  psh[2] = sh[p2];
  psh[3] = sh[p3];
  int_erep_v(flags|INT_NOPERM, psh, psizes);
  for (i[0]=0; i[0]<psizes[0]; i[0]++) {
      for (i[1]=0; i[1]<psizes[1]; i[1]++) {
          for (i[2]=0; i[2]<psizes[2]; i[2]++) {
              for (i[3]=0; i[3]<psizes[3]; i[3]++) {
                  if (fabs(buff[i[ip[0]]][i[ip[1]]][i[ip[2]]][i[ip[3]]]
                           - integrals[index]) > 1.0e-13) {
                      printf("perm %d %d %d %d failed:"
                             "((%d %d)(%d %d)|(%d %d)(%d %d)):"
                             " %18.14f != %18.14f "
                             "((%d %d)(%d %d)|(%d %d)(%d %d))\n",
                             p0, p1, p2, p3,
                             sh[0],i[0], sh[1],i[1], sh[2],i[2], sh[3],i[3],
                             buff[i[ip[0]]][i[ip[1]]][i[ip[2]]][i[ip[3]]],
                             integrals[index],
                             psh[0],i[p[0]], psh[1],i[p[1]],
                             psh[2],i[p[2]], psh[3],i[p[3]]);
                    }
                  index++;
                }
            }
        }
    }
}

void
do_shell_quartet_test(int flags, double* buffer,
                      int print, int bounds, int permute,
                      const RefKeyVal& keyval, centers_t& centers,
                      int i, int j, int k, int l)
{
  int sh[4], sizes[4];
  int ibuf;
  int ii, jj, kk, ll;
  sh[0] = i;
  sh[1] = j;
  sh[2] = k;
  sh[3] = l;
  double maxintegral, integralbound;
  if (bounds) {
      integralbound = int_bound_to_double(int_erep_4bound(i,j,k,l));
    }
  int_erep_v(flags,sh,sizes);
  ibuf = 0;
  maxintegral = 0.0;
  for (ii=0; ii<sizes[0]; ii++) {
      for (jj=0; jj<sizes[1]; jj++) {
          for (kk=0; kk<sizes[2]; kk++) {
              for (ll=0; ll<sizes[3]; ll++) {
                  if (fabs(buffer[ibuf]) > maxintegral) {
                      maxintegral = fabs(buffer[ibuf]);
                    }
                  if (bounds && fabs(buffer[ibuf]) > integralbound) {
                      printf("((%d %d)(%d %d)|(%d %d)(%d %d)) = %15.11f, "
                             "bound = %15.11f\n",
                             sh[0], ii, sh[1], jj, sh[2], kk, sh[3], ll,
                             buffer[ibuf], integralbound);
                    }
                  if (print && fabs(buffer[ibuf]) > 1.0e-9) {
                      printf(" ((%d %d)(%d %d)|(%d %d)(%d %d))"
                             " = %15.11f\n",
                             sh[0],ii,
                             sh[1],jj,
                             sh[2],kk,
                             sh[3],ll,
                             buffer[ibuf]);
                    }
                  ibuf++;
                }
            }
        }
    }

  if (permute) {
      double buff1[maxint][maxint][maxint][maxint];
      sh[0] = i;
      sh[1] = j;
      sh[2] = k;
      sh[3] = l;
      init_shell_perm(flags, buffer, buff1, sh, sizes);
      check_shell_perm(flags, buffer, buff1, sh, sizes, 0, 1, 2, 3);
      check_shell_perm(flags, buffer, buff1, sh, sizes, 1, 0, 2, 3);
      check_shell_perm(flags, buffer, buff1, sh, sizes, 0, 1, 3, 2);
      check_shell_perm(flags, buffer, buff1, sh, sizes, 1, 0, 3, 2);
      check_shell_perm(flags, buffer, buff1, sh, sizes, 2, 3, 0, 1);
      check_shell_perm(flags, buffer, buff1, sh, sizes, 2, 3, 1, 0);
      check_shell_perm(flags, buffer, buff1, sh, sizes, 3, 2, 0, 1);
      check_shell_perm(flags, buffer, buff1, sh, sizes, 3, 2, 1, 0);
    }

#if 0
  if (bounds) {
      printf("max(%d,%d,%d,%d)=%15.11f, bound=%15.11f, "
             "Q(%d,%d)=%3d, Q(%d,%d)=%3d\n",
             i, j, k, l, maxintegral, integralbound,
             i, j, int_erep_2bound(i,j), k, l, int_erep_2bound(k,l));
    }
#endif
}

void
do_4_center_test(int flags, double* buffer, int print, int bounds, int permute,
                 const RefKeyVal& keyval, centers_t& centers)
{
  int ii,jj,kk,ll, i,j,k,l, ibuf;

  for (i=0; i<centers.nshell; i++) {
      for (j=0; j<centers.nshell; j++) {
          for (k=0; k<centers.nshell; k++) {
              for (l=0; l<centers.nshell; l++) {
                  do_shell_quartet_test(flags, buffer, print, bounds, permute,
                                        keyval, centers, i, j, k, l);
                }
            }
        }
    }
}

void
test_4_center(const RefKeyVal& keyval, centers_t& centers)
{
  int i;

  int flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2|INT_NOPERM|INT_REDUND;
  double *buffer =
    int_initialize_erep(flags,0,&centers,&centers,&centers,&centers);

  int storage = keyval->intvalue("storage");
  int niter = keyval->intvalue("niter");
  int print = keyval->booleanvalue("print");
  int bounds = keyval->booleanvalue("bounds");
  int permute = keyval->booleanvalue("permute");

  printf("4 center test:\n");
  printf("  storage = %d\n", storage);
  printf("  niter   = %d\n", niter);
  printf("  print   = %d\n", print);
  printf("  bounds  = %d\n", bounds);
  printf("  permute = %d\n", permute);

  if (bounds) int_init_bounds();

  if (storage) int_storage(storage);

  for (i=0; i<niter; i++) {
      do_4_center_test(flags, buffer, print, bounds, permute, keyval, centers);
    }

  if (keyval->count("quartet") == 4) {
      do_shell_quartet_test(flags, buffer, print, bounds, permute,
                            keyval, centers,
                            keyval->intvalue("quartet", 0),
                            keyval->intvalue("quartet", 1),
                            keyval->intvalue("quartet", 2),
                            keyval->intvalue("quartet", 3));
    }

  if (storage) int_done_storage();

  if (bounds) int_done_bounds();

  int_done_erep();
}

void
do_shell_quartet_der_test(int flags,
                          double* buffer, int print, int bounds, int permute,
                          const RefKeyVal& keyval, centers_t& centers,
                          int i, int j, int k, int l)
{
  int ii,jj,kk,ll, ibuf, ider, xyz;
  der_centers_t dercenters;

  int sh[4], sizes[4];
  sh[0] = i;
  sh[1] = j;
  sh[2] = k;
  sh[3] = l;
  int_erep_all1der_v(flags,sh,sizes,&dercenters);
  ibuf = 0;
  for (ider=0; ider<dercenters.n; ider++) {
      for (xyz=0; xyz<3; xyz++) {
          for (ii=0; ii<sizes[0]; ii++) {
              for (jj=0; jj<sizes[1]; jj++) {
                  for (kk=0; kk<sizes[2]; kk++) {
                      for (ll=0; ll<sizes[3]; ll++) {
                          if (INT_NONZERO(buffer[ibuf])) {
                              printf(" ((%d %d)(%d %d)"
                                     "|(%d %d)(%d %d))(%d %d)"
                                     " = %15.11f\n",
                                     sh[0],ii,
                                     sh[1],jj,
                                     sh[2],kk,
                                     sh[3],ll,
                                     dercenters.num[ider], xyz,
                                     buffer[ibuf]
                                  );
                            }
                          ibuf++;
                        }
                    }
                }
            }
        }
    }
}

void
do_test_4der_center(int flags,
                    double* buffer, int print, int bounds, int permute,
                    const RefKeyVal& keyval, centers_t& centers)
{
  int i,j,k,l;
  for (i=0; i<centers.nshell; i++) {
      for (j=0; j<centers.nshell; j++) {
          for (k=0; k<centers.nshell; k++) {
              for (l=0; l<centers.nshell; l++) {
                  do_shell_quartet_der_test(flags, buffer,
                                            print, bounds, permute,
                                            keyval, centers,
                                            i, j, k, l);
                }
            }
        }
    }
}

void
test_4der_center(const RefKeyVal& keyval, centers_t& centers)
{
  int i;

  int flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2|INT_NOPERM|INT_REDUND;
  double *buffer =
    int_initialize_erep(flags,1,&centers,&centers,&centers,&centers);

  int niter = keyval->intvalue("niter");
  int print = keyval->booleanvalue("print");
  int bounds = keyval->booleanvalue("bounds");
  int permute = keyval->booleanvalue("permute");

  printf("4 center derivative test:\n");
  printf("  niter   = %d\n", niter);
  printf("  print   = %d\n", print);
  printf("  bounds  = %d\n", bounds);
  printf("  permute = %d\n", permute);

  if (permute || bounds) {
      printf("permute and bounds not implemented\n");
    }

  if (bounds) int_init_bounds();

  for (i=0; i<niter; i++) {
      do_test_4der_center(flags, buffer,
                          print, bounds, permute, keyval, centers);
    }

  if (keyval->count("quartet") == 4) {
      do_shell_quartet_der_test(flags, buffer, print, bounds, permute,
                                keyval, centers,
                                keyval->intvalue("quartet", 0),
                                keyval->intvalue("quartet", 1),
                                keyval->intvalue("quartet", 2),
                                keyval->intvalue("quartet", 3));
    }

  if (bounds) int_done_bounds();

  int_done_erep();
}
