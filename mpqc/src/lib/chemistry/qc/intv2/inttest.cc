
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <math/array/math_lib.h>
}

#include <util/keyval/keyval.h>
#include <chemistry/qc/intv2/int_libv2.h>

void test_3_center(const RefKeyVal&, centers_t&);
void test_4_center(const RefKeyVal&, centers_t&);
void test_4der_center(const RefKeyVal&, centers_t&);

main()
{
  int ii, i,j,k,l,m,n;

  char *infile = new char[strlen(SRCDIR)+strlen("/inttest.in")+1];
  sprintf(infile,SRCDIR "/inttest.in");

  RefKeyVal pkv(new ParsedKeyVal(infile));
  RefKeyVal keyval(new PrefixKeyVal(":centers :basis",pkv));

  centers_t centers;

  int_read_centers(keyval,centers);

  int_normalize_centers(&centers);
  int_initialize_offsets2(&centers,&centers,&centers,&centers);
  
  if (keyval->booleanvalue(":test:print")) print_centers(stdout,&centers);

  if (keyval->booleanvalue(":test:3")) test_3_center(keyval, centers);
  if (keyval->booleanvalue(":test:4")) test_4_center(keyval, centers);
  if (keyval->booleanvalue(":test:4der")) test_4der_center(keyval, centers);

  int_done_offsets2(&centers,&centers,&centers,&centers);

  pkv = keyval = 0;

  return 0;
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
test_4_center(const RefKeyVal& keyval, centers_t& centers)
{
  int ii,jj,kk,ll, i,j,k,l, ibuf;

  int flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2|INT_NOPERM|INT_REDUND;
  double *buffer =
    int_initialize_erep(flags,0,&centers,&centers,&centers,&centers);

  for (i=0; i<centers.nshell; i++) {
      for (j=0; j<centers.nshell; j++) {
          for (k=0; k<centers.nshell; k++) {
              for (l=0; l<centers.nshell; l++) {
                  int sh[4], sizes[4];
                  sh[0] = i;
                  sh[1] = j;
                  sh[2] = k;
                  sh[3] = l;
                  int_erep_v(flags,sh,sizes);
                  ibuf = 0;
                  for (ii=0; ii<sizes[0]; ii++) {
                      for (jj=0; jj<sizes[1]; jj++) {
                          for (kk=0; kk<sizes[2]; kk++) {
                              for (ll=0; ll<sizes[3]; ll++) {
                                  if (INT_NONZERO(buffer[ibuf])) {
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
                }
            }
        }
    }

  int_done_erep();
}

void
test_4der_center(const RefKeyVal& keyval, centers_t& centers)
{
  int ii,jj,kk,ll, i,j,k,l, ibuf, ider, xyz;
  der_centers_t dercenters;

  int flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2|INT_NOPERM|INT_REDUND;
  double *buffer =
    int_initialize_erep(flags,1,&centers,&centers,&centers,&centers);

  for (i=0; i<centers.nshell; i++) {
      for (j=0; j<centers.nshell; j++) {
          for (k=0; k<centers.nshell; k++) {
              for (l=0; l<centers.nshell; l++) {
                  int sh[4], sizes[4];
                  sh[0] = i;
                  sh[1] = j;
                  sh[2] = k;
                  sh[3] = l;
                  int_erep_all1der_v(flags,sh,sizes,&dercenters);
                  ibuf = 0;
                  for (ider=0; ider<dercenters.n; ider++) {
                      for (ii=0; ii<sizes[0]; ii++) {
                          for (jj=0; jj<sizes[1]; jj++) {
                              for (kk=0; kk<sizes[2]; kk++) {
                                  for (ll=0; ll<sizes[3]; ll++) {
                                      for (xyz=0; xyz<3; xyz++) {
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
            }
        }
    }

  int_done_erep();
}
