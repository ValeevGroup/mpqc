
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <math/array/math_lib.h>
}

#include <util/keyval/keyval.h>
#include <chemistry/qc/intv2/int_libv2.h>

main()
{
  int ii, i,j,k,l,m,n;

  char *infile = new char[strlen(SRCDIR)+strlen("/inttest.in")+1];
  sprintf(infile,SRCDIR "/inttest.in");

  RefKeyVal pkv(new ParsedKeyVal(infile));
  RefKeyVal keyval(new PrefixKeyVal(":centers :basis",*pkv.pointer()));

  centers_t centers;

  int_read_centers(*keyval.pointer(),centers);

  pkv = keyval = 0;

  int_normalize_centers(&centers);
  int_initialize_offsets2(&centers,&centers,&centers,&centers);

  int flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2|INT_NOPERM|INT_REDUND;
  double *buffer =
    int_initialize_erep(flags,0,&centers,&centers,&centers,&centers);
  
  print_centers(stdout,&centers);

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
  int_done_offsets2(&centers,&centers,&centers,&centers);

  exit(0);
}
