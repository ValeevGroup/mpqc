
#include <stdio.h>
#include <util/ipv2/ip_libv2.h>
#include "test.h"
#include "testallc.h"
#include "testasgn.h"
#include "testzero.h"
#include "testinit.h"
#include "testiseq.h"
#include "testprnt.h"
#include "testip.h"

main(argc,argv)
int argc;
char** argv;
{
  char *labs[] = {"thing1","thing2","thing3","thing4"};
  int i;
  tests_t x1,x2;
  FILE *in,*out;

  if (argc != 2) {
    fprintf(stderr,"usage: %s filename\n",argv[0]);
    exit(1);
    }

  in = fopen(argv[1],"r");
  out = fopen("out.dat","w");

  ip_initialize(in,out);
  ip_cwk_add(":default");

  init_tests(&x1);

  ip_read_tests("tests",&x1,0);
  if(iseq_tests(&x1,&x2)) printf("\ntrue\n");
  else printf("\nfalse\n");

  init_tests(&x2);
  assign_tests(&x2,&x1);

  if(iseq_tests(&x1,&x2)) printf("\ntrue\n");
  else printf("\nfalse\n");

  print_tests(out,&x1);
  print_tests(out,&x2);

  zero_tests(&x2);
  print_tests(out,&x2);
  }
