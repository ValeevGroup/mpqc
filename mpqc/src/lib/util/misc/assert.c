
#include <stdio.h>
#include <tmpl.h>
#include "assert.gbl"
#include "assert.lcl"

GLOBAL_FUNCTION VOID
util_assert(file,line)
char *file;
int line;
{
#ifdef NCUBE
  char filename[50];
  char message[256];
  FILE *fp;
  int mynode = mynode0();
  if (mynode==0) {
    printf("Assertion failed on node 0: in %s at %d.\n",file,line);
    }
# if 0 /* On a big subcube this will open too many files. */
  sprintf(filename,"node.%d",mynode);
  sprintf(message,"Assertion failed on node %d: in %s at %d.\n",
          mynode,file,line);
  fp = fopen(filename,"w");
  fprintf(fp,"%s",message);
  fclose(fp);
# endif
#else
  fprintf(stderr,"Assertion failed: in %s at %d.\n",file,line);
#endif
  exit(23);
  }
