#include <stdio.h>
#include <stdlib.h>
#include <comm/picl/picl.h>
#include "piclext.h"

int loop_neigh(int offset)
{
  int nproc,top,ord,dir,me,host;
  int result;
  getarc0(&nproc,&top,&ord,&dir);
  who0(&nproc,&me,&host);

  /* hypercube topology */
  if (top == 1) {
    while (offset < 0) offset += nproc;
    result = gray0((my_loop_index() + offset) % nproc);
    }
  /* all other topologies */
  else {
      result = (me + offset)%nproc;
      if (result < 0) result += nproc;
    }
  return result;
}

int loop_out_neigh() { return loop_neigh(1); }

int loop_in_neigh() { return loop_neigh(-1); }

int my_loop_index()
{
  int nproc,top,ord,dir,me,host;
  int result;
  getarc0(&nproc,&top,&ord,&dir);
  who0(&nproc,&me,&host);

  if (top == 1) {
      result = ginv0(me);
    }
  else {
      result = me;
    }
  return result;
}

