
 /* The high level PICL routines specific for the Paragon. */

#include <stdio.h>
#include <comm/picl/picl.h>
#include <nx.h>

void gopf2(data,ndata,tmp,msg_type,dest_node,comb)
     char *data;
     int ndata;
     char *tmp;
     int dest_node,msg_type;
     long (*comb)();
{
  gopf(data,ndata,tmp,comb);
}

void bcast0(char*buf,int bytes,int type,int root)
{

  /* the root sends the message */
  if (root == mynode()) {
      csend(type,buf,bytes,-1,0);
    }
  /* everybody else receives */
  else {
      crecv(type,buf,bytes);
    }

}
