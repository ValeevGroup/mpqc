
/* The fortran stubs for the high level PICL routines. */

void gsum0_(data, ndata, data_type, msgtype, root)
     void *data;
     int *ndata, *data_type,*msgtype, *root;
{
  gsum0(data, *ndata, *data_type, *msgtype, *root);
}

bcast0_(buf,bytes,type,root)
void *buf;
int *bytes,*type,*root;
{
  bcast0(buf,*bytes,*type,*root);
  }

