
#include <stdio.h>
#include <sys/types.h>
#include <util/group/picl.h>
#include <util/group/message.h>

#ifndef OLDCLOCK
#  include <sys/types.h>
#  include <sys/time.h>
#  include <sys/resource.h>
#  if defined(sun) || defined(AIX)
#  include <unistd.h>
#  endif
#endif

static RefMessageGrp global_messagegrp;

int
host0()
{
  return 0;
}

void
check0(int)
{
}

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC ((int) 1000000)
#endif
#ifndef OLDCLOCK
double clock0()
{
  double res;
  struct rusage r;

  getrusage(RUSAGE_SELF,&r);

  res = r.ru_utime.tv_sec + r.ru_stime.tv_sec;
  res += 0.000001 * ( r.ru_utime.tv_usec + r.ru_stime.tv_usec );
  return res;
}
#else
double clock0()
{
  return(((double)clock())/CLOCKS_PER_SEC);
}
#endif

void
close0(int)
{
  global_messagegrp = 0;
}

void
message0(char*message)
{
  printf ("%s\n", message);
}

void
load0(char* file, int node)
{
  fprintf(stderr,"load0 called but this is not a PICL host\n");
  abort();
}

void
open0(int*numproc, int*me, int*host)
{
  *numproc = 1;
  *me = 0;
  *host = -1;
  global_messagegrp = new ProcMessageGrp();
}

void
open0_messagegrp(int*numproc, int*me, int*host, const RefMessageGrp&grp)
{
  global_messagegrp = grp;
  *numproc = grp->n();
  *me = grp->me();
  *host = -1;
}

int
probe0(int msgtype)
{
  return global_messagegrp->probet(msgtype);
}

void
recv0(void* buf, int bytes, int msgtype)
{
  global_messagegrp->raw_recvt(msgtype, buf, bytes);
}

void
recvinfo0(int* bytes, int* type, int* source)
{
  *bytes = global_messagegrp->last_size();
  *type = global_messagegrp->last_type();
  *source = global_messagegrp->last_source();
}

void
send0(void*buf, int bytes, int msgtype, int dest)
{
  global_messagegrp->raw_sendt(dest, msgtype, buf, bytes);
}

void
sync0()
{
  global_messagegrp->sync();
}

void
who0(int*numproc, int*me, int*host)
{
  *numproc = global_messagegrp->n();
  *me = global_messagegrp->me();
  *host = -1;
}

void
setarc0(int*nprocs,int*,int*,int*)
{
  if (*nprocs != global_messagegrp->n()) {
      fprintf(stderr,"setarc0: ignored attempt to change nprocs\n");
    }
}

void
getarc0(int*nprocs, int*top, int*ord, int*dir)
{
  *nprocs = global_messagegrp->n();
  *top = -1;
  *ord = 0;
  *dir = 1;
}

void
bcast0(void*buf, int bytes, int msgtype, int root)
{
  global_messagegrp->raw_bcast(buf, bytes, root);
}

static void(*global_comb)(void*target,void*source,int ndata, int datatype);
static int global_datatype;

template <class T>
void tcomb(T*target, T*source, int ndata)
{
  (*global_comb)((void*)target,(void*)source,ndata,global_datatype);
}

void
gcomb0(void*buf, int items, int datatype, int msgtype, int root,
       void(*comb)(void*target,void*source,int ndata, int datatype))
{
  global_comb = comb;
  global_datatype = datatype;

  switch(msgtype) {
  case 0:
      {
        void (*tmptcomb)(char*,char*,int) = tcomb;
        GrpFunctionReduce<char> r(tmptcomb);
        global_messagegrp->reduce((char*)buf, items, r, (char*)0, root);
      }
      break;
  case 1:
      {
        void (*tmptcomb)(short*,short*,int) = tcomb;
        GrpFunctionReduce<short> r(tmptcomb);
        global_messagegrp->reduce((short*)buf, items, r, (short*)0, root);
      }
      break;
  case 2:
      {
        void (*tmptcomb)(int*,int*,int) = tcomb;
        GrpFunctionReduce<int> r(tmptcomb);
        global_messagegrp->reduce((int*)buf, items, r, (int*)0, root);
      }
      break;
  case 3:
      {
        void (*tmptcomb)(long*,long*,int) = tcomb;
        GrpFunctionReduce<long> r(tmptcomb);
        global_messagegrp->reduce((long*)buf, items, r, (long*)0, root);
      }
      break;
  case 4:
      {
        void (*tmptcomb)(float*,float*,int) = tcomb;
        GrpFunctionReduce<float> r(tmptcomb);
        global_messagegrp->reduce((float*)buf, items, r, (float*)0, root);
      }
      break;
  case 5:
      {
        void (*tmptcomb)(double*,double*,int) = tcomb;
        GrpFunctionReduce<double> r(tmptcomb);
        global_messagegrp->reduce((double*)buf, items, r, (double*)0, root);
      }
      break;
  default:
      fprintf(stderr,"gcomb0: bad data type\n");
      abort();
    }
      
}

void
gor0(void*buf, int items, int datatype, int msgtype, int root)
{
  switch(msgtype) {
  case 0:
      { GrpArithmeticOrReduce<char> r;
        global_messagegrp->reduce((char*)buf, items, r, (char*)0, root);}
      break;
  case 1:
      { GrpArithmeticOrReduce<short> r;
        global_messagegrp->reduce((short*)buf, items, r, (short*)0, root);}
      break;
  case 2:
      { GrpArithmeticOrReduce<int> r;
        global_messagegrp->reduce((int*)buf, items, r, (int*)0, root);}
      break;
  case 3:
      { GrpArithmeticOrReduce<long> r;
        global_messagegrp->reduce((long*)buf, items, r, (long*)0, root);}
      break;
  default:
      fprintf(stderr,"gor0: bad data type\n");
      abort();
    }
}

void
gxor0(void*buf, int items, int datatype, int msgtype, int root)
{
  switch(msgtype) {
  case 0:
      { GrpArithmeticXOrReduce<char> r;
        global_messagegrp->reduce((char*)buf, items, r, (char*)0, root);}
      break;
  case 1:
      { GrpArithmeticXOrReduce<short> r;
        global_messagegrp->reduce((short*)buf, items, r, (short*)0, root);}
      break;
  case 2:
      { GrpArithmeticXOrReduce<int> r;
        global_messagegrp->reduce((int*)buf, items, r, (int*)0, root);}
      break;
  case 3:
      { GrpArithmeticXOrReduce<long> r;
        global_messagegrp->reduce((long*)buf, items, r, (long*)0, root);}
      break;
  default:
      fprintf(stderr,"gxor0: bad data type\n");
      abort();
    }
}

void
gand0(void*buf, int items, int datatype, int msgtype, int root)
{
  switch(msgtype) {
  case 0:
      { GrpArithmeticAndReduce<char> r;
        global_messagegrp->reduce((char*)buf, items, r, (char*)0, root);}
      break;
  case 1:
      { GrpArithmeticAndReduce<short> r;
        global_messagegrp->reduce((short*)buf, items, r, (short*)0, root);}
      break;
  case 2:
      { GrpArithmeticAndReduce<int> r;
        global_messagegrp->reduce((int*)buf, items, r, (int*)0, root);}
      break;
  case 3:
      { GrpArithmeticAndReduce<long> r;
        global_messagegrp->reduce((long*)buf, items, r, (long*)0, root);}
      break;
  default:
      fprintf(stderr,"gand0: bad data type\n");
      abort();
    }
}

void
gmax0(void*buf, int items, int datatype, int msgtype, int root)
{
  switch(msgtype) {
  case 0:
      { GrpMaxReduce<char> r;
        global_messagegrp->reduce((char*)buf, items, r, (char*)0, root);}
      break;
  case 1:
      { GrpMaxReduce<short> r;
        global_messagegrp->reduce((short*)buf, items, r, (short*)0, root);}
      break;
  case 2:
      { GrpMaxReduce<int> r;
        global_messagegrp->reduce((int*)buf, items, r, (int*)0, root);}
      break;
  case 3:
      { GrpMaxReduce<long> r;
        global_messagegrp->reduce((long*)buf, items, r, (long*)0, root);}
      break;
  case 4:
      { GrpMaxReduce<float> r;
        global_messagegrp->reduce((float*)buf, items, r, (float*)0, root);}
      break;
  case 5:
      { GrpMaxReduce<double> r;
        global_messagegrp->reduce((double*)buf, items, r, (double*)0, root);}
      break;
  default:
      fprintf(stderr,"gmax0: bad data type\n");
      abort();
    }
}

void
gprod0(void*buf, int items, int datatype, int msgtype, int root)
{
  switch(msgtype) {
  case 0:
      { GrpProductReduce<char> r;
        global_messagegrp->reduce((char*)buf, items, r, (char*)0, root);}
      break;
  case 1:
      { GrpProductReduce<short> r;
        global_messagegrp->reduce((short*)buf, items, r, (short*)0, root);}
      break;
  case 2:
      { GrpProductReduce<int> r;
        global_messagegrp->reduce((int*)buf, items, r, (int*)0, root);}
      break;
  case 3:
      { GrpProductReduce<long> r;
        global_messagegrp->reduce((long*)buf, items, r, (long*)0, root);}
      break;
  case 4:
      { GrpProductReduce<float> r;
        global_messagegrp->reduce((float*)buf, items, r, (float*)0, root);}
      break;
  case 5:
      { GrpProductReduce<double> r;
        global_messagegrp->reduce((double*)buf, items, r, (double*)0, root);}
      break;
  default:
      fprintf(stderr,"gprod0: bad data type\n");
      abort();
    }
}

void
gsum0(void*buf, int items, int datatype, int msgtype, int root)
{
  switch(msgtype) {
  case 0:
      { GrpSumReduce<char> r;
        global_messagegrp->reduce((char*)buf, items, r, (char*)0, root);}
      break;
  case 1:
      { GrpSumReduce<short> r;
        global_messagegrp->reduce((short*)buf, items, r, (short*)0, root);}
      break;
  case 2:
      { GrpSumReduce<int> r;
        global_messagegrp->reduce((int*)buf, items, r, (int*)0, root);}
      break;
  case 3:
      { GrpSumReduce<long> r;
        global_messagegrp->reduce((long*)buf, items, r, (long*)0, root);}
      break;
  case 4:
      { GrpSumReduce<float> r;
        global_messagegrp->reduce((float*)buf, items, r, (float*)0, root);}
      break;
  case 5:
      { GrpSumReduce<double> r;
        global_messagegrp->reduce((double*)buf, items, r, (double*)0, root);}
      break;
  default:
      fprintf(stderr,"gsum0: bad data type\n");
      abort();
    }
}

void
gmin0(void*buf, int items, int datatype, int msgtype, int root)
{
  switch(msgtype) {
  case 0:
      { GrpMinReduce<char> r;
        global_messagegrp->reduce((char*)buf, items, r, (char*)0, root);}
      break;
  case 1:
      { GrpMinReduce<short> r;
        global_messagegrp->reduce((short*)buf, items, r, (short*)0, root);}
      break;
  case 2:
      { GrpMinReduce<int> r;
        global_messagegrp->reduce((int*)buf, items, r, (int*)0, root);}
      break;
  case 3:
      { GrpMinReduce<long> r;
        global_messagegrp->reduce((long*)buf, items, r, (long*)0, root);}
      break;
  case 4:
      { GrpMinReduce<float> r;
        global_messagegrp->reduce((float*)buf, items, r, (float*)0, root);}
      break;
  case 5:
      { GrpMinReduce<double> r;
        global_messagegrp->reduce((double*)buf, items, r, (double*)0, root);}
      break;
  default:
      fprintf(stderr,"gmin0: bad data type\n");
      abort();
    }
}

void
barrier0()
{
  fprintf(stderr,"barrier0: not implemented\n");
  abort();
}

int
gray0(int i)
{
  return( (i>>1)^i ) ;
}

int
ginv0(int i)
{
  int k ;
  k = i ;
  while ( k > 0 ) {
      k >>= 1 ;
      i ^= k ;
    }
  return (i) ;
}
