
#define IN_PICL_CC

#include <stdio.h>
#include <sys/types.h>
#include <util/group/picl.h>
#include <util/group/message.h>

#if defined(PARAGON)
#  include <nx.h>
#elif defined(HAVE_MPI)
#  include <mpi.h>
#elif defined(OLDCLOCK)
#else
#  include <sys/types.h>
#  include <sys/time.h>
#  include <sys/resource.h>
#  if defined(sun) || defined(AIX)
#  include <unistd.h>
#  endif
#endif

#if defined(AIX)
extern "C" {
int getrusage (
  int Who,
  struct rusage *RUsage); }
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

#if defined(PARAGON)
double clock0()
{
  static double first;
  static initp = 0;
  double t;
  if (!initp) {
    first = dclock();
    initp = 1;
  }
  t = dclock();
  return (t - first);
}
#elif defined(HAVE_MPI)
double clock0()
{
  static double first;
  static initp = 0;
  double t;
  if (!initp) {
    first = MPI_Wtime();
    initp = 1;
  }
  t = MPI_Wtime();
  return (t - first);
}
#elif defined(OLDCLOCK)
double clock0()
{
  return(((double)clock())/CLOCKS_PER_SEC);
}
#else
double clock0()
{
  double res;
  struct rusage r;

  getrusage(RUSAGE_SELF,&r);

  res = r.ru_utime.tv_sec + r.ru_stime.tv_sec;
  res += 0.000001 * ( r.ru_utime.tv_usec + r.ru_stime.tv_usec );
  return res;
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
  global_messagegrp = MessageGrp::initial_messagegrp();
  if (global_messagegrp.null()) {
      global_messagegrp = MessageGrp::get_default_messagegrp();
    }
  *numproc = global_messagegrp->n();
  *me = global_messagegrp->me();
  *host = -1;
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
      fprintf(stderr,"setarc0: ignored attempt to change nprocs:\n");
      fprintf(stderr," nprocs = %d, attempted to set to %d\n",
              global_messagegrp->n(), *nprocs);
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

typedef void (*combroutine)(void *, void *, int, int);

void
gcomb0(void*buf, int items, int datatype, int msgtype, int root,
       combroutine comb)
{
  global_comb = comb;
  global_datatype = datatype;

  switch(datatype) {
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
  switch(datatype) {
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
  switch(datatype) {
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
  switch(datatype) {
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
  switch(datatype) {
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
  switch(datatype) {
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
  switch(datatype) {
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
  switch(datatype) {
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

///////////////////////////////////////////////////////////////////////////
// These routines are not a part of picl, but their use crept into the code.

#if defined(PARAGON)
# include <nx.h>
  // things missing from nx.h:
  extern "C" {
    void gdlow(double x[], long n, double work[]);
    void gdhigh(double x[], long n, double work[]);
    void gdsum(double x[], long n, double work[]);
  }
#elif defined(I860)
# include <cube.h>
#endif

#ifdef I860
/* the following is a fast global sum for long vectors.
 * It was written by Stan Erwin, an Intel employee at the NIH.
 *
 * FORCE_TYPE is defined in <cube.h>
 *
 * modifications of van de Geijn's gdcomp to
 * 1. vectorize the add (via a daxpy call)
 * 2. change the three trip protocol to a 2 trip protocal using
 *    force types and control messages
 */

void
gdcomb(int n, double* x, double* y, int dim, int idim)
{
  int i,ibit, l1, l2 , me, tempdim,msgid,msgid2,dummy,ione;
  double done;
  ione = 1;
  done = 1.0;
  me = mynode();

  if (dim == idim) return;

  l1 = n/2;
  l2 = n-l1;
  ibit = 1<<dim;
  if ((me&ibit) == 0) {
    msgid = _irecv(FORCE_TYPE + me^ibit,&y[l1], l2*sizeof(double));
    _csend(me,dummy, 0, me^ibit,0);
    _crecv(me^ibit,dummy, 0);
    msgid2 = _isend(FORCE_TYPE + me,x, l1*sizeof(double), me^ibit,0);
    _msgwait(msgid);

    daxpy_(&l2, &done, &(y[l1]), &ione, &(x[l1]), &ione);

    tempdim = dim + 1;
    _msgwait(msgid2);

    gdcomb(l2, &x[l1], &y[l1],tempdim,idim);

    msgid = _irecv(FORCE_TYPE + me^ibit,x, l1*sizeof(double));
    _csend(me,dummy, 0, me^ibit,0);
    _crecv(me^ibit,dummy, 0);
    _csend(FORCE_TYPE + me,&x[l1], l2*sizeof(double), me^ibit,0);
    _msgwait(msgid);
    }
  else {
    msgid = _irecv(FORCE_TYPE + me^ibit,y,l1*sizeof(double));
    _csend(me,dummy, 0, me^ibit,0);
    _crecv(me^ibit,dummy, 0);
    msgid2 =_isend(FORCE_TYPE + me,&x[l1], l2*sizeof(double), me^ibit,0);
    _msgwait(msgid);

    daxpy_(&l1, &done, y, &ione, x, &ione);

    tempdim = dim + 1;
    _msgwait(msgid2);

    gdcomb(l1, x, y,tempdim,idim);

    msgid = _irecv(FORCE_TYPE + me^ibit,&x[l1], l2*sizeof(double));
    _csend(me,dummy, 0, me^ibit,0);
    _crecv(me^ibit,dummy, 0);
    _csend(FORCE_TYPE + me,x, l1*sizeof(double), me^ibit,0);
    _msgwait(msgid);
    }
  }

#endif // I860

void gop1(double* val, int len, double* tmp, int op, int type)
{
#if defined(PARAGON)
  switch (op) {
    case 'm':
      gdlow(val,len,tmp);
      break;
    case 'M':
      gdhigh(val,len,tmp);
      break;
    case '+':
      gdsum(val,len,tmp);
      break;
    }
#elif defined(I860)
  switch (op) {
    case 'm':
      gmin0(val,len,5,type,0);
      bcast0(val,len,type,0);
      break;
    case 'M':
      gmax0(val,len,5,type,0);
      bcast0(val,len,type,0);
      break;
    case '+':
      gdcomb(len,val,tmp,0,cubedim0());
      break;
    }
#else
  switch (op) {
    case 'm':
      gmin0(val,len,5,type,0);
      break;
    case 'M':
      gmax0(val,len,5,type,0);
      break;
    case '+':
      gsum0(val,len,5,type,0);
      break;
    }
  bcast0(val,len*sizeof(double),type,0);
#endif

}

void gop0(double *val,int len,int op,int type)
{
  double *tmp = (double*)malloc(sizeof(double)*len);
  if (!tmp) {
      fprintf(stderr,"gop0: couldn't allocate %d bytes\n",len);
      exit(1);
    }
  gop1(val,len,tmp,op,type);
  free(tmp);
}

/* These routines are used to implement a gop0 for signed chars. */
static void min_schar (signed char*x,signed char*t,int l,int type)
{
  int i;
  for (i=0; i<l; i++) { if (t[i] < x[i]) x[i] = t[i]; }
}
static void max_schar (signed char*x,signed char*t,int l,int type)
{
  int i;
  for (i=0; i<l; i++) { if (t[i] > x[i]) x[i] = t[i]; }
}
static void sum_schar (signed char*x,signed char*t,int l,int type)
{
  int i;
  for (i=0; i<l; i++) x[i] += t[i];
}

/* A gop0 for signed chars. */
void gop0_sc(signed char*val,int len,int op,int type)
{
  switch (op) {
    case 'm':
      gcomb0(val,len,0,type,0,(combroutine)min_schar);
      break;
    case 'M':
      gcomb0(val,len,0,type,0,(combroutine)max_schar);
      break;
    case '+':
      gcomb0(val,len,0,type,0,(combroutine)sum_schar);
      break;
    }
  bcast0(val,len,type,0);
}

#if !defined(I860)
int cubedim_() {
  int n,p,me, host,i=0;
  who0(&p,&me,&host);
  n=p;
  while (n > 0) {n >>= 1; i++;}
  return (i-1);
}
#ifndef NCUBE
int cubedim() {
  int n,p,me, host,i=0;
  who0(&p,&me,&host);
  n=p;
  while (n > 0) {n >>= 1; i++;}
  return (i-1);
}
int mynode() {
  int p, me, host;
  who0(&p, &me, &host);
  return (me);
}
int mynode_() {
  int p, me, host;
  who0(&p, &me, &host);
  return (me);
}
#endif /*NCUBE*/
int numnodes_() {
  int p, me, host;
  who0(&p, &me, &host);
  return (p);
}
int numnodes() {
  int p, me, host;
  who0(&p, &me, &host);
  return (p);
}

int infocount() {
  int bytes,type,source;
  recvinfo0(&bytes,&type,&source);
  return(bytes);
}
#endif /* !I860 */

int mynode0() {
  int p, me, host;
  who0(&p, &me, &host);
  return (me);
}
int mynode0_() {
  int p, me, host;
  who0(&p, &me, &host);
  return (me);
}

int cubedim0_() {
  int n,p,me, host,i=0;
  who0(&p,&me,&host);
  n=p;
  while (n > 0) {n >>= 1; i++;}
  return (i-1);
}
int cubedim0() {
  int n,p,me, host,i=0;
  who0(&p,&me,&host);
  n=p;
  while (n > 0) {n >>= 1; i++;}
  return (i-1);
}
int numnodes0_() {
  int p, me, host;
  who0(&p, &me, &host);
  return (p);
}
int numnodes0() {
  int p, me, host;
  who0(&p, &me, &host);
  return (p);
}

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

void
picl_prober()
{
  int bytes,type,source;

  if (!probe0(-1)) return;

  recvinfo0(&bytes,&type,&source);

  printf("On node %3d found a message of type %3d, bytes %3d, from %3d\n",
         mynode0(),type,bytes,source);
}
