/*      PICL to PVM interface library      
  All questions and criticisms should be directed to 
 Juan Meza  Center for Computational Engineering
	   Sandia National Laboratory
    		Livermore CA 94551                         */

#undef ARCH
#ifdef SGI
#define ARCH "SGI"
#endif

#ifdef SUN4
#define ARCH "SUN4"
#endif

#ifdef RS6000
#define ARCH "RIOS"
#endif

extern char* argv0;
static char *node_jobname;

#include <stdio.h>
#include <string.h>
#include "piclpvm.h"

#ifndef OLDCLOCK
#  include <sys/types.h>
#  include <sys/time.h>
#  include <sys/resource.h>
#  if defined(sun) || defined(AIX)
#  include <unistd.h>
#  endif
#endif

/* 
  Here we specify a HOSTID.  Its id must be different than all of
  the node ids.
 */

#define HOSTID (MAXPROCS-1)


/* Set the maximum number of processors (only used in inst_index) */

#define MAXPROCS 1025

/* Define some externals that will store info between calls */
static int last_count;   /* holds the size of the last msg sent or recvd */
static int last_node;    /* holds the id of where last msg came from */
static int last_type;    /* holds the type of where last msg recvd */
static int nprocs;       /* holds the number of processors allocated */

/* These are set by the setarc0 and getarc0 call. */
static int arc_nproc;
static int arc_top=-1;
static int arc_ord=0;
static int arc_dir=1;

/* There is never a host program in this version. */
static int hostavail = 0;

void send0();
void who0();

int host0() {return hostavail;}

void sync0()
{
  int i;
  int nproc,me,host,next;
  int sync;

  who0(&nproc,&me,&host);

  next = me + 1;
  if (next == nproc) next = 0;

  for (i=0; i<2; i++) {
    if (me==0) {
      send0(&sync,sizeof(sync),8881+i,next);
      }
    else {
      recv0(&sync,sizeof(sync),8881+i);
      send0(&sync,sizeof(sync),8881+i,next);
      }
    if (me==0) {
      recv0(&sync,sizeof(sync),8881+i);
      }
    }

}

void setarc0(np,top,ord,dir)
int *np;
int *top;
int *ord;
int *dir;
{
  arc_nproc = nprocs;
  *np = arc_nproc;
  *top = arc_top;
  *ord = arc_ord;
  *dir = arc_dir;
}

void getarc0(nprocs,top,ord,dir)
int *nprocs;
int *top;
int *ord;
int *dir;
{
  *nprocs = arc_nproc;
  *top = arc_top;
  *ord = arc_ord;
  *dir = arc_dir;
}

void message0(message)
char *message;
{
  printf ("%s\n", message);
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

void close0(release)
     int release;
{
  leave();
}

/* This is not available on the nodes. */
hidden_load0(file, node)
     char *file;
     int node;
{
  int i;
  int instance;

  if (node != -1) {
    printf("Error. This version will only allow node = -1.\n");
    return(-1);
  }

  for (i=0; i < (hostavail?nprocs:(nprocs-1)); i++) {
    instance = initiate(file,ARCH);
    }

}

void open0(numproc, me, host)
     int *numproc, *me, *host;
{
  int node;

  node_jobname = strrchr(argv0,'/');
  if (!node_jobname) node_jobname = argv0;
  else node_jobname++;

  node = enroll(node_jobname);

  if (node == 0) {
    nprocs = atoi(getenv("NUMPROC"));
    if (!nprocs) nprocs = 1;
    if (!hostavail) {
      int i;
      hidden_load0(node_jobname,-1);
      for (i=1; i<nprocs; i++) {
        initsend();
        putnint(&nprocs, 1);
        snd(node_jobname, i, 12);
        }
      }
    }
  else {
    /* the nodes other than zero get info from the host */
    rcv(12);
    getnint(&nprocs,1);
    }

#if 0
  { char comp[100]; int instance;
    whoami(comp,&instance);
    printf(" component = \"%s\", instance = %d, node = %d\n",
           comp,instance,node);
   }
#endif

  *me      = node;
  *numproc = nprocs;
  *host    = HOSTID;

}


int probe0(msgtype)
int msgtype;
{
  int probeval = probe(msgtype);
  if (probeval < 0) probeval = 0;
  return probeval;
}

recv0(buf,bytes,type)
     char *buf;
     int bytes, type;
{

  rcv(type);
  getbytes(buf, bytes);

}

void recvinfo0(bytes, msgtype, source)
int *bytes, *msgtype, *source;
{
  int instance;
  char component[100];

  rcvinfo(bytes, msgtype, component, &instance);
  *source  = instance;
}
 

void send0(buf, bytes, type, dest)
     char *buf;
     int bytes, type, dest;
{
  int  instance;

  initsend();
  putbytes(buf, bytes);

  /* printf("snd(\"%s\",%d,%d)\n",node_jobname,dest,type); */
  snd(node_jobname, dest, type);
}

void who0(numproc,me,host)
     int *numproc,*me,*host;
{
  char component[100];
  int  instance;

  whoami(component,&instance);

  *numproc = nprocs;
  *me      = instance;
  *host    = HOSTID;
}

