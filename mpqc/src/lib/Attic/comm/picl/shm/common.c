
int host0_() {return host0();}

double clock0_()
{
  return(clock0());
}

void close0_(release)
int *release;
{
  close0(*release);
}

void open0_(numproc,me,host)
int *numproc,*me,*host;
{
  open0(numproc,me,host);
}

void recv0_(buf,bytes,type)
int *type, *bytes;
char *buf;
{
  recv0(buf,*bytes,*type);
}

void recvinfo0_(bytes,msgtype,source)
int *bytes,*msgtype,*source;
{
  recvinfo0(bytes,msgtype,source);
}

void send0_(buf,bytes,msgtype,dest)
int *bytes,*msgtype,*dest;
char *buf;
{
  send0(buf,*bytes,*msgtype,*dest);
}

void who0_(numproc,me,host)
int *numproc,*me,*host;
{
  who0(numproc,me,host);
}

