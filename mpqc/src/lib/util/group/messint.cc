
#include <util/group/message.h>

#define CLASSNAME intMessageGrp
#define PARENTS public MessageGrp
#include <util/class/classia.h>
void *
intMessageGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  MessageGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

intMessageGrp::intMessageGrp()
{
}

intMessageGrp::intMessageGrp(const RefKeyVal& keyval):
  MessageGrp(keyval)
{
}

intMessageGrp::~intMessageGrp()
{
  delete[] source_seq;
  delete[] target_seq;
}

void
intMessageGrp::initialize(int me, int n, int nbits)
{
  int i;
  
  // Initialize the arrays storing the next sequence number.
  source_seq = new int[n];
  target_seq = new int[n];
  for (i=0; i<n; i++) {
      source_seq[i] = 0;
      target_seq[i] = 0;
    }

  // Find how many bits are needed to store the source information.
  int tmp = 1;
  for (i=1; tmp<n; i++) {
      tmp = tmp<<1;
    }
  src_nbit = i;

  // The remaining bits are for sequence information.
  seq_nbit = nbits - ctl_nbit - src_nbit;
  if (seq_nbit < 8) {
      fprintf(stderr,"intMessageGrp: not enough bits in underlying msgtype\n");
      abort();
    }

  // The bits available to user defined messages overlap
  // with the seq and the src bits.
  typ_nbit = seq_nbit + src_nbit;

  // Compute the shifts needed to construct a message type.
  seq_shift = 0;
  typ_shift = 0;
  src_shift = seq_shift + seq_nbit;
  ctl_shift = src_shift + src_nbit;

  // Compute the masks for each field.
  seq_mask = (1<<seq_nbit) - 1;
  typ_mask = (1<<typ_nbit) - 1;
  src_mask = (1<<src_nbit) - 1;

  // Complete initialization of the base class.
  MessageGrp::initialize(me, n);
}

int
intMessageGrp::msgtype_typ(int msgtype)
{
  return msgtype>>typ_shift & typ_mask;
}

int
intMessageGrp::typ_msgtype(int usrtype)
{
  return usrtype<<typ_shift | 1<<ctl_shift;
}

int
intMessageGrp::seq_msgtype(int source, int seq)
{
  return source<<src_shift | seq<<seq_shift | 2<<ctl_shift;
}

void
intMessageGrp::raw_send(int target, void* data, int nbyte)
{
  int& seq = target_seq[target];
  int msgtype = seq_msgtype(me(),seq);
#ifdef DEBUG
  printf("node %d sending to %d(%d) msgtype = %d\n",
         me(),target,seq,msgtype);
#endif
  basic_send(target, msgtype, data, nbyte);
  if (seq >= seq_mask) seq = 0;
  else seq++;
}

void
intMessageGrp::raw_recv(int sender, void* data, int nbyte)
{
  int& seq = source_seq[sender];
  int msgtype = seq_msgtype(sender,seq);
#ifdef DEBUG
  printf("node %d receiving from %d(%d) msgtype = %d\n",
         me(),sender,seq,msgtype);
#endif
  basic_recv(msgtype, data, nbyte);
#ifdef DEBUG
  printf("node %d received %d\n",me(),msgtype);
#endif
  if (seq >= seq_mask) seq = 0;
  else seq++;
}

void
intMessageGrp::raw_sendt(int target, int msgtype, void* data, int nbyte)
{
  basic_send(target, typ_msgtype(msgtype), data, nbyte);
}

void
intMessageGrp::raw_recvt(int type, void* data, int nbyte)
{
  basic_recv(typ_msgtype(type), data, nbyte);
}

int
intMessageGrp::probet(int type)
{
  return basic_probe(typ_msgtype(type));
}
