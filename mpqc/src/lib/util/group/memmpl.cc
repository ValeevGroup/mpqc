
#ifndef _util_group_memmpl_cc
#define _util_group_memmpl_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <unistd.h>
#include <util/group/memmpl.h>

#include <mpi.h>
#include <mpproto.h>

#define DISABLE do { fflush(stdout); } while(0)
#define ENABLE do { fflush(stdout); } while(0)


#define DEBUG 0

#if DEBUG
#  undef PRINTF
#  define PRINTF(args) do { DISABLE; \
                            printf args; \
                            fflush(stdout); \
                            ENABLE; \
                           } while(0)
#else
#  undef PRINTF
#  define PRINTF(args)
#endif

///////////////////////////////////////////////////////////////////////
// The handler function and its data

static volatile int global_source, global_type, global_mid;
static MPLMemoryGrp *global_mpl_mem = 0;

static void
mpl_memory_handler(int*msgid_arg)
{
  global_mpl_mem->active_ = 0;
  int i;
  int mid;
  int type;
  int source;
  int dsize;
  int junk;
  size_t count;
  double *source_data;
  double *data;

  //PRINTF(("entered mpl_memory_handler()\n"));

  if (msgid_arg) mpc_wait(msgid_arg, &count);

  const char *request_name
      = global_mpl_mem->data_request_buffer_.request_string();
  MemoryDataRequest::Request request
      = global_mpl_mem->data_request_buffer_.request();
  int offset = global_mpl_mem->data_request_buffer_.offset();
  int size = global_mpl_mem->data_request_buffer_.size();
  int node = global_mpl_mem->data_request_buffer_.node();
  int from_type = global_mpl_mem->data_type_from_handler_;
  int to_type = global_mpl_mem->data_type_to_handler_;

  global_mpl_mem->data_request_buffer_.print("handler:");

  switch (request) {
  case MemoryDataRequest::Deactivate:
      mpc_send(&junk, sizeof(junk), global_mpl_mem->me(), from_type, &mid);
      mpc_wait(&mid, &count);
      break;
  case MemoryDataRequest::Sync:
      global_mpl_mem->nsync_++;
      //if (global_mpl_mem->reactivate_
      //    || (global_mpl_mem->me() == 0
      //        &&(global_mpl_mem->nsync_ < global_mpl_mem->n() - 1))
      //    || (global_mpl_mem->me() != 0
      //        &&(global_mpl_mem->nsync_ < 1))) {
      //    global_mpl_mem->activate();
      //  }
      if (1||msgid_arg) global_mpl_mem->activate();
      break;
  case MemoryDataRequest::Retrieve:
      mpc_send(&global_mpl_mem->data_[offset], size, node, from_type, &mid);
      mpc_wait(&mid, &count);
      PRINTF(("%d:: send %d bytes at byte offset %d (mid = %d)\n",
              global_mpl_mem->me(), size, offset, mid));
      if (1||msgid_arg) global_mpl_mem->activate();
      break;
  case MemoryDataRequest::Replace:
      //PRINTF(("%d:: about to replace %d bytes at byte offset %d\n",
      //        global_mpl_mem->me(), size, offset));
      source = node;
      type = to_type;
      mpc_recv(&global_mpl_mem->data_[offset], size, &source, &type, &mid);
      mpc_wait(&mid, &count);
      PRINTF(("%d:: replaced %d bytes at byte offset %d (mid = %d)\n",
              global_mpl_mem->me(), size, offset, mid));
      mpc_send(&junk, sizeof(junk), node, from_type, &mid);
      mpc_wait(&mid, &count);
      PRINTF(("%d:: sent go ahead to %d (mid = %d)\n",
              global_mpl_mem->me(), node, mid));
      if (1||msgid_arg) global_mpl_mem->activate();
      break;
  case MemoryDataRequest::DoubleSum:
      dsize = size/sizeof(double);
      data = new double[dsize];
      source = node;
      type = to_type;
      mpc_recv(data, size, &source, &type, &mid);
      mpc_wait(&mid, &count);
      PRINTF(("%d:: summing %d bytes (%d doubles) (mid = %d)\n",
              global_mpl_mem->me(), size, size/sizeof(double), mid));
      source_data = (double*) &global_mpl_mem->data_[offset];
      for (i=0; i<dsize; i++) {
          source_data[i] += data[i];
        }
      delete[] data;
      mpc_send(&junk, sizeof(junk), node, from_type, &mid);
      mpc_wait(&mid, &count);
      PRINTF(("%d:: sent go ahead to %d (mid = %d)\n",
              global_mpl_mem->me(), node, mid));
      if (1||msgid_arg) global_mpl_mem->activate();
      break;
  default:
      fprintf(stderr, "mpl_memory_handler: bad request id\n");
      sleep(1);
      abort();
    }
}

///////////////////////////////////////////////////////////////////////
// The MPLMemoryGrp class

#define CLASSNAME MPLMemoryGrp
#define PARENTS public ActiveMsgMemoryGrp
#include <util/class/classi.h>
void *
MPLMemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  ActiveMsgMemoryGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

MPLMemoryGrp::MPLMemoryGrp(const RefMessageGrp& msg,
                                   int localsize):
  ActiveMsgMemoryGrp(msg, localsize)
{
  if (global_mpl_mem) {
      fprintf(stderr, "MPLMemoryGrp: only one allowed at a time\n");
      sleep(1);
      abort();
    }
  data_request_type_ = 13;
  data_type_to_handler_ = 14;
  data_type_from_handler_ = 15;
  global_mpl_mem = this;
  active_ = 0;
  nsync_ = 0;

  activate();

  PRINTF(("%d: data_ = 0x%x\n", me(), data_));
}

MPLMemoryGrp::~MPLMemoryGrp()
{
  PRINTF(("%d: ~MPLMemoryGrp\n", me()));
  global_mpl_mem = 0;
  deactivate();
}

void
MPLMemoryGrp::print_memreq(MemoryDataRequest &req,
                           const char *msg, int target)
{
#if DEBUG
  if (msg == 0) msg = "";

  char newmsg[80];
  if (target == -1) {
      sprintf(newmsg, "%d: %s", me(), msg);
    }
  else {
      sprintf(newmsg, "%d->%d: %s", me(), target, msg);
    }
  req.print(newmsg);
#endif
}

void
MPLMemoryGrp::activate()
{
  if (!active_) {
      global_type = data_request_type_;
      global_source = DONTCARE;
      mpc_rcvncall(data_request_buffer_.data(),
                   data_request_buffer_.nbytes(),
                   (int*)&global_source, (int*)&global_type, (int*)&global_mid,
                   mpl_memory_handler);
      PRINTF(("%d: activate: type = %d, buf = 0x%x, mid = %d\n", me(),
              data_request_type_, data_request_buffer_.data(),
              global_mid));
    }
  active_ = 1;
  reactivate_ = 1;
}

void
MPLMemoryGrp::deactivate()
{
  if (active_) {
#if 1
      reactivate_ = 0;
      sync();

      //global_type = data_request_type_;
      //global_source = DONTCARE;
      //mpc_rcvncall(data_request_buffer_.data(),
      //             data_request_buffer_.nbytes(),
      //             (int*)&global_source, (int*)&global_type, (int*)&global_mid,
      //             0);
#else
      PRINTF(("%d: MPLMemoryGrp::deactivate init\n", me()));
      MemoryDataRequest buf(MemoryDataRequest::Deactivate);
      int mid;
      mpc_send(buf.data(), buf.nbytes(), me(), data_request_type_, &mid);
      PRINTF(("%d: MPLMemoryGrp::deactivate pending\n", me()));
      size_t count;
      mpc_wait(&mid, &count); // noncritical
      int junk;
      int source = global_mpl_mem->me();
      int type = data_type_from_handler_;
      mpc_recv(&junk, sizeof(junk), &source, &type, &mid);
      mpc_wait(&mid, &count); // noncritical
      PRINTF(("%d: MPLMemoryGrp::deactivate complete\n", me()));
#endif
    }
  active_ = 0;
}

void
MPLMemoryGrp::retrieve_data(void *data, int node, int offset, int size)
{
  int needs_activation = 0;

  int oldlock;
  mpc_lockrnc(1, &oldlock);

  int tmpmid;

  PRINTF(("%d: retrieve_data: int offset = %d int size = %d\n",
          me(), offset/sizeof(int), size/sizeof(int)));

  MemoryDataRequest buf(MemoryDataRequest::Retrieve,
                         me(), offset, size);
  print_memreq(buf, "retrieve:", node);
  PRINTF(("%d: sent request to %d\n", me(), node));
  int mid;
  size_t count;
  mpc_send(buf.data(), buf.nbytes(), node, data_request_type_, &mid);
  do {
      tmpmid = DONTCARE;
      mpc_wait(&tmpmid, &count); // critical
      PRINTF(("retrieve_data: get tmpmid %d (global = %d, mid = %d)\n",
              tmpmid, global_mid, mid));
      if (tmpmid == global_mid) {
          mpl_memory_handler(0);
          needs_activation = 1;
        }
      else if (tmpmid != mid) {
          printf("MPLMemoryGrp: retrieve_data: stray message\n");
          sleep(1);
          abort();
        }
    } while(tmpmid != mid);

  PRINTF(("%d: waiting for data from %d\n", me(), node));
  int source = node;
  int type = data_type_from_handler_;
  mpc_recv(data, size, &source, &type, &mid);
  mpc_wait(&mid, &count); // noncritical
  PRINTF(("%d: got data from %d\n", me(), node));

  mpc_lockrnc(oldlock, &oldlock);

  if (0&&needs_activation) activate();
}

void
MPLMemoryGrp::replace_data(void *data, int node, int offset, int size)
{
  int needs_activation = 0;

  int oldlock;
  mpc_lockrnc(1, &oldlock);

  int tmpmid;

  PRINTF(("%d: replace_data: int offset = %d int size = %d\n",
          me(), offset/sizeof(int), size/sizeof(int)));

  MemoryDataRequest buf(MemoryDataRequest::Replace,
                        me(), offset, size);
  print_memreq(buf, "replace:", node);
  int mid;
  size_t count;
  mpc_send(buf.data(), buf.nbytes(), node, data_request_type_, &mid);
  do {
      tmpmid = DONTCARE;
      mpc_wait(&tmpmid, &count); // critical
      PRINTF(("replace_data: get tmpmid %d (global = %d, mid = %d)\n",
              tmpmid, global_mid, mid));
      if (tmpmid == global_mid) {
          mpl_memory_handler(0);
          needs_activation = 1;
        }
      else if (tmpmid != mid) {
          printf("MPLMemoryGrp: replace_data: stray message\n");
          sleep(1);
          abort();
        }
    } while(tmpmid != mid);

  mpc_send(data, size, node, data_type_to_handler_, &mid);
  mpc_wait(&mid, &count); // noncritical

  int junk;
  int source = node;
  int type = data_type_from_handler_;
  mpc_recv(&junk, sizeof(junk), &source, &type, &mid);
  mpc_wait(&mid, &count); // noncritical

  mpc_lockrnc(oldlock, &oldlock);

  if (0&&needs_activation) activate();
}

void
MPLMemoryGrp::sum_data(double *data, int node, int offset, int size)
{
  int needs_activation = 0;

  int oldlock;
  mpc_lockrnc(1, &oldlock);

  int tmpmid;
  int doffset = offset/sizeof(double);
  int dsize = size/sizeof(double);

  PRINTF(("%d: sum_data: doffset = %d dsize = %d node = %d\n",
          me(), doffset, dsize, node));

  MemoryDataRequest buf(MemoryDataRequest::DoubleSum,
                        me(), offset, size);
  print_memreq(buf, "sum:", node);
  int mid;
  size_t count;
  mpc_send(buf.data(), buf.nbytes(), node, data_request_type_, &mid);
  PRINTF(("MPLMemoryGrp: sum_data: mid = %d global_mid = %d\n",
          mid, global_mid));
  do {
      tmpmid = DONTCARE;
      mpc_wait(&tmpmid, &count); // critical
      PRINTF(("sum_data: get tmpmid %d (global = %d, mid = %d)\n",
              tmpmid, global_mid, mid));
      if (tmpmid == global_mid) {
          mpl_memory_handler(0);
          needs_activation = 1;
        }
      else if (tmpmid != mid) {
          printf("MPLMemoryGrp: sum_data: stray message id = %d size = %d\n",
                 tmpmid, count);
          sleep(1);
          abort();
        }
    } while (tmpmid != mid);

  PRINTF(("%d: sum_data: sent request, sending data\n", me()));

  mpc_send(data, size, node, data_type_to_handler_, &mid);
  mpc_wait(&mid, &count); // noncritical
  

  PRINTF(("%d: sum_data: sent data, waiting for ack\n", me()));

  int junk;
  int source = node;
  int type = data_type_from_handler_;
  mpc_recv(&junk, sizeof(junk), &source, &type, &mid);
  mpc_wait(&mid, &count); // noncritical

  PRINTF(("%d: sum_data: got ack, done\n", me()));

  mpc_lockrnc(oldlock, &oldlock);

  if (0&&needs_activation) activate();
}

void
MPLMemoryGrp::sync()
{
  if (!active_) {
      msg_->sync();
      return;
    }

  int needs_activation = 0;

  int i;
  int mid;
  size_t count;

  int oldlock;
  mpc_lockrnc(1, &oldlock);

  //sleep(1);

  PRINTF(("MPLMemoryGrp::sync() entered\n"));

  if (me() == 0) {
      // keep processing requests until all nodes have sent a
      // sync request
      PRINTF(("outside while loop nsync_ = %d, n() = %d\n", nsync_, n()));
      PRINTF(("outside while 2\n"));
      while (nsync_ < n() - 1) {
          PRINTF(("inside while\n"));
          mid = DONTCARE;
          PRINTF(("waiting for sync or other msg\n"));
          mpc_wait(&mid, &count);
          PRINTF(("got mid = %d (global_mid = %d)\n", mid, global_mid));
          if (mid == global_mid) {
              mpl_memory_handler(0);
              needs_activation = 1;
            }
          else {
              printf("WARNING: MPLMemoryGrp::sync: stray message\n");
            }
        }
      nsync_ = 0;

      PRINTF(("notifying nodes that sync is complete\n"));
      // tell all nodes that they can proceed
      MemoryDataRequest buf(MemoryDataRequest::Sync);
      for (i=1; i<n(); i++) {
          print_memreq(buf, "sync:", i);
          mpc_send(buf.data(), buf.nbytes(), i, data_request_type_, &mid);
          PRINTF(("node %d can proceed (mid = %d)\n", i, mid));
          mpc_wait(&mid, &count); // noncritical
          PRINTF(("node %d got proceed message\n", i, mid));
        }
    }
  else {
      // let node 0 know that i'm done
      MemoryDataRequest buf(MemoryDataRequest::Sync);
      print_memreq(buf, "sync:", 0);
      mpc_send(buf.data(), buf.nbytes(), 0, data_request_type_, &mid);
      PRINTF(("sending sync (mid = %d)\n", mid));
      mpc_wait(&mid, &count); // critical
      // watch for the done message from 0 or request messages
      while (!nsync_) {
          mid = DONTCARE;
          PRINTF(("waiting for sync\n"));
          mpc_wait(&mid, &count); // critical (handled correctly)
          PRINTF(("in sync got mid = %d\n", mid));
          if (mid == global_mid) {
              mpl_memory_handler(0);
              needs_activation = 1;
            }
          else {
              printf("WARNING: MPLMemoryGrp::sync: stray message\n");
            }
        }
      nsync_ = 0;
    }

  PRINTF(("MPLMemoryGrp::sync() done\n"));

  mpc_lockrnc(oldlock, &oldlock);

  //sleep(1);

  if (0&&needs_activation) activate();
}

#endif
