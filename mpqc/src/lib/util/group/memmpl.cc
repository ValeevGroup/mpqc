
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


#define ACK 0

#define DEBUG 0
#define DEBREQ 0

#if DEBUG
#  undef PRINTF
#  define PRINTF(args) do { DISABLE; \
                            printf args; \
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
static int ack_serial_number = 0;

static void
mpl_memory_handler_(MemoryDataRequest& buffer, int*msgid_arg)
{
  if (!global_mpl_mem) {
      fprintf(stderr, "mpl_memory_handler: called while inactive\n");
      sleep(1);
      abort();
    }

  //global_mpl_mem->active_ = 0;
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

  const char *request_name = buffer.request_string();
  MemoryDataRequest::Request request = buffer.request();
  int offset = buffer.offset();
  int size = buffer.size();
  int node = buffer.node();
  int from_type = global_mpl_mem->data_type_from_handler_;
  int to_type = global_mpl_mem->data_type_to_handler_;

  const char *handlerstr;
  if (msgid_arg) handlerstr = "====:";
  else handlerstr = "----:";

  if (DEBREQ) buffer.print(handlerstr);

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
      if (msgid_arg) global_mpl_mem->activate();
      break;
  case MemoryDataRequest::Retrieve:
      mpc_send(&global_mpl_mem->data_[offset], size, node, from_type, &mid);
      mpc_wait(&mid, &count);
      PRINTF(("%s sent %d bytes at byte offset %d (mid = %d)\n",
              handlerstr, size, offset, mid));
      if (msgid_arg) global_mpl_mem->activate();
      break;
  case MemoryDataRequest::Replace:
      //PRINTF(("%d:: about to replace %d bytes at byte offset %d\n",
      //        global_mpl_mem->me(), size, offset));
      source = node;
      type = to_type;
      mpc_recv(&global_mpl_mem->data_[offset], size, &source, &type, &mid);
      mpc_wait(&mid, &count);
      PRINTF(("%s replaced %d bytes at byte offset %d (mid = %d)\n",
              handlerstr, size, offset, mid));
#if ACK
      junk = (global_mpl_mem->me() & 0xff) + (ack_serial_number++ << 8);
      mpc_send(&junk, sizeof(junk), node, from_type, &mid);
      mpc_wait(&mid, &count);
      PRINTF(("%s sent ack %d:%d to %d\n",
              handlerstr, junk&0xff, junk>>8, node));
#endif
      if (msgid_arg) global_mpl_mem->activate();
      break;
  case MemoryDataRequest::DoubleSum:
      dsize = size/sizeof(double);
      data = new double[dsize];
      source = node;
      type = to_type;
      mpc_recv(data, size, &source, &type, &mid);
      mpc_wait(&mid, &count);
      if (count != size) {
          fprintf(stderr, "MPLMessageGrp: handler: DoubleSum: wrong size\n");
          sleep(1);
          abort();
        }
      PRINTF(("%s summing %d bytes (%d doubles) (mid = %d)\n",
              handlerstr, size, size/sizeof(double), mid));
      source_data = (double*) &global_mpl_mem->data_[offset];
      for (i=0; i<dsize; i++) {
          source_data[i] += data[i];
        }
      delete[] data;
#if ACK
      junk = (global_mpl_mem->me() & 0xff) + (ack_serial_number++ << 8);
      mpc_send(&junk, sizeof(junk), node, from_type, &mid);
      mpc_wait(&mid, &count);
      PRINTF(("%s sent ack %d:%d to %d\n",
              handlerstr, junk&0xff, junk>>8, node));
#endif
      if (msgid_arg) global_mpl_mem->activate();
      break;
  default:
      fprintf(stderr, "mpl_memory_handler: bad request id\n");
      sleep(1);
      abort();
    }
}

static void
mpl_memory_handler(int*msgid_arg)
{
  mpl_memory_handler_(global_mpl_mem->data_request_buffer_, msgid_arg);
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
  data_request_type_ = 113;
  data_type_to_handler_ = 114;
  data_type_from_handler_ = 115;
  global_mpl_mem = this;
  //active_ = 0;
  nsync_ = 0;

  activate();

  PRINTF(("%d: data_ = 0x%x\n", me(), data_));
}

MPLMemoryGrp::~MPLMemoryGrp()
{
  PRINTF(("%d: ~MPLMemoryGrp\n", me()));
  deactivate();
}

void
MPLMemoryGrp::print_memreq(MemoryDataRequest &req,
                           const char *msg, int target)
{
  if (msg == 0) msg = "";

  char newmsg[80];
  if (target == -1) {
      sprintf(newmsg, "%d: %s", me(), msg);
    }
  else {
      sprintf(newmsg, "%d->%d: %s", me(), target, msg);
    }
  req.print(newmsg);
}

void
MPLMemoryGrp::activate()
{
  //if (!active_) {
      global_type = data_request_type_;
      global_source = DONTCARE;
      mpc_rcvncall(data_request_buffer_.data(),
                   data_request_buffer_.nbytes(),
                   (int*)&global_source, (int*)&global_type, (int*)&global_mid,
                   mpl_memory_handler);
      //mpc_recv(data_request_buffer_.data(),
      //         data_request_buffer_.nbytes(),
      //         (int*)&global_source, (int*)&global_type, (int*)&global_mid);
      PRINTF(("activated memory request handler (mid = %d)\n",
              global_mid));
  //  }
  //active_ = 1;
  //reactivate_ = 1;
}

void
MPLMemoryGrp::deactivate()
{
  if (!global_mpl_mem) return;

  //if (active_) {
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
  //  }
  //active_ = 0;

  int oldlock;
  mpc_lockrnc(1, &oldlock);

  global_mpl_mem = 0;

  mpc_lockrnc(oldlock, &oldlock);
}

void
MPLMemoryGrp::retrieve_data(void *data, int node, int offset, int size)
{
  MemoryDataRequestQueue q;

  int oldlock;
  mpc_lockrnc(1, &oldlock);

  MemoryDataRequest buf(MemoryDataRequest::Retrieve,
                         me(), offset, size);
  if (DEBREQ) print_memreq(buf, "retrieve:", node);
  int mid;
  mpc_send(buf.data(), buf.nbytes(), node, data_request_type_, &mid);
  do_wait("retrieve: send req", mid, q, buf.nbytes());

  int source = node;
  int type = data_type_from_handler_;
  mpc_recv(data, size, &source, &type, &mid);
  do_wait("retrieve: recv dat", mid, q, size);

  flush_queue(q);

  mpc_lockrnc(oldlock, &oldlock);
}

void
MPLMemoryGrp::replace_data(void *data, int node, int offset, int size)
{
  MemoryDataRequestQueue q;

  int oldlock;
  mpc_lockrnc(1, &oldlock);

  MemoryDataRequest buf(MemoryDataRequest::Replace,
                        me(), offset, size);
  if (DEBREQ) print_memreq(buf, "replace:", node);
  int mid;
  mpc_send(buf.data(), buf.nbytes(), node, data_request_type_, &mid);
  do_wait("replace: send req", mid, q, buf.nbytes());

  mpc_send(data, size, node, data_type_to_handler_, &mid);
  do_wait("replace: send dat", mid, q, size);

#if ACK
  int junk;
  int source = node;
  int type = data_type_from_handler_;
  mpc_recv(&junk, sizeof(junk), &source, &type, &mid);
  do_wait("replace: recv ack", mid, q, sizeof(junk));
  PRINTF(("replace: got ack %d:%d\n", junk&0xff, junk>>8));
#endif

  flush_queue(q);

  mpc_lockrnc(oldlock, &oldlock);
}

void
MPLMemoryGrp::sum_data(double *data, int node, int offset, int size)
{
  MemoryDataRequestQueue q;

  int oldlock;
  mpc_lockrnc(1, &oldlock);

  int doffset = offset/sizeof(double);
  int dsize = size/sizeof(double);

  MemoryDataRequest buf(MemoryDataRequest::DoubleSum,
                        me(), offset, size);
  if (DEBREQ) print_memreq(buf, "sum:", node);
  int mid;
  mpc_send(buf.data(), buf.nbytes(), node, data_request_type_, &mid);
  do_wait("sum: send req", mid, q, buf.nbytes());

  mpc_send(data, size, node, data_type_to_handler_, &mid);
  do_wait("sum: send dat", mid, q, size);

#if ACK
  int junk;
  int source = node;
  int type = data_type_from_handler_;
  mpc_recv(&junk, sizeof(junk), &source, &type, &mid);
  do_wait("sum: recv ack", mid, q, sizeof(junk));
  PRINTF(("sum: got ack %d:%d\n", junk&0xff, junk>>8));
#endif

  flush_queue(q);

  mpc_lockrnc(oldlock, &oldlock);
}

void
MPLMemoryGrp::sync()
{
  //if (!active_) {
  //    msg_->sync();
  //    return;
  //  }

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
              activate();
            }
          else {
              printf("WARNING: MPLMemoryGrp::sync: stray message\n");
            }
        }
      nsync_ = 0;

      PRINTF(("notifying nodes that sync is complete\n"));
      // tell all nodes that they can proceed
      MemoryDataRequest buf(MemoryDataRequest::Sync, me());
      for (i=1; i<n(); i++) {
          if (DEBREQ) print_memreq(buf, "sync:", i);
          mpc_send(buf.data(), buf.nbytes(), i, data_request_type_, &mid);
          PRINTF(("node %d can proceed (mid = %d)\n", i, mid));
          mpc_wait(&mid, &count); // noncritical
          PRINTF(("node %d got proceed message\n", i, mid));
        }
    }
  else {
      // let node 0 know that i'm done
      MemoryDataRequest buf(MemoryDataRequest::Sync, me());
      if (DEBREQ) print_memreq(buf, "sync:", 0);
      mpc_send(buf.data(), buf.nbytes(), 0, data_request_type_, &mid);
      PRINTF(("sending sync (mid = %d)\n", mid));
      int tmpmid;
      do {
          tmpmid = DONTCARE;
          mpc_wait(&tmpmid, &count);
          if (tmpmid == global_mid) {
              PRINTF(("sync: wait: handling request %d-%d\n",
                      data_request_buffer_.node(),
                      data_request_buffer_.serial_number()));
              mpl_memory_handler(0);
              activate();
            }
          else if (tmpmid != mid) {
              printf("MPLMemoryGrp: sync: stray message id = %d size = %d\n",
                     tmpmid, count);
              sleep(1);
              abort();
            }
        } while (tmpmid != mid);
      // watch for the done message from 0 or request messages
      while (!nsync_) {
          mid = DONTCARE;
          PRINTF(("waiting for sync\n"));
          mpc_wait(&mid, &count); // critical (handled correctly)
          PRINTF(("in sync got mid = %d\n", mid));
          if (mid == global_mid) {
              mpl_memory_handler(0);
              activate();
            }
          else {
              printf("WARNING: MPLMemoryGrp::sync: stray message\n");
            }
        }
      nsync_ = 0;
    }

  PRINTF(("MPLMemoryGrp::sync() done\n"));

  mpc_lockrnc(oldlock, &oldlock);
}

void
MPLMemoryGrp::do_wait(const char *msg, int mid,
                      MemoryDataRequestQueue &q, size_t expectedsize)
{
  int oldlock;
  mpc_lockrnc(1, &oldlock);

  int tmpmid;
  size_t count;
  do {
      tmpmid = DONTCARE;
      PRINTF(("%s: wait: waiting\n", msg));
      mpc_wait(&tmpmid, &count);
      if (tmpmid == global_mid) {
          if (me() < data_request_buffer_.node()) {
              PRINTF(("%s: wait: queueing request %d-%d\n",
                      msg,
                      data_request_buffer_.node(),
                      data_request_buffer_.serial_number()));
              q.push(data_request_buffer_);
            }
          else {
              PRINTF(("%s: wait: handling request %d-%d\n",
                      msg,
                      data_request_buffer_.node(),
                      data_request_buffer_.serial_number()));
              mpl_memory_handler(0);
            }
          activate();
        }
      else if (tmpmid != mid) {
          printf("MPLMemoryGrp: %s: stray message id = %d size = %d\n",
                 msg, tmpmid, count);
          sleep(1);
          abort();
        }
      else {
          PRINTF(("%s: wait: wait complete\n", msg));
          if (count != expectedsize) {
              printf("MPLMemoryGrp: %s: got wrong size %d, expected %d\n",
                     msg, count, expectedsize);
              sleep(1);
              abort();
            }
        }
    } while (tmpmid != mid);

  mpc_lockrnc(oldlock, &oldlock);
}

void
MPLMemoryGrp::flush_queue(MemoryDataRequestQueue &q)
{
  int oldlock;
  mpc_lockrnc(1, &oldlock);

  for (int i=0; i<q.n(); i++) {
      PRINTF(("processing queued request %d-%d\n",
              data_request_buffer_.node(),
              data_request_buffer_.serial_number()));
      mpl_memory_handler_(q[i], 0);
    }
  q.clear();

  mpc_lockrnc(oldlock, &oldlock);
}

#endif
