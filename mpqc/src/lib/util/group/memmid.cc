
#ifndef _util_group_memmid_cc
#define _util_group_memmid_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <unistd.h>
#include <util/group/memmid.h>

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

static int ack_serial_number = 0;

void
MIDMemoryGrp::handler(long *mid)
{
  handler(data_request_buffer_, mid);
}

void
MIDMemoryGrp::handler(MemoryDataRequest& buffer, long *msgid_arg)
{
  int i;
  int mid;
  int dsize;
  int junk;
  double *source_data;
  double *data;

  if (msgid_arg) wait(*msgid_arg);

  const char *request_name = buffer.request_string();
  MemoryDataRequest::Request request = buffer.request();
  int offset = buffer.offset();
  int size = buffer.size();
  int node = buffer.node();
  int from_type = data_type_from_handler_;
  int to_type = data_type_to_handler_;

  const char *handlerstr;
  if (msgid_arg) handlerstr = "====:";
  else handlerstr = "----:";

  if (DEBREQ) buffer.print(handlerstr);

  switch (request) {
  case MemoryDataRequest::Deactivate:
      mid = send(&junk, sizeof(junk), me(), from_type);
      wait(mid);
      break;
  case MemoryDataRequest::Sync:
      nsync_++;
      if (msgid_arg) activate();
      break;
  case MemoryDataRequest::Retrieve:
      mid = send(&data_[offset], size, node, from_type);
      wait(mid);
      PRINTF(("%s sent %d bytes at byte offset %d (mid = %d)\n",
              handlerstr, size, offset, mid));
      if (msgid_arg) activate();
      break;
  case MemoryDataRequest::Replace:
      mid = recv(&data_[offset], size, node, to_type);
      wait(mid);
      PRINTF(("%s replaced %d bytes at byte offset %d (mid = %d)\n",
              handlerstr, size, offset, mid));
      if (use_acknowledgments_) {
          junk = (me() & 0xff) + (ack_serial_number++ << 8);
          mid = send(&junk, sizeof(junk), node, from_type);
          wait(mid);
          PRINTF(("%s sent ack %d:%d to %d\n",
                  handlerstr, junk&0xff, junk>>8, node));
        }
      if (msgid_arg) activate();
      break;
  case MemoryDataRequest::DoubleSum:
      dsize = size/sizeof(double);
      data = new double[dsize];
      mid = recv(data, size, node, to_type);
      wait(mid);
      PRINTF(("%s summing %d bytes (%d doubles) (mid = %d)\n",
              handlerstr, size, size/sizeof(double), mid));
      source_data = (double*) &data_[offset];
      for (i=0; i<dsize; i++) {
          source_data[i] += data[i];
        }
      delete[] data;
      if (use_acknowledgments_) {
          junk = (me() & 0xff) + (ack_serial_number++ << 8);
          mid = send(&junk, sizeof(junk), node, from_type);
          wait(mid);
          PRINTF(("%s sent ack %d:%d to %d\n",
                  handlerstr, junk&0xff, junk>>8, node));
        }
      if (msgid_arg) activate();
      break;
  default:
      fprintf(stderr, "mid_memory_handler: bad request id\n");
      sleep(1);
      abort();
    }
}

///////////////////////////////////////////////////////////////////////
// The MIDMemoryGrp class

#define CLASSNAME MIDMemoryGrp
#define PARENTS public ActiveMsgMemoryGrp
#include <util/class/classi.h>
void *
MIDMemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  ActiveMsgMemoryGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

MIDMemoryGrp::MIDMemoryGrp(const RefMessageGrp& msg):
  ActiveMsgMemoryGrp(msg)
{
  PRINTF(("MIDMemoryGrp CTOR\n", me()));
  data_request_type_ = 113;
  data_type_to_handler_ = 114;
  data_type_from_handler_ = 115;
  nsync_ = 0;
  use_acknowledgments_ = 0;
  use_active_messages_ = 1;
  active_ = 0;
}

MIDMemoryGrp::~MIDMemoryGrp()
{
  PRINTF(("%d: ~MIDMemoryGrp\n", me()));
}

void
MIDMemoryGrp::print_memreq(MemoryDataRequest &req,
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
MIDMemoryGrp::activate()
{
  if (!active_) {
      if (use_active_messages_) {
          data_request_mid_ = postrecv(data_request_buffer_.data(),
                                       data_request_buffer_.nbytes(),
                                       data_request_type_);
        }
      else {
          data_request_mid_ = recv(data_request_buffer_.data(),
                                   data_request_buffer_.nbytes(),
                                   -1,
                                   data_request_type_);
        }
      PRINTF(("activated memory request handler\n"));
      active_ = 1;
    }
}

void
MIDMemoryGrp::deactivate()
{
  if (active_) {
      sync();
      active_ = 0;
    }
}

void
MIDMemoryGrp::retrieve_data(void *data, int node, int offset, int size)
{
  MemoryDataRequestQueue q;

  long oldlock = lock();

  MemoryDataRequest buf(MemoryDataRequest::Retrieve,
                         me(), offset, size);
  if (DEBREQ) print_memreq(buf, "retrieve:", node);
  int mid = send(buf.data(), buf.nbytes(), node, data_request_type_);
  do_wait("retrieve: send req", mid, q, buf.nbytes());

  mid = recv(data, size, node, data_type_from_handler_);
  do_wait("retrieve: recv dat", mid, q, size);

  flush_queue(q);

  unlock(oldlock);
}

void
MIDMemoryGrp::replace_data(void *data, int node, int offset, int size)
{
  MemoryDataRequestQueue q;

  long oldlock = lock();

  MemoryDataRequest buf(MemoryDataRequest::Replace,
                        me(), offset, size);
  if (DEBREQ) print_memreq(buf, "replace:", node);
  int mid = send(buf.data(), buf.nbytes(), node, data_request_type_);
  do_wait("replace: send req", mid, q, buf.nbytes());

  mid = send(data, size, node, data_type_to_handler_);
  do_wait("replace: send dat", mid, q, size);

  if (use_acknowledgments_) {
      int junk;
      mid = recv(&junk, sizeof(junk), node, data_type_from_handler_);
      do_wait("replace: recv ack", mid, q, sizeof(junk));
      PRINTF(("replace: got ack %d:%d\n", junk&0xff, junk>>8));
    }

  flush_queue(q);

  unlock(oldlock);
}

void
MIDMemoryGrp::sum_data(double *data, int node, int offset, int size)
{
  MemoryDataRequestQueue q;

  long oldlock = lock();

  MemoryDataRequest buf(MemoryDataRequest::DoubleSum,
                        me(), offset, size);
  if (DEBREQ) print_memreq(buf, "sum:", node);
  int mid = send(buf.data(), buf.nbytes(), node, data_request_type_);
  do_wait("sum: send req", mid, q, buf.nbytes());

  mid = send(data, size, node, data_type_to_handler_);
  do_wait("sum: send dat", mid, q, size);

  if (use_acknowledgments_) {
      int junk;
      recv(&junk, sizeof(junk), node, data_type_from_handler_);
      do_wait("sum: recv ack", mid, q, sizeof(junk));
      PRINTF(("sum: got ack %d:%d\n", junk&0xff, junk>>8));
    }

  flush_queue(q);

  unlock(oldlock);
}

void
MIDMemoryGrp::sync()
{
  int i;
  int mid;

  long oldlock = lock();

  //sleep(1);

  PRINTF(("MIDMemoryGrp::sync() entered\n"));

  if (me() == 0) {
      // keep processing requests until all nodes have sent a
      // sync request
      PRINTF(("outside while loop nsync_ = %d, n() = %d\n", nsync_, n()));
      while (nsync_ < n() - 1) {
          PRINTF(("waiting for sync or other msg data_request_mid = %d\n",
                  data_request_mid_));
          mid = wait(data_request_mid_);
          PRINTF(("got mid = %d (data_request_mid_ = %d)\n", mid,
                  data_request_mid_));
          if (mid == data_request_mid_) {
              handler();
              activate();
            }
          else {
              printf("WARNING: MIDMemoryGrp::sync: stray message\n");
            }
        }
      nsync_ = 0;

      PRINTF(("notifying nodes that sync is complete\n"));
      // tell all nodes that they can proceed
      MemoryDataRequest buf(MemoryDataRequest::Sync, me());
      for (i=1; i<n(); i++) {
          if (DEBREQ) print_memreq(buf, "sync:", i);
          mid = send(buf.data(), buf.nbytes(), i, data_request_type_);
          PRINTF(("node %d can proceed (mid = %d)\n", i, mid));
          mid = wait(mid);
          PRINTF(("node %d got proceed message\n", i, mid));
        }
    }
  else {
      // let node 0 know that i'm done
      MemoryDataRequest buf(MemoryDataRequest::Sync, me());
      if (DEBREQ) print_memreq(buf, "sync:", 0);
      mid = send(buf.data(), buf.nbytes(), 0, data_request_type_);
      PRINTF(("sending sync (mid = %d)\n", mid));
      int tmpmid;
      do {
          tmpmid = wait(mid, data_request_mid_);
          if (tmpmid == data_request_mid_) {
              PRINTF(("sync: wait: handling request %d-%d\n",
                      data_request_buffer_.node(),
                      data_request_buffer_.serial_number()));
              handler();
              activate();
            }
          else if (tmpmid != mid) {
              printf("MIDMemoryGrp: sync: stray message id = %d\n",
                     tmpmid);
              sleep(1);
              abort();
            }
        } while (tmpmid != mid);
      // watch for the done message from 0 or request messages
      while (!nsync_) {
          PRINTF(("waiting for sync\n"));
          mid = wait(data_request_mid_);
          PRINTF(("in sync got mid = %d\n", mid));
          if (mid == data_request_mid_) {
              handler();
              activate();
            }
          else {
              printf("WARNING: MIDMemoryGrp::sync: stray message\n");
            }
        }
      nsync_ = 0;
    }

  PRINTF(("MIDMemoryGrp::sync() done\n"));

  unlock(oldlock);
}

void
MIDMemoryGrp::do_wait(const char *msg, int mid,
                      MemoryDataRequestQueue &q, size_t expectedsize)
{
  long oldlock = lock();

  int tmpmid;
  do {
      PRINTF(("%s: wait: waiting\n", msg));
      tmpmid = wait(data_request_mid_, mid);
      if (tmpmid == data_request_mid_) {
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
              handler();
            }
          activate();
        }
      else if (tmpmid != mid) {
          printf("MIDMemoryGrp: %s: stray message id = %d\n",
                 msg, tmpmid);
          sleep(1);
          abort();
        }
      else {
          PRINTF(("%s: wait: wait complete\n", msg));
        }
    } while (tmpmid != mid);

  unlock(oldlock);
}

void
MIDMemoryGrp::flush_queue(MemoryDataRequestQueue &q)
{
  long oldlock = lock();

  for (int i=0; i<q.n(); i++) {
      PRINTF(("processing queued request %d-%d\n",
              data_request_buffer_.node(),
              data_request_buffer_.serial_number()));
      handler(q[i]);
    }
  q.clear();

  unlock(oldlock);
}

#endif
