
#ifndef _util_group_mempgon_cc
#define _util_group_mempgon_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/group/mempgon.h>

extern "C" {
#include <nx.h>
}

typedef void (*handlertype)(...);

#define DISABLE do { masktrap(1); fflush(stdout); } while(0)
#define ENABLE do { fflush(stdout); masktrap(0); } while(0)

#define PRINTF(args) do { DISABLE; \
                          printf args; \
                          fflush(stdout); \
                          ENABLE; \
                         } while(0)

#undef PRINTF
#define PRINTF(args)

///////////////////////////////////////////////////////////////////////
// The handler function and its data

static ParagonMemoryGrp *global_pgon_mem = 0;

void
paragon_memory_handler(long ptype, long pcount, long pnode, long pptype)
{
  global_pgon_mem->active_ = 0;
  PRINTF(("%d:: entered paragon_memory_handler()\n", global_pgon_mem->me()));
  int offset = global_pgon_mem->data_request_buffer_.offset();
  int size = global_pgon_mem->data_request_buffer_.size();
  int node = global_pgon_mem->data_request_buffer_.node();
  int from_type = global_pgon_mem->data_type_from_handler_;
  int to_type = global_pgon_mem->data_type_to_handler_;
  PRINTF(("%d:: handler: request = \"%s\" int offset = %d, "
          "int size = %d, node = %d\n",
          global_pgon_mem->me(),
          global_pgon_mem->data_request_buffer_.request_string(),
          offset/sizeof(int), size/sizeof(int), node));
  switch (global_pgon_mem->data_request_buffer_.request()) {
  case MemoryDataRequest::Deactivate:
      int junk;
      csend(from_type, (char*) &junk, sizeof(junk), global_pgon_mem->me(), 0);
      break;
  case MemoryDataRequest::Retrieve:
      csend(from_type,
            &global_pgon_mem->data_[offset],
            size, node, 0);
      PRINTF(("%d:: send %d ints at int offset %d\n",
              global_pgon_mem->me(), size/sizeof(int), offset/sizeof(int)));
      global_pgon_mem->activate();
      break;
  case MemoryDataRequest::Replace:
      PRINTF(("%d:: about to replace %d ints at int offset %d\n",
              global_pgon_mem->me(), size/sizeof(int), offset/sizeof(int)));
      crecvx(to_type, &global_pgon_mem->data_[offset], size,
             node, 0, msginfo);
      PRINTF(("%d:: replaced %d ints at int offset %d\n",
              global_pgon_mem->me(), size/sizeof(int), offset/sizeof(int)));
      csend(from_type, (char*) &junk, sizeof(junk), node, 0);
      global_pgon_mem->activate();
      break;
  case MemoryDataRequest::DoubleSum:
  {
      int dsize = size/sizeof(double);
      double *data = new double[dsize];
      PRINTF(("%d:: receiving %d bytes (%d doubles)\n",
              global_pgon_mem->me(), size, size/sizeof(double)));
      crecvx(to_type, (char *) data, size, node, 0, msginfo);
      double *source_data = (double*) &global_pgon_mem->data_[offset];
      for (int i=0; i<dsize; i++) {
          source_data[i] += data[i];
        }
      delete[] data;
      csend(from_type, (char*) &junk, sizeof(junk), node, 0);
      global_pgon_mem->activate();
  }
      break;

  default:
      cerr << scprintf("paragon_memory_handler: bad request id\n");
      abort();
    }
}

///////////////////////////////////////////////////////////////////////
// The ParagonMemoryGrp class

#define CLASSNAME ParagonMemoryGrp
#define HAVE_KEYVAL_CTOR
#define PARENTS public ActiveMsgMemoryGrp
#include <util/class/classi.h>
void *
ParagonMemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  ActiveMsgMemoryGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

ParagonMemoryGrp::ParagonMemoryGrp(const RefMessageGrp& msg):
  ActiveMsgMemoryGrp(msg)
{
  if (global_pgon_mem) {
      cerr << scprintf("ParagonMemoryGrp: only one allowed at a time\n");
      abort();
    }
  data_request_type_ = 13;
  data_type_to_handler_ = 14;
  data_type_from_handler_ = 15;
  global_pgon_mem = this;
  active_ = 0;
}

ParagonMemoryGrp::ParagonMemoryGrp(const RefKeyVal& keyval):
  ActiveMsgMemoryGrp(keyval)
{
  if (global_pgon_mem) {
      cerr << scprintf("ParagonMemoryGrp: only one allowed at a time\n");
      abort();
    }
  data_request_type_ = 13;
  data_type_to_handler_ = 14;
  data_type_from_handler_ = 15;
  global_pgon_mem = this;
  active_ = 0;
}

ParagonMemoryGrp::~ParagonMemoryGrp()
{
  PRINTF(("%d: ~ParagonMemoryGrp\n", me()));
  global_pgon_mem = 0;
  deactivate();
}

void
ParagonMemoryGrp::activate()
{
  if (!active_) {
      PRINTF(("%d: activate: hrecv(%d, 0x%x, %d, handler)\n", me(),
              data_request_type_, (char *) data_request_buffer_.data(),
              data_request_buffer_.nbytes()));
      hrecv(data_request_type_, (char *) data_request_buffer_.data(),
            data_request_buffer_.nbytes(),
            (handlertype)paragon_memory_handler);
      active_ = 1;
    }
}

void
ParagonMemoryGrp::deactivate()
{
  if (active_) {
#if 1
      hrecv(data_request_type_, 0, 0, 0);
#else
      PRINTF(("%d: ParagonMemoryGrp::deactivate init\n", me()));
      MemoryDataRequest buf(MemoryDataRequest::Deactivate);
      int msgid = isend(data_request_type_,
                        (char *) buf.data(), buf.nbytes(), me(), 0);
      PRINTF(("%d: ParagonMemoryGrp::deactivate pending\n", me()));
      msgwait(msgid);
      int junk;
      crecv(data_type_from_handler_, (char*) &junk, sizeof(junk));
      PRINTF(("%d: ParagonMemoryGrp::deactivate complete\n", me()));
#endif
      active_ = 0;
    }
}

void
ParagonMemoryGrp::retrieve_data(void *data, int node, int offset, int size)
{
  PRINTF(("%d: retrieve_data: int offset = %d int size = %d\n",
          me(), offset/sizeof(int), size/sizeof(int)));

  MemoryDataRequest buf(MemoryDataRequest::Retrieve,
                         me(), offset, size);
  PRINTF(("%d: sent request to %d\n", me(), node));
  csend(data_request_type_, (char *) buf.data(),
        buf.nbytes(), node, 0);

  PRINTF(("%d: waiting for data from %d\n", me(), node));
  crecv(data_type_from_handler_, (char *) data, size);
  PRINTF(("%d: got data from %d\n", me(), node));
}

void
ParagonMemoryGrp::replace_data(void *data, int node, int offset, int size)
{
  PRINTF(("%d: replace_data: int offset = %d int size = %d\n",
          me(), offset/sizeof(int), size/sizeof(int)));

  MemoryDataRequest buf(MemoryDataRequest::Replace,
                        me(), offset, size);
  csend(data_request_type_, (char *) buf.data(),
        buf.nbytes(), node, 0);

  csend(data_type_to_handler_, (char *) data, size, node, 0);

  int junk;
  crecv(data_type_from_handler_, (char*) &junk, sizeof(junk));
}

void
ParagonMemoryGrp::sum_data(double *data, int node, int offset, int size)
{
  int doffset = offset/sizeof(double);
  int dsize = size/sizeof(double);

  PRINTF(("%d: sum_data: doffset = %d dsize = %d node = %d\n",
          me(), doffset, dsize, node));

  MemoryDataRequest buf(MemoryDataRequest::DoubleSum,
                        me(), offset, size);
  csend(data_request_type_, (char *) buf.data(),
        buf.nbytes(), node, 0);

  PRINTF(("%d: sum_data: sent request, sending data\n", me()));

  csend(data_type_to_handler_, (char *) data, size, node, 0);

  PRINTF(("%d: sum_data: sent data, waiting for ack\n", me()));

  int junk;
  crecv(data_type_from_handler_, (char *) &junk, sizeof(junk));

  PRINTF(("%d: sum_data: got ack, done\n", me()));
}

#endif
