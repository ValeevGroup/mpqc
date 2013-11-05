#ifndef MPQC_ARRAY_THREAD_HPP
#define MPQC_ARRAY_THREAD_HPP

#include "mpqc/mpi.hpp"
#include "mpqc/utility/timer.hpp"
#include "mpqc/array/forward.hpp"
#include "mpqc/array/socket.hpp"

#include <memory>

#include <boost/thread/thread.hpp>
#include "boost/thread/tss.hpp"

namespace mpqc {
namespace detail {
namespace ArrayServer {

    struct array_proxy {
        
        static const size_t BUFFER = (32<<20);

        struct Descriptor {
            ArrayBase *object;
            size_t rank;
            size_t count;         
        };

        struct Segment {
            size_t size;
            std::vector<range> extents;
        };

        array_proxy(ArrayBase *object,
                    const std::vector<range> &r) {
            this->object = object;
            this->rank = r.size();
            this->count = 0;

            size_t block = 1;
	    size_t N = BUFFER/sizeof(double);

            for (int i = 0; i < rank-1; ++i)
                block *= r[i].size();
            assert(block < N);
            
            foreach (range rj, split(r.back(), N/block)) {
                for (int i = 0; i < rank-1; ++i) {
                    data.push_back(*r[i].begin());
                    data.push_back(*r[i].end());
                }
                data.push_back(*rj.begin());
                data.push_back(*rj.end());
                ++count;
            }
        }

        array_proxy(Descriptor ds) {
            this->object = ds.object;
            this->rank = ds.rank;
            this->count = ds.count;
            this->data.resize(2*rank*count);
        }

        Descriptor descriptor() const {
            Descriptor ds;
            ds.object = this->object;
            ds.rank = this->rank;
            ds.count = this->count;
            return ds;
        }

        std::vector<Segment> segments() const {
            std::vector<Segment> segments;
            auto data = this->data.begin();
            for (int j = 0; j < count; ++j) {
                std::vector<range> r;
                size_t size = 1;
                for (int i = 0; i < rank; ++i) {
                    r.push_back(range(data[0], data[1]));
                    size *= r.back().size();
                    data += 2;
                }
                segments.push_back(Segment());
                segments.back().size = size;
                segments.back().extents = r;
            }
            return segments;
        }

    public:

        ArrayBase *object;
        int rank, count;
        std::vector<int> data;

    };


    struct Message {

        enum Request { INVALID = 0x0,
                       JOIN   = 0x100,
                       SYNC   = 0x101,
                       WRITE  = 0x200,
                       READ   = 0x201
        };

        Request request;
        array_proxy::Descriptor dataspace;
        int src, tag;

        explicit Message(int tag = 0, Request r = INVALID)
            : tag(tag), request(r) {}

        Message(int tag, Request op, array_proxy::Descriptor ds) {
            this->request = op;
            this->dataspace = ds;
            this->tag = tag;
        }

    };


    struct Thread : boost::noncopyable {

    private:
        explicit Thread(MPI::Comm comm)
	    : comm_(comm), tag_(1<<20)
	{
            buffer_ = malloc(array_proxy::BUFFER);
            this->socket_.start();
            this->servers_ = comm.allgather(this->socket_.address());
	    this->thread_ = new boost::thread(&Thread::run, this);
        }

    public:

	enum {
	    RECV_MASK = 1<<22,
	    SEND_MASK = 1<<23
	};

        ~Thread() {
	    //comm_.printf("~Thread\n");
	    sync();
	    join();
	    delete this->thread_;
            free(this->buffer_);
        }

	void join() {
	    int tag = next();
	    send(Message(tag, Message::JOIN), comm_.rank());
	    this->thread_->join();
	    //comm_.printf("thread joined\n");
	}

	void sync() const {
	    int tag = next();
            send(Message(tag, Message::SYNC), comm_.rank());
            comm_.recv<Message>(comm_.rank(), tag | SEND_MASK);
	}

	void send(Message msg, int proc) const {
            msg.src = comm_.rank();
            ArraySocket::send(&msg, this->servers_.at(proc));
	    //send(&msg, sizeof(Message), MPI_BYTE, proc, this->tag_);
	}

	/** send to server thread */
	void send(const void *data,
		  size_t count, MPI_Datatype type,
		  int proc, int tag) const {
	    assert(!(tag & SEND_MASK));
	    assert(!(tag & RECV_MASK));
	    comm_.send(data, count, type, proc, tag | RECV_MASK);
	}

	/** recv from server thread */
	void recv(void *data,
		  size_t count, MPI_Datatype type,
		  int proc, int tag) const {
	    assert(!(tag & SEND_MASK));
	    assert(!(tag & RECV_MASK));
	    comm_.recv(data, count, type, proc, tag | SEND_MASK);
	}

	static std::shared_ptr<Thread>& instance() {
	    static std::shared_ptr<Thread> thread;
	    if (!thread.get()) {
		MPI::initialize(MPI_THREAD_MULTIPLE);
		thread.reset(new Thread(MPI::Comm(MPI_COMM_WORLD)));
	    }
	    return thread;
	}

        static void run(Thread *thread) {
	    //mutex_.unlock();
	    thread->loop();
	    //mutex_.lock();
	}

    public:

        int next() const {
	    const unsigned int N = 1 << 21;
	    boost::mutex::scoped_lock lock(mutex_);
            return int(N + (next_++ % N));
        }

	int translate(MPI::Comm comm1, int rank1) const {
	    int rank2;;
	    MPI_Group group1, group2;
	    MPI_Comm_group(comm1, &group1);
	    MPI_Comm_group(this->comm_, &group2);
	    MPI_Group_translate_ranks(group1, 1, &rank1, group2, &rank2);
	    assert(rank2 != MPI_UNDEFINED);
	    return rank2;
	}	    

    private:

	void loop() {

            // std::cout << "thread/rank: "
            //           << boost::this_thread::get_id() << "/" << comm_.rank()
            //           << std::endl;

            while (1) {

		Message msg;

                this->socket_.wait(&msg);

		//std::cout << MPI::get_processor_name() << std::endl;

                // MPI_Status status;
                // status = comm_.recv(&msg, sizeof(Message), MPI_BYTE,
		// 		    MPI_ANY_SOURCE, this->tag_ | RECV_MASK);

                // MPI_Request request =
                //     comm_.irecv(&msg, sizeof(Message), MPI_BYTE,
		// 	       MPI_ANY_SOURCE, this->tag_);
                // status = MPI::wait(request, 10);

		//comm_.printf("Message received %i\n", msg.request);

                if (msg.request == Message::READ) {
                    //printf("Message::READ\n");
                    read(msg, msg.dataspace);
                    continue;
                }
                if (msg.request == Message::WRITE) {
                    //printf("Message::WRITE\n");
                    write(msg, msg.dataspace);
                    continue;
                }
                if (msg.request == Message::SYNC) {
                    sync(msg);
                    continue;
                }
                if (msg.request == Message::JOIN) {
                    //comm_.printf("Message::JOIN\n");
		    //comm_.ssend(msg, status.MPI_SOURCE, msg.tag);
                    break;
                }
		printf("invalid message request %i\n", msg.request);
		throw std::runtime_error("invalid message");
            }
        }

    private:

	MPI::Comm comm_;
        void *buffer_;
	boost::thread *thread_;

        ArraySocket socket_;
        std::vector<ArraySocket::Address> servers_;

        int tag_;
	mutable unsigned int next_;
	mutable boost::mutex mutex_;

        void sync(Message msg) {
	    //comm_.printf("thread message/sync dst=%i tag=%i\n", status.MPI_SOURCE, tag);
            comm_.send(Message(Message::SYNC), msg.src, msg.tag | SEND_MASK);
        }

        void read(Message msg, array_proxy::Descriptor ds) { 
            io<Message::READ>(array_proxy(ds), msg.src, msg.tag); 
        }

        void write(Message msg, array_proxy::Descriptor ds) {
            io<Message::WRITE>(array_proxy(ds), msg.src, msg.tag);
        }

    private:

        template<Message::Request OP>
        void io(array_proxy ds, int proc, int tag) {
	    // comm_.printf("thread recv descriptor bytes=%lu src=%i tag=%i\n",
	    // 		 ds.data.size()*sizeof(int), proc, tag);
            comm_.recv(&ds.data[0], ds.data.size(), MPI_INT, proc, tag | RECV_MASK);
            const auto &segments = ds.segments();
            double* buffer = static_cast<double*>(this->buffer_);
            for (int i = 0; i < segments.size(); ++i) {
                const auto &extents = segments[i].extents;
                if (OP == Message::WRITE) {
		    // comm_.printf("thread recv/write bytes=%lu, dst=%i tag=%i\n",
		    // 		 segments[i].size*8, proc, tag);
                    comm_.recv(buffer, segments[i].size, MPI_DOUBLE,
			       proc, tag | RECV_MASK);
		    //std::cout << "write " << o.extents() << std::endl;
                    ds.object->put(extents, buffer);
                }
                if (OP == Message::READ) {
		    // comm_.printf("thread read/send bytes=%lu, dst=%i tag=%i\n",
		    // 		 segments[i].size*8, proc, tag);
		    //std::cout << "read " << o.extents() << std::endl;
                    ds.object->get(extents, buffer);
                    comm_.send(buffer, segments[i].size, MPI_DOUBLE,
			       proc, tag | SEND_MASK);
                }
            }
        }

    };


}
} // namespace detail
} // namespace mpqc

namespace mpqc {
namespace detail {

    struct array_thread_comm {

        typedef ArrayServer::Message Message;
	typedef ArrayServer::Thread Thread;
        typedef ArrayServer::array_proxy array_proxy;

        explicit array_thread_comm(MPI_Comm comm)
            : comm_(comm)
        {
	    thread_ = Thread::instance();
            sync();
            //printf("mpqc::Comm initialized\n");
        }

	Thread& thread() {
	    return *this->thread_;
	}

        void sync() const {
            //comm_.printf("sync\n");
            comm_.barrier();
	    thread_->sync();
            comm_.barrier();
        }
        
        template<typename T>
        void write(const T *data, ArrayBase *object,
                   const std::vector<range> &r, int rank) const {
            //printf("Comm::write\n");
            io<Message::WRITE>((T*)data, object, r, rank);
        }

        template<typename T>
        void read(T *data, ArrayBase *object,
                  const std::vector<range> &r, int rank) const {
            //printf("Comm::read\n");
            io<Message::READ>(data, object, r, rank);
        }

        MPI::Comm comm() const {
          return comm_;
        }

    private:

        MPI::Comm comm_;
	std::shared_ptr<Thread> thread_;

	int tag() const {
	    static boost::thread_specific_ptr<int> tag;
	    if (!tag.get()) tag.reset(new int(thread_->next()));
	    return *tag;
	}

        template<Message::Request OP, typename T>
        void io(T* buffer, ArrayBase *object,
		const std::vector<range> &r,
                int proc) const
        {

            static_assert(OP == Message::WRITE ||
                          OP == Message::READ,
                          "invalid OP");

	    //proc = thread_->translate(this->comm_, proc);
            int tag = this->tag();
            array_proxy ds(object, r);
	    
	    {
		//MPQC_PROFILE_LINE;
		//printf("message dst=%i tag=%i\n", proc, thread_->tag());
		thread_->send(Message(tag, OP, ds.descriptor()), proc);
	    }

	    {
		//MPQC_PROFILE_LINE;
		//printf("descriptor dst=%i tag=%i\n", proc, tag);
		thread_->send(&ds.data[0], ds.data.size(), MPI_INT, proc, tag);
	    }

            auto segments = ds.segments();
            std::vector<MPI_Request> requests(ds.count);

	    mpqc::timer t;
	    size_t total = 0;
            for (int i = 0; i < segments.size(); ++i) {
                size_t size = segments[i].size;
                if (OP == Message::READ) {
		    //MPQC_PROFILE_LINE;
		    // printf("recv segment %i bytes=%lu proc=%i tag=%i\n",
		    // 	   i, size*sizeof(T), proc, tag);
                    /*requests[i] = i*/
                    thread_->recv(buffer, size*sizeof(T), MPI_BYTE, proc, tag);
		}
                if (OP == Message::WRITE) {
		    //MPQC_PROFILE_LINE;
 		    // printf("send segment %i bytes=%lu proc=%i tag=%i\n",
		    // 	   i, size*sizeof(T), proc, tag);
                    /*requests[i] = i*/
                    thread_->send(buffer, size*sizeof(T), MPI_BYTE, proc, tag);
		}
                buffer += size;
		total += size;
            }
            //MPI_Waitall(n, &requests[0], MPI_STATUSES_IGNORE);
	    // double mb = (total*sizeof(double))/1e6;
	    // printf("I/O %s: %f Mbytes, %f Mbytes/s\n",
	    // 	   (OP == Message::WRITE ? "WRITE" : "READ"), mb, mb/t);
        }

    };


}
}

#endif // MPQC_ARRAY_THREAD_HPP
