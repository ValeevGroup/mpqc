
/*
 * Copyright 2009 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This file is a part of the MPQC LMP2 library.
 *
 * The MPQC LMP2 library is free software: you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 */

#include <string.h>

#include <list>
#include <iostream>

#include <mpqc_config.h>
#include <chemistry/qc/lmp2/collective.h>

#ifdef HAVE_MPI
#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#endif

namespace sc {

#ifdef HAVE_MPI
void custom_alltoallv(void *sendbuf,
                      int *sendcnts,
                      int *sdispls,
		      MPI_Datatype sendtype,
                      void *recvbuf,
                      int *recvcnts,
                      int *rdispls,
		      MPI_Datatype recvtype,
		      MPI_Comm comm)
{
#if 1
  int me, nproc;
  MPI_Aint sendextent, recvextent, sendlb, recvlb;
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Type_get_extent(sendtype,&sendlb,&sendextent);
  MPI_Type_get_extent(recvtype,&recvlb,&recvextent);
  int source = (me + nproc - 1)%nproc; //  me - 1 mod nproc
  for (int next_node = (me+1)%nproc;
       next_node != me;
       next_node = (next_node+1)%nproc) {
      MPI_Request send_req, recv_req;
      MPI_Irecv(&((char*)recvbuf)[recvextent*rdispls[source]],recvcnts[source],recvtype,source,
                100+source,comm,&recv_req);
      MPI_Isend(&((char*)sendbuf)[sendextent*sdispls[next_node]],sendcnts[next_node],sendtype,next_node,
                100+me,comm,&send_req);
      MPI_Status status;
      MPI_Wait(&recv_req,&status); // wait for the recv to complete
      MPI_Wait(&send_req,&status); // wait for the send to complete
      source = (source + nproc - 1)%nproc; // source - 1 mod nproc
    }
  int selfbytes;
  if (sendextent*sendcnts[me] > recvextent*recvcnts[me]) selfbytes = recvextent*recvcnts[me];
  else selfbytes = sendextent*sendcnts[me];
  memcpy(&((char*)recvbuf)[recvextent*rdispls[me]],&((char*)sendbuf)[sendextent*sdispls[me]],selfbytes);
#else
  MPI_Alltoallv(sendbuf, sendcnts, sdispls, sendtype,
                                recvbuf, recvcnts, rdispls, recvtype,
                                comm);
#endif
}
#endif

#ifdef HAVE_MPI
class AllSumPipeline {
    double *data_;
    double *scratch_data_;
    int n_;
    int batch_;
    int up_;
    int dn_;
    int stage_;
    enum Stages { NoOp, RvUp, RvDn, SdUp, SdDn, Wait, Test, DSum, Done };
    static const char *stagenames[];
    static Stages stages_one1[];
    static Stages stages_top1[];
    static Stages stages_bot1[];
    static Stages stages_mid1[];
    static Stages stages_one2[];
    static Stages stages_top2[];
    static Stages stages_bot2[];
    static Stages stages_mid2[];
    Stages *stages_;
    MPI_Request req_;
    MPI_Comm comm_;
  public:
    AllSumPipeline(double *data, int n, int upstream, int downstream,
                   int batch, MPI_Comm comm, int half) {
      data_ = data;
      n_ = n;
      up_ = upstream;
      dn_ = downstream;
      if (half == 1) {
          if (up_ == -1 && dn_ == -1) stages_ = stages_one1;
          if (up_ == -1 && dn_ != -1) stages_ = stages_top1;
          if (up_ != -1 && dn_ == -1) stages_ = stages_bot1;
          if (up_ != -1 && dn_ != -1) stages_ = stages_mid1;
        }
      else {
          if (up_ == -1 && dn_ == -1) stages_ = stages_one2;
          if (up_ == -1 && dn_ != -1) stages_ = stages_top2;
          if (up_ != -1 && dn_ == -1) stages_ = stages_bot2;
          if (up_ != -1 && dn_ != -1) stages_ = stages_mid2;
        }
      stage_ = 0;
      batch_ = batch;
      comm_ = comm;
      scratch_data_ = new double[n];
    }
    AllSumPipeline(const AllSumPipeline &p) {
      data_ = p.data_;
      n_ = p.n_;
      up_ = p.up_;
      dn_ = p.dn_;
      stage_ = 0;
      stages_ = p.stages_;
      batch_ = p.batch_;
      comm_ = p.comm_;
      scratch_data_ = new double[n_];
    }
    ~AllSumPipeline() {
      delete[] scratch_data_;
    }
    const char *stagename(Stages s) { return stagenames[int(s)]; }
    void next() {
      bool failed;
      int count, rc = MPI_SUCCESS;
      MPI_Status status;
//       std::cout << "stage = " << stagename(stages_[stage_])
//                 << " batch = " << batch_
//                 << std::endl;
      switch (stages_[stage_]) {
      case Done:
          break;
      case NoOp:
          stage_++;
          break;
      case RvUp:
//           std::cout << "irecv(" << n_ << "," << up_ << ")" << std::endl;
          rc = MPI_Irecv(scratch_data_, n_, MPI_DOUBLE, up_, batch_, comm_, &req_);
          stage_++;
          break;
      case RvDn:
//           std::cout << "irecv(" << n_ << "," << dn_ << ")" << std::endl;
          rc = MPI_Irecv(data_, n_, MPI_DOUBLE, dn_, batch_, comm_, &req_);
          stage_++;
          break;
      case SdUp:
//           std::cout << "isend(" << n_ << "," << up_ << ")" << std::endl;
          rc = MPI_Isend(data_, n_, MPI_DOUBLE, up_, batch_, comm_, &req_);
          stage_++;
          break;
      case SdDn:
//           std::cout << "isend(" << n_ << "," << dn_ << ")" << std::endl;
          rc = MPI_Isend(data_, n_, MPI_DOUBLE, dn_, batch_, comm_, &req_);
          stage_++;
          break;
      case Wait:
          failed = false;
          rc = MPI_Wait(&req_,&status);
          if (stages_[stage_-1] == RvUp || stages_[stage_-1] == RvDn) {
              int count, crc;
              crc = MPI_Get_count(&status,MPI_DOUBLE,&count);
              if (crc != MPI_SUCCESS || count!= n_) {
                  std::cerr << "AllSumPipeline: failed: bad count" << std::endl;
                  std::cerr << "count   = " << count << std::endl;
                  std::cerr << "n       = " << n_ << std::endl;
                  std::cerr << "countrc = " << crc << std::endl;
                  MPI_Abort(MPI_COMM_WORLD,1);
                }
            }
          stage_++;
          break;
      case Test:
          failed = false;
          int flag;
          rc = MPI_Test(&req_,&flag,&status);
          if (flag) stage_++;
          break;
      case DSum:
          for (int i=0; i<n_; i++) {
              data_[i] += scratch_data_[i];
            }
          stage_++;
          break;
      default:
          break;
        }
      if (rc != MPI_SUCCESS) {
          std::cerr << "AllSumPipeline: MPI operation failed" << std::endl;
          std::cerr << "rc     = " << rc << std::endl;
          std::cerr << "stage_ = " << stage_ << std::endl;
        }
    }
    bool done() { return stages_[stage_] == Done; }
};

AllSumPipeline::Stages AllSumPipeline::stages_one1[] = {AllSumPipeline::Done};
AllSumPipeline::Stages AllSumPipeline::stages_one2[] = {AllSumPipeline::Done};
AllSumPipeline::Stages AllSumPipeline::stages_top1[] = {AllSumPipeline::SdDn,
                                                        AllSumPipeline::Wait,
                                                        AllSumPipeline::Done};
AllSumPipeline::Stages AllSumPipeline::stages_top2[] = {AllSumPipeline::RvDn,
                                                        AllSumPipeline::Wait,
                                                        AllSumPipeline::Done};
AllSumPipeline::Stages AllSumPipeline::stages_bot1[] = {AllSumPipeline::RvUp,
                                                        AllSumPipeline::Wait,
                                                        AllSumPipeline::DSum,
                                                        AllSumPipeline::Done};
AllSumPipeline::Stages AllSumPipeline::stages_bot2[] = {AllSumPipeline::SdUp,
                                                        AllSumPipeline::Wait,
                                                        AllSumPipeline::Done};
AllSumPipeline::Stages AllSumPipeline::stages_mid1[] = {AllSumPipeline::RvUp,
                                                        AllSumPipeline::Wait,
                                                        AllSumPipeline::DSum,
                                                        AllSumPipeline::SdDn,
                                                        AllSumPipeline::Wait,
                                                        AllSumPipeline::Done};
AllSumPipeline::Stages AllSumPipeline::stages_mid2[] = {AllSumPipeline::RvDn,
                                                        AllSumPipeline::Wait,
                                                        AllSumPipeline::SdUp,
                                                        AllSumPipeline::Wait,
                                                        AllSumPipeline::Done};
const char *AllSumPipeline::stagenames[] =  { "NoOp",
                                              "RvUp",
                                              "RvDn",
                                              "SdUp",
                                              "SdDn",
                                              "Wait",
                                              "Test",
                                              "DSum",
                                              "Done" };

void
step_pipes(std::list<AllSumPipeline>&pipes, int max_steppers = -1)
{
  if (pipes.front().done()) {
      pipes.pop_front();
    }
  int ipipe = 0;
  for (std::list<AllSumPipeline>::iterator i = pipes.begin();
      i != pipes.end() && ((max_steppers==-1)?true:(ipipe < max_steppers));
      i++, ipipe++) {
      i->next();
    }
}

void
allsum(double *data, long n)
{
  int me, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  int dn = ((me == 0)? -1: me-1);
  int up = ((me == nproc-1)? -1: me+1);

  std::list<AllSumPipeline> pipes;

  const int max_chunk = 3000;

  long nchunk = n/max_chunk;
  if (n%max_chunk) nchunk++;

  const bool always_dup = false;
  static MPI_Comm comm;
  static bool inited = false;

  if (always_dup || !inited) {
      MPI_Comm_dup(MPI_COMM_WORLD, &comm);
      inited = true;
    }

  // i goes over the 1st half: the reduce.
  // j goes over the 2nd half: the bcast.

  long jbatch=0, j=0; 

  for (long ibatch=0, i=0; i<n; ibatch++,i+=max_chunk) {
      long nichunk = ((n-i>max_chunk)?max_chunk:(n-i));

      // add the current chunk to the list of pipes
      pipes.push_back(AllSumPipeline(&data[i], nichunk, up, dn, ibatch, comm, 1));

      if (0 && j<n) {
          long njchunk = ((n-j>max_chunk)?max_chunk:(n-j));
          pipes.push_back(AllSumPipeline(&data[j], njchunk, up, dn, jbatch, comm, 2));
          jbatch++;
          j+=max_chunk;
        }

      // Advance through the pipes.
      step_pipes(pipes);
    }

  // finish outstanding pipeline steps
  while (pipes.begin() != pipes.end()) {
      if (0 && j<n) {
          long njchunk = ((n-j>max_chunk)?max_chunk:(n-j));
          pipes.push_back(AllSumPipeline(&data[j], njchunk, up, dn, jbatch, comm, 2));
          jbatch++;
          j+=max_chunk;
        }
      step_pipes(pipes);
    }

  for (; j<n; jbatch++,j+=max_chunk) {
      long njchunk = ((n-j>max_chunk)?max_chunk:(n-j));

      // add the current chunk to the list of pipes
      pipes.push_back(AllSumPipeline(&data[j], njchunk, up, dn, jbatch, comm, 2));

      // Advance through the pipes.
      step_pipes(pipes);
    }

  // finish outstanding pipeline steps
  while (pipes.begin() != pipes.end()) {
      step_pipes(pipes);
    }

  if (always_dup) MPI_Comm_free(&comm);

}
#endif

}
