
#include <stdio.h>
#include <util/group/message.h>
#include <util/group/message.cc>

/////////////////////////////////////////////////////////////////////////
// instantiate templates

#ifdef __GNUG__
template class GrpReduce<double>;
template class GrpReduce<int>;
template class GrpReduce<long>;
template class GrpReduce<float>;
template class GrpReduce<short>;
template class GrpReduce<char>;
template class GrpReduce<unsigned char>;

template class GrpFunctionReduce<double>;
template class GrpFunctionReduce<int>;
template class GrpFunctionReduce<long>;
template class GrpFunctionReduce<float>;
template class GrpFunctionReduce<short>;
template class GrpFunctionReduce<char>;
template class GrpFunctionReduce<unsigned char>;

template class GrpMinReduce<double>;
template class GrpMinReduce<int>;
template class GrpMinReduce<long>;
template class GrpMinReduce<float>;
template class GrpMinReduce<short>;
template class GrpMinReduce<char>;
template class GrpMinReduce<unsigned char>;

template class GrpMaxReduce<double>;
template class GrpMaxReduce<int>;
template class GrpMaxReduce<long>;
template class GrpMaxReduce<float>;
template class GrpMaxReduce<short>;
template class GrpMaxReduce<char>;
template class GrpMaxReduce<unsigned char>;

template class GrpSumReduce<double>;
template class GrpSumReduce<int>;
template class GrpSumReduce<long>;
template class GrpSumReduce<float>;
template class GrpSumReduce<short>;
template class GrpSumReduce<char>;
template class GrpSumReduce<unsigned char>;

template class GrpProductReduce<double>;
template class GrpProductReduce<int>;
template class GrpProductReduce<long>;
template class GrpProductReduce<float>;
template class GrpProductReduce<short>;
template class GrpProductReduce<char>;
template class GrpProductReduce<unsigned char>;

template class GrpArithmeticOrReduce<int>;
template class GrpArithmeticOrReduce<long>;
template class GrpArithmeticOrReduce<short>;
template class GrpArithmeticOrReduce<char>;
template class GrpArithmeticOrReduce<unsigned char>;

template class GrpArithmeticAndReduce<int>;
template class GrpArithmeticAndReduce<long>;
template class GrpArithmeticAndReduce<short>;
template class GrpArithmeticAndReduce<char>;
template class GrpArithmeticAndReduce<unsigned char>;

template class GrpArithmeticXOrReduce<int>;
template class GrpArithmeticXOrReduce<long>;
template class GrpArithmeticXOrReduce<short>;
template class GrpArithmeticXOrReduce<char>;
template class GrpArithmeticXOrReduce<unsigned char>;
#endif

/////////////////////////////////////////////////////////////////////////
// sum reduction members

template <class T>
void
do_sum(MessageGrp* grp, T* data, int n, T* tmp, int target)
{
  GrpSumReduce<T> gred;
  grp->reduce(data, n, gred, tmp, target);
}

void
MessageGrp::sum(double* data, int n, double* tmp, int target)
{
  do_sum(this, data, n, tmp, target);
}

void
MessageGrp::sum(int* data, int n, int* tmp, int target)
{
  do_sum(this, data, n, tmp, target);
}

void
MessageGrp::sum(char* data, int n, char* tmp, int target)
{
  do_sum(this, data, n, tmp, target);
}

void
MessageGrp::sum(unsigned char* data, int n, unsigned char* tmp, int target)
{
  do_sum(this, data, n, tmp, target);
}

/////////////////////////////////////////////////////////////////////////
// min reduction members

template <class T>
void
do_max(MessageGrp* grp, T* data, int n, T* tmp, int target)
{
  GrpMaxReduce<T> gred;
  grp->reduce(data, n, gred, tmp, target);
}

void
MessageGrp::max(double* data, int n, double* tmp, int target)
{
  do_max(this, data, n, tmp, target);
}

void
MessageGrp::max(int* data, int n, int* tmp, int target)
{
  do_max(this, data, n, tmp, target);
}

void
MessageGrp::max(char* data, int n, char* tmp, int target)
{
  do_max(this, data, n, tmp, target);
}

void
MessageGrp::max(unsigned char* data, int n, unsigned char* tmp, int target)
{
  do_max(this, data, n, tmp, target);
}

/////////////////////////////////////////////////////////////////////////
// max reduction members

template <class T>
void
do_min(MessageGrp* grp, T* data, int n, T* tmp, int target)
{
  GrpMinReduce<T> gred;
  grp->reduce(data, n, gred, tmp, target);
}

void
MessageGrp::min(double* data, int n, double* tmp, int target)
{
  do_min(this, data, n, tmp, target);
}

void
MessageGrp::min(int* data, int n, int* tmp, int target)
{
  do_min(this, data, n, tmp, target);
}

void
MessageGrp::min(char* data, int n, char* tmp, int target)
{
  do_min(this, data, n, tmp, target);
}

void
MessageGrp::min(unsigned char* data, int n, unsigned char* tmp, int target)
{
  do_min(this, data, n, tmp, target);
}

/////////////////////////////////////////////////////////////////////////
// generic reduction

void
MessageGrp::reduce(double* data, int n, GrpReduce<double>& red,
                   double* scratch, int target)
{
  int tgop_max = gop_max_/sizeof(double);
  if (tgop_max == 0) tgop_max = gop_max_?1:n;

  int passed_scratch;
  if (!scratch) {
      scratch = new double[n>tgop_max?tgop_max:n];
      passed_scratch = 0;
    }
  else passed_scratch = 1;

  for (int idat=0; idat<n; idat+=tgop_max) {
      int ndat = (idat+tgop_max>n)?(n-idat):tgop_max;
      RefGlobalMsgIter i(topology_->global_msg_iter(this,
                                                    (target== -1?0:target)));
      for (i->backwards(); !i->done(); i->next()) {
          if (i->send()) {
              send(i->sendto(), &data[idat], ndat);
            }
          if (i->recv()) {
              recv(i->recvfrom(), scratch, ndat);
              red.reduce(&data[idat], scratch, ndat);
            }
        }
    }

  if (target == -1) {
      bcast(data, n, 0);
    }

  if (!passed_scratch) delete[] scratch;
}

void
MessageGrp::reduce(int* data, int n, GrpReduce<int>& red,
                   int* scratch, int target)
{
  int tgop_max = gop_max_/sizeof(int);
  if (tgop_max == 0) tgop_max = gop_max_?1:n;

  int passed_scratch;
  if (!scratch) {
      scratch = new int[n>tgop_max?tgop_max:n];
      passed_scratch = 0;
    }
  else passed_scratch = 1;

  for (int idat=0; idat<n; idat+=tgop_max) {
      int ndat = (idat+tgop_max>n)?(n-idat):tgop_max;
      RefGlobalMsgIter i(topology_->global_msg_iter(this,
                                                    (target== -1?0:target)));
      for (i->backwards(); !i->done(); i->next()) {
          if (i->send()) {
              send(i->sendto(), &data[idat], ndat);
            }
          if (i->recv()) {
              recv(i->recvfrom(), scratch, ndat);
              red.reduce(&data[idat], scratch, ndat);
            }
        }
    }

  if (target == -1) {
      bcast(data, n, 0);
    }

  if (!passed_scratch) delete[] scratch;
}

void
MessageGrp::reduce(char* data, int n, GrpReduce<char>& red,
                   char* scratch, int target)
{
  int tgop_max = gop_max_/sizeof(char);
  if (tgop_max == 0) tgop_max = gop_max_?1:n;

  int passed_scratch;
  if (!scratch) {
      scratch = new char[n>tgop_max?tgop_max:n];
      passed_scratch = 0;
    }
  else passed_scratch = 1;

  for (int idat=0; idat<n; idat+=tgop_max) {
      int ndat = (idat+tgop_max>n)?(n-idat):tgop_max;
      RefGlobalMsgIter i(topology_->global_msg_iter(this,
                                                    (target== -1?0:target)));
      for (i->backwards(); !i->done(); i->next()) {
          if (i->send()) {
              send(i->sendto(), &data[idat], ndat);
            }
          if (i->recv()) {
              recv(i->recvfrom(), scratch, ndat);
              red.reduce(&data[idat], scratch, ndat);
            }
        }
    }

  if (target == -1) {
      bcast(data, n, 0);
    }

  if (!passed_scratch) delete[] scratch;
}

void
MessageGrp::reduce(unsigned char* data, int n, GrpReduce<unsigned char>& red,
                   unsigned char* scratch, int target)
{
  int tgop_max = gop_max_/sizeof(unsigned char);
  if (tgop_max == 0) tgop_max = gop_max_?1:n;

  int passed_scratch;
  if (!scratch) {
      scratch = new unsigned char[n>tgop_max?tgop_max:n];
      passed_scratch = 0;
    }
  else passed_scratch = 1;

  for (int idat=0; idat<n; idat+=tgop_max) {
      int ndat = (idat+tgop_max>n)?(n-idat):tgop_max;
      RefGlobalMsgIter i(topology_->global_msg_iter(this,
                                                    (target== -1?0:target)));
      for (i->backwards(); !i->done(); i->next()) {
          if (i->send()) {
              send(i->sendto(), &data[idat], ndat);
            }
          if (i->recv()) {
              recv(i->recvfrom(), scratch, ndat);
              red.reduce(&data[idat], scratch, ndat);
            }
        }
    }

  if (target == -1) {
      bcast(data, n, 0);
    }

  if (!passed_scratch) delete[] scratch;
}

void
MessageGrp::reduce(short* data, int n, GrpReduce<short>& red,
                   short* scratch, int target)
{
  int tgop_max = gop_max_/sizeof(short);
  if (tgop_max == 0) tgop_max = gop_max_?1:n;

  int passed_scratch;
  if (!scratch) {
      scratch = new short[n>tgop_max?tgop_max:n];
      passed_scratch = 0;
    }
  else passed_scratch = 1;

  for (int idat=0; idat<n; idat+=tgop_max) {
      int ndat = (idat+tgop_max>n)?(n-idat):tgop_max;
      RefGlobalMsgIter i(topology_->global_msg_iter(this,
                                                    (target== -1?0:target)));
      for (i->backwards(); !i->done(); i->next()) {
          if (i->send()) {
              send(i->sendto(), &data[idat], ndat);
            }
          if (i->recv()) {
              recv(i->recvfrom(), scratch, ndat);
              red.reduce(&data[idat], scratch, ndat);
            }
        }
    }

  if (target == -1) {
      bcast(data, n, 0);
    }

  if (!passed_scratch) delete[] scratch;
}

void
MessageGrp::reduce(float* data, int n, GrpReduce<float>& red,
                   float* scratch, int target)
{
  int tgop_max = gop_max_/sizeof(float);
  if (tgop_max == 0) tgop_max = gop_max_?1:n;

  int passed_scratch;
  if (!scratch) {
      scratch = new float[n>tgop_max?tgop_max:n];
      passed_scratch = 0;
    }
  else passed_scratch = 1;

  for (int idat=0; idat<n; idat+=tgop_max) {
      int ndat = (idat+tgop_max>n)?(n-idat):tgop_max;
      RefGlobalMsgIter i(topology_->global_msg_iter(this,
                                                    (target== -1?0:target)));
      for (i->backwards(); !i->done(); i->next()) {
          if (i->send()) {
              send(i->sendto(), &data[idat], ndat);
            }
          if (i->recv()) {
              recv(i->recvfrom(), scratch, ndat);
              red.reduce(&data[idat], scratch, ndat);
            }
        }
    }

  if (target == -1) {
      bcast(data, n, 0);
    }

  if (!passed_scratch) delete[] scratch;
}

void
MessageGrp::reduce(long* data, int n, GrpReduce<long>& red,
                   long* scratch, int target)
{
  int tgop_max = gop_max_/sizeof(long);
  if (tgop_max == 0) tgop_max = gop_max_?1:n;

  int passed_scratch;
  if (!scratch) {
      scratch = new long[n>tgop_max?tgop_max:n];
      passed_scratch = 0;
    }
  else passed_scratch = 1;

  for (int idat=0; idat<n; idat+=tgop_max) {
      int ndat = (idat+tgop_max>n)?(n-idat):tgop_max;
      RefGlobalMsgIter i(topology_->global_msg_iter(this,
                                                    (target== -1?0:target)));
      for (i->backwards(); !i->done(); i->next()) {
          if (i->send()) {
              send(i->sendto(), &data[idat], ndat);
            }
          if (i->recv()) {
              recv(i->recvfrom(), scratch, ndat);
              red.reduce(&data[idat], scratch, ndat);
            }
        }
    }

  if (target == -1) {
      bcast(data, n, 0);
    }

  if (!passed_scratch) delete[] scratch;
}

#ifdef __GNUG__
#define INSTANTIATE_DO_X(func,type) \
    template void func(MessageGrp*, type *, int, type *, int)

INSTANTIATE_DO_X(do_sum,int);
INSTANTIATE_DO_X(do_sum,double);
INSTANTIATE_DO_X(do_sum,char);
INSTANTIATE_DO_X(do_sum,unsigned char);

INSTANTIATE_DO_X(do_max,int);
INSTANTIATE_DO_X(do_max,double);
INSTANTIATE_DO_X(do_max,char);
INSTANTIATE_DO_X(do_max,unsigned char);

INSTANTIATE_DO_X(do_min,int);
INSTANTIATE_DO_X(do_min,double);
INSTANTIATE_DO_X(do_min,char);
INSTANTIATE_DO_X(do_min,unsigned char);

#endif
