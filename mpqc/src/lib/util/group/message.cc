
#ifndef _util_group_message_cc
#define _util_group_message_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/group/message.h>

template <class T>
void
GrpSumReduce<T>::reduce(T*target, T*data, int nelement)
{
  for (int i=0; i<nelement; i++) {
      target[i] += data[i];
    }
}

template <class T>
void
GrpMinReduce<T>::reduce(T*target, T*data, int nelement)
{
  for (int i=0; i<nelement; i++) {
      if (target[i] > data[i]) target[i] = data[i];
    }
}

template <class T>
void
GrpMaxReduce<T>::reduce(T*target, T*data, int nelement)
{
  for (int i=0; i<nelement; i++) {
      if (target[i] < data[i]) target[i] = data[i];
    }
}

template <class T>
void
GrpArithmeticAndReduce<T>::reduce(T*target, T*data, int nelement)
{
  for (int i=0; i<nelement; i++) {
      target[i] = target[i] & data[i];
    }
}

template <class T>
void
GrpArithmeticOrReduce<T>::reduce(T*target, T*data, int nelement)
{
  for (int i=0; i<nelement; i++) {
      target[i] = target[i] | data[i];
    }
}

template <class T>
void
GrpArithmeticXOrReduce<T>::reduce(T*target, T*data, int nelement)
{
  for (int i=0; i<nelement; i++) {
      target[i] = target[i] ^ data[i];
    }
}

template <class T>
void
GrpProductReduce<T>::reduce(T*target, T*data, int nelement)
{
  for (int i=0; i<nelement; i++) {
      target[i] *= data[i];
    }
}

template <class T>
void
GrpFunctionReduce<T>::reduce(T*target, T*data, int nelement)
{
  (*func_)(target,data,nelement);
}

#endif
