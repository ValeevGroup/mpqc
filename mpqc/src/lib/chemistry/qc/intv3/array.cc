//
// array.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdlib.h>
#include <util/misc/exenv.h>
#include <chemistry/qc/intv3/array.h>

using namespace std;
using namespace sc;

static void
no_storage(const char *msg)
{
  ExEnv::errn() << msg << ": ran out of memory" << endl;
  abort();
}

////////////////////////////////////////////////////////////////////////////

IntV3Arraydouble2::IntV3Arraydouble2()
{
  n1_ = n2_ = 0;
  data_ = 0;
}

IntV3Arraydouble2::~IntV3Arraydouble2()
{
  for (int i=0; i<n1_; i++) {
      delete[] data_[i];
    }
  delete[] data_;
}

void
IntV3Arraydouble2::set_dim(int n1, int n2)
{
  n1_ = n1;
  n2_ = n2;
  data_ = new double*[n1_];
  if (data_ == 0) no_storage("IntV3Arraydouble2");
  for (int i=0; i<n1_; i++) {
      data_[i] = new double[n2_];
      if (data_[i] == 0) no_storage("IntV3Arraydouble2");
    }
}

void
IntV3Arraydouble2::print(ostream &o)
{
  for (int i=0; i<n1_; i++) {
      o << "i = " << i << endl;
      for (int j=0; j<n2_; j++) {
          o << data_[i][j];
        }
      o << endl;
    }
}

int
IntV3Arraydouble2::nbyte() const
{
  return n1_ * (sizeof(double*) + n2_ * sizeof(double));
}

////////////////////////////////////////////////////////////////////////////

IntV3Arraydouble3::IntV3Arraydouble3()
{
  n1_ = n2_ = n3_ = 0;
  data_ = 0;
}

IntV3Arraydouble3::~IntV3Arraydouble3()
{
  for (int i=0; i<n1_; i++) {
      for (int j=0; j<n2_; j++) {
          delete[] data_[i][j];
        }
      delete[] data_[i];
    }
  delete[] data_;
}

void
IntV3Arraydouble3::set_dim(int n1, int n2, int n3)
{
  n1_ = n1;
  n2_ = n2;
  n3_ = n3;
  data_ = new double**[n1_];
  if (data_ == 0) no_storage("IntV3Arraydouble3");
  for (int i=0; i<n1_; i++) {
      data_[i] = new double*[n2_];
      if (data_[i] == 0) no_storage("IntV3Arraydouble3");
      for (int j=0; j<n2_; j++) {
          data_[i][j] = new double[n3_];
          if (data_[i][j] == 0) no_storage("IntV3Arraydouble3");
        }
    }
}

void
IntV3Arraydouble3::print(ostream &o)
{
  for (int i=0; i<n1_; i++) {
      for (int j=0; j<n2_; j++) {
          o << "i, j = " << i << j << endl;
          for (int k=0; k<n3_; k++) {
              o << data_[i][j][k];
            }
        }
      o << endl;
    }
}

int
IntV3Arraydouble3::nbyte() const
{
  return n1_*(sizeof(double**) + n2_*(sizeof(double*) + n3_*sizeof(double)));
}

////////////////////////////////////////////////////////////////////////////

IntV3Arraydoublep2::IntV3Arraydoublep2()
{
  n1_ = n2_ = 0;
  data_ = 0;
}

IntV3Arraydoublep2::~IntV3Arraydoublep2()
{
  for (int i=0; i<n1_; i++) {
      delete[] data_[i];
    }
  delete[] data_;
}

void
IntV3Arraydoublep2::set_dim(int n1, int n2)
{
  n1_ = n1;
  n2_ = n2;
  data_ = new double**[n1_];
  if (data_ == 0) no_storage("IntV3Arraydoublep2");
  for (int i=0; i<n1_; i++) {
      data_[i] = new double*[n2_];
      if (data_[i] == 0) no_storage("IntV3Arraydoublep2");
    }
}

void
IntV3Arraydoublep2::print(ostream &o)
{
  for (int i=0; i<n1_; i++) {
      o << "i = " << i << endl;
      for (int j=0; j<n2_; j++) {
          o << data_[i][j];
        }
      o << endl;
    }
}

int
IntV3Arraydoublep2::nbyte() const
{
  return n1_*(sizeof(double**)
              + n2_*(sizeof(double*)));
}

////////////////////////////////////////////////////////////////////////////

IntV3Arraydoublep3::IntV3Arraydoublep3()
{
  n1_ = n2_ = n3_ = 0;
  data_ = 0;
}

IntV3Arraydoublep3::~IntV3Arraydoublep3()
{
  for (int i=0; i<n1_; i++) {
      for (int j=0; j<n2_; j++) {
          delete[] data_[i][j];
        }
      delete[] data_[i];
    }
  delete[] data_;
}

void
IntV3Arraydoublep3::set_dim(int n1, int n2, int n3)
{
  n1_ = n1;
  n2_ = n2;
  n3_ = n3;
  data_ = new double***[n1_];
  if (data_ == 0) no_storage("IntV3Arraydoublep3");
  for (int i=0; i<n1_; i++) {
      data_[i] = new double**[n2_];
      if (data_[i] == 0) no_storage("IntV3Arraydoublep3");
      for (int j=0; j<n2_; j++) {
          data_[i][j] = new double*[n3_];
          if (data_[i][j] == 0) no_storage("IntV3Arraydoublep3");
        }
    }
}

void
IntV3Arraydoublep3::delete_data()
{
  for (int i=0; i<n1_; i++) {
      double ***datai = data_[i];
      for (int j=0; j<n2_; j++) {
          double **dataj = datai[j];
          for (int k=0; k<n3_; k++) {
              delete[] dataj[k];
            }
        }
    }
}

void
IntV3Arraydoublep3::print(ostream &o)
{
  for (int i=0; i<n1_; i++) {
      for (int j=0; j<n2_; j++) {
          o << "i, j = " << i << j << endl;
          for (int k=0; k<n3_; k++) {
              o << data_[i][j][k];
            }
          o << endl;
        }
    }
}

int
IntV3Arraydoublep3::nbyte() const
{
  return n1_*(sizeof(double***)
              + n2_*(sizeof(double**)
                     + n3_*sizeof(double*)));
}

////////////////////////////////////////////////////////////////////////////

IntV3Arraydoublep4::IntV3Arraydoublep4()
{
  n1_ = n2_ = n3_ = n4_ = 0;
  data_ = 0;
}

IntV3Arraydoublep4::~IntV3Arraydoublep4()
{
  for (int i=0; i<n1_; i++) {
      for (int j=0; j<n2_; j++) {
          for (int k=0; k<n3_; k++) {
              delete[] data_[i][j][k];
            }
          delete[] data_[i][j];
        }
      delete[] data_[i];
    }
  delete[] data_;
  data_=0;
}

void
IntV3Arraydoublep4::set_dim(int n1, int n2, int n3, int n4)
{
  if (data_) {
      for (int i=0; i<n1_; i++) {
          for (int j=0; j<n2_; j++) {
              for (int k=0; k<n3_; k++) {
                  delete[] data_[i][j][k];
                }
              delete[] data_[i][j];
            }
          delete[] data_[i];
        }
      delete[] data_;
    }
  n1_ = n1;
  n2_ = n2;
  n3_ = n3;
  n4_ = n4;
  data_ = new double****[n1_];
  if (data_ == 0) no_storage("IntV3Arraydoublep4");
  for (int i=0; i<n1_; i++) {
      data_[i] = new double***[n2_];
      if (data_[i] == 0) no_storage("IntV3Arraydoublep4");
      for (int j=0; j<n2_; j++) {
          data_[i][j] = new double**[n3_];
          if (data_[i][j] == 0) no_storage("IntV3Arraydoublep4");
          for (int k=0; k<n3_ ;k++) {
              data_[i][j][k] = new double*[n4_];
              if (data_[i][j][k] == 0) no_storage("IntV3Arraydoublep4");
            }
        }
    }
}

void
IntV3Arraydoublep4::print(ostream &o)
{
  for (int i=0; i<n1_; i++) {
      for (int j=0; j<n2_; j++) {
          for (int k=0; k<n3_; k++) {
              o << "i, j, k = " << i << j << k << endl;
              for (int l=0; l<n4_; l++) {
                  o << data_[i][j][k][l];
                }
              o << endl;
            }
        }
    }
}

int
IntV3Arraydoublep4::nbyte() const
{
  return n1_*(sizeof(double****)
              + n2_*(sizeof(double***)
                     + n3_*(sizeof(double**)
                            + n4_*sizeof(double*))));
}

////////////////////////////////////////////////////////////////////////////

IntV3Arrayint3::IntV3Arrayint3()
{
  n1_ = n2_ = n3_ = 0;
  data_ = 0;
}

IntV3Arrayint3::~IntV3Arrayint3()
{
  for (int i=0; i<n1_; i++) {
      for (int j=0; j<n2_; j++) {
          delete[] data_[i][j];
        }
      delete[] data_[i];
    }
  delete[] data_;
}

void
IntV3Arrayint3::set_dim(int n1, int n2, int n3)
{
  n1_ = n1;
  n2_ = n2;
  n3_ = n3;
  data_ = new int**[n1_];
  if (data_ == 0) no_storage("IntV3Arrayint3");
  for (int i=0; i<n1_; i++) {
      data_[i] = new int*[n2_];
      if (data_[i] == 0) no_storage("IntV3Arrayint3");
      for (int j=0; j<n2_; j++) {
          data_[i][j] = new int[n3_];
          if (data_[i][j] == 0) no_storage("IntV3Arrayint3");
        }
    }
}

void
IntV3Arrayint3::print(ostream &o)
{
  for (int i=0; i<n1_; i++) {
      for (int j=0; j<n2_; j++) {
          o << "i, j = " << i << j << endl;
          for (int k=0; k<n3_; k++) {
              o << data_[i][j][k];
            }
          o << endl;
        }
    }
}

int
IntV3Arrayint3::nbyte() const
{
  return n1_*(sizeof(int**)
              + n2_*(sizeof(int*)
                     + n3_*sizeof(int)));
}

////////////////////////////////////////////////////////////////////////////

IntV3Arrayint4::IntV3Arrayint4()
{
  n1_ = n2_ = n3_ = n4_ = 0;
  data_ = 0;
}

IntV3Arrayint4::~IntV3Arrayint4()
{
  for (int i=0; i<n1_; i++) {
      for (int j=0; j<n2_; j++) {
          for (int k=0; k<n3_; k++) {
              delete[] data_[i][j][k];
            }
          delete[] data_[i][j];
        }
      delete[] data_[i];
    }
  delete[] data_;
}

void
IntV3Arrayint4::set_dim(int n1, int n2, int n3, int n4)
{
  n1_ = n1;
  n2_ = n2;
  n3_ = n3;
  n4_ = n4;
  data_ = new int***[n1_];
  if (data_ == 0) no_storage("IntV3Arrayint4");
  for (int i=0; i<n1_; i++) {
      data_[i] = new int**[n2_];
      if (data_[i] == 0) no_storage("IntV3Arrayint4");
      for (int j=0; j<n2_; j++) {
          data_[i][j] = new int*[n3_];
          if (data_[i][j] == 0) no_storage("IntV3Arrayint4");
          for (int k=0; k<n3_ ;k++) {
              data_[i][j][k] = new int[n4_];
              if (data_[i][j][k] == 0) no_storage("IntV3Arrayint4");
            }
        }
    }
}

void
IntV3Arrayint4::print(ostream &o)
{
  for (int i=0; i<n1_; i++) {
      for (int j=0; j<n2_; j++) {
          for (int k=0; k<n3_; k++) {
              o << "i, j, k = " << i << j << k << endl;
              for (int l=0; l<n4_; l++) {
                  o << data_[i][j][k][l];
                }
              o << endl;
            }
        }
    }
}

int
IntV3Arrayint4::nbyte() const
{
  return n1_*(sizeof(int***)
              + n2_*(sizeof(int**)
                     + n3_*(sizeof(int*)
                            + n4_*sizeof(int))));
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
