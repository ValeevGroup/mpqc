
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdlib.h>
#include <chemistry/qc/intv3/array.h>

static void
no_storage(const char *msg)
{
  cerr << msg << ": ran out of memory" << endl;
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
}

void
IntV3Arraydoublep4::set_dim(int n1, int n2, int n3, int n4)
{
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
