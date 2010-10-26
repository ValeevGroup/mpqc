//
// carray.h --- C Style Arrays
//

#ifndef _util_container_carray_h
#define _util_container_carray_h

namespace sc {
    
template <class T>
T **
new_c_array2(int l, int m, T)
{
  T *a = 0;
  T **b = 0;
  if (l*m) a = new T[l*m];
  if (l) b = new T*[l];
  for (int i=0; i<l; i++) b[i] = &a[i*m];
  return b;
}

template <class T>
T **
new_zero_c_array2(int l, int m, T)
{
  T *a = 0;
  T **b = 0;
  if (l*m) a = new T[l*m];
  if (l) b = new T*[l];
  for (int i=0; i<l; i++) {
      b[i] = &a[i*m];
      for (int j=0; j<m; j++) {
          b[i][j] = 0;
        }
    }
  return b;
}

template <class T>
void
delete_c_array2(T**b)
{
  if (b) delete[] b[0];
  delete[] b;
}

template <class T>
T ***
new_c_array3(int l, int m, int n, T)
{
  T *a = 0;
  T **b = 0;
  T ***c = 0;
  if (l*m*n) a = new T[l*m*n];
  if (l*m) b = new T*[l*m];
  if (l) c = new T**[l];
  for (int i=0,ij=0; i<l; i++) {
      c[i] = &b[i*m];
      for (int j=0; j<m; j++,ij++) {
          c[i][j] = &a[ij*n];
        }
    }
  return c;
}

template <class T>
T ***
new_zero_c_array3(int l, int m, int n, T)
{
  T *a = 0;
  T **b = 0;
  T ***c = 0;
  if (l*m*n) a = new T[l*m*n];
  if (l*m) b = new T*[l*m];
  if (l) c = new T**[l];
  for (int i=0,ij=0; i<l; i++) {
      c[i] = &b[i*m];
      for (int j=0; j<m; j++,ij++) {
          c[i][j] = &a[ij*n];
          for (int k=0; k<n; k++) {
              c[i][j][k] = 0;
            }
        }
    }
  return c;
}

template <class T>
void
delete_c_array3(T***b)
{
  if (b && b[0]) delete[] b[0][0];
  if (b) delete[] b[0];
  delete[] b;
}

}

#endif

// ///////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
