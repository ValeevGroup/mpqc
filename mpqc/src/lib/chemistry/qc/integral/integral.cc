
#include <stdio.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/integral/integral.h>
#include <math/scmat/block.h>

///////////////////////////////////////////////////////////////////////

OneBodyInt::OneBodyInt(const RefGaussianBasisSet&b):
  bs(b)
{
  // allocate a buffer
  int biggest_shell = b->max_nfunction_in_shell();
  if (biggest_shell) {
      buffer_ = new double[biggest_shell * biggest_shell];
    }
  else {
      buffer_ = 0;
    }
}

OneBodyInt::~OneBodyInt()
{
  if (buffer_) delete[] buffer_;
}

int OneBodyInt::nbasis()
{
  return bs->nbasis();
}

int OneBodyInt::nshell()
{
  return bs->nshell();
}

void
OneBodyInt::process(SCMatrixBlockIter&)
{
  fprintf(stderr,"OneBodyInt::process(SCMatrixBlockIter&):"
          " cannot handle generic case\n");
  abort();
}

void
OneBodyInt::process(SCMatrixRectBlock* b)
{
  // convert basis function indices into shell indices
  int ishstart = bs->function_to_shell(b->istart);
  int jshstart = bs->function_to_shell(b->jstart);
  int ishend = bs->function_to_shell(b->iend);
  int jshend = bs->function_to_shell(b->jend);

  int njdata = b->jend - b->jstart;

  for (int ish=ishstart;ish<ishend;ish++) {
      int nish = bs->operator[](ish).nfunction();
      for (int jsh=jshstart;jsh<jshend;jsh++) {
          int njsh = bs->operator[](jsh).nfunction();

          // compute a set of shell integrals
          compute_shell(ish,jsh,buffer_);

          double*tmp = buffer_;
          int ifn = bs->shell_to_function(ish);
          int jfnsave = bs->shell_to_function(jsh);
          for (int i=0; i<nish; i++,ifn++) {
              if (ifn < b->istart || ifn >= b->iend) {
                  tmp += njsh;
                }
              else {
                  int jfn = jfnsave;
                  int data_index = (ifn - b->istart)*njdata + jfn - b->jstart;
                  for (int j=0; j<njsh; j++,jfn++) {
                      if (jfn >= b->jstart && jfn < b->jend) {
                          b->data[data_index] = *tmp;
                          data_index++;
                        }
                      tmp++;
                    }
                }
            }
          
        }
    }
}

void
OneBodyInt::process(SCMatrixLTriBlock* b)
{
  // convert basis function indices into shell indices
  int fnstart = b->start;
  int fnend = b->end;
  int shstart = bs->function_to_shell(fnstart);
  int shend = (fnend?bs->function_to_shell(fnend - 1) + 1 : 0);

  // loop over all needed shells
  for (int ish=shstart;ish<shend;ish++) {
      int nish = bs->operator[](ish).nfunction();
      for (int jsh=shstart;jsh<=ish;jsh++) {
          int njsh = bs->operator[](jsh).nfunction();

          // compute a set of shell integrals
          compute_shell(ish,jsh,buffer_);

          // take the integrals from buffer and put them into the LTri block
          double*tmp = buffer_;
          int ifn = bs->shell_to_function(ish);
          int jfnsave = bs->shell_to_function(jsh);
          for (int i=0; i<nish; i++,ifn++) {
              // skip over basis functions that are not needed
              if (ifn < fnstart || ifn >= fnend) {
                  tmp += njsh;
                }
              else {
                  int jfn = jfnsave;
                  int irelfn = ifn - fnstart;
                  int data_index = ((irelfn+1)*irelfn>>1) + jfn - fnstart;
                  for (int j=0; j<njsh; j++,jfn++) {
                      // skip over basis functions that are not needed
                      if (jfn <= ifn && jfn >= fnstart) {
                          b->data[data_index] = *tmp;
                          data_index++;
                        }
                      tmp++;
                    }
                }
            }
          
        }
    }
}

///////////////////////////////////////////////////////////////////////

OneBody3Int::OneBody3Int(const RefGaussianBasisSet&b):
  bs(b)
{
  // allocate a buffer
  int biggest_shell = b->max_nfunction_in_shell();
  if (biggest_shell) {
      buffer_ = new double[biggest_shell * biggest_shell * 3];
    }
  else {
      buffer_ = 0;
    }
}

OneBody3Int::~OneBody3Int()
{
  if (buffer_) delete[] buffer_;
}

int OneBody3Int::nbasis()
{
  return bs->nbasis();
}

int OneBody3Int::nshell()
{
  return bs->nshell();
}

void
OneBody3Int::process(SCMatrixBlockIter&,
                     SCMatrixBlockIter&,
                     SCMatrixBlockIter&)
{
  fprintf(stderr,"OneBody3Int::process(SCMatrixBlockIter&):"
          " cannot handle generic case\n");
  abort();
}

void
OneBody3Int::process(SCMatrixRectBlock* a,
                     SCMatrixRectBlock* b,
                     SCMatrixRectBlock* c)
{
  // convert basis function indices into shell indices
  int ishstart = bs->function_to_shell(b->istart);
  int jshstart = bs->function_to_shell(b->jstart);
  int ishend = bs->function_to_shell(b->iend);
  int jshend = bs->function_to_shell(b->jend);

  int njdata = b->jend - b->jstart;

  for (int ish=ishstart;ish<ishend;ish++) {
      int nish = bs->operator[](ish).nfunction();
      for (int jsh=jshstart;jsh<jshend;jsh++) {
          int njsh = bs->operator[](jsh).nfunction();

          // compute a set of shell integrals
          compute_shell(ish,jsh,buffer_);

          double*tmp = buffer_;
          int ifn = bs->shell_to_function(ish);
          int jfnsave = bs->shell_to_function(jsh);
          for (int i=0; i<nish; i++,ifn++) {
              if (ifn < b->istart || ifn >= b->iend) {
                  tmp += njsh * 3;
                }
              else {
                  int jfn = jfnsave;
                  int data_index = (ifn - b->istart)*njdata + jfn - b->jstart;
                  for (int j=0; j<njsh; j++,jfn++) {
                      if (jfn >= b->jstart && jfn < b->jend) {
                          a->data[data_index] = tmp[0];
                          b->data[data_index] = tmp[1];
                          b->data[data_index] = tmp[2];
                          data_index++;
                        }
                      tmp += 3;
                    }
                }
            }
          
        }
    }
}

void
OneBody3Int::process(SCMatrixLTriBlock* a,
                     SCMatrixLTriBlock* b,
                     SCMatrixLTriBlock* c)
{
  // convert basis function indices into shell indices
  int fnstart = b->start;
  int fnend = b->end;
  int shstart = bs->function_to_shell(fnstart);
  int shend = (fnend?bs->function_to_shell(fnend - 1) + 1 : 0);

  // loop over all needed shells
  for (int ish=shstart;ish<shend;ish++) {
      int nish = bs->operator[](ish).nfunction();
      for (int jsh=shstart;jsh<=ish;jsh++) {
          int njsh = bs->operator[](jsh).nfunction();

          // compute a set of shell integrals
          compute_shell(ish,jsh,buffer_);

          // take the integrals from buffer and put them into the LTri block
          double*tmp = buffer_;
          int ifn = bs->shell_to_function(ish);
          int jfnsave = bs->shell_to_function(jsh);
          for (int i=0; i<nish; i++,ifn++) {
              // skip over basis functions that are not needed
              if (ifn < fnstart || ifn >= fnend) {
                  tmp += njsh * 3;
                }
              else {
                  int jfn = jfnsave;
                  int irelfn = ifn - fnstart;
                  int data_index = ((irelfn+1)*irelfn>>1) + jfn - fnstart;
                  for (int j=0; j<njsh; j++,jfn++) {
                      // skip over basis functions that are not needed
                      if (jfn <= ifn && jfn >= fnstart) {
                          a->data[data_index] = tmp[0];
                          b->data[data_index] = tmp[1];
                          c->data[data_index] = tmp[2];
                          data_index++;
                        }
                      tmp += 3;
                    }
                }
            }
          
        }
    }
}
