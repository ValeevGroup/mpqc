
#include <stdio.h>
#include <chemistry/qc/basis/basis.h>
#include "integral.h"
#include <math/scmat/block.h>

///////////////////////////////////////////////////////////////////////

OneBodyInt::OneBodyInt(const GaussianBasisSet*b):
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

const GaussianBasisSet& OneBodyInt::basis()
{
  return *bs;
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
  int shend = bs->function_to_shell(fnend);

  // loop over all needed shells
  int ioff = 0;
  for (int ish=shstart;ish<shend;ish++) {
      int nish = bs->operator[](ish).nfunction();
      for (int jsh=shstart;jsh<=shend;jsh++) {
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
                  int data_index = ((irelfn+1)*irelfn)>>1 + jfn - fnstart;
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
