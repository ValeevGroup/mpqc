
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdio.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/integral/integral.h>
#include <math/scmat/block.h>

///////////////////////////////////////////////////////////////////////

OneBodyInt::OneBodyInt(const RefGaussianBasisSet&b, OneBodyIntIter *it):
  bs1(b), bs2(b), iter(it)
{
  // allocate a buffer
  int biggest_shell = b->max_nfunction_in_shell();
  if (biggest_shell) {
      buffer_ = new double[biggest_shell * biggest_shell];
    }
  else {
      buffer_ = 0;
    }

  if (!iter)
    iter = new OneBodyIntIter;
}

OneBodyInt::OneBodyInt(const RefGaussianBasisSet&b1,
                       const RefGaussianBasisSet&b2,
                       OneBodyIntIter *it) :
  bs1(b1), bs2(b2), iter(it)
{
  // allocate a buffer
  int biggest_shell = b1->max_nfunction_in_shell() *
                      b2->max_nfunction_in_shell();
  if (biggest_shell) {
      buffer_ = new double[biggest_shell];
    }
  else {
      buffer_ = 0;
    }

  if (!iter)
    iter = new OneBodyIntIter;
}

OneBodyInt::~OneBodyInt()
{
  if (buffer_) {
    delete[] buffer_;
    buffer_=0;
  }
  if (iter) {
    delete iter;
    iter=0;
  }
}

int OneBodyInt::nbasis()
{
  return bs1->nbasis();
}

int OneBodyInt::nbasis1()
{
  return bs1->nbasis();
}

int OneBodyInt::nbasis2()
{
  return bs2->nbasis();
}

int OneBodyInt::nshell()
{
  return bs1->nshell();
}

int OneBodyInt::nshell1()
{
  return bs1->nshell();
}

int OneBodyInt::nshell2()
{
  return bs2->nshell();
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
  int ishstart = bs1->function_to_shell(b->istart);
  int jshstart = bs2->function_to_shell(b->jstart);
  int b1end = b->iend;
  int ishend = (b1end?bs1->function_to_shell(b1end-1) + 1 : 0);
  int b2end = b->jend;
  int jshend = (b2end?bs2->function_to_shell(b2end-1) + 1 : 0);

  int njdata = b->jend - b->jstart;

  iter->reset(ishstart, ishend, jshstart, jshend);
  for (iter->start(); iter->ready(); iter->next()) {
      int ish=iter->ishell();
      int jsh=iter->jshell();

      int nish = bs1->operator[](ish).nfunction();
      int njsh = bs2->operator[](jsh).nfunction();

      double scale = iter->scale();
    
      // compute a set of shell integrals
      compute_shell(ish,jsh,buffer_);

      double *tmp = buffer_;
      int ifn = bs1->shell_to_function(ish);
      int jfnsave = bs2->shell_to_function(jsh);
      for (int i=0; i<nish; i++,ifn++) {
          if (ifn < b->istart || ifn >= b->iend) {
              tmp += njsh;
            }
          else {
              int jfn = jfnsave;
              int data_index = (ifn - b->istart)*njdata + jfn - b->jstart;
              for (int j=0; j<njsh; j++,jfn++) {
                  if (jfn >= b->jstart && jfn < b->jend) {
                      b->data[data_index] += *tmp * scale;
                      data_index++;
                    }
                  tmp++;
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
  int shstart = bs1->function_to_shell(fnstart);
  int shend = (fnend?bs1->function_to_shell(fnend - 1) + 1 : 0);

  // loop over all needed shells
  iter->reset(shstart, shend, 0, 0);
  for (iter->start_ltri(); iter->ready_ltri(); iter->next_ltri()) {
      int ish=iter->ishell();
      int jsh=iter->jshell();

      int nish = bs1->operator[](ish).nfunction();
      int njsh = bs1->operator[](jsh).nfunction();

      double scale = iter->scale();

      // compute a set of shell integrals
      compute_shell(ish,jsh,buffer_);

      // take the integrals from buffer and put them into the LTri block
      double*tmp = buffer_;
      int ifn = bs1->shell_to_function(ish);
      int jfnsave = bs1->shell_to_function(jsh);
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
                      b->data[data_index] += *tmp * scale;
                      data_index++;
                    }
                  tmp++;
                }
            }
        }
    }
}

int
OneBodyInt::has_side_effects()
{
  return 1;
}

///////////////////////////////////////////////////////////////////////

OneBody3Int::OneBody3Int(const RefGaussianBasisSet&b, OneBodyIntIter *it) :
  bs1(b), bs2(b), iter(it)
{
  // allocate a buffer
  int biggest_shell = b->max_nfunction_in_shell();
  if (biggest_shell) {
      buffer_ = new double[biggest_shell * biggest_shell * 3];
    }
  else {
      buffer_ = 0;
    }

  if (!iter)
    iter = new OneBodyIntIter;
}

OneBody3Int::OneBody3Int(const RefGaussianBasisSet&b1,
                         const RefGaussianBasisSet&b2,
                         OneBodyIntIter *it) :
  bs1(b1), bs2(b2), iter(it)
{
  // allocate a buffer
  int biggest_shell = b1->max_nfunction_in_shell() *
                      b2->max_nfunction_in_shell();
  if (biggest_shell) {
      buffer_ = new double[biggest_shell * 3];
    }
  else {
      buffer_ = 0;
   }

  if (!iter)
    iter = new OneBodyIntIter;
}

OneBody3Int::~OneBody3Int()
{
  if (buffer_) {
    delete[] buffer_;
    buffer_=0;
  }
  if (iter) {
    delete iter;
    iter=0;
  }
}

int OneBody3Int::nbasis()
{
  return bs1->nbasis();
}

int OneBody3Int::nbasis1()
{
  return bs1->nbasis();
}

int OneBody3Int::nbasis2()
{
  return bs2->nbasis();
}

int OneBody3Int::nshell()
{
  return bs1->nshell();
}

int OneBody3Int::nshell1()
{
  return bs1->nshell();
}

int OneBody3Int::nshell2()
{
  return bs2->nshell();
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
  int ishstart = bs1->function_to_shell(b->istart);
  int jshstart = bs2->function_to_shell(b->jstart);
  int ishend = bs1->function_to_shell(b->iend);
  int jshend = bs2->function_to_shell(b->jend);

  int njdata = b->jend - b->jstart;

  iter->reset(ishstart, ishend, jshstart, jshend);
  for (iter->start(); iter->ready(); iter->next()) {
      int ish=iter->ishell();
      int jsh=iter->jshell();

      int nish = bs1->operator[](ish).nfunction();
      int njsh = bs2->operator[](jsh).nfunction();

      double scale = iter->scale();

      // compute a set of shell integrals
      compute_shell(ish,jsh,buffer_);

      double*tmp = buffer_;
      int ifn = bs1->shell_to_function(ish);
      int jfnsave = bs2->shell_to_function(jsh);
      for (int i=0; i<nish; i++,ifn++) {
          if (ifn < b->istart || ifn >= b->iend) {
              tmp += njsh * 3;
            }
          else {
              int jfn = jfnsave;
              int data_index = (ifn - b->istart)*njdata + jfn - b->jstart;
              for (int j=0; j<njsh; j++,jfn++) {
                  if (jfn >= b->jstart && jfn < b->jend) {
                      a->data[data_index] += tmp[0] * scale;
                      b->data[data_index] += tmp[1] * scale;
                      c->data[data_index] += tmp[2] * scale;
                      data_index++;
                    }
                  tmp += 3;
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
  int shstart = bs1->function_to_shell(fnstart);
  int shend = (fnend?bs1->function_to_shell(fnend - 1) + 1 : 0);

  // loop over all needed shells
  iter->reset(shstart, shend, 0, 0);
  for (iter->start_ltri(); iter->ready_ltri(); iter->next_ltri()) {
      int ish=iter->ishell();
      int jsh=iter->jshell();

      int nish = bs1->operator[](ish).nfunction();
      int njsh = bs1->operator[](jsh).nfunction();

      double scale = iter->scale();
    
      // compute a set of shell integrals
      compute_shell(ish,jsh,buffer_);

      // take the integrals from buffer and put them into the LTri block
      double*tmp = buffer_;
      int ifn = bs1->shell_to_function(ish);
      int jfnsave = bs1->shell_to_function(jsh);
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
                      a->data[data_index] += tmp[0] * scale;
                      b->data[data_index] += tmp[1] * scale;
                      c->data[data_index] += tmp[2] * scale;
                      data_index++;
                    }
                  tmp += 3;
                }
            }
        }
    }
}

int
OneBody3Int::has_side_effects()
{
  return 1;
}

int
OneBody3Int::has_side_effects_in_arg1()
{
  return 1;
}

int
OneBody3Int::has_side_effects_in_arg2()
{
  return 1;
}
