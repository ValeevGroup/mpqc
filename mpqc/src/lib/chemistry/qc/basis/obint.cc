
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdio.h>

#include <math/scmat/block.h>

#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/basis.h>

///////////////////////////////////////////////////////////////////////

OneBodyInt::OneBodyInt(const RefGaussianBasisSet&b) :
  bs1(b), bs2(b)
{
  // allocate a buffer
  int biggest_shell = b->max_nfunction_in_shell();
  biggest_shell *= biggest_shell;

  if (biggest_shell) {
    buffer_ = new double[biggest_shell];
  } else {
    buffer_ = 0;
  }
}

OneBodyInt::OneBodyInt(const RefGaussianBasisSet&b1,
                       const RefGaussianBasisSet&b2) :
  bs1(b1), bs2(b2)
{
  // allocate a buffer
  int biggest_shell = b1->max_nfunction_in_shell() *
                      b2->max_nfunction_in_shell();
    
  if (biggest_shell) {
    buffer_ = new double[biggest_shell];
  } else {
    buffer_ = 0;
  }
}

OneBodyInt::~OneBodyInt()
{
  if (buffer_) {
    delete[] buffer_;
    buffer_=0;
  }
}

int
OneBodyInt::nbasis() const
{
  return bs1->nbasis();
}

int
OneBodyInt::nbasis1() const
{
  return bs1->nbasis();
}

int
OneBodyInt::nbasis2() const
{
  return bs2->nbasis();
}

int
OneBodyInt::nshell() const
{
  return bs1->nshell();
}

int
OneBodyInt::nshell1() const
{
  return bs1->nshell();
}

int
OneBodyInt::nshell2() const
{
  return bs2->nshell();
}

RefGaussianBasisSet
OneBodyInt::basis()
{
  return bs1;
}

RefGaussianBasisSet
OneBodyInt::basis1()
{
  return bs1;
}

RefGaussianBasisSet
OneBodyInt::basis2()
{
  return bs2;
}

const double *
OneBodyInt::buffer() const
{
  return buffer_;
}

///////////////////////////////////////////////////////////////////////

ShellPairIter::ShellPairIter()
{
}

ShellPairIter::~ShellPairIter()
{
}

void
ShellPairIter::init(const double * b, int ishell, int jshell,
                    int fi, int fj, int ni, int nj,
                    int red, double scl)
{
  e12 = ((ishell==jshell) && red);
  
  ioffset=fi;
  joffset=fj;

  iend=ni;
  jend=nj;

  buf=b;
  scale_=scl;
}

///////////////////////////////////////////////////////////////////////

OneBodyIntIter::OneBodyIntIter()
{
}

OneBodyIntIter::OneBodyIntIter(const RefOneBodyInt& o) :
  obi(o)
{
}

OneBodyIntIter::~OneBodyIntIter()
{
}

void
OneBodyIntIter::start(int ist, int jst, int ien, int jen)
{
  istart=ist;
  jstart=jst;
  iend=ien;
  jend=jen;
  
  icur=istart;
  jcur=jstart;

  if (!iend) {
    iend=obi->nshell1();
    jend=obi->nshell2();
  }

  ij = (icur*(icur+1)>>1) + jcur;
}

static inline int
min(int i, int j)
{
  return (i<j) ? i : j;
}

void
OneBodyIntIter::next()
{
  int jlast = (redund) ? min(icur,jend-1) : jend-1;
  
  if (jcur < jlast) {
    jcur++;
    ij++;
    return;
  }

  jcur=jstart;
  icur++;

  ij = (icur*(icur+1)>>1) + jcur;
}

double
OneBodyIntIter::scale() const
{
  return 1.0;
}

ShellPairIter&
OneBodyIntIter::current_pair()
{
  obi->compute_shell(icur,jcur);
  spi.init(obi->buffer(), icur, jcur,
           obi->basis1()->shell_to_function(icur),
           obi->basis2()->shell_to_function(jcur),
           obi->basis1()->operator()(icur).nfunction(),
           obi->basis2()->operator()(jcur).nfunction(),
           redund, scale()
           );

  return spi;
}

///////////////////////////////////////////////////////////////////////

OneBodyIntOp::OneBodyIntOp(const RefOneBodyInt& it)
{
  iter = new OneBodyIntIter(it);
}

OneBodyIntOp::OneBodyIntOp(const RefOneBodyIntIter& it) :
  iter(it)
{
}

OneBodyIntOp::~OneBodyIntOp()
{
}

void
OneBodyIntOp::process(SCMatrixBlockIter&)
{
  fprintf(stderr,"OneBodyIntOp::process(SCMatrixBlockIter&):"
          " cannot handle generic case\n");
  abort();
}

void
OneBodyIntOp::process(SCMatrixRectBlock* b)
{
  RefGaussianBasisSet bs1 = iter->one_body_int()->basis1();
  RefGaussianBasisSet bs2 = iter->one_body_int()->basis2();
  
  // convert basis function indices into shell indices
  int ishstart = bs1->function_to_shell(b->istart);
  int jshstart = bs2->function_to_shell(b->jstart);

  int b1end = b->iend;
  int ishend = (b1end?bs1->function_to_shell(b1end-1) + 1 : 0);

  int b2end = b->jend;
  int jshend = (b2end?bs2->function_to_shell(b2end-1) + 1 : 0);

  int njdata = b->jend - b->jstart;

  iter->redundant(0);

  for (iter->start(ishstart,jshstart,ishend,jshend);
       iter->ready(); iter->next()) {
    int ish=iter->ishell();
    int jsh=iter->jshell();

    ShellPairIter& spi = iter->current_pair();

    for (spi.start(); spi.ready(); spi.next()) {
      int ifn = spi.i();
      int jfn = spi.j();
      
      if (ifn < b->istart || ifn >= b->iend ||
          jfn < b->jstart || jfn >= b->jend)
        continue;
      
      int data_index = (ifn - b->istart)*njdata + jfn - b->jstart;
      b->data[data_index] += spi.val();
    }
  }
}

void
OneBodyIntOp::process(SCMatrixLTriBlock* b)
{
  RefGaussianBasisSet bs1 = iter->one_body_int()->basis1();

  // convert basis function indices into shell indices
  int fnstart = b->start;
  int fnend = b->end;
  int shstart = bs1->function_to_shell(fnstart);
  int shend = (fnend?bs1->function_to_shell(fnend - 1) + 1 : 0);

  iter->redundant(1);

  // loop over all needed shells
  for (iter->start(shstart,shstart,shend,shend); iter->ready(); iter->next()) {
    int ish=iter->ishell();
    int jsh=iter->jshell();

    ShellPairIter& spi = iter->current_pair();

    // compute a set of shell integrals
    for (spi.start(); spi.ready(); spi.next()) {
      int ifn = spi.i();
      int jfn = spi.j();
      
      if (ifn < fnstart || ifn >= fnend)
        continue;
      
      int ioff = ifn-fnstart;
      int joff = jfn-fnstart;
      
      int data_index = (ioff*(ioff+1)>>1)+joff;
      
      //printf("%d %d %f\n",ifn,jfn,spi.val());
      b->data[data_index] += spi.val();
      data_index++;
    }
  }
}

int
OneBodyIntOp::has_side_effects()
{
  return 1;
}

///////////////////////////////////////////////////////////////////////

OneBody3IntOp::OneBody3IntOp(const RefOneBodyInt& it)
{
  iter = new OneBodyIntIter(it);
}

OneBody3IntOp::OneBody3IntOp(const RefOneBodyIntIter& it) :
  iter(it)
{
}

OneBody3IntOp::~OneBody3IntOp()
{
}

void
OneBody3IntOp::process(SCMatrixBlockIter&,
                       SCMatrixBlockIter&,
                       SCMatrixBlockIter&)
{
  fprintf(stderr,"OneBody3IntOp::process(SCMatrixBlockIter&):"
          " cannot handle generic case\n");
  abort();
}

void
OneBody3IntOp::process(SCMatrixRectBlock* a,
                       SCMatrixRectBlock* b,
                       SCMatrixRectBlock* c)
{
  RefGaussianBasisSet bs1 = iter->one_body_int()->basis1();
  RefGaussianBasisSet bs2 = iter->one_body_int()->basis2();

  // convert basis function indices into shell indices
  int ishstart = bs1->function_to_shell(b->istart);
  int jshstart = bs2->function_to_shell(b->jstart);

  int ishend = bs1->function_to_shell(b->iend);
  int jshend = bs2->function_to_shell(b->jend);

  int njdata = b->jend - b->jstart;

  iter->redundant(0);

  for (iter->start(ishstart,jshstart,ishend,jshend);
       iter->ready(); iter->next()) {
    int ish=iter->ishell();
    int jsh=iter->jshell();

    // compute a set of shell integrals
    ShellPairIter& spi = iter->current_pair();

    for (spi.start(); spi.ready(); spi.next()) {
      int ifn = spi.i();
      int jfn = spi.j();
        
      if (ifn < b->istart || ifn >= b->iend ||
          jfn < b->jstart || jfn >= b->jend)
        continue;
      
      int data_index = (ifn - b->istart)*njdata + jfn - b->jstart;

#if 0
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
#endif
    }
  }
}

void
OneBody3IntOp::process(SCMatrixLTriBlock* a,
                       SCMatrixLTriBlock* b,
                       SCMatrixLTriBlock* c)
{
#if 0
  RefGaussianBasisSet bs1 = iter->one_body_int()->basis1();

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
#endif
}

int
OneBody3IntOp::has_side_effects()
{
  return 1;
}

int
OneBody3IntOp::has_side_effects_in_arg1()
{
  return 1;
}

int
OneBody3IntOp::has_side_effects_in_arg2()
{
  return 1;
}
