
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdio.h>

#include <math/scmat/offset.h>

#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/basis/basis.h>

///////////////////////////////////////////////////////////////////////

TwoBodyInt::TwoBodyInt(const RefGaussianBasisSet&b) :
  bs1(b), bs2(b), bs3(b), bs4(b)
{
  // allocate a buffer
  int biggest_shell = b->max_nfunction_in_shell();
  biggest_shell *= biggest_shell;
  biggest_shell *= biggest_shell;

  if (biggest_shell) {
    buffer_ = new double[biggest_shell];
  } else {
    buffer_ = 0;
  }

  store1_=store2_=1;
  int_store_=0;
}

TwoBodyInt::TwoBodyInt(const RefGaussianBasisSet&b1,
                       const RefGaussianBasisSet&b2,
                       const RefGaussianBasisSet&b3,
                       const RefGaussianBasisSet&b4) :
  bs1(b1), bs2(b2), bs3(b3), bs4(b4)
{
  // allocate a buffer
  int biggest_shell = b1->max_nfunction_in_shell() *
                      b2->max_nfunction_in_shell() *
                      b3->max_nfunction_in_shell() *
                      b4->max_nfunction_in_shell();
    
  if (biggest_shell) {
    buffer_ = new double[biggest_shell];
  } else {
    buffer_ = 0;
  }

  store1_=store2_=1;
  int_store_=0;
}

TwoBodyInt::~TwoBodyInt()
{
  if (buffer_) {
    delete[] buffer_;
    buffer_=0;
  }
}

int
TwoBodyInt::nbasis() const
{
  return bs1->nbasis();
}

int
TwoBodyInt::nbasis1() const
{
  return bs1->nbasis();
}

int
TwoBodyInt::nbasis2() const
{
  return bs2->nbasis();
}

int
TwoBodyInt::nbasis3() const
{
  return bs3->nbasis();
}

int
TwoBodyInt::nbasis4() const
{
  return bs4->nbasis();
}

int
TwoBodyInt::nshell() const
{
  return bs1->nshell();
}

int
TwoBodyInt::nshell1() const
{
  return bs1->nshell();
}

int
TwoBodyInt::nshell2() const
{
  return bs2->nshell();
}

int
TwoBodyInt::nshell3() const
{
  return bs3->nshell();
}

int
TwoBodyInt::nshell4() const
{
  return bs4->nshell();
}

RefGaussianBasisSet
TwoBodyInt::basis()
{
  return bs1;
}

RefGaussianBasisSet
TwoBodyInt::basis1()
{
  return bs1;
}

RefGaussianBasisSet
TwoBodyInt::basis2()
{
  return bs2;
}

RefGaussianBasisSet
TwoBodyInt::basis3()
{
  return bs3;
}

RefGaussianBasisSet
TwoBodyInt::basis4()
{
  return bs4;
}

const double *
TwoBodyInt::buffer() const
{
  return buffer_;
}

///////////////////////////////////////////////////////////////////////

ShellQuartetIter::ShellQuartetIter()
{
}

ShellQuartetIter::~ShellQuartetIter()
{
}

void
ShellQuartetIter::init(const double * b,
                       int is, int js, int ks, int ls,
                       int fi, int fj, int fk, int fl,
                       int ni, int nj, int nk, int nl,
                       double scl)
{
  e12 = (is==js);
  e34 = (ks==ls);
  e13e24 = (is==ks) && (js==ls);
  
  istart=fi;
  jstart=fj;
  kstart=fk;
  lstart=fl;

  index=0;
  
  iend=ni;
  jend=nj;
  kend=nk;
  lend=nl;

  buf=b;
  scale_=scl;
}

void
ShellQuartetIter::start()
{
  icur=0; i_ = istart;
  jcur=0; j_ = jstart;
  kcur=0; k_ = kstart;
  lcur=0; l_ = lstart;
}

void
ShellQuartetIter::next()
{
  index++;
  
  if (lcur < ((e34) ? (((e13e24)&&((kcur)==(icur)))?(jcur):(kcur))
              : ((e13e24)&&((kcur)==(icur)))?(jcur):(lend)-1)) {
    lcur++;
    l_++;
    return;
  }

  lcur=0;
  l_=lstart;

  if (kcur < ((e13e24)?(icur):((kend)-1))) {
    kcur++;
    k_++;
    return;
  }

  kcur=0;
  k_=kstart;

  if (jcur < ((e12)?(icur):((jend)-1))) {
    jcur++;
    j_++;
    return;
  }

  jcur=0;
  j_=jstart;
  
  icur++;
  i_++;
}

///////////////////////////////////////////////////////////////////////

TwoBodyIntIter::TwoBodyIntIter()
{
}

TwoBodyIntIter::TwoBodyIntIter(const RefTwoBodyInt& t) :
  tbi(t)
{
}

TwoBodyIntIter::~TwoBodyIntIter()
{
}

void
TwoBodyIntIter::start()
{
  icur=0;
  jcur=0;
  kcur=0;
  lcur=0;

  iend = tbi->nshell();
}

void
TwoBodyIntIter::next()
{
  if (lcur < ((icur==kcur) ? jcur : kcur)) { // increment l loop?
    lcur++;
    return;
  }

  // restart l loop
  lcur=0;

  if (kcur < icur) { // increment k loop?
    kcur++;
    return;
  }

  // restart k loop
  kcur=0;

  if (jcur < icur) { // increment j loop?
    jcur++;
    return;
  }

  // restart j loop
  jcur=0;

  // increment i loop
  icur++;
}

double
TwoBodyIntIter::scale() const
{
  return 1.0;
}

ShellQuartetIter&
TwoBodyIntIter::current_quartet()
{
  tbi->compute_shell(icur,jcur,kcur,lcur);
  
  sqi.init(tbi->buffer(),
           icur, jcur, kcur, lcur,
           tbi->basis()->shell_to_function(icur),
           tbi->basis()->shell_to_function(jcur),
           tbi->basis()->shell_to_function(kcur),
           tbi->basis()->shell_to_function(lcur),
           tbi->basis()->operator()(icur).nfunction(),
           tbi->basis()->operator()(jcur).nfunction(),
           tbi->basis()->operator()(kcur).nfunction(),
           tbi->basis()->operator()(lcur).nfunction(),
           scale()
    );

  return sqi;
}

///////////////////////////////////////////////////////////////////////

TwoBodyDerivInt::TwoBodyDerivInt(const RefGaussianBasisSet&b) :
  bs1(b), bs2(b), bs3(b), bs4(b)
{
  // allocate a buffer
  int biggest_shell = b->max_nfunction_in_shell();
  biggest_shell *= biggest_shell;
  biggest_shell *= biggest_shell;
  biggest_shell *= 12;

  if (biggest_shell) {
    buffer_ = new double[biggest_shell];
  } else {
    buffer_ = 0;
  }

  store1_=store2_=1;
  int_store_=0;
}

TwoBodyDerivInt::TwoBodyDerivInt(const RefGaussianBasisSet&b1,
                                 const RefGaussianBasisSet&b2,
                                 const RefGaussianBasisSet&b3,
                                 const RefGaussianBasisSet&b4) :
  bs1(b1), bs2(b2), bs3(b3), bs4(b4)
{
  // allocate a buffer
  int biggest_shell = b1->max_nfunction_in_shell() *
                      b2->max_nfunction_in_shell() *
                      b3->max_nfunction_in_shell() *
                      b4->max_nfunction_in_shell() * 12;
    
  if (biggest_shell) {
    buffer_ = new double[biggest_shell];
  } else {
    buffer_ = 0;
  }

  store1_=store2_=1;
  int_store_=0;
}

TwoBodyDerivInt::~TwoBodyDerivInt()
{
  if (buffer_) {
    delete[] buffer_;
    buffer_=0;
  }
}

int
TwoBodyDerivInt::nbasis() const
{
  return bs1->nbasis();
}

int
TwoBodyDerivInt::nbasis1() const
{
  return bs1->nbasis();
}

int
TwoBodyDerivInt::nbasis2() const
{
  return bs2->nbasis();
}

int
TwoBodyDerivInt::nbasis3() const
{
  return bs3->nbasis();
}

int
TwoBodyDerivInt::nbasis4() const
{
  return bs4->nbasis();
}

int
TwoBodyDerivInt::nshell() const
{
  return bs1->nshell();
}

int
TwoBodyDerivInt::nshell1() const
{
  return bs1->nshell();
}

int
TwoBodyDerivInt::nshell2() const
{
  return bs2->nshell();
}

int
TwoBodyDerivInt::nshell3() const
{
  return bs3->nshell();
}

int
TwoBodyDerivInt::nshell4() const
{
  return bs4->nshell();
}

RefGaussianBasisSet
TwoBodyDerivInt::basis()
{
  return bs1;
}

RefGaussianBasisSet
TwoBodyDerivInt::basis1()
{
  return bs1;
}

RefGaussianBasisSet
TwoBodyDerivInt::basis2()
{
  return bs2;
}

RefGaussianBasisSet
TwoBodyDerivInt::basis3()
{
  return bs3;
}

RefGaussianBasisSet
TwoBodyDerivInt::basis4()
{
  return bs4;
}

const double *
TwoBodyDerivInt::buffer() const
{
  return buffer_;
}

