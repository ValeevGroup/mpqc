
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdio.h>

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
  icur=0;
  jcur=0;
  kcur=0;
  lcur=0;
}

void
ShellQuartetIter::next()
{
  index++;
  
  if (lcur < ((e34) ? (((e13e24)&&((kcur)==(icur)))?(jcur):(kcur))
              : ((e13e24)&&((kcur)==(icur)))?(jcur):(lend)-1)) {
    lcur++;
    return;
  }

  lcur=0;

  if (kcur < ((e13e24)?(icur):((kend)-1))) {
    kcur++;
    return;
  }

  kcur=0;

  if (jcur < ((e12)?(icur):((jend)-1))) {
    jcur++;
    return;
  }

  jcur=0;
  icur++;
}

ShellQuartetIter::operator int()
{
  return (icur < iend);
}

int
ShellQuartetIter::i() const
{
  return icur+istart;
}

int
ShellQuartetIter::j() const
{
  return jcur+jstart;
}

int
ShellQuartetIter::k() const
{
  return kcur+kstart;
}

int
ShellQuartetIter::l() const
{
  return lcur+lstart;
}

int
ShellQuartetIter::ij() const
{
  return i()*(i()+1)>>1 + j();
}

int
ShellQuartetIter::ik() const
{
  return i()*(i()+1)>>1 + k();
}

int
ShellQuartetIter::il() const
{
  return i()*(i()+1)>>1 + l();
}

int
ShellQuartetIter::kl() const
{
  return k()*(k()+1)>>1 + l();
}

int
ShellQuartetIter::jl() const
{
  return j()*(j()+1)>>1 + l();
}

int
ShellQuartetIter::jk() const
{
  return j()*(j()+1)>>1 + k();
}

int
ShellQuartetIter::ijkl() const
{
  return ij()*(ij()+1)>>1 + kl();
}

double
ShellQuartetIter::val() const
{
  return buf[index]*scale_;
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

TwoBodyIntIter::operator int()
{
  return (icur < iend);
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

int
TwoBodyIntIter::ishell() const
{
  return icur;
}

int
TwoBodyIntIter::jshell() const
{
  return jcur;
}

int
TwoBodyIntIter::kshell() const
{
  return kcur;
}

int
TwoBodyIntIter::lshell() const
{
  return lcur;
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
           tbi->basis()->function_to_shell(icur),
           tbi->basis()->function_to_shell(jcur),
           tbi->basis()->function_to_shell(kcur),
           tbi->basis()->function_to_shell(lcur),
           tbi->basis()->operator()(icur).nfunction(),
           tbi->basis()->operator()(jcur).nfunction(),
           tbi->basis()->operator()(kcur).nfunction(),
           tbi->basis()->operator()(lcur).nfunction(),
           scale()
    );

  return sqi;
}
