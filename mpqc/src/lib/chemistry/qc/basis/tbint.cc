
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/scmat/offset.h>

#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/basis/basis.h>

///////////////////////////////////////////////////////////////////////

TwoBodyInt::TwoBodyInt(const RefGaussianBasisSet&b1,
                       const RefGaussianBasisSet&b2,
                       const RefGaussianBasisSet&b3,
                       const RefGaussianBasisSet&b4) :
  bs1_(b1), bs2_(b2), bs3_(b3), bs4_(b4), redundant_(1)
{
  buffer_ = 0;
}

TwoBodyInt::~TwoBodyInt()
{
}

int
TwoBodyInt::nbasis() const
{
  return bs1_->nbasis();
}

int
TwoBodyInt::nbasis1() const
{
  return bs1_->nbasis();
}

int
TwoBodyInt::nbasis2() const
{
  return bs2_->nbasis();
}

int
TwoBodyInt::nbasis3() const
{
  return bs3_->nbasis();
}

int
TwoBodyInt::nbasis4() const
{
  return bs4_->nbasis();
}

int
TwoBodyInt::nshell() const
{
  return bs1_->nshell();
}

int
TwoBodyInt::nshell1() const
{
  return bs1_->nshell();
}

int
TwoBodyInt::nshell2() const
{
  return bs2_->nshell();
}

int
TwoBodyInt::nshell3() const
{
  return bs3_->nshell();
}

int
TwoBodyInt::nshell4() const
{
  return bs4_->nshell();
}

RefGaussianBasisSet
TwoBodyInt::basis()
{
  return bs1_;
}

RefGaussianBasisSet
TwoBodyInt::basis1()
{
  return bs1_;
}

RefGaussianBasisSet
TwoBodyInt::basis2()
{
  return bs2_;
}

RefGaussianBasisSet
TwoBodyInt::basis3()
{
  return bs3_;
}

RefGaussianBasisSet
TwoBodyInt::basis4()
{
  return bs4_;
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
                       double scl, int redund)
{
  redund_ = redund;
  
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

  if (redund_) {
    if (lcur < lend-1) {
      lcur++;
      l_++;
      return;
    }

    lcur=0;
    l_=lstart;
    
    if (kcur < kend-1) {
      kcur++;
      k_++;
      return;
    }

    kcur=0;
    k_=kstart;

    if (jcur < jend-1) {
      jcur++;
      j_++;
      return;
    }

    jcur=0;
    j_=jstart;
  
    icur++;
    i_++;

  } else {
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
           scale(),
           tbi->redundant()
    );

  return sqi;
}

///////////////////////////////////////////////////////////////////////

TwoBodyDerivInt::TwoBodyDerivInt(const RefGaussianBasisSet&b1,
                                 const RefGaussianBasisSet&b2,
                                 const RefGaussianBasisSet&b3,
                                 const RefGaussianBasisSet&b4):
  bs1_(b1), bs2_(b2), bs3_(b3), bs4_(b4)
{
  buffer_ = 0;
}

TwoBodyDerivInt::~TwoBodyDerivInt()
{
}

int
TwoBodyDerivInt::nbasis() const
{
  return bs1_->nbasis();
}

int
TwoBodyDerivInt::nbasis1() const
{
  return bs1_->nbasis();
}

int
TwoBodyDerivInt::nbasis2() const
{
  return bs2_->nbasis();
}

int
TwoBodyDerivInt::nbasis3() const
{
  return bs3_->nbasis();
}

int
TwoBodyDerivInt::nbasis4() const
{
  return bs4_->nbasis();
}

int
TwoBodyDerivInt::nshell() const
{
  return bs1_->nshell();
}

int
TwoBodyDerivInt::nshell1() const
{
  return bs1_->nshell();
}

int
TwoBodyDerivInt::nshell2() const
{
  return bs2_->nshell();
}

int
TwoBodyDerivInt::nshell3() const
{
  return bs3_->nshell();
}

int
TwoBodyDerivInt::nshell4() const
{
  return bs4_->nshell();
}

RefGaussianBasisSet
TwoBodyDerivInt::basis()
{
  return bs1_;
}

RefGaussianBasisSet
TwoBodyDerivInt::basis1()
{
  return bs1_;
}

RefGaussianBasisSet
TwoBodyDerivInt::basis2()
{
  return bs2_;
}

RefGaussianBasisSet
TwoBodyDerivInt::basis3()
{
  return bs3_;
}

RefGaussianBasisSet
TwoBodyDerivInt::basis4()
{
  return bs4_;
}

const double *
TwoBodyDerivInt::buffer() const
{
  return buffer_;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
