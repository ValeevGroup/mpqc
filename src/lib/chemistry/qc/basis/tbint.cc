//
// tbint.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

#include <cassert>
#include <limits.h>

#include <util/misc/scexception.h>
#include <math/scmat/offset.h>

#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/basis.h>

using namespace sc;

double* init_log2_to_double();

///////////////////////////////////////////////////////////////////////

TwoBodyInt::TwoBodyInt(Integral *integral,
                       const Ref<GaussianBasisSet>&b1,
                       const Ref<GaussianBasisSet>&b2,
                       const Ref<GaussianBasisSet>&b3,
                       const Ref<GaussianBasisSet>&b4) :
  integral_(integral),
  bs1_(b1), bs2_(b2), bs3_(b3), bs4_(b4), redundant_(1)
{
  MPQC_ASSERT(bs1_.nonnull() && bs2_.nonnull() && bs3_.nonnull() && bs4_.nonnull());
  integral_->reference();
  buffer_ = 0;
  log2_to_double_ = init_log2_to_double();
}

TwoBodyInt::~TwoBodyInt()
{
  integral_->dereference();
  if (integral_->nreference() == 0) delete integral_;
  delete[] log2_to_double_;
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

Ref<GaussianBasisSet>
TwoBodyInt::basis(size_t c)
{
  if (c >= 4)
    throw ProgrammingError("TwoBodyInt::basis(c): c >= 4",
                           __FILE__, __LINE__);
  switch (c) {
    case 0: return bs1_; break;
    case 1: return bs2_; break;
    case 2: return bs3_; break;
    case 3: return bs4_; break;
    default: MPQC_ASSERT(false); // unreachable
  }
  return 0; // unreachable
}

Ref<GaussianBasisSet>
TwoBodyInt::basis1()
{
  return bs1_;
}

Ref<GaussianBasisSet>
TwoBodyInt::basis2()
{
  return bs2_;
}

Ref<GaussianBasisSet>
TwoBodyInt::basis3()
{
  return bs3_;
}

Ref<GaussianBasisSet>
TwoBodyInt::basis4()
{
  return bs4_;
}

const double *
TwoBodyInt::buffer(TwoBodyOper::type i) const
{
  if (i == TwoBodyOper::eri) return buffer_;
  return 0;
}

std::pair<std::map<TwoBodyOper::type,const double*>,std::array<unsigned long, 4> >
TwoBodyInt::compute_shell_arrays(int i,int j, int k, int l)
{
  int saved_redundant = redundant();
  set_redundant(1);

  compute_shell(i,j,k,l);

  set_redundant(saved_redundant);

  std::pair<std::map<TwoBodyOper::type,const double*>,std::array<unsigned long, 4> > r;

  const Ref<TwoBodyOperSetDescr>& descr = this->descr();
  for (unsigned int it=0; it<descr->size(); it++) {
      TwoBodyOper::type tbit = descr->opertype(it);
      r.first[tbit] = buffer(tbit);
  }

  r.second[0] = basis1()->shell(i).nfunction();
  r.second[1] = basis2()->shell(j).nfunction();
  r.second[2] = basis3()->shell(k).nfunction();
  r.second[3] = basis4()->shell(l).nfunction();

  return r;
}

void
TwoBodyInt::set_integral_storage(size_t storage)
{
}

double
TwoBodyInt::shell_bound(int s1, int s2, int s3, int s4)
{
  int ibound = log2_shell_bound(s1,s2,s3,s4);
  if( ibound < SCHAR_MIN )
    return log2_to_double_[0];
  else if( ibound > SCHAR_MAX )
    return log2_to_double_[ SCHAR_MAX - SCHAR_MIN ];
  return log2_to_double_[ ibound - SCHAR_MIN ];
}

bool
TwoBodyInt::cloneable() const
{
  return false;
}

Ref<TwoBodyInt>
TwoBodyInt::clone()
{
  throw FeatureNotImplemented("TwoBodyInt::clone() not implemented",__FILE__,__LINE__);
}

///////////////////////////////////////////////////////////////////////

TwoBodyThreeCenterInt::TwoBodyThreeCenterInt(Integral *integral,
                       const Ref<GaussianBasisSet>&b1,
                       const Ref<GaussianBasisSet>&b2,
                       const Ref<GaussianBasisSet>&b3) :
  integral_(integral),
  bs1_(b1), bs2_(b2), bs3_(b3), redundant_(1)
{
  integral_->reference();
  buffer_ = 0;
  log2_to_double_ = init_log2_to_double();
}

TwoBodyThreeCenterInt::~TwoBodyThreeCenterInt()
{
  integral_->dereference();
  if (integral_->nreference() == 0) delete integral_;
  delete[] log2_to_double_;
}

int
TwoBodyThreeCenterInt::nbasis() const
{
  return bs1_->nbasis();
}

int
TwoBodyThreeCenterInt::nbasis1() const
{
  return bs1_->nbasis();
}

int
TwoBodyThreeCenterInt::nbasis2() const
{
  return bs2_->nbasis();
}

int
TwoBodyThreeCenterInt::nbasis3() const
{
  return bs3_->nbasis();
}

int
TwoBodyThreeCenterInt::nshell() const
{
  return bs1_->nshell();
}

int
TwoBodyThreeCenterInt::nshell1() const
{
  return bs1_->nshell();
}

int
TwoBodyThreeCenterInt::nshell2() const
{
  return bs2_->nshell();
}

int
TwoBodyThreeCenterInt::nshell3() const
{
  return bs3_->nshell();
}

Ref<GaussianBasisSet>
TwoBodyThreeCenterInt::basis(size_t c)
{
  if (c >= 3)
    throw ProgrammingError("TwoBodyThreeCenterInt::basis(c): c >= 4",
                           __FILE__, __LINE__);
  switch (c) {
    case 0: return bs1_; break;
    case 1: return bs2_; break;
    case 2: return bs3_; break;
    default: MPQC_ASSERT(false); // unreachable
  }
  return 0; // unreachable
}

Ref<GaussianBasisSet>
TwoBodyThreeCenterInt::basis1()
{
  return bs1_;
}

Ref<GaussianBasisSet>
TwoBodyThreeCenterInt::basis2()
{
  return bs2_;
}

Ref<GaussianBasisSet>
TwoBodyThreeCenterInt::basis3()
{
  return bs3_;
}

const double *
TwoBodyThreeCenterInt::buffer(TwoBodyOper::type i) const
{
  if (i==TwoBodyOper::eri) return buffer_;
  return 0;
}

void
TwoBodyThreeCenterInt::set_integral_storage(size_t storage)
{
}

double
TwoBodyThreeCenterInt::shell_bound(int s1, int s2, int s3)
{
  int ibound = log2_shell_bound(s1,s2,s3);
  if( ibound < SCHAR_MIN )
    return log2_to_double_[0];
  else if( ibound > SCHAR_MAX )
    return log2_to_double_[ SCHAR_MAX - SCHAR_MIN ];
  return log2_to_double_[ ibound - SCHAR_MIN ];
}

bool
TwoBodyThreeCenterInt::cloneable() const
{
  return false;
}

Ref<TwoBodyThreeCenterInt>
TwoBodyThreeCenterInt::clone()
{
  throw FeatureNotImplemented("TwoBodyThreeCenterInt::clone() not implemented",__FILE__,__LINE__);
}

///////////////////////////////////////////////////////////////////////

TwoBodyTwoCenterInt::TwoBodyTwoCenterInt(Integral *integral,
                       const Ref<GaussianBasisSet>&b1,
                       const Ref<GaussianBasisSet>&b2) :
  integral_(integral),
  bs1_(b1), bs2_(b2), redundant_(1)
{
  integral_->reference();
  buffer_ = 0;
  log2_to_double_ = init_log2_to_double();
}

TwoBodyTwoCenterInt::~TwoBodyTwoCenterInt()
{
  integral_->dereference();
  if (integral_->nreference() == 0) delete integral_;
  delete[] log2_to_double_;
}

int
TwoBodyTwoCenterInt::nbasis() const
{
  return bs1_->nbasis();
}

int
TwoBodyTwoCenterInt::nbasis1() const
{
  return bs1_->nbasis();
}

int
TwoBodyTwoCenterInt::nbasis2() const
{
  return bs2_->nbasis();
}

int
TwoBodyTwoCenterInt::nshell() const
{
  return bs1_->nshell();
}

int
TwoBodyTwoCenterInt::nshell1() const
{
  return bs1_->nshell();
}

int
TwoBodyTwoCenterInt::nshell2() const
{
  return bs2_->nshell();
}

Ref<GaussianBasisSet>
TwoBodyTwoCenterInt::basis(size_t c)
{
  if (c >= 2)
    throw ProgrammingError("TwoBodyTwoCenterInt::basis(c): c >= 4",
                           __FILE__, __LINE__);
  switch (c) {
    case 0: return bs1_; break;
    case 1: return bs2_; break;
    default: MPQC_ASSERT(false); // unreachable
  }
  return 0; // unreachable
}

Ref<GaussianBasisSet>
TwoBodyTwoCenterInt::basis1()
{
  return bs1_;
}

Ref<GaussianBasisSet>
TwoBodyTwoCenterInt::basis2()
{
  return bs2_;
}

const double *
TwoBodyTwoCenterInt::buffer(TwoBodyOper::type i) const
{
  if (i==TwoBodyOper::eri) return buffer_;
  return 0;
}

void
TwoBodyTwoCenterInt::set_integral_storage(size_t storage)
{
}

double
TwoBodyTwoCenterInt::shell_bound(int s1, int s2)
{
  int ibound = log2_shell_bound(s1,s2);
  if( ibound < SCHAR_MIN )
    return log2_to_double_[0];
  else if( ibound > SCHAR_MAX )
    return log2_to_double_[ SCHAR_MAX - SCHAR_MIN ];
  return log2_to_double_[ ibound - SCHAR_MIN ];
}

bool
TwoBodyTwoCenterInt::cloneable() const
{
  return false;
}

Ref<TwoBodyTwoCenterInt>
TwoBodyTwoCenterInt::clone()
{
  throw FeatureNotImplemented("TwoBodyTwoCenterInt::clone() not implemented",__FILE__,__LINE__);
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

TwoBodyIntIter::TwoBodyIntIter(const Ref<TwoBodyInt>& t) :
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

TwoBodyDerivInt::TwoBodyDerivInt(Integral *integral,
                                 const Ref<GaussianBasisSet>&b1,
                                 const Ref<GaussianBasisSet>&b2,
                                 const Ref<GaussianBasisSet>&b3,
                                 const Ref<GaussianBasisSet>&b4):
  integral_(integral),
  bs1_(b1), bs2_(b2), bs3_(b3), bs4_(b4)
{
  integral_->reference();
  buffer_ = 0;
  log2_to_double_ = init_log2_to_double();
}

TwoBodyDerivInt::~TwoBodyDerivInt()
{
  integral_->dereference();
  if (integral_->nreference() == 0) delete integral_;
  delete[] log2_to_double_;
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

Ref<GaussianBasisSet>
TwoBodyDerivInt::basis()
{
  return bs1_;
}

Ref<GaussianBasisSet>
TwoBodyDerivInt::basis1()
{
  return bs1_;
}

Ref<GaussianBasisSet>
TwoBodyDerivInt::basis2()
{
  return bs2_;
}

Ref<GaussianBasisSet>
TwoBodyDerivInt::basis3()
{
  return bs3_;
}

Ref<GaussianBasisSet>
TwoBodyDerivInt::basis4()
{
  return bs4_;
}

const double *
TwoBodyDerivInt::buffer() const
{
  return buffer_;
}

double
TwoBodyDerivInt::shell_bound(int s1, int s2, int s3, int s4)
{
  int ibound = log2_shell_bound(s1,s2,s3,s4);
  if( ibound < SCHAR_MIN )
    return log2_to_double_[0];
  else if( ibound > SCHAR_MAX )
    return log2_to_double_[ SCHAR_MAX - SCHAR_MIN ];
  return log2_to_double_[ ibound - SCHAR_MIN ];
}

bool
TwoBodyDerivInt::cloneable() const
{
  return false;
}

Ref<TwoBodyDerivInt>
TwoBodyDerivInt::clone()
{
  throw FeatureNotImplemented("TwoBodyDerivInt::clone() not implemented",__FILE__,__LINE__);
}

///////////////////////////////////////////////////////////////////////

TwoBodyThreeCenterDerivInt::TwoBodyThreeCenterDerivInt(Integral *integral,
                                 const Ref<GaussianBasisSet>&b1,
                                 const Ref<GaussianBasisSet>&b2,
                                 const Ref<GaussianBasisSet>&b3):
  integral_(integral),
  bs1_(b1), bs2_(b2), bs3_(b3)
{
  integral_->reference();
  buffer_ = 0;
  log2_to_double_ = init_log2_to_double();
}

TwoBodyThreeCenterDerivInt::~TwoBodyThreeCenterDerivInt()
{
  integral_->dereference();
  if (integral_->nreference() == 0) delete integral_;
  delete[] log2_to_double_;
}

int
TwoBodyThreeCenterDerivInt::nbasis() const
{
  return bs1_->nbasis();
}

int
TwoBodyThreeCenterDerivInt::nbasis1() const
{
  return bs1_->nbasis();
}

int
TwoBodyThreeCenterDerivInt::nbasis2() const
{
  return bs2_->nbasis();
}

int
TwoBodyThreeCenterDerivInt::nbasis3() const
{
  return bs3_->nbasis();
}

int
TwoBodyThreeCenterDerivInt::nshell() const
{
  return bs1_->nshell();
}

int
TwoBodyThreeCenterDerivInt::nshell1() const
{
  return bs1_->nshell();
}

int
TwoBodyThreeCenterDerivInt::nshell2() const
{
  return bs2_->nshell();
}

int
TwoBodyThreeCenterDerivInt::nshell3() const
{
  return bs3_->nshell();
}

Ref<GaussianBasisSet>
TwoBodyThreeCenterDerivInt::basis()
{
  return bs1_;
}

Ref<GaussianBasisSet>
TwoBodyThreeCenterDerivInt::basis1()
{
  return bs1_;
}

Ref<GaussianBasisSet>
TwoBodyThreeCenterDerivInt::basis2()
{
  return bs2_;
}

Ref<GaussianBasisSet>
TwoBodyThreeCenterDerivInt::basis3()
{
  return bs3_;
}

const double *
TwoBodyThreeCenterDerivInt::buffer() const
{
  return buffer_;
}

double
TwoBodyThreeCenterDerivInt::shell_bound(int s1, int s2, int s3)
{
  int ibound = log2_shell_bound(s1,s2,s3);
  if( ibound < SCHAR_MIN )
    return log2_to_double_[0];
  else if( ibound > SCHAR_MAX )
    return log2_to_double_[ SCHAR_MAX - SCHAR_MIN ];
  return log2_to_double_[ ibound - SCHAR_MIN ];
}

///////////////////////////////////////////////////////////////////////

TwoBodyTwoCenterDerivInt::TwoBodyTwoCenterDerivInt(Integral *integral,
                                 const Ref<GaussianBasisSet>&b1,
                                 const Ref<GaussianBasisSet>&b2):
  integral_(integral),
  bs1_(b1), bs2_(b2)
{
  integral_->reference();
  buffer_ = 0;
  log2_to_double_ = init_log2_to_double();
}

TwoBodyTwoCenterDerivInt::~TwoBodyTwoCenterDerivInt()
{
  integral_->dereference();
  if (integral_->nreference() == 0) delete integral_;
  delete[] log2_to_double_;
}

int
TwoBodyTwoCenterDerivInt::nbasis() const
{
  return bs1_->nbasis();
}

int
TwoBodyTwoCenterDerivInt::nbasis1() const
{
  return bs1_->nbasis();
}

int
TwoBodyTwoCenterDerivInt::nbasis2() const
{
  return bs2_->nbasis();
}

int
TwoBodyTwoCenterDerivInt::nshell() const
{
  return bs1_->nshell();
}

int
TwoBodyTwoCenterDerivInt::nshell1() const
{
  return bs1_->nshell();
}

int
TwoBodyTwoCenterDerivInt::nshell2() const
{
  return bs2_->nshell();
}

Ref<GaussianBasisSet>
TwoBodyTwoCenterDerivInt::basis()
{
  return bs1_;
}

Ref<GaussianBasisSet>
TwoBodyTwoCenterDerivInt::basis1()
{
  return bs1_;
}

Ref<GaussianBasisSet>
TwoBodyTwoCenterDerivInt::basis2()
{
  return bs2_;
}

const double *
TwoBodyTwoCenterDerivInt::buffer() const
{
  return buffer_;
}

double
TwoBodyTwoCenterDerivInt::shell_bound(int s1, int s2)
{
  int ibound = log2_shell_bound(s1,s2);
  if( ibound < SCHAR_MIN )
    return log2_to_double_[0];
  else if( ibound > SCHAR_MAX )
    return log2_to_double_[ SCHAR_MAX - SCHAR_MIN ];
  return log2_to_double_[ ibound - SCHAR_MIN ];
}

/////////////////////////////////////////////////////////////////////////////

double*
init_log2_to_double()
{
  int n = -1 * SCHAR_MIN + SCHAR_MAX + 1;
  double* ptr = new double[n];

  int i=-1;
  for(int log2=SCHAR_MIN; log2<=SCHAR_MAX; ++log2)
    ptr[++i] = pow(2.0,log2);

  return ptr;
}

///////////////////////////////////////////////////////////////////////

TwoBodyTwoCenterIntIter::TwoBodyTwoCenterIntIter()
{
}

TwoBodyTwoCenterIntIter::TwoBodyTwoCenterIntIter(const Ref<TwoBodyTwoCenterInt>& o,
                                                 TwoBodyOper::type t) :
  tbi(o), type(t)
{
}

TwoBodyTwoCenterIntIter::~TwoBodyTwoCenterIntIter()
{
}

void
TwoBodyTwoCenterIntIter::start(int ist, int jst, int ien, int jen)
{
  istart=ist;
  jstart=jst;
  iend=ien;
  jend=jen;

  icur=istart;
  jcur=jstart;

  if (!iend) {
    iend=tbi->nshell1();
    jend=tbi->nshell2();
  }

  ij = (icur*(icur+1)>>1) + jcur;
}

void
TwoBodyTwoCenterIntIter::next()
{
  int jlast = (redund) ? std::min(icur,jend-1) : jend-1;

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
TwoBodyTwoCenterIntIter::scale() const
{
  return 1.0;
}

ShellPairIter&
TwoBodyTwoCenterIntIter::current_pair()
{
  tbi->compute_shell(icur,jcur);
  spi.init(tbi->buffer(type), icur, jcur,
           tbi->basis1()->shell_to_function(icur),
           tbi->basis2()->shell_to_function(jcur),
           tbi->basis1()->operator()(icur).nfunction(),
           tbi->basis2()->operator()(jcur).nfunction(),
           redund, scale()
           );

  return spi;
}

bool
TwoBodyTwoCenterIntIter::cloneable() const
{
  return tbi->cloneable();
}

Ref<TwoBodyTwoCenterIntIter>
TwoBodyTwoCenterIntIter::clone()
{
  return new TwoBodyTwoCenterIntIter(tbi->clone(),
                                     type);
}

///////////////////////////////////////////////////////////////////////

TwoBodyTwoCenterIntOp::TwoBodyTwoCenterIntOp(const Ref<TwoBodyTwoCenterInt>& it,
                                             TwoBodyOper::type type)
{
  iter = new TwoBodyTwoCenterIntIter(it, type);
}

TwoBodyTwoCenterIntOp::TwoBodyTwoCenterIntOp(const Ref<TwoBodyTwoCenterIntIter>& it) :
  iter(it)
{
}

TwoBodyTwoCenterIntOp::~TwoBodyTwoCenterIntOp()
{
}

bool
TwoBodyTwoCenterIntOp::cloneable() const
{
  return iter->cloneable();
}

Ref<SCElementOp>
TwoBodyTwoCenterIntOp::clone()
{
  return new TwoBodyTwoCenterIntOp(iter->clone());
}

void
TwoBodyTwoCenterIntOp::process(SCMatrixBlockIter& b)
{
  ExEnv::err0() << indent
       << "TwoBodyTwoCenterIntOp::process: cannot handle generic case\n";
  abort();
}

void
TwoBodyTwoCenterIntOp::process_spec_rect(SCMatrixRectBlock* b)
{
  Ref<GaussianBasisSet> bs1 = iter->two_body_int()->basis1();
  Ref<GaussianBasisSet> bs2 = iter->two_body_int()->basis2();

  // convert basis function indices into shell indices
  int ishstart = bs1->function_to_shell(b->istart);
  int jshstart = bs2->function_to_shell(b->jstart);

  int b1end = b->iend;
  int ishend = (b1end?bs1->function_to_shell(b1end-1) + 1 : 0);

  int b2end = b->jend;
  int jshend = (b2end?bs2->function_to_shell(b2end-1) + 1 : 0);

  int njdata = b->jend - b->jstart;

  iter->set_redundant(0);

  for (iter->start(ishstart,jshstart,ishend,jshend);
       iter->ready(); iter->next()) {
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
TwoBodyTwoCenterIntOp::process_spec_ltri(SCMatrixLTriBlock* b)
{
  Ref<GaussianBasisSet> bs1 = iter->two_body_int()->basis1();

  // convert basis function indices into shell indices
  int fnstart = b->start;
  int fnend = b->end;
  int shstart = bs1->function_to_shell(fnstart);
  int shend = (fnend?bs1->function_to_shell(fnend - 1) + 1 : 0);

  iter->set_redundant(1);

  // loop over all needed shells
  for (iter->start(shstart,shstart,shend,shend); iter->ready(); iter->next()) {
    ShellPairIter& spi = iter->current_pair();

    // compute a set of shell integrals
    for (spi.start(); spi.ready(); spi.next()) {
      int ifn = spi.i();
      int jfn = spi.j();

      if (ifn < fnstart || ifn >= fnend)
        continue;

      int ioff = ifn-fnstart;
      int joff = jfn-fnstart;

      int data_index = i_offset(ioff)+joff;

      b->data[data_index] += spi.val();
    }
  }
}

void
TwoBodyTwoCenterIntOp::process_spec_rectsub(SCMatrixRectSubBlock* b)
{
  Ref<GaussianBasisSet> bs1 = iter->two_body_int()->basis1();
  Ref<GaussianBasisSet> bs2 = iter->two_body_int()->basis2();

  // convert basis function indices into shell indices
  int istart = b->istart;
  int jstart = b->jstart;
  int iend = b->iend;
  int jend = b->jend;

  int ishstart = bs1->function_to_shell(istart);
  int jshstart = bs2->function_to_shell(jstart);

  int ishend = (iend ? bs1->function_to_shell(iend-1) + 1 : 0);
  int jshend = (jend ? bs2->function_to_shell(jend-1) + 1 : 0);

  int njdata = b->istride;

  iter->set_redundant(0);

  for (iter->start(ishstart,jshstart,ishend,jshend);
       iter->ready(); iter->next()) {
    ShellPairIter& spi = iter->current_pair();

    for (spi.start(); spi.ready(); spi.next()) {
      int ifn = spi.i();
      int jfn = spi.j();

      if (ifn < istart || ifn >= iend || jfn < jstart || jfn >= jend)
        continue;

      int data_index = ifn*njdata + jfn;
      b->data[data_index] += spi.val();
    }
  }
}

void
TwoBodyTwoCenterIntOp::process_spec_ltrisub(SCMatrixLTriSubBlock* b)
{
  Ref<GaussianBasisSet> bs1 = iter->two_body_int()->basis1();

  // convert basis function indices into shell indices
  int istart = b->istart;
  int iend = b->iend;

  int jstart = b->jstart;
  int jend = b->jend;

  int ishstart = bs1->function_to_shell(istart);
  int jshstart = bs1->function_to_shell(jstart);

  int ishend = (iend ? bs1->function_to_shell(iend-1) + 1 : 0);
  int jshend = (jend ? bs1->function_to_shell(jend-1) + 1 : 0);

  iter->set_redundant(1);

  // loop over all needed shells
  for (iter->start(ishstart,jshstart,ishend,jshend);
       iter->ready(); iter->next()) {
    ShellPairIter& spi = iter->current_pair();

    // compute a set of shell integrals
    for (spi.start(); spi.ready(); spi.next()) {
      int ifn = spi.i();
      int jfn = spi.j();

      if (ifn < istart || ifn >= iend || jfn < jstart || jfn >= jend)
        continue;

      int data_index = i_offset(ifn)+jfn;
      b->data[data_index] += spi.val();
    }
  }
}

int
TwoBodyTwoCenterIntOp::has_side_effects()
{
  return 1;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
