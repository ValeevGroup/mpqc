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

#ifdef __GNUC__
#pragma implementation
#endif

#include <limits.h>

#include <util/class/scexception.h>
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
TwoBodyInt::basis()
{
  return bs1_;
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
TwoBodyInt::buffer(tbint_type i) const
{
  if (i==eri) return buffer_;
  return 0;
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

const Ref<TwoBodyIntTypeDescr>&
TwoBodyInt::inttypedescr(TwoBodyInt::tbint_type type)
{
    static Ref<TwoBodyIntTypeDescr> symm_type = new TwoBodyIntTypeDescr(2,+1,+1,+1);
    static Ref<TwoBodyIntTypeDescr> t1r12_inttype = new TwoBodyIntTypeDescr(2,-1,+1,0);
    static Ref<TwoBodyIntTypeDescr> t2r12_inttype = new TwoBodyIntTypeDescr(2,+1,-1,0);
    switch(type) {
    case TwoBodyInt::r12t1:
    case TwoBodyInt::t1g12:
	return t1r12_inttype;

    case TwoBodyInt::r12t2:
    case TwoBodyInt::t2g12:
	return t2r12_inttype;

    case TwoBodyInt::eri:
    case TwoBodyInt::r12:
    case TwoBodyInt::r12_0_g12:
    case TwoBodyInt::r12_m1_g12:
    case TwoBodyInt::g12t1g12:
    case TwoBodyInt::anti_g12g12:
    case TwoBodyInt::r12_0_gg12:
    case TwoBodyInt::r12_m1_gg12:
    case TwoBodyInt::gg12t1gg12:
	return symm_type;
    }
    throw ProgrammingError("TwoBodyInt::inttype() -- incorrect type");
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
TwoBodyThreeCenterInt::basis()
{
  return bs1_;
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
TwoBodyThreeCenterInt::buffer(tbint_type i) const
{
  if (i==eri) return buffer_;
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
TwoBodyTwoCenterInt::basis()
{
  return bs1_;
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
TwoBodyTwoCenterInt::buffer(tbint_type i) const
{
  if (i==eri) return buffer_;
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

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
