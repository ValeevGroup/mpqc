//
// sointegral.cc --- implementation of the Integral class
//
// Copyright (C) 1998 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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

#include <util/misc/formio.h>
#include <util/class/scexception.h>
#include <chemistry/qc/basis/sointegral.h>

using namespace std;
using namespace sc;

#define DEBUG 0

/////////////////////////////////////////////////////////////////////////////
// OneBodySOInt

OneBodySOInt::OneBodySOInt(const Ref<OneBodyInt> &ob)
{
  ob_ = ob;

  b1_ = new SOBasis(ob->basis1(), ob->integral());

  if (ob->basis2() == ob->basis1()) b2_ = b1_;
  else b2_ = new SOBasis(ob->basis2(), ob->integral());

  only_totally_symmetric_ = 0;

  buffer_ = new double[b1_->max_nfunction_in_shell()
                      *b2_->max_nfunction_in_shell()];
}

OneBodySOInt::~OneBodySOInt()
{
  delete[] buffer_;
}

Ref<SOBasis>
OneBodySOInt::basis() const
{
  return b1_;
}

Ref<SOBasis>
OneBodySOInt::basis1() const
{
  return b1_;
}

Ref<SOBasis>
OneBodySOInt::basis2() const
{
  return b2_;
}

void
OneBodySOInt::compute_shell(int ish, int jsh)
{
  const double *aobuf = ob_->buffer();

  const SOTransform &t1 = b1_->trans(ish);
  const SOTransform &t2 = b2_->trans(jsh);

  int nso1 = b1_->nfunction(ish);
  int nso2 = b2_->nfunction(jsh);

  memset(buffer_, 0, nso1*nso2*sizeof(double));

  int nao2 = b2_->naofunction(jsh);

  // loop through the AO shells that make up this SO shell
  for (int i=0; i<t1.naoshell; i++) {
    const SOTransformShell &s1 = t1.aoshell[i];
    for (int j=0; j<t2.naoshell; j++) {
      const SOTransformShell &s2 = t2.aoshell[j];
      ob_->compute_shell(s1.aoshell, s2.aoshell);
      for (int itr=0; itr<s1.nfunc; itr++) {
        const SOTransformFunction &ifunc = s1.func[itr];
        double icoef = ifunc.coef;
        int iaofunc = ifunc.aofunc;
        int isofunc = b1_->function_offset_within_shell(ish, ifunc.irrep)
          + ifunc.sofunc;
        int iaooff = iaofunc;
        int isooff = isofunc;
        for (int jtr=0; jtr<s2.nfunc; jtr++) {
          const SOTransformFunction &jfunc = s2.func[jtr];
          double jcoef = jfunc.coef * icoef;
          int jaofunc = jfunc.aofunc;
          int jsofunc = b2_->function_offset_within_shell(jsh, jfunc.irrep)
            + jfunc.sofunc;
          int jaooff = iaooff*nao2 + jaofunc;
          int jsooff = isooff*nso2 + jsofunc;
          buffer_[jsooff] += jcoef * aobuf[jaooff];
#if DEBUG
          if (fabs(aobuf[jaooff]*jcoef) > 1.0e-10) {
            ExEnv::outn() <<"("<<isofunc
                 <<"|"<<jsofunc
                 <<") += "<<scprintf("%8.5f",jcoef)<<" * "
                 <<"("<<iaofunc
                 <<"|"<<jaofunc
                 <<"): "
                 <<scprintf("%8.5f -> %8.5f",
                            aobuf[jaooff], buffer_[jsooff])
                 << endl;
          }
#endif
        }
      }
    }
  }
}

void
OneBodySOInt::reinitialize()
{
  ob_->reinitialize();
}

/////////////////////////////////////////////////////////////////////////////
// OneBodySODerivInt

OneBodySODerivInt::OneBodySODerivInt(const Ref<OneBodyDerivInt> &obd)
{
  throw FeatureNotImplemented("Justin, get on it!", __FILE__, __LINE__);
  obd_ = obd;

  b1_ = new SOBasis(obd->basis1(), obd->integral());

  if (obd->basis2() == obd->basis1()) b2_ = b1_;
  else b2_ = new SOBasis(obd->basis2(), obd->integral());

  only_totally_symmetric_ = 0;

  const int num_salc_coords = 12 * obd->basis1()->molecule()->point_group()->order();

  buffer_ = new double[b1_->max_nfunction_in_shell()
                      *b2_->max_nfunction_in_shell()
                      *num_salc_coords];
}

OneBodySODerivInt::~OneBodySODerivInt()
{
  delete[] buffer_;
}

Ref<SOBasis>
OneBodySODerivInt::basis() const
{
  return b1_;
}

Ref<SOBasis>
OneBodySODerivInt::basis1() const
{
  return b1_;
}

Ref<SOBasis>
OneBodySODerivInt::basis2() const
{
  return b2_;
}

void
OneBodySODerivInt::compute_shell(int ish, int jsh)
{
  const double *aobuf = obd_->buffer();

  const SOTransform &t1 = b1_->trans(ish);
  const SOTransform &t2 = b2_->trans(jsh);

  int nso1 = b1_->nfunction(ish);
  int nso2 = b2_->nfunction(jsh);

  //TODO: wrong size here, and in buffer_ allocation
  memset(buffer_, 0, nso1*nso2*sizeof(double));

  int nao2 = b2_->naofunction(jsh);



  // loop through the AO shells that make up this SO shell
  for (int i=0; i<t1.naoshell; i++) {
    const SOTransformShell &s1 = t1.aoshell[i];
    for (int j=0; j<t2.naoshell; j++) {
      const SOTransformShell &s2 = t2.aoshell[j];
      /** Compute the derivative integrals and place the result in the buffer
              returned by buffer(). */
       // TODO: last arg is center on which to compute derivative
       obd_->compute_shell(s1.aoshell, s2.aoshell, 0);
      // obd_->compute_shell(s1.aoshell, s2.aoshell);
      for (int itr=0; itr<s1.nfunc; itr++) {
        const SOTransformFunction &ifunc = s1.func[itr];
        double icoef = ifunc.coef;
        int iaofunc = ifunc.aofunc;
        int isofunc = b1_->function_offset_within_shell(ish, ifunc.irrep)
          + ifunc.sofunc;
        int iaooff = iaofunc;
        int isooff = isofunc;
        for (int jtr=0; jtr<s2.nfunc; jtr++) {
          const SOTransformFunction &jfunc = s2.func[jtr];
          double jcoef = jfunc.coef * icoef;
          int jaofunc = jfunc.aofunc;
          int jsofunc = b2_->function_offset_within_shell(jsh, jfunc.irrep)
            + jfunc.sofunc;
          int jaooff = iaooff*nao2 + jaofunc;
          int jsooff = isooff*nso2 + jsofunc;
          buffer_[jsooff] += jcoef * aobuf[jaooff];
#if DEBUG
          if (fabs(aobuf[jaooff]*jcoef) > 1.0e-10) {
            ExEnv::outn() <<"("<<isofunc
                 <<"|"<<jsofunc
                 <<") += "<<scprintf("%8.5f",jcoef)<<" * "
                 <<"("<<iaofunc
                 <<"|"<<jaofunc
                 <<"): "
                 <<scprintf("%8.5f -> %8.5f",
                            aobuf[jaooff], buffer_[jsooff])
                 << endl;
          }
#endif
        }
      }
    }
  }
}

void
OneBodySODerivInt::reinitialize()
{
  obd_->reinitialize();
}

/////////////////////////////////////////////////////////////////////////////
// TwoBodySOInt

TwoBodySOInt::TwoBodySOInt(const Ref<TwoBodyInt> &tb)
{
  tb_ = tb;

  b1_ = new SOBasis(tb->basis1(), tb->integral());

  if (tb->basis2() == tb->basis1()) b2_ = b1_;
  else b2_ = new SOBasis(tb->basis2(), tb->integral());

  if (tb->basis3() == tb->basis1()) b3_ = b1_;
  else if (tb->basis3() == tb->basis2()) b3_ = b2_;
  else b3_ = new SOBasis(tb->basis3(), tb->integral());

  if (tb->basis4() == tb->basis1()) b4_ = b1_;
  else if (tb->basis4() == tb->basis2()) b4_ = b2_;
  else if (tb->basis4() == tb->basis3()) b4_ = b3_;
  else b4_ = new SOBasis(tb->basis4(), tb->integral());

  redundant_ = 1;
  only_totally_symmetric_ = 0;

  buffer_ = new double[b1_->max_nfunction_in_shell()
                      *b2_->max_nfunction_in_shell()
                      *b3_->max_nfunction_in_shell()
                      *b4_->max_nfunction_in_shell()];
}

TwoBodySOInt::~TwoBodySOInt()
{
  delete[] buffer_;
}

Ref<SOBasis>
TwoBodySOInt::basis() const
{
  return b1_;
}

Ref<SOBasis>
TwoBodySOInt::basis1() const
{
  return b1_;
}

Ref<SOBasis>
TwoBodySOInt::basis2() const
{
  return b2_;
}

Ref<SOBasis>
TwoBodySOInt::basis3() const
{
  return b3_;
}

Ref<SOBasis>
TwoBodySOInt::basis4() const
{
  return b4_;
}

void
TwoBodySOInt::compute_shell(int ish, int jsh, int ksh, int lsh)
{
  tb_->set_redundant(1);
  const double *aobuf = tb_->buffer();

  const SOTransform &t1 = b1_->trans(ish);
  const SOTransform &t2 = b2_->trans(jsh);
  const SOTransform &t3 = b3_->trans(ksh);
  const SOTransform &t4 = b4_->trans(lsh);

  int nso1 = b1_->nfunction(ish);
  int nso2 = b2_->nfunction(jsh);
  int nso3 = b3_->nfunction(ksh);
  int nso4 = b4_->nfunction(lsh);

  memset(buffer_, 0, nso1*nso2*nso3*nso4*sizeof(double));

  int nao2 = b2_->naofunction(jsh);
  int nao3 = b3_->naofunction(ksh);
  int nao4 = b4_->naofunction(lsh);

  // loop through the ao shells that make up this so shell
  for (int i=0; i<t1.naoshell; i++) {
    const SOTransformShell &s1 = t1.aoshell[i];
    for (int j=0; j<t2.naoshell; j++) {
      const SOTransformShell &s2 = t2.aoshell[j];
      for (int k=0; k<t3.naoshell; k++) {
        const SOTransformShell &s3 = t3.aoshell[k];
        for (int l=0; l<t4.naoshell; l++) {
          const SOTransformShell &s4 = t4.aoshell[l];
          tb_->compute_shell(s1.aoshell, s2.aoshell, s3.aoshell, s4.aoshell);
          for (int itr=0; itr<s1.nfunc; itr++) {
            const SOTransformFunction &ifunc = s1.func[itr];
            double icoef = ifunc.coef;
            int iaofunc = ifunc.aofunc;
            int isofunc = b1_->function_offset_within_shell(ish,
                                                            ifunc.irrep)
                          + ifunc.sofunc;
            int iaooff = iaofunc;
            int isooff = isofunc;
            for (int jtr=0; jtr<s2.nfunc; jtr++) {
              const SOTransformFunction &jfunc = s2.func[jtr];
              double jcoef = jfunc.coef * icoef;
              int jaofunc = jfunc.aofunc;
              int jsofunc = b2_->function_offset_within_shell(jsh,
                                                              jfunc.irrep)
                            + jfunc.sofunc;
              int jaooff = iaooff*nao2 + jaofunc;
              int jsooff = isooff*nso2 + jsofunc;
              for (int ktr=0; ktr<s3.nfunc; ktr++) {
                const SOTransformFunction &kfunc = s3.func[ktr];
                double kcoef = kfunc.coef * jcoef;
                int kaofunc = kfunc.aofunc;
                int ksofunc = b3_->function_offset_within_shell(ksh,
                                                                kfunc.irrep)
                              + kfunc.sofunc;
                int kaooff = jaooff*nao3 + kaofunc;
                int ksooff = jsooff*nso3 + ksofunc;
                for (int ltr=0; ltr<s4.nfunc; ltr++) {
                  const SOTransformFunction &lfunc = s4.func[ltr];
                  double lcoef = lfunc.coef * kcoef;
                  int laofunc = lfunc.aofunc;
                  int lsofunc = b4_->function_offset_within_shell(lsh,
                                                                  lfunc.irrep)
                                + lfunc.sofunc;
                  int laooff = kaooff*nao4 + laofunc;
                  int lsooff = ksooff*nso4 + lsofunc;
                  buffer_[lsooff] += lcoef * aobuf[laooff];
#if DEBUG
                  if (fabs(aobuf[laooff]*lcoef) > 1.0e-10) {
                      ExEnv::outn() <<"("<<isofunc<<jsofunc
                           <<"|"<<ksofunc<<lsofunc
                           <<") += "<<scprintf("%8.5f",lcoef)<<" * "
                           <<"("<<iaofunc<<jaofunc
                           <<"|"<<kaofunc<<laofunc
                           <<"): "
                           <<scprintf("%8.5f -> %8.5f",
                                      aobuf[laooff], buffer_[lsooff])
                           << endl;
                    }
#endif
                  }
                }
              }
            }
          }
        }
      }
    }
}

/////////////////////////////////////////////////////////////////////////////
// TwoBodySODerivInt

TwoBodySODerivInt::TwoBodySODerivInt(const Ref<TwoBodyDerivInt> &tb)
{
  throw FeatureNotImplemented("Justin, get on it!", __FILE__, __LINE__);
  tb_ = tb;

  b1_ = new SOBasis(tb->basis1(), tb->integral());

  if (tb->basis2() == tb->basis1()) b2_ = b1_;
  else b2_ = new SOBasis(tb->basis2(), tb->integral());

  if (tb->basis3() == tb->basis1()) b3_ = b1_;
  else if (tb->basis3() == tb->basis2()) b3_ = b2_;
  else b3_ = new SOBasis(tb->basis3(), tb->integral());

  if (tb->basis4() == tb->basis1()) b4_ = b1_;
  else if (tb->basis4() == tb->basis2()) b4_ = b2_;
  else if (tb->basis4() == tb->basis3()) b4_ = b3_;
  else b4_ = new SOBasis(tb->basis4(), tb->integral());

  redundant_ = 1;
  only_totally_symmetric_ = 0;

  const int num_salc_coords = 12 * tb->basis1()->molecule()->point_group()->order();

  buffer_ = new double[b1_->max_nfunction_in_shell()
                      *b2_->max_nfunction_in_shell()
                      *b3_->max_nfunction_in_shell()
                      *b4_->max_nfunction_in_shell()
                      * num_salc_coords ];
}

TwoBodySODerivInt::~TwoBodySODerivInt()
{
  delete[] buffer_;
}

Ref<SOBasis>
TwoBodySODerivInt::basis() const
{
  return b1_;
}

Ref<SOBasis>
TwoBodySODerivInt::basis1() const
{
  return b1_;
}

Ref<SOBasis>
TwoBodySODerivInt::basis2() const
{
  return b2_;
}

Ref<SOBasis>
TwoBodySODerivInt::basis3() const
{
  return b3_;
}

Ref<SOBasis>
TwoBodySODerivInt::basis4() const
{
  return b4_;
}

void
TwoBodySODerivInt::compute_shell(int ish, int jsh, int ksh, int lsh)
{
  // TODO figure out if we need this here tb_->set_redundant(1);
  const double *aobuf = tb_->buffer();

  const SOTransform &t1 = b1_->trans(ish);
  const SOTransform &t2 = b2_->trans(jsh);
  const SOTransform &t3 = b3_->trans(ksh);
  const SOTransform &t4 = b4_->trans(lsh);

  int nso1 = b1_->nfunction(ish);
  int nso2 = b2_->nfunction(jsh);
  int nso3 = b3_->nfunction(ksh);
  int nso4 = b4_->nfunction(lsh);

  memset(buffer_, 0, nso1*nso2*nso3*nso4*sizeof(double));

  int nao2 = b2_->naofunction(jsh);
  int nao3 = b3_->naofunction(ksh);
  int nao4 = b4_->naofunction(lsh);

  // loop through the ao shells that make up this so shell
  for (int i=0; i<t1.naoshell; i++) {
    const SOTransformShell &s1 = t1.aoshell[i];
    for (int j=0; j<t2.naoshell; j++) {
      const SOTransformShell &s2 = t2.aoshell[j];
      for (int k=0; k<t3.naoshell; k++) {
        const SOTransformShell &s3 = t3.aoshell[k];
        for (int l=0; l<t4.naoshell; l++) {
          const SOTransformShell &s4 = t4.aoshell[l];
          //tb_->compute_shell(s1.aoshell, s2.aoshell, s3.aoshell, s4.aoshell);
          for (int itr=0; itr<s1.nfunc; itr++) {
            const SOTransformFunction &ifunc = s1.func[itr];
            double icoef = ifunc.coef;
            int iaofunc = ifunc.aofunc;
            int isofunc = b1_->function_offset_within_shell(ish,
                                                            ifunc.irrep)
                          + ifunc.sofunc;
            int iaooff = iaofunc;
            int isooff = isofunc;
            for (int jtr=0; jtr<s2.nfunc; jtr++) {
              const SOTransformFunction &jfunc = s2.func[jtr];
              double jcoef = jfunc.coef * icoef;
              int jaofunc = jfunc.aofunc;
              int jsofunc = b2_->function_offset_within_shell(jsh,
                                                              jfunc.irrep)
                            + jfunc.sofunc;
              int jaooff = iaooff*nao2 + jaofunc;
              int jsooff = isooff*nso2 + jsofunc;
              for (int ktr=0; ktr<s3.nfunc; ktr++) {
                const SOTransformFunction &kfunc = s3.func[ktr];
                double kcoef = kfunc.coef * jcoef;
                int kaofunc = kfunc.aofunc;
                int ksofunc = b3_->function_offset_within_shell(ksh,
                                                                kfunc.irrep)
                              + kfunc.sofunc;
                int kaooff = jaooff*nao3 + kaofunc;
                int ksooff = jsooff*nso3 + ksofunc;
                for (int ltr=0; ltr<s4.nfunc; ltr++) {
                  const SOTransformFunction &lfunc = s4.func[ltr];
                  double lcoef = lfunc.coef * kcoef;
                  int laofunc = lfunc.aofunc;
                  int lsofunc = b4_->function_offset_within_shell(lsh,
                                                                  lfunc.irrep)
                                + lfunc.sofunc;
                  int laooff = kaooff*nao4 + laofunc;
                  int lsooff = ksooff*nso4 + lsofunc;
                  buffer_[lsooff] += lcoef * aobuf[laooff];
#if DEBUG
                  if (fabs(aobuf[laooff]*lcoef) > 1.0e-10) {
                      ExEnv::outn() <<"("<<isofunc<<jsofunc
                           <<"|"<<ksofunc<<lsofunc
                           <<") += "<<scprintf("%8.5f",lcoef)<<" * "
                           <<"("<<iaofunc<<jaofunc
                           <<"|"<<kaofunc<<laofunc
                           <<"): "
                           <<scprintf("%8.5f -> %8.5f",
                                      aobuf[laooff], buffer_[lsooff])
                           << endl;
                    }
#endif
                  }
                }
              }
            }
          }
        }
      }
    }
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
