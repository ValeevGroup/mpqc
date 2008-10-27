//
// tbint.h
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

#ifndef _chemistry_qc_basis_tbint_h
#define _chemistry_qc_basis_tbint_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/ref/ref.h>
#include <util/group/message.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/dercent.h>
#include <chemistry/qc/basis/inttypedescr.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

class Integral;

/** This is an abstract base type for classes that
    compute integrals involving two electrons and 2 functions per electron.
 */
class TwoBodyInt : public RefCount {

  private:
    double *log2_to_double_;
 
  protected:
    // this is who created me
    Integral *integral_;

    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;
    Ref<GaussianBasisSet> bs3_;
    Ref<GaussianBasisSet> bs4_;

    double *buffer_;

    int redundant_;
    
    TwoBodyInt(Integral *integral,
               const Ref<GaussianBasisSet>&bs1,
               const Ref<GaussianBasisSet>&bs2,
               const Ref<GaussianBasisSet>&bs3,
               const Ref<GaussianBasisSet>&bs4);
  public:
    virtual ~TwoBodyInt();
  
    /// Return the number of basis functions on center one.
    int nbasis() const;
    
    /// Return the number of basis functions on center one.
    int nbasis1() const;
    /// Return the number of basis functions on center two.
    int nbasis2() const;
    /// Return the number of basis functions on center three.
    int nbasis3() const;
    /// Return the number of basis functions on center four.
    int nbasis4() const;

    /// Return the number of shells on center one.
    int nshell() const;
    
    /// Return the number of shells on center one.
    int nshell1() const;
    /// Return the number of shells on center two.
    int nshell2() const;
    /// Return the number of shells on center three.
    int nshell3() const;
    /// Return the number of shells on center four.
    int nshell4() const;

    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis();

    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis1();
    /// Return the basis set on center two.
    Ref<GaussianBasisSet> basis2();
    /// Return the basis set on center three.
    Ref<GaussianBasisSet> basis3();
    /// Return the basis set on center four.
    Ref<GaussianBasisSet> basis4();

  /** Types of two-body integrals that TwoBodyInt understands:
      eri -- integral of \f$r_{12}^{-1}\f$,
      r12 -- integral of \f$r_{12}\f$,
      r12ti -- integral of \f$[r_{12},\hat{T}_i]\f$,
      r12_0_g12 -- integral of \f$ g_{12}=\exp(-\gamma r_{12}^2) \f$,
      r12_m1_g12 -- integral of \f$g_{12}/r_{12}\f$,
      tig12 -- integral of \f$[\hat{T}_i,g_{12}]\f$,
      g12t1g12 -- integral of \f$[g_{12},[\hat{T}_1,g_{12}]]\f$,
      g12p4g12_m_g12t1g12t1 -- integral of
         \f$[g_{12}, [\hat{p}^4_1 + \hat{p}^4_2, g_{12}]] - 2 [g_{12}, [\hat{T}_1 + \hat{T}_2, g_{12}]](\hat{T}_1 + \hat{T}_2)\f$,
      anti_g12g12 -- integral of 
       */
    enum tbint_type { eri =0, r12 =1, r12t1 =2, r12t2 =3,
                      r12_0_g12 =4, r12_m1_g12 =5, t1g12 =6, t2g12 =7,
                      g12t1g12 =8, g12p4g12_m_g12t1g12t1 =9, anti_g12g12 =10,
                      r12_0_gg12 =11, r12_m1_gg12 =12, gg12t1gg12 =13};
    /// The max number of such types
    static const int max_num_tbint_types = 14;
    /// The number of types supported by this TwoBodyInt instance
    virtual unsigned int num_tbint_types() const =0;
    /// Maps integral type t to its index
    virtual unsigned int inttype(tbint_type t) const =0;
    /// Maps integral type index t to type
    virtual TwoBodyInt::tbint_type inttype(unsigned int t) const =0;
    /// Returns a descriptor for integral type t
    static const Ref<TwoBodyIntTypeDescr>& inttypedescr(tbint_type t);

    /** The computed shell integrals will be put in the buffer returned
        by this member.  Some TwoBodyInt specializations have more than
	one buffer:  The type arguments selects which buffer is returned.
	If the requested type is not supported, then 0 is returned. */
    virtual const double * buffer(tbint_type type = eri) const;

    /** Given four shell indices, integrals will be computed and placed in
        the buffer.  The first two indices correspond to electron 1 and the
        second two indices correspond to electron 2.*/
    virtual void compute_shell(int,int,int,int) = 0;

    /** Given four shell indices, supported two body integral types
        are computed and returned.  The first two indices correspond
        to electron 1 and the second two indices correspond to
        electron 2. This is used in the python interface where the
        return type is automatically converted to a map of numpy
        arrays. */
    std::pair<std::map<tbint_type,const double*>,unsigned long[4]>
    compute_shell_arrays(int,int,int,int);

    /** Return log base 2 of the maximum magnitude of any integral in a
        shell block obtained from compute_shell.  An index of -1 for any
        argument indicates any shell.  */
    virtual int log2_shell_bound(int= -1,int= -1,int= -1,int= -1) = 0;

    /** Return the maximum magnitude of any integral in a
        shell block obtained from compute_shell.  An index of -1 for any
        argument indicates any shell.  */
    double shell_bound(int= -1,int= -1,int= -1,int= -1);

    /** If redundant is true, then keep redundant integrals in the buffer.
        The default is true. */
    virtual int redundant() const { return redundant_; }
    /// See redundant().
    virtual void set_redundant(int i) { redundant_ = i; }

    /// This storage is used to cache computed integrals.
    virtual void set_integral_storage(size_t storage);

    /// Return the integral factory that was used to create this object.
    Integral *integral() const { return integral_; }

};

// //////////////////////////////////////////////////////////////////////////

/** This is an abstract base type for classes that compute integrals
    involving two electrons in three Gaussian functions.
 */
class TwoBodyThreeCenterInt : public RefCount {

  private:
    double *log2_to_double_;

  protected:
    // this is who created me
    Integral *integral_;

    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;
    Ref<GaussianBasisSet> bs3_;

    double *buffer_;

    int redundant_;
    
    TwoBodyThreeCenterInt(Integral *integral,
                          const Ref<GaussianBasisSet>&bs1,
                          const Ref<GaussianBasisSet>&bs2,
                          const Ref<GaussianBasisSet>&bs3);
  public:
    virtual ~TwoBodyThreeCenterInt();
  
    /// Return the number of basis functions on center one.
    int nbasis() const;
    
    /// Return the number of basis functions on center one.
    int nbasis1() const;
    /// Return the number of basis functions on center two.
    int nbasis2() const;
    /// Return the number of basis functions on center three.
    int nbasis3() const;

    /// Return the number of shells on center one.
    int nshell() const;
    
    /// Return the number of shells on center one.
    int nshell1() const;
    /// Return the number of shells on center two.
    int nshell2() const;
    /// Return the number of shells on center three.
    int nshell3() const;

    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis();

    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis1();
    /// Return the basis set on center two.
    Ref<GaussianBasisSet> basis2();
    /// Return the basis set on center three.
    Ref<GaussianBasisSet> basis3();

    /// Types of two-body integrals that this evaluator understands
    typedef TwoBodyInt::tbint_type tbint_type;
    /// The total number of such types
    static const int max_num_tbint_types = TwoBodyInt::max_num_tbint_types;

    /** The computed shell integrals will be put in the buffer returned
        by this member.  Some TwoBodyInt specializations have more than
	one buffer:  The type arguments selects which buffer is returned.
	If the requested type is not supported, then 0 is returned. */
    virtual const double * buffer(tbint_type type = TwoBodyInt::eri) const;

    /** Given three shell indices, integrals will be computed and placed in
        the buffer.  The first two indices correspond to electron 1 and the
        second index corresponds to electron 2.*/
    virtual void compute_shell(int,int,int) = 0;

    /** Return log base 2 of the maximum magnitude of any integral in a
        shell block obtained from compute_shell.  An index of -1 for any
        argument indicates any shell.  */
    virtual int log2_shell_bound(int= -1,int= -1,int= -1) = 0;

    /** Return the maximum magnitude of any integral in a
        shell block obtained from compute_shell.  An index of -1 for any
        argument indicates any shell.  */
    double shell_bound(int= -1,int= -1,int= -1);

    /** If redundant is true, then keep redundant integrals in the buffer.
        The default is true. */
    int redundant() const { return redundant_; }
    /// See redundant().
    void set_redundant(int i) { redundant_ = i; }

    /// This storage is used to cache computed integrals.
    virtual void set_integral_storage(size_t storage);

    /// Return the integral factory that was used to create this object.
    Integral *integral() const { return integral_; }

};

// //////////////////////////////////////////////////////////////////////////

/** This is an abstract base type for classes that
    compute integrals involving two electrons in two
    Gaussian functions.
 */
class TwoBodyTwoCenterInt : public RefCount {

  private:
    double *log2_to_double_;

  protected:
    // this is who created me
    Integral *integral_;

    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;

    double *buffer_;

    int redundant_;
    
    TwoBodyTwoCenterInt(Integral *integral,
                        const Ref<GaussianBasisSet>&bs1,
                        const Ref<GaussianBasisSet>&bs2);
  public:
    virtual ~TwoBodyTwoCenterInt();
  
    /// Return the number of basis functions on center one.
    int nbasis() const;
    
    /// Return the number of basis functions on center one.
    int nbasis1() const;
    /// Return the number of basis functions on center two.
    int nbasis2() const;

    /// Return the number of shells on center one.
    int nshell() const;
    
    /// Return the number of shells on center one.
    int nshell1() const;
    /// Return the number of shells on center two.
    int nshell2() const;

    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis();

    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis1();
    /// Return the basis set on center two.
    Ref<GaussianBasisSet> basis2();

    /// Types of two-body integrals that this evaluator understands
    typedef TwoBodyInt::tbint_type tbint_type;
    /// The total number of such types
    static const int max_num_tbint_types = TwoBodyInt::max_num_tbint_types;

    /** The computed shell integrals will be put in the buffer returned
        by this member.  Some TwoBodyInt specializations have more than
	one buffer:  The type arguments selects which buffer is returned.
	If the requested type is not supported, then 0 is returned. */
    virtual const double * buffer(tbint_type type = TwoBodyInt::eri) const;

    /** Given four shell indices, integrals will be computed and placed in
        the buffer.  The first index corresponds to electron 1 and the
        second index corresponds to electron 2.*/
    virtual void compute_shell(int,int) = 0;

    /** Return log base 2 of the maximum magnitude of any integral in a
        shell block obtained from compute_shell.  An index of -1 for any
        argument indicates any shell.  */
    virtual int log2_shell_bound(int= -1,int= -1) = 0;

    /** Return the maximum magnitude (as a double) of any integral in a
        shell block obtained from compute_shell.  An index of -1 for any
        argument indicates any shell.  */
    double shell_bound(int= -1,int= -1);

    /** If redundant is true, then keep redundant integrals in the buffer.
        The default is true. */
    int redundant() const { return redundant_; }
    /// See redundant().
    void set_redundant(int i) { redundant_ = i; }

    /// This storage is used to cache computed integrals.
    virtual void set_integral_storage(size_t storage);

    /// Return the integral factory that was used to create this object.
    Integral *integral() const { return integral_; }

};

// //////////////////////////////////////////////////////////////////////////

class ShellQuartetIter {
  protected:
    const double * buf;
    double scale_;

    int redund_;
    
    int e12;
    int e34;
    int e13e24;

    int index;
    
    int istart;
    int jstart;
    int kstart;
    int lstart;

    int iend;
    int jend;
    int kend;
    int lend;

    int icur;
    int jcur;
    int kcur;
    int lcur;

    int i_;
    int j_;
    int k_;
    int l_;
    
  public:
    ShellQuartetIter();
    virtual ~ShellQuartetIter();

    virtual void init(const double *,
                      int, int, int, int,
                      int, int, int, int,
                      int, int, int, int,
                      double, int);

    virtual void start();
    virtual void next();

    int ready() const { return icur < iend; }

    int i() const { return i_; }
    int j() const { return j_; }
    int k() const { return k_; }
    int l() const { return l_; }

    int nint() const { return iend*jend*kend*lend; }
    
    double val() const { return buf[index]*scale_; }
};

class TwoBodyIntIter {
  protected:
    Ref<TwoBodyInt> tbi;
    ShellQuartetIter sqi;
    
    int iend;
    
    int icur;
    int jcur;
    int kcur;
    int lcur;
    
  public:
    TwoBodyIntIter();
    TwoBodyIntIter(const Ref<TwoBodyInt>&);

    virtual ~TwoBodyIntIter();
    
    virtual void start();
    virtual void next();

    int ready() const { return (icur < iend); }

    int ishell() const { return icur; }
    int jshell() const { return jcur; }
    int kshell() const { return kcur; }
    int lshell() const { return lcur; }

    virtual double scale() const;

    ShellQuartetIter& current_quartet();
};

// //////////////////////////////////////////////////////////////////////////

/** This is an abstract base type for classes that
    compute integrals involving two electrons.
 */
class TwoBodyDerivInt : public RefCount {

  private:
    double *log2_to_double_;

  protected:
    // this is who created me
    Integral *integral_;

    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;
    Ref<GaussianBasisSet> bs3_;
    Ref<GaussianBasisSet> bs4_;

    double *buffer_;

    TwoBodyDerivInt(Integral* integral,
                    const Ref<GaussianBasisSet>&b1,
                    const Ref<GaussianBasisSet>&b2,
                    const Ref<GaussianBasisSet>&b3,
                    const Ref<GaussianBasisSet>&b4);
  public:
    virtual ~TwoBodyDerivInt();
  
    /// Return the number of basis functions on center one.
    int nbasis() const;

    /// Return the number of basis functions on center one.
    int nbasis1() const;
    /// Return the number of basis functions on center two.
    int nbasis2() const;
    /// Return the number of basis functions on center three.
    int nbasis3() const;
    /// Return the number of basis functions on center four.
    int nbasis4() const;

    /// Return the number of shells on center one.
    int nshell() const;

    /// Return the number of shells on center one.
    int nshell1() const;
    /// Return the number of shells on center two.
    int nshell2() const;
    /// Return the number of shells on center three.
    int nshell3() const;
    /// Return the number of shells on center four.
    int nshell4() const;

    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis();

    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis1();
    /// Return the basis set on center two.
    Ref<GaussianBasisSet> basis2();
    /// Return the basis set on center three.
    Ref<GaussianBasisSet> basis3();
    /// Return the basis set on center four.
    Ref<GaussianBasisSet> basis4();

    /** The computed shell integrals will be put in the buffer returned
        by this member.
    */
    const double * buffer() const;
    
    /** Given for shell indices, this will cause the integral buffer
        to be filled in. */
    virtual void compute_shell(int,int,int,int,DerivCenters&) = 0;

    /** Return log base 2 of the maximum magnitude of any integral in a
        shell block.  An index of -1 for any argument indicates any shell.  */
    virtual int log2_shell_bound(int= -1,int= -1,int= -1,int= -1) = 0;

    /** Return the maximum magnitude of any integral in a
        shell block.  An index of -1 for any argument indicates any shell.  */
    double shell_bound(int= -1,int= -1,int= -1,int= -1);
};

// //////////////////////////////////////////////////////////////////////////

/** This is an abstract base type for classes that
    compute three centers integrals involving two electrons.
 */
class TwoBodyThreeCenterDerivInt : public RefCount {

  private:
    double *log2_to_double_;

  protected:
    // this is who created me
    Integral *integral_;

    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;
    Ref<GaussianBasisSet> bs3_;

    double *buffer_;

    TwoBodyThreeCenterDerivInt(Integral* integral,
                    const Ref<GaussianBasisSet>&b1,
                    const Ref<GaussianBasisSet>&b2,
                    const Ref<GaussianBasisSet>&b3);
  public:
    virtual ~TwoBodyThreeCenterDerivInt();
  
    /// Return the number of basis functions on center one.
    int nbasis() const;

    /// Return the number of basis functions on center one.
    int nbasis1() const;
    /// Return the number of basis functions on center two.
    int nbasis2() const;
    /// Return the number of basis functions on center three.
    int nbasis3() const;

    /// Return the number of shells on center one.
    int nshell() const;

    /// Return the number of shells on center one.
    int nshell1() const;
    /// Return the number of shells on center two.
    int nshell2() const;
    /// Return the number of shells on center three.
    int nshell3() const;

    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis();

    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis1();
    /// Return the basis set on center two.
    Ref<GaussianBasisSet> basis2();
    /// Return the basis set on center three.
    Ref<GaussianBasisSet> basis3();

    /** The computed shell integrals will be put in the buffer returned
        by this member.
    */
    const double * buffer() const;
    
    /** Given for shell indices, this will cause the integral buffer
        to be filled in. */
    virtual void compute_shell(int,int,int,DerivCenters&) = 0;

    /** Return log base 2 of the maximum magnitude of any integral in a
        shell block.  An index of -1 for any argument indicates any shell.  */
    virtual int log2_shell_bound(int= -1,int= -1,int= -1) = 0;

    /** Return the maximum magnitude of any integral in a
        shell block.  An index of -1 for any argument indicates any shell.  */
    double shell_bound(int= -1,int= -1,int= -1);

};

// //////////////////////////////////////////////////////////////////////////

/** This is an abstract base type for classes that
    compute two centers integrals involving two electrons.
 */
class TwoBodyTwoCenterDerivInt : public RefCount {

  private:
    double *log2_to_double_;

  protected:
    // this is who created me
    Integral *integral_;

    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;

    double *buffer_;

    TwoBodyTwoCenterDerivInt(Integral* integral,
                    const Ref<GaussianBasisSet>&b1,
                    const Ref<GaussianBasisSet>&b2);
  public:
    virtual ~TwoBodyTwoCenterDerivInt();
  
    /// Return the number of basis functions on center one.
    int nbasis() const;

    /// Return the number of basis functions on center one.
    int nbasis1() const;
    /// Return the number of basis functions on center two.
    int nbasis2() const;

    /// Return the number of shells on center one.
    int nshell() const;

    /// Return the number of shells on center one.
    int nshell1() const;
    /// Return the number of shells on center two.
    int nshell2() const;

    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis();

    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis1();
    /// Return the basis set on center two.
    Ref<GaussianBasisSet> basis2();

    /** The computed shell integrals will be put in the buffer returned
        by this member.
    */
    const double * buffer() const;
    
    /** Given for shell indices, this will cause the integral buffer
        to be filled in. */
    virtual void compute_shell(int,int,DerivCenters&) = 0;

    /** Return log base 2 of the maximum magnitude of any integral in a
        shell block.  An index of -1 for any argument indicates any shell.  */
    virtual int log2_shell_bound(int= -1,int= -1) = 0;

    /** Return the maximum magnitude of any integral in a
        shell block.  An index of -1 for any argument indicates any shell.  */
    double shell_bound(int= -1,int= -1);

};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
