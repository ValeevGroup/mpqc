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

// //////////////////////////////////////////////////////////////////////////

class Integral;

/** This is an abstract base type for classes that
    compute integrals involving two electrons.
 */
class TwoBodyInt : public RefCount {
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
    
    /// Return the number of basis functions on the given center.
    //@{
    int nbasis1() const;
    int nbasis2() const;
    int nbasis3() const;
    int nbasis4() const;
    //@}

    /// Return the number of shells on center one.
    int nshell() const;
    
    /// Return the number of shells on the given center.
    //@{
    int nshell1() const;
    int nshell2() const;
    int nshell3() const;
    int nshell4() const;
    //@}

    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis();

    /// Return the basis set on the given center.
    //@{
    Ref<GaussianBasisSet> basis1();
    Ref<GaussianBasisSet> basis2();
    Ref<GaussianBasisSet> basis3();
    Ref<GaussianBasisSet> basis4();
    //@}

    /** The computed shell integrals will be put in the buffer returned
        by this member.
    */
    const double * buffer() const;
    
    /** Given for shell indices, this will cause the integral buffer
        to be filled in. */
    virtual void compute_shell(int,int,int,int) = 0;

    /** Return log base 2 of the maximum magnitude of any integral in a
        shell block.  An index of -1 for any argument indicates any shell.  */
    virtual int log2_shell_bound(int= -1,int= -1,int= -1,int= -1) = 0;

    /** If redundant is true, then keep redundant integrals in the buffer.
        The default is true. */
    //@{
    int redundant() const { return redundant_; }
    void set_redundant(int i) { redundant_ = i; }
    //@}

    /// This storage is used to cache computed integrals.
    virtual void set_integral_storage(int storage);

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

    /// Return the number of basis functions on the given center.
    //@{
    int nbasis1() const;
    int nbasis2() const;
    int nbasis3() const;
    int nbasis4() const;
    //@}

    /// Return the number of shells on center one.
    int nshell() const;

    /// Return the number of shells on the given center.
    //@{
    int nshell1() const;
    int nshell2() const;
    int nshell3() const;
    int nshell4() const;
    //@}

    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis();

    /// Return the basis set on the given center.
    //@{
    Ref<GaussianBasisSet> basis1();
    Ref<GaussianBasisSet> basis2();
    Ref<GaussianBasisSet> basis3();
    Ref<GaussianBasisSet> basis4();
    //@}

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
};



#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
