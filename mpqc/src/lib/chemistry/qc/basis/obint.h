//
// obint.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
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

#ifndef _chemistry_qc_basis_obint_h
#define _chemistry_qc_basis_obint_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/ref/ref.h>
#include <util/state/state.h>
#include <math/scmat/matrix.h>
#include <math/scmat/elemop.h>

#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/dercent.h>

namespace sc {

class Integral;

// //////////////////////////////////////////////////////////////////////////

class EfieldDotVectorData: public RefCount
{
  public:
    EfieldDotVectorData() {};
    ~EfieldDotVectorData();

    double position[3];
    double vector[3];

    void set_position(double*);
    void set_vector(double*);
};


class DipoleData: public RefCount
{
  public:
    double origin[3];

    DipoleData(double *d) {origin[0]=d[0]; origin[1]=d[1]; origin[2]=d[2];}
    DipoleData() {origin[0]=origin[1]=origin[2]=0.0;}
    ~DipoleData();
    void set_origin(double*);
};


class PointChargeData: public RefCount
{
  private:
    int ncharges_;
    const double *charges_;
    const double *const*positions_;
    double *alloced_charges_;
    double **alloced_positions_;

  public:
    // If copy_data is 0, the passed positions and charges will
    // be stored (but not freed).
    PointChargeData(int ncharge,
                    const double *const*positions, const double *charges,
                    int copy_data = 0);
    ~PointChargeData();

    int ncharges() const { return ncharges_; }
    const double *charges() const { return charges_; }
    const double *const*positions() const { return positions_; }
};


/** OneBodyInt is an abstract base class for objects that
    compute integrals between two basis functions. */
class OneBodyInt : public RefCount {
  protected:
    // this is who created me
    Integral *integral_;

    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;

    double *buffer_;

    OneBodyInt(Integral *integral,
               const Ref<GaussianBasisSet>&b1,
               const Ref<GaussianBasisSet>&b2 = 0);

  public:
    virtual ~OneBodyInt();
  
    /// Returns the number of basis functions on center one.
    int nbasis() const;

    /// Returns the number of basis functions on the center one.
    int nbasis1() const;
    /// Returns the number of basis functions on the center two.
    int nbasis2() const;

    /// Return the number of shells on center one.
    int nshell() const;

    /// Return the number of shells on the center one.
    int nshell1() const;
    /// Return the number of shells on the center two.
    int nshell2() const;

    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis();

    /// Return the basis set on the center one.
    Ref<GaussianBasisSet> basis1();
    /// Return the basis set on the center two.
    Ref<GaussianBasisSet> basis2();

    /// Returns the buffer where the integrals are placed.
    const double * buffer() const;
    
    /** Computes the integrals between basis functions in the
        given shell pair. */
    virtual void compute_shell(int,int) = 0;

    /** This is called for one body integrals that take data to let
        them know that the data they reference has changed. */
    virtual void reinitialize();

    /** Return true if the clone member can be called.  The default
     * implementation returns false. */
    virtual bool cloneable();

    /** Returns a clone of this.  The default implementation throws an
     * exception. */
    virtual Ref<OneBodyInt> clone();

    Integral *integral() const { return integral_; }
};

// //////////////////////////////////////////////////////////////////////////

/** OneBodyOneCenterInt is an abstract base class for objects that
    compute integrals between two basis functions. */
class OneBodyOneCenterInt : public RefCount {
  protected:
    // this is who created me
    Integral *integral_;

    Ref<GaussianBasisSet> bs1_;

    double *buffer_;

    OneBodyOneCenterInt(Integral *integral,
               const Ref<GaussianBasisSet>&b1);

  public:
    virtual ~OneBodyOneCenterInt();
  
    /// Returns the number of basis functions on center one.
    int nbasis() const;

    /// Returns the number of basis functions on the center one.
    int nbasis1() const;

    /// Return the number of shells on center one.
    int nshell() const;

    /// Return the number of shells on the center one.
    int nshell1() const;

    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis();

    /// Return the basis set on the center one.
    Ref<GaussianBasisSet> basis1();

    /// Returns the buffer where the integrals are placed.
    const double * buffer() const;
    
    /** Computes the integrals for basis functions on the
        given shell. */
    virtual void compute_shell(int) = 0;

    /** This is called for one body integrals that take data to let
        them know that the data they reference has changed. */
    virtual void reinitialize();

    /** Return true if the clone member can be called.  The default
     * implementation returns false. */
    virtual bool cloneable();

    /** Returns a clone of this.  The default implementation throws an
     * exception. */
    virtual Ref<OneBodyOneCenterInt> clone();

    Integral *integral() const { return integral_; }
};

// //////////////////////////////////////////////////////////////////////////

class OneBodyOneCenterWrapper : public OneBodyOneCenterInt {
    Ref<OneBodyInt> ob_;
    int jsh_;
  public:
    OneBodyOneCenterWrapper(const Ref<OneBodyInt>& ob,
                            int sh2 = 0);
    void compute_shell(int);
};

// //////////////////////////////////////////////////////////////////////////

class ShellPairIter {
  private:
    const double * buf;
    double scale_;

    int e12;

    int index;
    
    int ioffset;
    int joffset;

    int iend;
    int jend;

    int icur;
    int jcur;
    
  public:
    ShellPairIter();
    ~ShellPairIter();

    void init(const double * buffer, int ishell, int jshell,
              int ioff, int joff, int nfunci, int nfuncj, int redund=0,
              double scale=1.0);

    void start() { icur=jcur=index=0; }
    int ready() const { return (icur < iend); }

    void next() {
      if (jcur < ((e12)?(icur):((jend)-1))) {
        index++;
        jcur++;
        return;
      }

      jcur=0;
      icur++;

      index = icur*jend;
    }

    int current_i() const { return icur; }
    int current_j() const { return jcur; }

    int i() const { return icur+ioffset; }
    int j() const { return jcur+joffset; }

    int nint() const { return iend*jend; }
    
    double val() const { return buf[index]*scale_; }
};

// //////////////////////////////////////////////////////////////////////////

class OneBodyIntIter : public RefCount {
  protected:
    Ref<OneBodyInt> obi; // help me obi wan
    ShellPairIter spi;
    
    int redund;
    
    int istart;
    int jstart;
    
    int iend;
    int jend;

    int icur;
    int jcur;

    int ij;
    
  public:
    OneBodyIntIter();
    OneBodyIntIter(const Ref<OneBodyInt>&);
    virtual ~OneBodyIntIter();
    
    virtual void start(int ist=0, int jst=0, int ien=0, int jen=0);
    virtual void next();

    int ready() const { return (icur < iend); }

    int ishell() const { return icur; }
    int jshell() const { return jcur; }

    int ijshell() const { return ij; }

    int redundant() const { return redund; }
    void set_redundant(int i) { redund=i; }
    
    virtual double scale() const;

    Ref<OneBodyInt> one_body_int() { return obi; }

    ShellPairIter& current_pair();

    virtual bool cloneable();
    virtual Ref<OneBodyIntIter> clone();
};



// //////////////////////////////////////////////////////////////////////////

class OneBodyIntOp: public SCElementOp {
  protected:
    Ref<OneBodyIntIter> iter;

  public:
    OneBodyIntOp(const Ref<OneBodyInt>&);
    OneBodyIntOp(const Ref<OneBodyIntIter>&);
    virtual ~OneBodyIntOp();
  
    void process(SCMatrixBlockIter&);
    void process_spec_rect(SCMatrixRectBlock*);
    void process_spec_ltri(SCMatrixLTriBlock*);
    void process_spec_rectsub(SCMatrixRectSubBlock*);
    void process_spec_ltrisub(SCMatrixLTriSubBlock*);

    bool cloneable();
    Ref<SCElementOp> clone();

    int has_side_effects();
};

class OneBody3IntOp: public SCElementOp3 {
  private:
    Ref<OneBodyIntIter> iter;

  public:
    OneBody3IntOp(const Ref<OneBodyInt>&b);
    OneBody3IntOp(const Ref<OneBodyIntIter>&);
    virtual ~OneBody3IntOp();
  
    void process(SCMatrixBlockIter&,
                 SCMatrixBlockIter&,
                 SCMatrixBlockIter&);
    void process_spec_rect(SCMatrixRectBlock*,
                           SCMatrixRectBlock*,
                           SCMatrixRectBlock*);
    void process_spec_ltri(SCMatrixLTriBlock*,
                           SCMatrixLTriBlock*,
                           SCMatrixLTriBlock*);

    int has_side_effects();
    int has_side_effects_in_arg1();
    int has_side_effects_in_arg2();

};

// //////////////////////////////////////////////////////////////////////////

/** OneBodyDerivInt is an abstract base class for objects that
    compute one body derivative integrals. */
class OneBodyDerivInt : public RefCount {
  protected:
    // this is who created me
    Integral *integral_;

    Ref<GaussianBasisSet> bs1;
    Ref<GaussianBasisSet> bs2;

    double *buffer_;

  public:
    OneBodyDerivInt(Integral *, const Ref<GaussianBasisSet>&b);
    OneBodyDerivInt(Integral *,
                    const Ref<GaussianBasisSet>&b1,
                    const Ref<GaussianBasisSet>&b2);
    virtual ~OneBodyDerivInt();
  
    /// Return the number of basis functions on center one.
    int nbasis() const;
    /// Return the number of basis functions on the center one.
    int nbasis1() const;
    /// Return the number of basis functions on the center two.
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

    /** The computed shell integrals will be put in the buffer returned by
        this member.  */
    const double * buffer() const;
    
    /** Compute the derivative integrals and place the result in the buffer
        returned by buffer(). */
    virtual void compute_shell(int ish, int jsh, DerivCenters&) = 0;
    /** Compute the derivative integrals with respect to the given center
        and place the result in the buffer returned by buffer(). */
    virtual void compute_shell(int ish, int jsh, int center) = 0;
};

// //////////////////////////////////////////////////////////////////////////

/** OneBodyOneCenterDerivInt is an abstract base class for objects that
    compute one body derivative integrals on a single center. */
class OneBodyOneCenterDerivInt : public RefCount {
  protected:
    // this is who created me
    Integral *integral_;

    Ref<GaussianBasisSet> bs1;

    double *buffer_;

  public:
    OneBodyOneCenterDerivInt(Integral *, const Ref<GaussianBasisSet>&b);
    virtual ~OneBodyOneCenterDerivInt();
  
    /// Return the number of basis functions on center one.
    int nbasis() const;
    /// Return the number of basis functions on center one.
    int nbasis1() const;

    /// Return the number of shells on center one.
    int nshell() const;
    /// Return the number of shells on center one.
    int nshell1() const;

    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis();
    /// Return the basis set on center one.
    Ref<GaussianBasisSet> basis1();

    /** The computed shell integrals will be put in the buffer returned by
        this member.  */
    const double * buffer() const;
    
    /** Compute the derivative integrals and place the result in the buffer
        returned by buffer(). */
    virtual void compute_shell(int ish, DerivCenters&) = 0;
    /** Compute the derivative integrals with respect to the given center
        and place the result in the buffer returned by buffer(). */
    virtual void compute_shell(int ish, int center) = 0;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
