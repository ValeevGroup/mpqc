//
// transform.h
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

#ifndef _chemistry_qc_basis_transform_h
#define _chemistry_qc_basis_transform_h

#include <vector>

namespace sc {

// ///////////////////////////////////////////////////////////////////////////

/** This is a base class for a container for a component of a sparse
    Cartesian to solid harmonic basis function transformation.  */
class SphericalTransformComponent {
  protected:
    double coef_;
    int a_, b_, c_, cartindex_, pureindex_;

  public:
    virtual ~SphericalTransformComponent();

    /// Returns the exponent of x.
    int a() const { return a_; }
    /// Returns the exponent of y.
    int b() const { return b_; }
    /// Returns the exponent of z.
    int c() const { return c_; }
    /// Returns the index of the Cartesian basis function.
    int cartindex() const { return cartindex_; }
    /// Returns the index solid harmonic basis function.
    int pureindex() const { return pureindex_; }
    /// Returns the coefficient of this component of the transformation.
    double coef() const { return coef_; }

    /** Initialize this object.  This must be provided in all
        specializations of this class to establish the ordering between a,
        b and c and the index of the Cartesian basis function.  Other
        things such as adjustment of the coefficient to account for
        normalization differences can be done as well.  The default
        SphericalTransform::init() implementation requires that only the
        x<sup>l</sup>, y<sup>l</sup> and z<sup>l</sup> basis functions are
        normalized to unity. */
    virtual void init(int a, int b, int c, double coef, int pureindex) =0;
};

// ///////////////////////////////////////////////////////////////////////////

/** This is a base class for a container for a sparse Cartesian to solid
    harmonic basis function transformation.  */
class SphericalTransform {
  protected:
    int n_;
    int l_;
    int subl_;
    std::vector<SphericalTransformComponent*> components_;

    SphericalTransform();

    /** This constructs the SphericalTransform for the given Cartesian
        angular momentum l and solid harmonic angular momentum subl.
        Usually, l and subl will be the same.  They would differ when the S
        component of a D Cartesian shell or the P component of an F
        Cartesian shell is desired, for example (see the natural atomic
        orbital code for an example of such use).  The init member must be
        called to complete initialization.  */
    SphericalTransform(int l, int subl = -1);

    /** This determines all of the components of the transformation.  It
        should be possible to implement the
        SphericalTransformComponent::init specialization in such a way that
        the default SphericalTransform::init can be used.  */
    virtual void init();
    
  public:
    virtual ~SphericalTransform();

    /** Adds another SphericalTransformComponent */
    void add(int a, int b, int c, double coef, int pureindex);

    /// Returns the Cartesian basis function index of component i.
    int cartindex(int i) const { return components_[i]->cartindex(); }
    /// Returns the solid harmonic basis function index of component i.
    int pureindex(int i) const { return components_[i]->pureindex(); }
    /// Returns the transform coefficient of component i.
    double coef(int i) const { return components_[i]->coef(); }
    /// Returns the Cartesian basis function's x exponent of component i.
    int a(int i) const { return components_[i]->a(); }
    /// Returns the Cartesian basis function's y exponent of component i.
    int b(int i) const { return components_[i]->b(); }
    /// Returns the Cartesian basis function's z exponent of component i.
    int c(int i) const { return components_[i]->c(); }
    /// Returns the angular momentum.
    int l() const { return l_; }
    /// Returns the number of components in the transformation.
    int n() const { return n_; }

    /** This must create a SphericalTransformComponent of the
        appropriate specialization. */
    virtual SphericalTransformComponent * new_component() = 0;
};

/// This describes a solid harmonic to Cartesian transform.
class ISphericalTransform: public SphericalTransform {
  protected:
    ISphericalTransform();
    ISphericalTransform(int l,int subl=-1);
    void init();
};

// ///////////////////////////////////////////////////////////////////////////

/// This iterates through the components of a SphericalTransform.
class SphericalTransformIter {
  private:
    int i_;

  protected:
    const SphericalTransform *transform_;
    
  public:
    SphericalTransformIter();
    SphericalTransformIter(const SphericalTransform*);

    void begin() { i_ = 0; }
    void start() { begin(); }
    void next() { i_++; }
    int ready() { return i_ < transform_->n(); }
    operator int() { return ready(); }
    int l() { return transform_->l(); }
    int cartindex() { return transform_->cartindex(i_); }
    int pureindex() { return transform_->pureindex(i_); }
    int bfn() { return pureindex(); }
    double coef() { return transform_->coef(i_); }
    int a() { return transform_->a(i_); }
    int b() { return transform_->b(i_); }
    int c() { return transform_->c(i_); }
    int l(int i) { return i?(i==1?b():c()):a(); }
    int n() { return 2*l() + 1; }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
