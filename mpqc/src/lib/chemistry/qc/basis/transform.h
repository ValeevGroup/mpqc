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

#if defined(__GNUC__)
#pragma interface
#endif

#ifndef _chemistry_qc_basis_tranform_h
#define _chemistry_qc_basis_tranform_h

/////////////////////////////////////////////////////////////////////////////

class SphericalTransformComponent {
  protected:
    double coef_;
    int a_, b_, c_, cartindex_, pureindex_;

  public:
    int a() const { return a_; }
    int b() const { return b_; }
    int c() const { return c_; }
    int cartindex() const { return cartindex_; }
    int pureindex() const { return pureindex_; }
    double coef() const { return coef_; }

    virtual void init(int a, int b, int c, double coef, int pureindex) =0;
};

/////////////////////////////////////////////////////////////////////////////

class SphericalTransform {
  protected:
    int n_;
    int l_;
    SphericalTransformComponent *components_;

    SphericalTransform();
    SphericalTransform(int l);

    virtual void init();
    
  public:
    virtual ~SphericalTransform();

    void add(int a, int b, int c, double coef, int pureindex);

    int cartindex(int i) const { return components_[i].cartindex(); }
    int pureindex(int i) const { return components_[i].pureindex(); }
    double coef(int i) const { return components_[i].coef(); }
    int a(int i) const { return components_[i].a(); }
    int b(int i) const { return components_[i].b(); }
    int c(int i) const { return components_[i].c(); }
    int l() const { return l_; }
    int n() const { return n_; }

    virtual SphericalTransformComponent * new_components() = 0;
};

// The inverse transforms
class ISphericalTransform: public SphericalTransform {
  protected:
    ISphericalTransform();
    ISphericalTransform(int l);
    void init();
};

/////////////////////////////////////////////////////////////////////////////

class SphericalTransformIter {
  private:
    int i_;

  protected:
    const SphericalTransform *transform_;
    
  public:
    SphericalTransformIter();
    SphericalTransformIter(SphericalTransform*);

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

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
