//
// sointegral.h --- definition of the Integral class
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

#ifndef _chemistry_qc_basis_sointegral_h
#define _chemistry_qc_basis_sointegral_h

#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/basis/sobasis.h>

namespace sc {

/** OneBodySOInt computes two-center one-electron integrals in a symmetry-adapted basis
 */
class OneBodySOInt : public RefCount {
  protected:
    Ref<OneBodyInt> ob_;

    Ref<SOBasis> b1_;
    Ref<SOBasis> b2_;

    double *buffer_;

    int only_totally_symmetric_;
  public:
    OneBodySOInt(const Ref<OneBodyInt> &);
    virtual ~OneBodySOInt();

    Ref<SOBasis> basis() const;
    Ref<SOBasis> basis1() const;
    Ref<SOBasis> basis2() const;

    /// TODO document OneBodySOInt::buffer()
    const double * buffer() const { return buffer_; }
    
    /// computes an SO shell doublet of integrals
    /// @param so_shell1 the SO shell index for center 1, so_shell1 \f$ \in \f$ [ 0, this->basis1()->nshell() )
    /// @param so_shell2 the SO shell index for center 2, so_shell1 \f$ \in \f$ [ 0, this->basis2()->nshell() )
    virtual void compute_shell(int so_shell1, int so_shell2);

    // an index of -1 for a shell indicates any shell
    //virtual int log2_shell_bound(int= -1,int= -1) = 0;

    virtual void reinitialize();

    int only_totally_symmetric() const { return only_totally_symmetric_; }
    void set_only_totally_symmetric(int i) { only_totally_symmetric_ = i; }
};


/** OneBodySODerivInt computes two-center one-electron integrals in a symmetry-adapted basis
 */
class OneBodySODerivInt : public RefCount {
  protected:
    Ref<OneBodyDerivInt> obd_;

    Ref<SOBasis> b1_;
    Ref<SOBasis> b2_;

    double *buffer_;

    int only_totally_symmetric_;
  public:
    OneBodySODerivInt(const Ref<OneBodyDerivInt> &);
    virtual ~OneBodySODerivInt();

    Ref<SOBasis> basis() const;
    Ref<SOBasis> basis1() const;
    Ref<SOBasis> basis2() const;

    /// TODO document OneBodySODerivInt::buffer()
    const double * buffer() const { return buffer_; }

    /// computes an SO shell doublet of integrals
    /// @param so_shell1 the SO shell index for center 1, so_shell1 \f$ \in \f$ [ 0, this->basis1()->nshell() )
    /// @param so_shell2 the SO shell index for center 2, so_shell1 \f$ \in \f$ [ 0, this->basis2()->nshell() )
    virtual void compute_shell(int so_shell1, int so_shell2);

    // an index of -1 for a shell indicates any shell
    //virtual int log2_shell_bound(int= -1,int= -1) = 0;

    virtual void reinitialize();

    int only_totally_symmetric() const { return only_totally_symmetric_; }
    void set_only_totally_symmetric(int i) { only_totally_symmetric_ = i; }
};

/** TwoBodySOInt computes four-center two-electron integrals in a symmetry-adapted basis
 */
class TwoBodySOInt : public RefCount {
  protected:
    Ref<TwoBodyInt> tb_;

    Ref<SOBasis> b1_;
    Ref<SOBasis> b2_;
    Ref<SOBasis> b3_;
    Ref<SOBasis> b4_;

    double *buffer_;

    int redundant_;
    int only_totally_symmetric_;
  public:
    TwoBodySOInt(const Ref<TwoBodyInt> &);
    virtual ~TwoBodySOInt();

    Ref<SOBasis> basis() const;
    Ref<SOBasis> basis1() const;
    Ref<SOBasis> basis2() const;
    Ref<SOBasis> basis3() const;
    Ref<SOBasis> basis4() const;

    /// TODO document TwoBodySOInt::buffer()
    const double * buffer() const { return buffer_; }
    
    /// computes an SO shell quartet of integrals
    /// @param so_shell1 the SO shell index for center 1, so_shell1 \f$ \in \f$ [ 0, this->basis1()->nshell() )
    /// @param so_shell2 the SO shell index for center 2, so_shell1 \f$ \in \f$ [ 0, this->basis2()->nshell() )
    /// @param so_shell3 the SO shell index for center 3, so_shell1 \f$ \in \f$ [ 0, this->basis3()->nshell() )
    /// @param so_shell4 the SO shell index for center 4, so_shell1 \f$ \in \f$ [ 0, this->basis4()->nshell() )
    virtual void compute_shell(int so_shell1, int so_shell2, int so_shell3, int so_shell4);

    // an index of -1 for a shell indicates any shell
    //virtual int log2_shell_bound(int= -1,int= -1,int= -1,int= -1) = 0;

    // if redundant is true, then keep redundant integrals in buffer_.  The
    // default is true.
    int redundant() const { return redundant_; }
    // cannot do nonredundant at the moment
    //void set_redundant(int i) { redundant_ = i; }

    int only_totally_symmetric() const { return only_totally_symmetric_; }
    void set_only_totally_symmetric(int i) { only_totally_symmetric_ = i; }
};

/** TwoBodySODerivInt computes four-center two-electron derivative integrals in a symmetry-adapted basis
 */
class TwoBodySODerivInt : public RefCount {
  protected:
    Ref<TwoBodyDerivInt> tb_;

    Ref<SOBasis> b1_;
    Ref<SOBasis> b2_;
    Ref<SOBasis> b3_;
    Ref<SOBasis> b4_;

    double *buffer_;

    int redundant_;
    int only_totally_symmetric_;
  public:
    TwoBodySODerivInt(const Ref<TwoBodyDerivInt> &);
    virtual ~TwoBodySODerivInt();

    Ref<SOBasis> basis() const;
    Ref<SOBasis> basis1() const;
    Ref<SOBasis> basis2() const;
    Ref<SOBasis> basis3() const;
    Ref<SOBasis> basis4() const;

    /// TODO document TwoBodySODerivInt::buffer()
    const double * buffer() const { return buffer_; }

    /// computes an SO shell quartet of derivative integrals
    /// @param so_shell1 the SO shell index for center 1, so_shell1 \f$ \in \f$ [ 0, this->basis1()->nshell() )
    /// @param so_shell2 the SO shell index for center 2, so_shell1 \f$ \in \f$ [ 0, this->basis2()->nshell() )
    /// @param so_shell3 the SO shell index for center 3, so_shell1 \f$ \in \f$ [ 0, this->basis3()->nshell() )
    /// @param so_shell4 the SO shell index for center 4, so_shell1 \f$ \in \f$ [ 0, this->basis4()->nshell() )
    virtual void compute_shell(int so_shell1, int so_shell2, int so_shell3, int so_shell4);

    // an index of -1 for a shell indicates any shell
    //virtual int log2_shell_bound(int= -1,int= -1,int= -1,int= -1) = 0;

    // if redundant is true, then keep redundant integrals in buffer_.  The
    // default is true.
    int redundant() const { return redundant_; }
    // cannot do nonredundant at the moment
    //void set_redundant(int i) { redundant_ = i; }

    int only_totally_symmetric() const { return only_totally_symmetric_; }
    void set_only_totally_symmetric(int i) { only_totally_symmetric_ = i; }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
