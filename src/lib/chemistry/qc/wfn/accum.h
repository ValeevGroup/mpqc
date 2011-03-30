//
// accum.h
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

#ifndef _chemistry_qc_wfn_accum_h
#define _chemistry_qc_wfn_accum_h

#include <chemistry/qc/wfn/wfn.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

/** AccumH computes additions to the one body Hamiltonian.  Specializations
    of AccumH can be given to derivatives of the SCF class, and there they
    will be used to modify the one body Hamiltonian.
 */
class AccumH: virtual public SavableState {
  protected:
    Ref<Wavefunction> wfn_;

  public:
    AccumH();
    AccumH(StateIn&);
    /** The KeyVal constructor.

        <dl>

        <dt><tt>wavefunction</tt><dd> This gives a Wavefunction object.
        Sometimes additional contributions to the Hamiltonian depend on the
        current Wavefunction, and, in these cases, the contribution must be
        obtained iteratively.  It is usually not necessary to specify this,
        since it is given in the init call below.

        </dl>
    */
    AccumH(const Ref<KeyVal>&);
    virtual ~AccumH();

    void save_data_state(StateOut&);
    
    /** Sets the current Wavefunction.  This is needed if the contribution
        depends on the current Wavefunction.  This would override a
        Wavefunction givin in the KeyVal CTOR.
    */
    virtual void init(const Ref<Wavefunction>&);
    /** Sum the contribution from this object into h. */
    virtual void accum(const RefSymmSCMatrix& h) =0;
    /** Print information about the contribution. */
    virtual void print_summary();
    /** Should be called after we are finished with this AccumH.  The
        reference to current Wavefunction object will be removed, and accum
        cannot be called until another init call is made.
    */
    virtual void done();

    /** Returns the scalar contribution to the energy.
        Available only after accum is called.
    */
    virtual double e();
};


/** This specialization of AccumH does nothing.
 */
class AccumHNull: public AccumH {
  public:
    AccumHNull();
    AccumHNull(StateIn&);
    AccumHNull(const Ref<KeyVal>&);
    ~AccumHNull();

    void save_data_state(StateOut&);
    
    void accum(const RefSymmSCMatrix& h);
};

/** This specialization of AccumHNull does nothing.  Sums the
    results of several AccumH objects.
 */
class SumAccumH: public AccumH {
  protected:
    int n_;
    Ref<AccumH> *accums_;

  public:
    SumAccumH(StateIn&);
    /** The KeyVal constructor.

        <dl>

        <dt><tt>accums</tt><dd>This gives an array of
        AccumH objects that are each called to obtain the total
        contribution.

        </dl>
    */
    SumAccumH(const Ref<KeyVal>&);
    ~SumAccumH();

    void save_data_state(StateOut&);

    void init(const Ref<Wavefunction>&);
    void accum(const RefSymmSCMatrix& h);
    void done();

    double e();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
