//
// xscf.h --- definition of the excited-state open shell singlet SCF class
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

#ifndef _chemistry_qc_scf_xscf_h
#define _chemistry_qc_scf_xscf_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/scmat/elemop.h>
#include <math/scmat/block.h>
#include <math/scmat/blkiter.h>
#include <math/optimize/scextrap.h>
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/wfn/effh.h>

#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>

////////////////////////////////////////////////////////////////////////////

class XSCF: public OneBodyWavefunction
{
#   define CLASSNAME XSCF
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
 protected:
    RefSelfConsistentExtrapolation _extrap;
    RefSCExtrapData _data;
    RefSCExtrapError _error;

    RefAccumDIH _accumdih;
    RefAccumDDH _accumddh;
    RefAccumEffectiveH _accumeffh;

    RefSymmSCMatrix _fockc;
    RefSymmSCMatrix _focka;
    RefSymmSCMatrix _fockb;
    RefSymmSCMatrix _fockab;

    RefSymmSCMatrix _ka;
    RefSymmSCMatrix _kb;

    RefDiagSCMatrix _fock_evalsc;
    RefDiagSCMatrix _fock_evalsa;
    RefDiagSCMatrix _fock_evalsb;
    
    int _ndocc;
    int aorb;
    int borb;

    double occa;
    double occb;

    double ci1;
    double ci2;

    double sab;
    double eop;
    
    int _density_reset_freq;

    int _maxiter;
    int _eliminate;

    char *ckptdir;
    char *fname;

    // these are temporary data, so they should not be checkpointed
    RefSymmSCMatrix _densc;
    RefSymmSCMatrix _densa;
    RefSymmSCMatrix _densb;
    RefSymmSCMatrix _densab;
    RefSymmSCMatrix _densab2;
    RefSCVector _ca;
    RefSCVector _cb;

    RefSymmSCMatrix _gr_hcore;
    RefSCMatrix _gr_vector;
    
    void init();
    void compute();
    void do_vector(double&,double&);
    void form_ao_fock(centers_t *, double*, double&);
    void do_gradient(const RefSCVector&);
    
  public:
    XSCF(StateIn&);
    XSCF(const XSCF&);
    XSCF(const RefKeyVal&);
    XSCF(const OneBodyWavefunction&);
    ~XSCF();

    XSCF& operator=(const XSCF&);
    
    void save_data_state(StateOut&);

    void print(ostream&o=cout);

    RefSCMatrix eigenvectors();

    double occupation(int vectornum);

    int value_implemented();
    int gradient_implemented();
    int hessian_implemented();
};
SavableState_REF_dec(XSCF);

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
// End:
