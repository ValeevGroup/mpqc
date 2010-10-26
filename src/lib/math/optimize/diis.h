//
// diis.h
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

#ifndef _math_optimize_diis_h
#define _math_optimize_diis_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/optimize/scextrap.h>

namespace sc {

/** The DIIS class provides DIIS extrapolation. */
class DIIS: public SelfConsistentExtrapolation {
  protected:
    int start;
    int ndiis;
    int iter;
    int ngroup;
    int ngroupdiis;
    double damping_factor;
    double mixing_fraction;

    double * btemp;
    double ** bold;
    double ** bmat;

    Ref<SCExtrapData> Ldata;

    Ref<SCExtrapData> *diism_data;
    Ref<SCExtrapData> *diism_datain;
    Ref<SCExtrapError> *diism_error;

    void init();
  public:
    DIIS(int strt=1, int ndi=5, double dmp =0, int ngr=1, int ngrdiis=1,
         double mf=0);
    DIIS(StateIn&);
    /** The DIIS KeyVal constructor recognizes the following keywords:

        <dl>

        <dt><tt>n</tt><dd> This integer maximum number of data sets to
        retain.  The default is 5.

        <dt><tt>start</tt><dd> The DIIS extrapolation will begin on the
        iteration given by this integer.  The default is 1.

        <dt><tt>damping_factor</tt><dd> This nonnegative floating point
        number is used to dampen the DIIS extrapolation.  The default is
        0.0.

        <dt><tt>mixing_fraction</tt><dd> This floating point number in
        [0,1] is used to dampen the DIIS extrapolation by mixing the input
        data with the output data for each iteration. The default is 0.0,
        which performs no mixing. The approach described in
        Kerker, Phys. Rev. B, 23, p3082, 1981.

        <dt><tt>ngroup</tt><dd> The number of iterations in a DIIS group.
        DIIS extrapolation is only used for the first ngroupdiis of these
        interations.  The default is 1.  If ngroup is 1 and ngroupdiis is
        greater than 0, then DIIS will be used on all iterations after and
        including the start iteration.

        <dt><tt>ngroupdiis</tt><dd> The number of DIIS extrapolations to do
        at the beginning of an iteration group.  See the documentation for
        ngroup.

        </dl> */
    DIIS(const Ref<KeyVal>&);
    ~DIIS();

    void save_data_state(StateOut&);
    
    int extrapolate(const Ref<SCExtrapData>& data,
                    const Ref<SCExtrapError>& error);

    void start_extrapolation();

    void reinitialize(Ref<SCExtrapData> data=0);

    /// Override DescribedClass::print.
    void print(std::ostream& = ExEnv::out0()) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
