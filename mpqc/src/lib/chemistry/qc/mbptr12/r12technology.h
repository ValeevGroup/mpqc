//
// r12technology.h
//
// Copyright (C) 2007 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#ifndef _chemistry_qc_mbptr12_r12technology_h
#define _chemistry_qc_mbptr12_r12technology_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/mbptr12/ansatz.h>
#include <chemistry/qc/mbptr12/linearr12.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

/** R12Technology describes technical features of the R12 approach. */
class R12Technology: virtual public SavableState {

    bool abs_eq_obs_;
    bool vbs_eq_obs_;

    Ref<LinearR12::CorrelationFactor> corrfactor_;
    LinearR12::StandardApproximation stdapprox_;
    Ref<LinearR12Ansatz> ansatz_;
    LinearR12::ABSMethod abs_method_;
    unsigned int maxnabs_;
    bool gbc_;
    bool ebc_;
    bool omit_P_;
    bool safety_check_;
    LinearR12::PositiveDefiniteB posdef_B_;

    // for debugging purposes only
    bool omit_B_;

  public:
    R12Technology(StateIn&);
    /** The KeyVal constructor.
        <dl>

        <dt><tt>corr_factor</tt><dd> This string specifies which correlation factor to use.
        Allowed values are "r12", "g12", "geng12", and "none". The default is "r12".
      
        <dt><tt>corr_param</tt><dd> This keyword specifies optional parameters
        of the correlation factor. <tt>corr_param</tt> can be a single floating-point value
        an array of floating-point values, or an array of arrays of 2-element arrays of
        floating-point values. Single value specifies the parameter of the single
        correlation function. The 1-d array form specifies a set of primitive correlation functions
        characterized by the corresponding parameters. The 3-d array form specifies
        a set of contracted correlation functions. For example,
        <tt>corr_param = 3.0</tt> specifies a single correlation function
        with parameter 3.0. <tt>corr_param = [ 1.0 3.0 10.0 ]</tt> specifies
        3 correlation functions with parameters 1.0, 3.0 and 10.0.
        <tt>corr_param = [ [[1.0 0.35][3.0 0.65]]  [[10.0 1.0]] ]</tt>
        specifies 2 correlation functions, first composed of 2 primitive functions
        with parameters 1.0 and 3.0 combined linearly with coefficients
        0.35 and 0.65, and second primitive function with parameter 10.0 .

        This keyword has no meaning for some correlation factors, e.g., "r12" and "none",
        and is not used. There is no default.
      
        <dt><tt>stdapprox</tt><dd> This gives a string that must take on one
        of the values below.  The default is A'. WARNING: standard approximation A is now obsolete.

        <dl>

          <dt><tt>A'</tt><dd> Use second order M&oslash;ller-Plesset perturbation theory
	  with linear R12 terms in standard approximation A' (MP2-R12/A').
          This will cause MP2-R12/A energies to be computed also.
          Only energies can be computed with the MP2-R12/A' method.

          <dt><tt>A''</tt><dd> Use second order M&oslash;ller-Plesset perturbation theory
	  with linear R12 terms in standard approximation A'' (MP2-R12/A'').
          Only energies can be computed with the MP2-R12/A'' method.

          <dt><tt>B</tt><dd> Use second order M&oslash;ller-Plesset perturbation theory
	  with linear R12 terms in standard approximation B.
          This will cause A and A' energies to be computed also.
          Only energies can be computed with the MP2-R12/B method.

          <dt><tt>C</tt><dd> Use second order M&oslash;ller-Plesset perturbation theory
	  with linear R12 terms in standard approximation C.
          Only energies can be computed with the MP2-R12/C method.

        </dl>
        
        <dt><tt>ansatz</tt><dd> This object specifies the ansatz (see LinearR12Ansatz).
        
        <dt><tt>gbc</tt><dd> This boolean specifies whether Generalized Brillouin
        Condition (GBC) is assumed to hold. The default is "true". This keyword is
        only valid if stdapprox=A'.
        The effect of setting this keyword to true is very small --
        hence it is not recommended to use this keyword.

        <dt><tt>ebc</tt><dd> This boolean specifies whether Extended Brillouin
        Condition (EBC) is assumed to hold. The default is "true". This keyword
        is only valid if stdapprox=A'.
        The effect of setting this keyword to true is small --
        hence it is not recommended to use this keyword.
       
        <dt><tt>maxnabs</tt><dd> This integer specifies the maximum number of ABS indices per integral.
        Valid values are between 1 and 2. The default is to include all terms necessary for a given method.
        For example, MP2-F12/B energy involves integrals with 2 ABS indices. Setting maxnabs to 1
        will leave out such terms.
       
        <dt><tt>abs_method</tt><dd> This string specifies whether the old ABS method, introduced
        by Klopper and Samson, or the new ABS variant, CABS, introduced by Valeev, should be used.
	Valid values are "ABS" (Klopper and Samson), "ABS+", "CABS", and "CABS+", where the "+" labels
	a method where the union of OBS and ABS is used to construct the RI basis. The default is "ABS".
        The default in 2.3.0 and later will be "CABS+".

	<dt><tt>safety_check</tt><dd> Set to true if you want to perform safety checks, e.g., for completeness
        of the RI basis, linear independence of the geminal basis, positive definiteness of B matrix, etc.
	The default is true (to perform the checks).

	<dt><tt>posdef_B</tt><dd> This keyword specifies whether and how to enforce the positive definiteness of
	matrix B. Valid choices are <tt>no</tt>, <tt>yes</tt> (enforce positive definite matrix B and its pair-dependent
	counterpart, tilde-B), <tt>weak</tt> (same as <tt>yes</tt>, except the positive-definiteness of tilde-B
	is not enforced). If this keyword is set to <tt>no</tt> then sometimes nonphysical results can be obtained, e.g.,
	positive pair energy corrections can result from using too many correlation functions.
	<tt>posdef_B = yes</tt> offers the best protection against nonphysical results.
	The default is <tt>weak</tt>, which is cheaper <tt>yes</tt> and is definitely safer than <tt>no</tt>.

        </dl> */
    R12Technology(const Ref<KeyVal>&,
		  const Ref<GaussianBasisSet>& bs,
		  const Ref<GaussianBasisSet>& vbs,
		  const Ref<GaussianBasisSet>& abs
	);
    ~R12Technology();

    void save_data_state(StateOut&);

    const Ref<LinearR12::CorrelationFactor>& corrfactor() const;
    /// this changes the correlation factor
    void corrfactor(const Ref<LinearR12::CorrelationFactor>&);
    unsigned int maxnabs() const;
    bool gbc() const;
    bool ebc() const;
    LinearR12::ABSMethod abs_method() const;
    LinearR12::StandardApproximation stdapprox() const;
    const Ref<LinearR12Ansatz>& ansatz() const;
    bool spinadapted() const;
    bool omit_P() const;
    bool include_mp1() const;
    bool safety_check() const;
    const LinearR12::PositiveDefiniteB& posdef_B() const;

    //
    // these are for debugging only
    //
    // omit expensive parts of B
    bool omit_B() const;
    
    // This checks if ints is suitable for R12 calculations. Throws, if false.
    void check_integral_factory(const Ref<Integral>& ints);

    void print(std::ostream&o=ExEnv::out0()) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
