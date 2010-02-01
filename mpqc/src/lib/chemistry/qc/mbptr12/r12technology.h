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

#include <chemistry/qc/mbptr12/linearr12.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

/** R12Technology describes technical features of the R12 approach. */
class R12Technology: virtual public SavableState {
  public:

  /**
    Projector of R12 methods:
    0: Q_{12} = 1
    1: Q_{12} = (1 - P_1)(1 - P_2)
    2: Q_{12} = (1 - V_1 V_2)(1 - O_1)(1 - O_2)
    3: Q_{12} = 1 - P_1 P_2
  */
  enum Projector {Projector_0 = 0,
    Projector_1 = 1,
    Projector_2 = 2,
    Projector_3 = 3};
  enum StandardApproximation {
    //StdApprox_A = 0, // is now obsolete
    StdApprox_Ap = 1,
    StdApprox_App = 2,
    StdApprox_B = 3,
    StdApprox_C = 4,
    StdApprox_Cp = 5
    };
  enum ABSMethod {ABS_ABS = 0,
    ABS_ABSPlus = 1,
    ABS_CABS = 2,
    ABS_CABSPlus = 3};

  /// geminal generating space
  enum OrbitalProduct_GG {
    OrbProdGG_ij = 0,
    OrbProdGG_pq = 1,
  };

  /// space of orbital products from which geminal substitutions are allowed
  enum OrbitalProduct_gg {
    OrbProdgg_ij = 0,
    OrbProdgg_pq = 1,
  };

  enum PositiveDefiniteB {
    PositiveDefiniteB_no = 0,
    PositiveDefiniteB_yes = 1,
    PositiveDefiniteB_weak = 2
  };

  enum GeminalAmplitudeAnsatz {
    GeminalAmplitudeAnsatz_fullopt = 0,
    GeminalAmplitudeAnsatz_fixed = 1,
    GeminalAmplitudeAnsatz_scaledfixed = 2
  };

  enum H0_dk_approx_pauli {
    H0_dk_approx_pauli_true = 0,
    H0_dk_approx_pauli_fHf = 1,
    H0_dk_approx_pauli_fHf_Q = 2,
    H0_dk_approx_pauli_false = 3
  };

  /**
   * R12Ansatz specifies the manner in which the R12 geminals are constructed.
   */
  class R12Ansatz : virtual public SavableState {
    public:
    /** The KeyVal constructor.
    <dl>

    <dt><tt>orbital_product_GG</tt><dd> This specifies how the geminal space is generated.
    Geminal functions are products of the correlation factor and 2 orbitals.
    This keyword specifies which orbital products are allowed.
    Valid choices are:
      <dl>
        <dt><tt>ij</tt><dd> Biproducts of occupied orbitals. This is the default.
        <dt><tt>pq</tt><dd> Biproducts of any Hartree-Fock orbitals. This has not been implemented yet.
      </dl>

    <dt><tt>orbital_product_gg</tt><dd> Space of orbital products from which geminal substitutions are allowed.
    Specified in the same way as orbital_product_GG.

    <dt><tt>projector</tt><dd> This specifies the form of the orthogonal projector.
    Valid values are:
      <dl>
        <dt><tt>0</tt><dd> 1. Should be used ONLY for testing. This implies the Weak Orthogonality Functional (WOF). Not implemented yet.
        <dt><tt>1</tt><dd> (1-P1)(1-P2). Not implemented yet.
        <dt><tt>2</tt><dd> (1-O1)(1-O2)(1-V1V2). This is the default.
        <dt><tt>3</tt><dd> 1-P1P2. Should be used ONLY for testing.
      </dl>

    <dt><tt>wof</tt><dd> Setting this to <tt>true</tt> will cause the Weak Orthogonality Functional to be used. The default is <tt>false</tt>,
    unless <tt>projector=0</tt>.

    <dt><tt>diag</tt><dd> Setting this to <tt>true</tt> will only keep the diagonal terms,
    which is equivalent to the "old" (pre-1992) form of R12 theory. The default is <tt>false</tt>,
    which corresponds to the orbital invariant ansatz of Klopper.

    <dt><tt>amplitudes</tt><dd> This keyword specifies how the geminal amplitudes are determined.
    Permitted values are <tt>optimized</tt> (for fully optimized amplitudes) and <tt>fixed</tt>
    (fixed using first-order cusp-conditions, a la Ten-no). The default is <tt>fixed</tt>
    if the diagonal ansatz is used with an appropriate correlation factor (either Slater-type geminal
    or linear), otherwise <tt>optimized</tt>.

    </dl>
    */
    R12Ansatz(const Ref<KeyVal>&);
    /// The StateIn constructor
    R12Ansatz(StateIn&);
    /// The default constructor creates orbital-invariant ansatz with projector 2
    R12Ansatz();
    ~R12Ansatz();

    void save_data_state(StateOut&);
    void print(std::ostream& o =ExEnv::out0()) const;

    R12Technology::Projector projector() const;
    bool diag() const;
    R12Technology::GeminalAmplitudeAnsatz amplitudes() const;
    bool wof() const;
    R12Technology::OrbitalProduct_GG orbital_product_GG() const;
    R12Technology::OrbitalProduct_gg orbital_product_gg() const;

    private:
    R12Technology::Projector projector_;
    bool diag_;
    R12Technology::GeminalAmplitudeAnsatz amplitudes_;
    bool scaled_;     //<
    bool wof_;
    R12Technology::OrbitalProduct_GG orbital_product_GG_;
    R12Technology::OrbitalProduct_gg orbital_product_gg_;
  };

  private:
    bool abs_eq_obs_;
    bool vbs_eq_obs_;

    Ref<LinearR12::CorrelationFactor> corrfactor_;
    StandardApproximation stdapprox_;
    Ref<R12Ansatz> ansatz_;
    ABSMethod abs_method_;
    unsigned int maxnabs_;
    bool gbc_;
    bool ebc_;
    bool coupling_;
    bool omit_P_;
    H0_dk_approx_pauli H0_dk_approx_pauli_;
    bool H0_dk_keep_;
    bool safety_check_;
    PositiveDefiniteB posdef_B_;

    // for debugging purposes only
    bool omit_B_;

    // no need to store this guy
#if 0
    // determines the weight function used to fit the correlation factor
    struct GTGFitWeight {
      typedef enum {TewKlopper, Cusp} Type;
    };
    GTGFitWeight::Type gtg_fit_weight_;
#endif

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

          <dt><tt>C</tt><dd> Use second order M&oslash;ller-Plesset perturbation theory
      with linear R12 terms in standard approximation C' (simplified variant of approximation C,
      in which only integrals with 1 ABS index are used. Refer to: Valeev, to be published).
          Only energies can be computed with the MP2-R12/C' method.

        </dl>

        <dt><tt>ansatz</tt><dd> This object specifies the ansatz (see R12Ansatz).

        <dt><tt>gbc</tt><dd> This boolean specifies whether Generalized Brillouin
        Condition (GBC) is assumed to hold. The default is "true". This keyword is
        only valid if stdapprox=A'.
        The effect of setting this keyword to false is very small --
        hence it is not recommended to use this keyword.

        <dt><tt>ebc</tt><dd> This boolean specifies whether Extended Brillouin
        Condition (EBC) is assumed to hold. The default is "true". This keyword
        is only valid if stdapprox=A'.
        The effect of setting this keyword to false is small --
        hence it is not recommended to use this keyword.

        <dt><tt>coupling</tt><dd> This boolean specifies whether the doubles-geminal block of the zeroth-order
        operator is included. The default is "false". This keyword
        is not valid if stdapprox=A'.
        The effect of setting this keyword to true is small --
        hence it is not recommended to use this keyword.

        <dt><tt>maxnabs</tt><dd> This integer specifies the maximum number of ABS indices per integral.
        Valid values are 1 or 2. The default is 2 except for R12/A'' method.

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

	<dt><tt>gtg_fit_weight</tt><dd> This keyword determines how the correlation factor is fit to Gaussians (hence
	only valid when <tt>corr_factor</tt> is set to <tt>stg-ng</tt>)
	The choices are <tt>tewklopper</tt>, which is appropriate for energy computations, and <tt>cusp</tt>, which is appropriate
	for accurate cusp region description. The default is <tt>tewklopper</tt>. Choosing <tt>cusp</tt> is probably only appropriate
	when many (9 or more) Gaussians are used for the fit.

    <dt><tt>H0_dk_approx_pauli</tt><dd> This string keyword determines how H0 DK Hamiltonian
    is approximated by Pauli Hamiltonian. The allowed values are
    <dl>

          <dt><tt>true/yes</tt><dd> Use Pauli everywhere. Valid in approximations A'' and C.
          In approximation C this choice involves the appearance of M intermediate
          (double commutator of mass-velocity term with 2 correlation factors)
          and use of Pauli Hamiltonian everywhere in P intermediate.
          In approximations A'' this involves the use of Pauli Hamiltonian in
          the double commutator, in single commutator, in Q intermediate.

          <dt><tt>fHf</tt><dd> Use Pauli in the "diagonal" term, i.e. in f12 H f12
          which involves the use of Pauli Hamiltonian in the double commutator
          and in Q intermediate. Only valid in approximation C.

          <dt><tt>fHf_Q</tt><dd> Same as <tt>fHf</tt> but use the full DK Hamiltonian
          in Q intermediate. This is equivalent to assuming that terms of higher order
          than Pauli commute with the correlation factor, hence they should be kept in Q.
          Only valid in approximation C.

          <dt><tt>false/no</tt><dd> Use full DKH operator (this is the default).
          Valid in approximations A'' and C.
          In approximation C this means treat relativistic terms like exchange,
          this affects the Q intermediate and fKf part of P intermediate.
          In approximation A'' this means that all relativistic terms are dropped from H0 --
          this is chosen to be consistent with the nonrelativistic A'' method
          where exchange operator is dropped completely because its
          commutators cannot be evaluated analytically.

        </dl>

    <dt><tt>H0_dk_keep</tt><dd> This boolean keyword specifies whether to keep
    relativistic terms or drop them. This is only considered in approximation A'' if <tt>H0_dk_approx_pauli=false</tt>.
    The default is <tt>false</tt>. Setting to <tt>true</tt> will keep the
    relativistic terms in Q intermediate of the A'' approximation.

    */
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
    bool coupling() const;
    ABSMethod abs_method() const;
    StandardApproximation stdapprox() const;
    const Ref<R12Ansatz>& ansatz() const;
    bool spinadapted() const;
    bool omit_P() const;
    H0_dk_approx_pauli H0_dk_approx() const;
    bool H0_dk_keep() const;
    bool safety_check() const;
    PositiveDefiniteB posdef_B() const;

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
