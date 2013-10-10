//
// rdm.h
//
// Copyright (C) 2009 Edward Valeev
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

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_rdm_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_rdm_h

#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/wfn/spin.h>
#include <math/distarray4/distarray4.h>

namespace sc {

  /// Rank of the RDM
  typedef enum {Zero=0, One = 1, Two = 2, Three = 3, Four = 4} Rank;

  /// @addtogroup ChemistryElectronicStructure
  /// @{

  namespace {
    template <Rank R> struct __spincase;
    template <> struct __spincase<One> {
      typedef SpinCase1 type;
    };
    template <> struct __spincase<Two> {
      typedef SpinCase2 type;
    };

    template <Rank R> struct __nspincases;
    template <> struct __nspincases<One> {
      static const int value = NSpinCases1;
    };
    template <> struct __nspincases<Two> {
        static const int value = NSpinCases2;
    };

  }

  template <Rank R> class RDMCumulant;
  class OrbitalSpace;

  /// RDM<R> is a reduced density matrix of rank R
  /// @tparam R Rank of the density
  template <Rank R>
  class RDM : public Compute, virtual public SavableState {
      typedef RDM<R> this_type;
      typedef typename __spincase<R>::type spincase;
    public:
      /** A KeyVal constructor is used to generate a RDM<R>
          object from the input. The full list of keywords
          that are accepted is below.

          <table border="1">

          <tr><td>%Keyword<td>Type<td>Default<td>Description

          <tr><td><tt>wfn</tt><td>Wavefunction<td>none<td>the Wavefunction object

          </table>
       */
      RDM(const Ref<KeyVal>& kv) {
        wfn_ = require_dynamic_cast<Wavefunction*>(
              kv->describedclassvalue("wfn").pointer(),
              "RDM<R>::RDM\n"
              );
      }
      RDM(StateIn& si) : SavableState(si) {
        wfn_ << SavableState::restore_state(si);
      }
      RDM(const Ref<Wavefunction>& wfn) : wfn_(wfn) {
      }
      ~RDM() {
      }
      void save_data_state(StateOut& so) {
        SavableState::save_state(wfn_.pointer(), so);
      }

      virtual void obsolete() {
        wfn_->obsolete();
        for(int s=0; s<__nspincases<R>::value; ++s)
          scmat_[s] = 0;
      }

      /// the corresponding Wavefunction
      Ref<Wavefunction> wfn() const { return wfn_; }
      virtual void compute() {
        const double energy = wfn_->value();
      }
      /// the orbital space of spincase s in which the density is reported
      virtual Ref<OrbitalSpace> orbs(SpinCase1 s) const =0;
      /// bra/ket dimension
      virtual size_t ndim(spincase spincase) const;
      /// returns the ket block for the given bra index
      virtual const double* obtain_block(spincase spin,  size_t bra) const {
        throw ProgrammingError("RDM<R>::obtain_block() is not yet implemented",
                               __FILE__,
                               __LINE__);
      }
      /// releases the ket block
      virtual void release_block(spincase spin, size_t bra, double*) const {
        throw ProgrammingError("RDM<R>::release_block() is not yet implemented",
                               __FILE__,
                               __LINE__);
      }
      /// full density matrix
      virtual RefSymmSCMatrix scmat(spincase spin) const;

      /// cumulant of rank R
      virtual Ref< RDMCumulant<R> > cumulant() const;
      /// RDM of rank decreased by 1
      virtual Ref< RDM< static_cast<Rank>(R-1) > > rdm_m_1() const;

    private:
      static ClassDesc class_desc_;
      Ref<Wavefunction> wfn_;

    protected:
      mutable RefSymmSCMatrix scmat_[__nspincases<R>::value];
  };

  template <Rank R>
  ClassDesc
  RDM<R>::class_desc_(typeid(this_type),
                      (std::string("RDM<") +
                       char('1' + R - 1) +
                       std::string(">")
                      ).c_str(), 1,
                      "virtual public SavableState", 0, 0, 0 );

  /// this specialization is needed to make RDM<R>::rdm_m_1() work
  template <> class RDM<Zero> : public RefCount {};

  ///////////////

  /// RDMCumulant<R> is a reduced density matrix cumulant of rank R
  /// @param R Rank of the density
  template <Rank R>
  class RDMCumulant : public Compute, virtual public SavableState {
      typedef RDMCumulant<R> this_type;
      typedef RDM<R> density_type;
      typedef typename __spincase<R>::type spincase;
    public:
      RDMCumulant(const Ref<density_type>& density) : density_(density) {
      }
      RDMCumulant(StateIn& si) : SavableState(si) {
        density_ << SavableState::restore_state(si);
      }
      ~RDMCumulant() {
      }
      void save_data_state(StateOut& so) {
        SavableState::save_state(density_.pointer(), so);
      }
      void obsolete() {
        density_->obsolete();
        for(int s=0; s<__nspincases<R>::value; ++s)
          scmat_[s] = 0;
      }

      /// the corresponding Wavefunction
      Ref<Wavefunction> wfn() const { return density_->wfn(); }
      void compute() { density_->compute(); }
      /// the corresponding Density
      Ref<density_type> density() const { return density_; }
      /// bra/ket dimension
      size_t ndim(spincase spincase) const { return density_->ndim(); }
      /// returns the ket block for the given bra index
      virtual const double* obtain_block(spincase spin,  size_t bra) const {
        throw ProgrammingError("RDMCumulant<R>::obtain_block() is not yet implemented",
                               __FILE__,
                               __LINE__);
      }
      /// releases the ket block
      virtual void release_block(spincase spin, size_t bra, double*) const {
        throw ProgrammingError("RDMCumulant<R>::release_block() is not yet implemented",
                               __FILE__,
                               __LINE__);
      }
      /// full cumulant matrix
      virtual RefSymmSCMatrix scmat(spincase spin) const;

    private:
      static ClassDesc class_desc_;
      Ref<density_type> density_;

    protected:
      mutable RefSymmSCMatrix scmat_[__nspincases<R>::value];
  };

  template <Rank R>
  ClassDesc
  RDMCumulant<R>::class_desc_(typeid(this_type),
			      (std::string("RDMCumulant<") +
			       char('1' + R - 1) +
			       std::string(">")
			      ).c_str(), 1,
			      "virtual public SavableState", 0, 0, 0 );


  /// SpinFreeRDM<R> is a spin-free reduced density matrix of rank R
  /// @tparam R Rank of the density
  template <Rank R>
  class SpinFreeRDM : public Compute, virtual public SavableState {
      typedef SpinFreeRDM<R> this_type;
      typedef typename __spincase<R>::type spincase;
    public:
      /** A KeyVal constructor is used to generate a SpinFreeRDM<R>
          object from the input. The full list of keywords
          that are accepted is below.

          <table border="1">

          <tr><td>%Keyword<td>Type<td>Default<td>Description

          <tr><td><tt>wfn</tt><td>Wavefunction<td>none<td>the Wavefunction object

          </table>
       */
      SpinFreeRDM(const Ref<KeyVal>& kv) {
        wfn_ = require_dynamic_cast<Wavefunction*>(
              kv->describedclassvalue("wfn").pointer(),
              "RDM<R>::RDM\n"
              );
      }
      SpinFreeRDM(StateIn& si) : SavableState(si) {
        wfn_ << SavableState::restore_state(si);
      }
      SpinFreeRDM(const Ref<Wavefunction>& wfn) : wfn_(wfn) {
      }
      ~SpinFreeRDM() {
      }
      void save_data_state(StateOut& so) {
        SavableState::save_state(wfn_.pointer(), so);
      }

      virtual void obsolete() {
        wfn_->obsolete();
        scmat_ = 0;
        da4_ = 0;
      }

      /// the corresponding Wavefunction
      Ref<Wavefunction> wfn() const { return wfn_; }
      virtual void compute() {
        const double energy = wfn_->value();
      }
      /// the orbital space of spincase s in which the density is reported
      virtual Ref<OrbitalSpace> orbs() const =0;
      /// bra/ket dimension
      virtual size_t ndim() const;
      /// returns the ket block for the given bra index
      virtual const double* obtain_block(spincase spin,  size_t bra) const {
        throw ProgrammingError("SpinFreeRDM<R>::obtain_block() is not yet implemented",
                               __FILE__,
                               __LINE__);
      }
      /// releases the ket block
      virtual void release_block(spincase spin, size_t bra, double*) const {
        throw ProgrammingError("SpinFreeRDM<R>::release_block() is not yet implemented",
                               __FILE__,
                               __LINE__);
      }
      /// full density matrix, can be used for RDM of any rank
      virtual RefSymmSCMatrix scmat() const;
      /** should only be used for R=2
          @return DistArray4 object that contains RDM2
        */
      virtual const Ref<DistArray4>& da4() const;

      /// RDM of rank decreased by 1
      virtual Ref< SpinFreeRDM< static_cast<Rank>(R-1) > > rdm_m_1() const;

    private:
      static ClassDesc class_desc_;
      Ref<Wavefunction> wfn_;

    protected:
      // used for R=1 (and for now for R=2)
      mutable RefSymmSCMatrix scmat_;
      // (will be) used for R=2
      mutable Ref<DistArray4> da4_;
  };

  template <Rank R>
  ClassDesc
  SpinFreeRDM<R>::class_desc_(typeid(this_type),
                      (std::string("SpinFreeRDM<") +
                       char('1' + R - 1) +
                       std::string(">")
                      ).c_str(), 1,
                      "virtual public SavableState", 0, 0, 0 );

  /// this specialization is needed to make SpinFreeRDM<R>::rdm_m_1() work
  template <> class SpinFreeRDM<Zero> : public RefCount {};

  ///////////////

  class OrbitalSpace;

  /// OBWfnRDMTwo is a 2-RDM from a OneBodyWavefunction
  class OBWfnRDMTwo : public RDM<Two> {
      typedef RDMCumulant<Two> cumulant_type;
      typedef RDM<One> rdm_m_1_type;
    public:
    /** A KeyVal constructor is used to generate a OBWfnRDMTwo
        object from the input. The full list of keywords
        that are accepted is below.

        <table border="1">

        <tr><td>%Keyword<td>Type<td>Default<td>Description

        <tr><td><tt>wfn</tt><td>OneBodyWavefunction<td>none<td>the OneBodyWavefunction object

        </table>
     */
      OBWfnRDMTwo(const Ref<KeyVal>& kv);
      OBWfnRDMTwo(StateIn& si);
      OBWfnRDMTwo(const Ref<OneBodyWavefunction>& wfn);
      virtual ~OBWfnRDMTwo();
      void save_data_state(StateOut& so);

      Ref<OneBodyWavefunction> wfn() const { return wfn_; }
      Ref<OrbitalSpace> orbs(SpinCase1 s) const;
      RefSymmSCMatrix scmat(SpinCase2 spincase) const;
      Ref<cumulant_type> cumulant() const;
      Ref<rdm_m_1_type> rdm_m_1() const;

    private:
      Ref<OneBodyWavefunction> wfn_;
      mutable RefSymmSCMatrix scmat_[NSpinCases2];
      mutable Ref<OrbitalSpace> orbs_[NSpinCases1];

      static ClassDesc class_desc_;
  };

  /// OBWfnRDMCumulantTwo is the cumulant of OBWfnRDMTwo
  class OBWfnRDMCumulantTwo : public RDMCumulant<Two> {
    public:
      OBWfnRDMCumulantTwo(const Ref<OBWfnRDMTwo>& density);
      OBWfnRDMCumulantTwo(StateIn& si);
      virtual ~OBWfnRDMCumulantTwo();
      void save_data_state(StateOut& so);

      void compute();
      const double* obtain_block(SpinCase2 spin, size_t bra) const;
      void release_block(SpinCase2, size_t bra, double*) const;
      RefSymmSCMatrix scmat(SpinCase2 spincase) const;

    private:
      Ref<OBWfnRDMTwo> density_;
      mutable RefSymmSCMatrix scmat_[NSpinCases2];

      static ClassDesc class_desc_;
  };

#if 0
  /// WfnRDMOne is a 1-RDM from a Wavefunction
  /// This implementation assumes that the AO-basis density matrix is available for this Wavefunction
  class WfnRDMOne : public RDM<One> {
      typedef RDMCumulant<One> cumulant_type;
    public:
    /** A KeyVal constructor is used to generate a WfnRDMOne
        object from the input. The full list of keywords
        that are accepted is below.

        <table border="1">

        <tr><td>%Keyword<td>Type<td>Default<td>Description

        <tr><td><tt>wfn</tt><td>Wavefunction<td>none<td>the Wavefunction object

        </table>
     */
      WfnRDMOne(const Ref<KeyVal>& kv);
      WfnRDMOne(StateIn& si);
      WfnRDMOne(const Ref<Wavefunction>& wfn);
      ~WfnRDMOne();
      void save_data_state(StateOut& so);

      Ref<Wavefunction> wfn() const { return wfn_; }
      void compute();
      size_t ndim(SpinCase1 spincase) const;
      const double* obtain_block(SpinCase1 spin, size_t bra) const;
      void release_block(SpinCase1, size_t bra, double*) const;
      RefSymmSCMatrix scmat(SpinCase1 spincase) const;
      Ref<cumulant_type> cumulant() const;
      Ref< RDM<Zero> > rdm_m_1() const;

    private:
      Ref<Wavefunction> wfn_;
      mutable RefSymmSCMatrix scmat_[NSpinCases1];

      static ClassDesc class_desc_;
  };

  /// WfnRDMCumulantOne is the cumulant of WfnRDMOne
  class WfnRDMCumulantOne : public RDMCumulant<One> {
    public:
      WfnRDMCumulantOne(const Ref<WfnRDMOne>& density);
      WfnRDMCumulantOne(StateIn& si);
      ~WfnRDMCumulantOne();
      void save_data_state(StateOut& so);

      void compute();
      const double* obtain_block(SpinCase1 spin, size_t bra) const;
      void release_block(SpinCase1, size_t bra, double*) const;
      RefSymmSCMatrix scmat(SpinCase1 spincase) const;

    private:
      Ref<WfnRDMOne> density_;
      mutable RefSymmSCMatrix scmat_[NSpinCases1];

      static ClassDesc class_desc_;
  };
#endif

  /// OBWfnRDMOne is a 1-RDM from a OneBodyWavefunction
  class OBWfnRDMOne : public RDM<One> {
    public:
    /** A KeyVal constructor is used to generate a OBWfnRDMOne
        object from the input. The full list of keywords
        that are accepted is below.

        <table border="1">

        <tr><td>%Keyword<td>Type<td>Default<td>Description

        <tr><td><tt>wfn</tt><td>OneBodyWavefunction<td>none<td>the OneBodyWavefunction object

        </table>
     */
      OBWfnRDMOne(const Ref<KeyVal>& kv);
      OBWfnRDMOne(StateIn& si);
      OBWfnRDMOne(const Ref<OneBodyWavefunction>& wfn);
      virtual ~OBWfnRDMOne();
      void save_data_state(StateOut& so);

      Ref<OneBodyWavefunction> wfn() const { return wfn_; }
      Ref<OrbitalSpace> orbs(SpinCase1 s) const;

    private:
      Ref<OneBodyWavefunction> wfn_;
      mutable Ref<OrbitalSpace> orbs_[NSpinCases1];

      static ClassDesc class_desc_;
  };

  /// @}
  // end of addtogroup ChemistryElectronicStructure

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
