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

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_rdm_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_rdm_h

#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/mbptr12/spin.h>

namespace sc {

  /// Rank of the RDM
  typedef enum {One = 1, Two = 2, Three = 3, Four = 4} Rank;

  namespace {
    template <Rank R> struct __spincase;
    template <> struct __spincase<One> {
      typedef SpinCase1 type;
    };
    template <> struct __spincase<Two> {
      typedef SpinCase2 type;
    };
  }

  template <Rank R> class RDMCumulant;

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

          <tr><td><tt>psiwfn</tt><td>Wavefunction<td>none<td>the Wavefunction object

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
      ~RDM() {
      }
      void save_data_state(StateOut& so) {
        SavableState::save_state(wfn_.pointer(), so);
      }

      /// the corresponding Wavefunction
      Ref<Wavefunction> wfn() const { return wfn_; }
      /// bra/ket dimension
      virtual size_t ndim(spincase spincase) const =0;
      /// returns the ket block for the given bra index
      virtual const double* obtain_block(spincase spin,  size_t bra) const =0;
      /// releases the ket block
      virtual void release_block(spincase spin, size_t bra, double*) const =0;
      /// full density matrix
      virtual RefSymmSCMatrix scmat(spincase spin) const =0;

      /// cumulant of rank R
      virtual Ref< RDMCumulant<R> > cumulant() const =0;
      /// RDM of rank decreased by 1
      //virtual Ref< RDM< static_cast<Rank>(R-1) > > rdm_m_1() const =0;

    private:
      static ClassDesc class_desc_;
      Ref<Wavefunction> wfn_;
  };

  template <Rank R>
  ClassDesc
  RDM<R>::class_desc_(typeid(this_type),
                      (std::string("RDM<") +
                       char("1" + R - 1) +
                       std::string(">")
                      ).c_str(), 1,
                      "virtual public SavableState", 0, 0, 0 );

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

      /// the corresponding Wavefunction
      Ref<Wavefunction> wfn() const { return density_->wfn(); }
      /// the corresponding Density
      Ref<density_type> density() const { return density_; }
      /// bra/ket dimension
      size_t ndim(spincase spincase) const { return density->ndim(); }
      /// returns the ket block for the given bra index
      virtual const double* obtain_block(spincase spin,  size_t bra) const =0;
      /// releases the ket block
      virtual void release_block(spincase spin, size_t bra, double*) const =0;
      /// full cumulant matrix
      virtual RefSymmSCMatrix scmat(spincase spin) const =0;

    private:
      static ClassDesc class_desc_;
      Ref<density_type> density_;
  };

  template <Rank R>
  ClassDesc
  RDMCumulant<R>::class_desc_(typeid(this_type),
			      (std::string("RDMCumulant<") +
			       char("1" + R - 1) +
			       std::string(">")
			      ).c_str(), 1,
			      "virtual public SavableState", 0, 0, 0 );


} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
