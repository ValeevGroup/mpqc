//
// creator.h
//
// Copyright (C) 2005 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#ifndef _chemistry_qc_mbptr12_creator_h
#define _chemistry_qc_mbptr12_creator_h

#include <chemistry/qc/mbptr12/transform_tbint.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/linearr12.h>

namespace sc {
  
  /** RangeCreator<T> is Functor which can be used up to n times to create objects
      of type T. operator() returns the objects, or 0 when done. Thus T must be Comparable to
      an int or Constructable from an int. */
  template <typename T>
    class RangeCreator {
      public:
      RangeCreator(unsigned int n) :
        n_(n), ncreated_(0) {}
      /// returns a new object T, or T(0) if done
      virtual T operator()() =0;
      
      protected:
      bool can_create() const { return ncreated_ < n_; }
      void next() { ++ncreated_; }
      
      private:
      unsigned int n_;
      unsigned int ncreated_;
    };

  /** Creates TwoBodyMOIntsTransforms using transforms known by name to R12IntEval */
  class NamedTransformCreator : public RangeCreator< Ref<TwoBodyMOIntsTransform> >
  {
    public:
    typedef Ref<TwoBodyMOIntsTransform> ObjT;
    
    NamedTransformCreator(Ref<R12IntEval>& r12eval,
                          const Ref<MOIndexSpace>& space1,
                          const Ref<MOIndexSpace>& space2,
                          const Ref<MOIndexSpace>& space3,
                          const Ref<MOIndexSpace>& space4,
                          bool CorrFunctionInBra = false,
                          bool CorrFunctionInKet = false);
    /// Implementation of RangeCreator::operator()
    ObjT operator()();
    
    private:
    Ref<R12IntEval> r12eval_;
    Ref<MOIndexSpace> space1_;
    Ref<MOIndexSpace> space2_;
    Ref<MOIndexSpace> space3_;
    Ref<MOIndexSpace> space4_;
    bool CorrFunctionInBraKet_;
    unsigned int nf12bra_;
    unsigned int nf12ket_;
    unsigned int braindex_;
    unsigned int ketindex_;
    
    void increment_indices();
  };
  
  /** Creates new TwoBodyMOIntsTransforms and adds them to the transform map (R12IntEval, at the moment) */
  class NewTransformCreator : public RangeCreator< Ref<TwoBodyMOIntsTransform> >
  {
    public:
    typedef Ref<TwoBodyMOIntsTransform> ObjT;
    
    NewTransformCreator(Ref<R12IntEval>& r12eval,
                        const Ref<MOIndexSpace>& space1,
                        const Ref<MOIndexSpace>& space2,
                        const Ref<MOIndexSpace>& space3,
                        const Ref<MOIndexSpace>& space4,
                        bool CorrFunctionInBra = false,
                        bool CorrFunctionInKet = false,
                        MOIntsTransformFactory::StorageType storage =
                          MOIntsTransformFactory::StorageType_13);
    /// Implementation of RangeCreator::operator()
    ObjT operator()();
    
    private:
    Ref<R12IntEval> r12eval_;
    Ref<MOIntsTransformFactory> tfactory_;
    Ref<MOIndexSpace> space1_;
    Ref<MOIndexSpace> space2_;
    Ref<MOIndexSpace> space3_;
    Ref<MOIndexSpace> space4_;
    MOIntsTransformFactory::StorageType storage_;
    bool CorrFunctionInBraKet_;
    bool CorrFunction_;
    unsigned int nf12bra_;
    unsigned int nf12ket_;
    unsigned int braindex_;
    unsigned int ketindex_;
    
    void increment_indices();
  };
  
  using LinearR12::CorrelationFactor;
  /** Creates TwoBodyIntDescr for correlation factpr C */
  class TwoBodyIntDescrCreator : public RangeCreator< Ref<TwoBodyIntDescr> >
  {
    public:
    typedef Ref<TwoBodyIntDescr> ObjT;
    
    TwoBodyIntDescrCreator(const Ref<CorrelationFactor>& corrfactor,
                           const Ref<Integral>& integral,
                           bool CorrFunctionInBra = false,
                           bool CorrFunctionInKet = false);
    /// Implementation of RangeCreator::operator()
    ObjT operator()();
    
    private:
    Ref<CorrelationFactor> corrfactor_;
    Ref<Integral> integral_;
    bool CorrFunctionInBraKet_;
    unsigned int nf12bra_;
    unsigned int nf12ket_;
    unsigned int braindex_;
    unsigned int ketindex_;
    
    void increment_indices();
  };
  
  
}

#endif

