//
// twoparticlecontraction.h
//
// Copyright (C) 2005 Edward Valeev
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

#include <util/state/state.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>

#ifndef _chemistry_qc_mbptr12_twoparticlecontraction_h
#define _chemistry_qc_mbptr12_twoparticlecontraction_h

namespace sc {
  
  namespace mbptr12 {
    
    /** TwoParticleContraction contracts nrow-by-ncol bra- or ket-blocks of two 2-particle tensors
        i.e. each ij element of the first tensor is multiplied with the ij element of the ij element
        of the second tensor and all ij-ij products are summed.
    */
    class TwoParticleContraction : virtual public SavableState {
      public:
      TwoParticleContraction(unsigned int nrow, unsigned int ncol);
      TwoParticleContraction(StateIn& si);
      virtual ~TwoParticleContraction() {}
      
      void save_data_state(StateOut& so);
      
      unsigned int nrow() const;
      unsigned int ncol() const;
      
      /// Computes contraction of blocks A and B
      virtual double contract(const double* A, const double* B) const =0;
      
      protected:
      double dot_prod(const double* A, const double* B) const;
      
      private:
      unsigned int nrow_;
      unsigned int ncol_;
    };
    
    /** Direct_Contraction is a straight scalar (dot) product of 2 rectangular blocks, scaled by scale.
        Does the same as TwoParticleContraction but also scales the result.
    */
    class Direct_Contraction : public TwoParticleContraction {
      public:
      Direct_Contraction(unsigned int nrow, unsigned int ncol, double scale);
      Direct_Contraction(StateIn& si);
      ~Direct_Contraction() {}
      void save_data_state(StateOut& so);
      
      /// Computes contraction of blocks A and B
      double contract(const double* A, const double* B) const;
      
      private:
      double scale_;
    };
    
    /** ABS_OBS_Contraction contracts 2 square nobs-by-nobs blocks
        for the ABS approach.
        Contracts only the occ-occ and vir-vir subblocks of the
        obs-obs blocks.
    */
    class ABS_OBS_Contraction : public TwoParticleContraction {
      public:
      /// nobs = rank of OBS (number of all orbitals), nocc1 and nocc2 are number of occupied orbitals for spaces of electron 1 and 2
      ABS_OBS_Contraction(unsigned int nobs, unsigned int nocc1, unsigned int nocc2);
      ABS_OBS_Contraction(StateIn&);
      ~ABS_OBS_Contraction() {}
      void save_data_state(StateOut&);
      /// Computes contraction of blocks A and B
      double contract(const double* A, const double* B) const;
      
      private:
      unsigned int nocc1_;
      unsigned int nocc2_;
    };
    
    /** CABS_OBS_Contraction contracts 2 square nobs-by-nobs blocks
        for the CABS approach. Does effectively the same as
        TwoParticleContraction, but multiplies the result by -1.
    */
    class CABS_OBS_Contraction : public TwoParticleContraction {
      public:
      /// CABS OBS contraction does not depend on the occupied spaces
      CABS_OBS_Contraction(unsigned int nobs);
      CABS_OBS_Contraction(StateIn&);
      ~CABS_OBS_Contraction() {}
      void save_data_state(StateOut&);
      /// Computes contraction of blocks A and B
      double contract(const double* A, const double* B) const;
    };
  }
}

#endif

