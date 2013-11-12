//
// femo.h
//
// Copyright (C) 2008 Edward Valeev
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

#ifndef _chemistry_qc_wfn_femo_h
#define _chemistry_qc_wfn_femo_h

#include <vector>
#include <math/scmat/matrix.h>

namespace sc {

  /// @addtogroup ChemistryElectronicStructureOneBody
  /// @{

  /** Describes a simple the free-electron molecular orbital model that can be
      used to guess the lowest-energy orbital configuration. */
  class FEMO : public RefCount {
    public:
      /** Construct a FEMO state using the eigenvalues. If beta eigenvalues are not specified,
          alpha eigenvalues will be used for both spincases. The eigenvalues must be provided in blocked form. */
      FEMO(int nalpha, int nbeta,
           const RefDiagSCMatrix& evalsa, const RefDiagSCMatrix& evalsb = 0);
      ~FEMO() {}

      /// returns the energy
      double E() const;
      /// returns the number of alpha electrons
      int nalpha() const;
      /// returns the number of beta electrons
      int nbeta() const;
      /// returns the number of alpha electrons in irrep i
      int nalpha(int i) const;
      /// returns the number of beta electrons in irrep i
      int nbeta(int i) const;
      
    private:
      FEMO() {}
      
      double E_;
      int nalpha_;
      int nbeta_;
      // # of electrons per irrep
      std::vector<int> alphapi_;
      std::vector<int> betapi_;
      
  };
  
  /** Finds the FEMO configuration that corresponds to the maximum multiplicity.
      Configuration A is preferred to B if it has higher muliplicity AND
      its energy is greater than that of B by at most Etol.
      if allow_closed_shell is set to false, only try open-shell configurations
       */
  class HundsFEMOSeeker {
    public:
      HundsFEMOSeeker(int nelectron, double Etol, bool allow_closed_shell,
                      const RefDiagSCMatrix& evalsa, const RefDiagSCMatrix& evalsb = 0);
      ~HundsFEMOSeeker() {}
      const Ref<FEMO>& result() const { return result_; }

      // the default tolerance
      static double tolerance;
    private:
      Ref<FEMO> result_;
  };
  
  /// @}
  // end of addtogroup ChemistryElectronicStructureOneBody

}

#endif /* header guard */
