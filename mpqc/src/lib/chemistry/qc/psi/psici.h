//
// psici.h
//
// Copyright (C) 2008 Martin Torheyden
//
// Author: Martin Torheyden <mtorhey@vt.edu>
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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _chemistry_qc_psi_psici_h
#define _chemistry_qc_psi_psici_h

#include <chemistry/qc/psi/psiwfn.h>
#include <chemistry/qc/mbptr12/spin.h>

namespace sc {

  ///////////////////////////////////////////////////////////////////
  /// PsiCI is a Psi configuration interaction wavefunction (detci)
  
  class PsiCI : public PsiCorrWavefunction {
    private:
      /// Psi keywords
      bool eval_opdm_;
      bool eval_tpdm_;
      bool opdm_print_;
      bool tpdm_print_;
      int root_;          /// compute a specific root of the wave function
      int num_roots_;
      int ex_lvl_;        /// CI excitation level
      bool repl_otf_;     /// do CI string replacements on the fly. saves memory, but is slower.
      int maxiter_;       /// maximum number of CI iterations.
    protected:
      void write_input(int convergence);
    public:
      PsiCI(const Ref<KeyVal>&);
      PsiCI(StateIn&);
      ~PsiCI();
      void save_data_state(StateOut&);
      void compute();
  };
  
}

#endif /*_chemistry_qc_psi_psici_h*/
