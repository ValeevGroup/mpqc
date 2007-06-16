//
// fjt.h
//
// Copyright (C) 2001 Edward Valeev
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

#ifndef _chemistry_qc_libint2_fjt_h
#define _chemistry_qc_libint2_fjt_h

#include <util/ref/ref.h>

class Taylor_Fjt_Eval : public RefCount {
  private:
    double **grid;            /* Table of "exact" Fm(T) values. Row index corresponds to
				 values of T (max_T+1 rows), column index to values
				 of m (max_m+1 columns) */
    double delT;              /* The step size for T, depends on cutoff */
    double cutoff;            /* Tolerance cutoff used in all computations of Fm(T) */
    int order_interp;         /* Order of (Taylor) interpolation */
    int max_m;                /* Maximum value of m in the table, depends on cutoff
				 and the number of terms in Taylor interpolation */
    int max_T;                /* Maximum index of T in the table, depends on cutoff
			         and m */
    double *T_crit;           /* Maximum T for each row, depends on cutoff;
				 for a given m and T_idx <= max_T_idx[m] use Taylor interpolation,
				 for a given m and T_idx > max_T_idx[m] use the asymptotic formula */
    double *Fjt_buffer;       /* Here computed values of Fj(T) are stored */

  public:
    Taylor_Fjt_Eval(unsigned int mmax, double accuracy);
    ~Taylor_Fjt_Eval();
    double *compute_Fjt(double T, unsigned int l);  /* The function which computes a set of Fm(T), 0<=m<=l
						       for given T and l */
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
