//
// hcore.cc
//
// Copyright (C) 2001 Edward Valeev
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

#include <util/misc/math.h>

#include <chemistry/qc/libint2/int1e.h>
#include <chemistry/qc/libint2/macros.h>

using namespace sc;

void Int1eLibint2::hcore(int sh1, int sh2)
{
  zero_buffers_();
  compute_doublet_info_(sh1, sh2);

  int maxam1 = int_shell1_->max_am();
  int minam1 = int_shell1_->min_am();
  int maxam2 = int_shell2_->max_am();
  int minam2 = int_shell2_->min_am();
  
  if (maxam1 != minam1 || maxam2 != minam2) {
    // fail();
    hcore_full_general_();
  }
  else {
    hcore_full_general_();
  }
}


void Int1eLibint2::hcore_full_general_()
{
  int maxam1 = int_shell1_->max_am();
  int maxam2 = int_shell2_->max_am();
  int z1weight = 1;
  int y1weight = maxam1 + 1;
  int x1weight = y1weight * y1weight; 
  int z2weight = 1;
  int y2weight = maxam2 + 1;
  int x2weight = y2weight * y2weight; 

  /* See if need to transform to spherical harmonics */
  bool need_cart2sph_transform = false;
  if (int_shell1_->has_pure() ||
      int_shell2_->has_pure())
    need_cart2sph_transform = true;

  /* See if contraction quartets need to be resorted into a shell quartet */
  bool need_sort_to_shell_doublet = false;
  int num_gen_shells = 0;
  if (int_shell1_->ncontraction() > 1)
    num_gen_shells++;
  if (int_shell2_->ncontraction() > 1)
    num_gen_shells++;
  if (maxam1 + maxam2 && num_gen_shells >= 1)
    need_sort_to_shell_doublet = true;

  /* Determine where integrals need to go at each stage */
  if (need_sort_to_shell_doublet) {
      prim_ints_ = cart_ints_;
      if (need_cart2sph_transform)
	contr_doublets_ = sphharm_ints_;
      else
	contr_doublets_ = cart_ints_;
      shell_doublet_ = target_ints_buffer_;
    }
    else {
      if (need_cart2sph_transform) {
	prim_ints_ = cart_ints_;
	contr_doublets_ = target_ints_buffer_;
	shell_doublet_ = target_ints_buffer_;
      }
      else {
	prim_ints_ = target_ints_buffer_;
	shell_doublet_ = target_ints_buffer_;
      }
    }

  /* Begin loops over primitives. */
  for (int p1=0; p1<int_shell1_->nprimitive(); p1++) {
    double a1 = int_shell1_->exponent(p1);
    for (int p2=0; p2<int_shell2_->nprimitive(); p2++) {
      double a2 = int_shell2_->exponent(p2);
      
      double gamma = a1+a2;
      double oog = 1.0/gamma;
      double over_pf = exp(-a1*a2*doublet_info_.AB2*oog)*sqrt(M_PI*oog)*M_PI*oog;
      
      double P[3], PA[3], PB[3], PC[3];
      for(int xyz=0; xyz<3; xyz++) {
	P[xyz] = (a1*doublet_info_.A[xyz] + a2*doublet_info_.B[xyz])*oog;
	PA[xyz] = P[xyz] - doublet_info_.A[xyz];
	PB[xyz] = P[xyz] - doublet_info_.B[xyz];
      }

      OI_OSrecurs_(OIX_,OIY_,OIZ_,PA,PB,gamma,maxam1+2,maxam2+2);

      /*--- contract each buffer into appropriate location ---*/
      double *ints_buf = prim_ints_;
      for (int gc1=0; gc1<int_shell1_->ncontraction(); gc1++) {
	double norm1 = int_shell1_->coefficient_unnorm(gc1,p1);
	int am1 = int_shell1_->am(gc1);
	for (int gc2=0; gc2<int_shell2_->ncontraction(); gc2++) {
	  double norm2 = int_shell2_->coefficient_unnorm(gc2,p2);
	  int am2 = int_shell2_->am(gc2);
	  double total_pf = over_pf * norm1 * norm2;
	  
	  int k1,l1,m1,k2,l2,m2;
	  FOR_CART(k1,l1,m1,am1)
	    FOR_CART(k2,l2,m2,am2)
	      double x0 = OIX_[k1][k2];
	      double y0 = OIY_[l1][l2];
	      double z0 = OIZ_[m1][m2];
	      double tx = a2*(2*k2+1)*OIX_[k1][k2] - 2*a2*a2*OIX_[k1][k2+2];
	      if (k2 >= 2)
		tx -= 0.5*k2*(k2-1)*OIX_[k1][k2-2];
	      double ty = a2*(2*l2+1)*OIY_[l1][l2] - 2*a2*a2*OIY_[l1][l2+2];
	      if (l2 >= 2)
		ty -= 0.5*l2*(l2-1)*OIY_[l1][l2-2];
	      double tz = a2*(2*m2+1)*OIZ_[m1][m2] - 2*a2*a2*OIZ_[m1][m2+2];
	      if (m2 >= 2)
		tz -= 0.5*m2*(m2-1)*OIZ_[m1][m2-2];
	      *(ints_buf++) += total_pf*(tx*y0*z0 + x0*ty*z0 + x0*y0*tz);
	    END_FOR_CART
	  END_FOR_CART

	}
      }

      if (bs1_->molecule() != bs2_->molecule()) {
	//	    fail();
      }

      int natom = bs1_->ncenter();
      for(int atom=0; atom<natom; atom++) {
	// if charge is 0 - skip to the next one
        double Z = bs1_->molecule()->charge(atom);
	if (Z == 0.0)
	  continue;
	PC[0] = P[0] - bs1_->r(atom,0);
	PC[1] = P[1] - bs1_->r(atom,1);
	PC[2] = P[2] - bs1_->r(atom,2);
	AI_OSrecurs_<0>(PA,PB,PC,gamma,maxam1,maxam2);

	/*--- contract each buffer into appropriate location ---*/
	double *ints_buf = prim_ints_;
	for (int gc1=0; gc1<int_shell1_->ncontraction(); gc1++) {
	  double norm1 = int_shell1_->coefficient_unnorm(gc1,p1);
	  int am1 = int_shell1_->am(gc1);
	  for (int gc2=0; gc2<int_shell2_->ncontraction(); gc2++) {
	    double norm2 = int_shell2_->coefficient_unnorm(gc2,p2);
	    int am2 = int_shell2_->am(gc2);
	    double total_pf = over_pf * norm1 * norm2 * Z;
	    
	    int k1,l1,m1,k2,l2,m2;
	    FOR_CART(k1,l1,m1,am1)
	      int ind1 = k1*x1weight + l1*y1weight + m1*z1weight;
	      FOR_CART(k2,l2,m2,am2)
		int ind2 = k2*x2weight + l2*y2weight + m2*z2weight;
	        *ints_buf -= AI0_[ind1][ind2][0] * total_pf;
		ints_buf++;
	      END_FOR_CART
	    END_FOR_CART
	  
	  }
	}
      }
      
    }
  }

  if (need_cart2sph_transform)
    transform_contrquartets_(prim_ints_,contr_doublets_);

  // If not CCA-compliant normalization -- re-normalize all integrals
#if INTEGRALLIBINT2_NORMCONV != INTEGRALLIBINT2_NORMCONV_CCA
  norm_contrcart1_(need_cart2sph_transform ? contr_doublets_ : prim_ints_);
#endif

  if (need_sort_to_shell_doublet)
    sort_contrdoublets_to_shelldoublet_(contr_doublets_,shell_doublet_);
}


void Int1eLibint2::hcore_sameam_general_()
{
  int tam1 = int_shell1_->am(0);
  int tam2 = int_shell2_->am(0);
  int z1weight = 1;
  int y1weight = tam1 + 1;
  int x1weight = y1weight * y1weight; 
  int z2weight = 1;
  int y2weight = tam2 + 1;
  int x2weight = y2weight * y2weight; 

  /* Begin loops over primitives. */
  for (int p1=0; p1<int_shell1_->nprimitive(); p1++) {
    double a1 = int_shell1_->exponent(p1);
    for (int p2=0; p2<int_shell2_->nprimitive(); p2++) {
      double a2 = int_shell2_->exponent(p2);
      
      double gamma = a1+a2;
      double oog = 1.0/gamma;
      double over_pf = exp(-a1*a2*doublet_info_.AB2*oog)*sqrt(M_PI*oog)*M_PI*oog;
      
      double P[3], PA[3], PB[3], PC[3];
      for(int xyz=0; xyz<3; xyz++) {
	P[xyz] = (a1*doublet_info_.A[xyz] + a2*doublet_info_.B[xyz])*oog;
	PA[xyz] = P[xyz] - doublet_info_.A[xyz];
	PB[xyz] = P[xyz] - doublet_info_.B[xyz];
      }

      OI_OSrecurs_(OIX_,OIY_,OIZ_,PA,PB,gamma,tam1+2,tam2+2);

      /*--- contract each buffer into appropriate location ---*/
      double *ints_buf = cart_ints_;
      for (int gc1=0; gc1<int_shell1_->ncontraction(); gc1++) {
	double norm1 = int_shell1_->coefficient_unnorm(gc1,p1);
	for (int gc2=0; gc2<int_shell2_->ncontraction(); gc2++) {
	  double norm2 = int_shell2_->coefficient_unnorm(gc2,p2);
	  double total_pf = over_pf * norm1 * norm2;
	  
	  int k1,l1,m1,k2,l2,m2;
	  FOR_CART(k1,l1,m1,tam1)
	    FOR_CART(k2,l2,m2,tam2)
	      double x0 = OIX_[k1][k2];
	      double y0 = OIY_[l1][l2];
	      double z0 = OIZ_[m1][m2];
	      double tx = a2*(2*k2+1)*OIX_[k1][k2] - 2*a2*a2*OIX_[k1][k2+2];
	      if (k2 >= 2)
		tx -= 0.5*k2*(k2-1)*OIX_[k1][k2-2];
	      double ty = a2*(2*l2+1)*OIY_[l1][l2] - 2*a2*a2*OIY_[l1][l2+2];
	      if (l2 >= 2)
		ty -= 0.5*l2*(l2-1)*OIY_[l1][l2-2];
	      double tz = a2*(2*m2+1)*OIZ_[m1][m2] - 2*a2*a2*OIZ_[m1][m2+2];
	      if (m2 >= 2)
		tz -= 0.5*m2*(m2-1)*OIZ_[m1][m2-2];
	      *(ints_buf++) += total_pf*(tx*y0*z0 + x0*ty*z0 + x0*y0*tz);
	    END_FOR_CART
	  END_FOR_CART

	}
      }

      if (bs1_->molecule() != bs2_->molecule()) {
	//	    fail();
      }

      int natom = bs1_->ncenter();
      for(int atom=0; atom<natom; atom++) {
	// if charge is 0 - skip to the next one
	double Z = bs1_->molecule()->charge(atom);
	if (Z == 0.0)
	  continue;
	PC[0] = P[0] - bs1_->r(atom,0);
	PC[1] = P[1] - bs1_->r(atom,1);
	PC[2] = P[2] - bs1_->r(atom,2);
	AI_OSrecurs_<0>(PA,PB,PC,gamma,tam1,tam2);

	/*--- contract each buffer into appropriate location ---*/
	double *ints_buf = cart_ints_;
	for (int gc1=0; gc1<int_shell1_->ncontraction(); gc1++) {
	  double norm1 = int_shell1_->coefficient_unnorm(gc1,p1);
	  for (int gc2=0; gc2<int_shell2_->ncontraction(); gc2++) {
	    double norm2 = int_shell2_->coefficient_unnorm(gc2,p2);
	    double total_pf = over_pf * norm1 * norm2 * Z;
	    
	    int k1,l1,m1,k2,l2,m2;
	    FOR_CART(k1,l1,m1,tam1)
	      int ind1 = k1*x1weight + l1*y1weight + m1*z1weight;
	      FOR_CART(k2,l2,m2,tam2)
		int ind2 = k2*x2weight + l2*y2weight + m2*z2weight;
	        *ints_buf -= AI0_[ind1][ind2][0] * total_pf;
		ints_buf++;
	      END_FOR_CART
	    END_FOR_CART
	  
	  }
	}
      }
      
    }
  }

  /*----------------------------------------------------------------------
    transform to spherical harmonics and/or resort to the target ordering
   ----------------------------------------------------------------------*/

  /*--- sort to the target ordering ---*/
  double *source_ints_buf = cart_ints_;
  double *target_ints_buf = target_ints_buffer_;
  int target_bf1_offset = 0;
  int target_bf2_offset = 0;
  int nbf2 = int_shell2_->nfunction();
  for (int gc1=0; gc1<int_shell1_->ncontraction(); gc1++) {
    int tsize1 = INT_NCART_NN(tam1);
    for (int gc2=0; gc2<int_shell2_->ncontraction(); gc2++) {
      int tsize2 = INT_NCART_NN(tam2);
      
      int k1,l1,m1,k2,l2,m2;
      int bf1 = 0;
      FOR_CART(k1,l1,m1,tam1)
	double *target_ints_buf = target_ints_buffer_ + (target_bf1_offset+bf1)*nbf2 +
	                          target_bf2_offset;
        FOR_CART(k2,l2,m2,tam2)
  	  *(target_ints_buf++) = *(source_ints_buf++);
        END_FOR_CART
	bf1++;
      END_FOR_CART
      target_bf2_offset += tsize2;
    }
    target_bf1_offset += tsize1;
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
