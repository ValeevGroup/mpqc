//
// permute2e.cc
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <chemistry/qc/cints/int2e.h>
#include <chemistry/qc/cints/macros.h>
#ifdef DMALLOC
#include <dmalloc.h>
#endif

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}
inline void fail()
{
  ExEnv::errn() << scprintf("failing module:\n%s",__FILE__) << endl;
  abort();
}

void Int2eCints::sort_contrquartets_to_shellquartet_(double *source_ints_buf, double *target_ints_buf)
{
  int target_bf1_offset = 0;
  int nbf2 = int_shell2_->nfunction();
  int nbf3 = int_shell3_->nfunction();
  int nbf4 = int_shell4_->nfunction();
  int nbf234 = nbf2*nbf3*nbf4;
  int nbf34 = nbf3*nbf4;
  for (int gc1=0; gc1<int_shell1_->ncontraction(); gc1++) {
    int am1 = int_shell1_->am(gc1);
    int tsize1 = int_shell1_->nfunction(gc1);

    int target_bf2_offset = 0;
    for (int gc2=0; gc2<int_shell2_->ncontraction(); gc2++) {
      int am2 = int_shell2_->am(gc2);
      int tsize2 = int_shell2_->nfunction(gc2);

      int target_bf3_offset = 0;
      for (int gc3=0; gc3<int_shell3_->ncontraction(); gc3++) {
	int am3 = int_shell3_->am(gc3);
	int tsize3 = int_shell3_->nfunction(gc3);
	
	int target_bf4_offset = 0;
	for (int gc4=0; gc4<int_shell4_->ncontraction(); gc4++) {
	  int am4 = int_shell4_->am(gc4);
	  int tsize4 = int_shell4_->nfunction(gc4);

	  double *target_offset = target_ints_buf +
	    ((target_bf1_offset*nbf2 +
	      target_bf2_offset)*nbf3 +
	     target_bf3_offset)*nbf4 +
	    target_bf4_offset;
	  for (int bf1 = 0; bf1 < tsize1; bf1++, target_offset+=(nbf234-tsize2*nbf34)) {
	    for (int bf2 = 0; bf2 < tsize2; bf2++, target_offset+=(nbf34-tsize3*nbf4)) {
	      for (int bf3 = 0; bf3 < tsize3; bf3++, target_offset+=(nbf4-tsize4)) {
		for (int bf4 = 0; bf4 < tsize4; bf4++) {
	  
		  *(target_offset++) = *(source_ints_buf++);

		}
	      }
	    }
	  }

	  target_bf4_offset += tsize4;
	}
	target_bf3_offset += tsize3;
      }
      target_bf2_offset += tsize2;
    }
    target_bf1_offset += tsize1;
  }
}


//void Int2eCints::get_nonredundant_ints_(double *source, double *target, int e13e24, int e12, int e34)
//{
//  if (e13e24 && e12 && e34)
//    get_nonredundant_1111_(source,target);
//  else if (e13e24)
//    get_nonredundant_1212_(source,target);
//  else
//    get_nonredundant_1234_(source,target);
//}


// source can be safely same as target

void Int2eCints::get_nonredundant_ints_(double *source, double *target, int e13e24, int e12, int e34)
{
  int i,j,k,l;

  int nbf1 = int_shell1_->nfunction();
  int nbf2 = int_shell2_->nfunction();
  int nbf3 = int_shell3_->nfunction();
  int nbf4 = int_shell4_->nfunction();

  double *redundant_ptr = source;
  double *nonredundant_ptr = target;

  int nbf34 = nbf3*nbf4;
  for (i=0; i<nbf1; i++) {
    int jmax = e12 ? i : nbf2-1;
    for (j=0; j<=jmax; j++) {
      int kmax = e13e24 ? i : nbf3-1;
      for (k=0; k<=kmax; k++) {
	int lmax = e34 ? ( (e13e24&&(i==k)) ? j : k) : ( (e13e24&&(i==k)) ? j : nbf4-1);
        for (l=0; l<=lmax; l++) {
          *(nonredundant_ptr++) = redundant_ptr[l];
	}
        redundant_ptr += nbf4;
      }
      redundant_ptr += (nbf3-(kmax+1))*nbf4;
    }
    redundant_ptr += (nbf2-(jmax+1))*nbf34;
  }
}

void Int2eCints::permute_target_(double *s, double *t, int p13p24, int p12, int p34)
{
  if (!p13p24) {
    if (p12)
      if (p34)
	permute_1234_to_2143_(s,t);
      else
	permute_1234_to_2134_(s,t);
    else
      permute_1234_to_1243_(s,t);
  }
  else {
    if (p12)
      if (p34)
	permute_1234_to_4321_(s,t);
      else
	permute_1234_to_4312_(s,t);
    else
      if (p34)
	permute_1234_to_3421_(s,t);
      else
	permute_1234_to_3412_(s,t);
  }
}

void Int2eCints::permute_1234_to_1243_(double *source_ints_buf, double *target_ints_buf)
{
  int nbf1 = int_shell1_->nfunction();
  int nbf2 = int_shell2_->nfunction();
  int nbf3 = int_shell4_->nfunction();
  int nbf4 = int_shell3_->nfunction();
  for (int bf1=0; bf1<nbf1; bf1++) {
    for (int bf2=0; bf2<nbf2; bf2++) {
      for (int bf4=0; bf4<nbf4; bf4++) {
	for (int bf3=0; bf3<nbf3; bf3++) {
	  double *target_ints_ptr = target_ints_buf + ((bf1*nbf2 + bf2)*nbf3 + bf3)*nbf4 + bf4;
  	  *(target_ints_ptr) = *(source_ints_buf++);
	}
      }
    }
  }
}

void Int2eCints::permute_1234_to_2134_(double *source_ints_buf, double *target_ints_buf)
{
  int nbf1 = int_shell2_->nfunction();
  int nbf2 = int_shell1_->nfunction();
  int nbf3 = int_shell3_->nfunction();
  int nbf4 = int_shell4_->nfunction();
  for (int bf2=0; bf2<nbf2; bf2++) {
    for (int bf1=0; bf1<nbf1; bf1++) {
      for (int bf3=0; bf3<nbf3; bf3++) {
	for (int bf4=0; bf4<nbf4; bf4++) {
	  double *target_ints_ptr = target_ints_buf + ((bf1*nbf2 + bf2)*nbf3 + bf3)*nbf4 + bf4;
  	  *(target_ints_ptr) = *(source_ints_buf++);
	}
      }
    }
  }
}

void Int2eCints::permute_1234_to_2143_(double *source_ints_buf, double *target_ints_buf)
{
  int nbf1 = int_shell2_->nfunction();
  int nbf2 = int_shell1_->nfunction();
  int nbf3 = int_shell4_->nfunction();
  int nbf4 = int_shell3_->nfunction();
  for (int bf2=0; bf2<nbf2; bf2++) {
    for (int bf1=0; bf1<nbf1; bf1++) {
      for (int bf4=0; bf4<nbf4; bf4++) {
	for (int bf3=0; bf3<nbf3; bf3++) {
	  double *target_ints_ptr = target_ints_buf + ((bf1*nbf2 + bf2)*nbf3 + bf3)*nbf4 + bf4;
  	  *(target_ints_ptr) = *(source_ints_buf++);
	}
      }
    }
  }
}

void Int2eCints::permute_1234_to_3412_(double *source_ints_buf, double *target_ints_buf)
{
  int nbf1 = int_shell3_->nfunction();
  int nbf2 = int_shell4_->nfunction();
  int nbf3 = int_shell1_->nfunction();
  int nbf4 = int_shell2_->nfunction();
  for (int bf3=0; bf3<nbf3; bf3++) {
    for (int bf4=0; bf4<nbf4; bf4++) {
      for (int bf1=0; bf1<nbf1; bf1++) {
	for (int bf2=0; bf2<nbf2; bf2++) {
	  double *target_ints_ptr = target_ints_buf + ((bf1*nbf2 + bf2)*nbf3 + bf3)*nbf4 + bf4;
  	  *(target_ints_ptr) = *(source_ints_buf++);
	}
      }
    }
  }
}

void Int2eCints::permute_1234_to_4312_(double *source_ints_buf, double *target_ints_buf)
{
  int nbf1 = int_shell4_->nfunction();
  int nbf2 = int_shell3_->nfunction();
  int nbf3 = int_shell1_->nfunction();
  int nbf4 = int_shell2_->nfunction();
  for (int bf3=0; bf3<nbf3; bf3++) {
    for (int bf4=0; bf4<nbf4; bf4++) {
      for (int bf2=0; bf2<nbf2; bf2++) {
	for (int bf1=0; bf1<nbf1; bf1++) {
	  double *target_ints_ptr = target_ints_buf + ((bf1*nbf2 + bf2)*nbf3 + bf3)*nbf4 + bf4;
  	  *(target_ints_ptr) = *(source_ints_buf++);
	}
      }
    }
  }
}

void Int2eCints::permute_1234_to_3421_(double *source_ints_buf, double *target_ints_buf)
{
  int nbf1 = int_shell3_->nfunction();
  int nbf2 = int_shell4_->nfunction();
  int nbf3 = int_shell2_->nfunction();
  int nbf4 = int_shell1_->nfunction();
  for (int bf4=0; bf4<nbf4; bf4++) {
    for (int bf3=0; bf3<nbf3; bf3++) {
      for (int bf1=0; bf1<nbf1; bf1++) {
	for (int bf2=0; bf2<nbf2; bf2++) {
	  double *target_ints_ptr = target_ints_buf + ((bf1*nbf2 + bf2)*nbf3 + bf3)*nbf4 + bf4;
  	  *(target_ints_ptr) = *(source_ints_buf++);
	}
      }
    }
  }
}

void Int2eCints::permute_1234_to_4321_(double *source_ints_buf, double *target_ints_buf)
{
  int nbf1 = int_shell4_->nfunction();
  int nbf2 = int_shell3_->nfunction();
  int nbf3 = int_shell2_->nfunction();
  int nbf4 = int_shell1_->nfunction();
  for (int bf4=0; bf4<nbf4; bf4++) {
    for (int bf3=0; bf3<nbf3; bf3++) {
      for (int bf2=0; bf2<nbf2; bf2++) {
	for (int bf1=0; bf1<nbf1; bf1++) {
	  double *target_ints_ptr = target_ints_buf + ((bf1*nbf2 + bf2)*nbf3 + bf3)*nbf4 + bf4;
  	  *(target_ints_ptr) = *(source_ints_buf++);
	}
      }
    }
  }
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
