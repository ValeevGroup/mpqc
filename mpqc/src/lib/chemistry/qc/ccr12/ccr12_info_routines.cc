//
// ccr12_info_routine.cc
//
// Copyright (C) 2008 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@theochem.uni-stuttgart.de>
// Maintainer: TS
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

#include <string>
#include <algorithm>
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/ccr12/ccr12_info.h>


using namespace sc;
using namespace std;

void CCR12_Info::smith_dgemm(const long i,const long j,const long k,
                             const double a,const double* da,const long l,
                             const double* db,const long m,const double b,
                             double* dc,const long n){
  const char chart='t';
  const char charn='n';
  const double one=1.0;
  const int i_=(const int)i;
  const int j_=(const int)j;
  const int k_=(const int)k;
  const int l_=(const int)l;
  const int m_=(const int)m;
  const int n_=(const int)n;
  F77_DGEMM(&chart,&charn,&i_,&j_,&k_,&a,da,&l_,db,&m_,&one,dc,&n_);
}


void CCR12_Info::restricted_2(const long inp1,const long inp2, 
                              long& out1,long& out2){
  
  if( restricted() && get_spin(inp1)+get_spin(inp2)==4L ) {
    out1=get_alpha(inp1);
    out2=get_alpha(inp2); 
  } else {
    out1=inp1;
    out2=inp2;
  }
} 

void CCR12_Info::restricted_4(const long inp1,const long inp2,const long inp3,const long inp4, 
                              long& out1,long& out2,long& out3,long& out4){
  
  if( restricted() && get_spin(inp1)+get_spin(inp2)+get_spin(inp3)+get_spin(inp4)==8L ) {
    out1=get_alpha(inp1);
    out2=get_alpha(inp2); 
    out3=get_alpha(inp3); 
    out4=get_alpha(inp4); 
  } else {
    out1=inp1;
    out2=inp2;
    out3=inp3;
    out4=inp4;
  }
} 


void CCR12_Info::restricted_6(const long inp1,const long inp2,const long inp3,
                              const long inp4,const long inp5,const long inp6, 
                              long& out1,long& out2,long& out3,
                              long& out4,long& out5,long& out6){
  
  if( restricted() && get_spin(inp1)+get_spin(inp2)+get_spin(inp3)
                     +get_spin(inp4)+get_spin(inp5)+get_spin(inp6)==12L ) {
    out1=get_alpha(inp1);
    out2=get_alpha(inp2); 
    out3=get_alpha(inp3); 
    out4=get_alpha(inp4); 
    out5=get_alpha(inp5); 
    out6=get_alpha(inp6); 
  } else {
    out1=inp1;
    out2=inp2;
    out3=inp3;
    out4=inp4;
    out5=inp5;
    out6=inp6;
  }
} 


void CCR12_Info::restricted_8(const long inp1,const long inp2,const long inp3,const long inp4,
                              const long inp5,const long inp6,const long inp7,const long inp8, 
                              long& out1,long& out2,long& out3,long& out4,
                              long& out5,long& out6,long& out7,long& out8){
  
  if( restricted() && get_spin(inp1)+get_spin(inp2)+get_spin(inp3)+get_spin(inp4)
                     +get_spin(inp5)+get_spin(inp6)+get_spin(inp7)+get_spin(inp8)==16L ) {
    out1=get_alpha(inp1);
    out2=get_alpha(inp2); 
    out3=get_alpha(inp3); 
    out4=get_alpha(inp4); 
    out5=get_alpha(inp5); 
    out6=get_alpha(inp6); 
    out7=get_alpha(inp7); 
    out8=get_alpha(inp8); 
  } else {
    out1=inp1;
    out2=inp2;
    out3=inp3;
    out4=inp4;
    out5=inp5;
    out6=inp6;
    out7=inp7;
    out8=inp8;
  }
} 


void CCR12_Info::sort_indices2(const double* unsorted,double* sorted,
                               const long a,const long b,
                               const int i,const int j,const double factor)
{ 
  // prototype
  if (i==0) {
    const int j0max=(const int)(a*b);
    for (int j0=0; j0<a*b; ++j0) sorted[j0]=unsorted[j0]*factor;
  } else { 
    int id[2];
    int jd[2] = {a, b};
    long iall=0;
    for(int j0=0;j0<(int)a;++j0){
      id[0]=j0;
      for(int j1=0;j1<(int)b;++j1,++iall){
        id[1]=j1;
        long ib=id[j]+jd[j]*id[i];
        sorted[ib]=unsorted[iall]*factor;
      }
    }
  }
  
  // is there any good code for matrix transposition?
} 


void CCR12_Info::sort_indices4(const double* unsorted,double* sorted,
                               const long a,const long b,const long c,const long d,
                               const int i,const int j,const int k,const int l,
                               const double factor)
{
  // prototype
  int id[4];
  int jd[4] = {a, b, c, d};

  long iall=0;
  for(int j0=0;j0<(int)a;++j0){
    id[0]=j0;
    for(int j1=0;j1<(int)b;++j1){
      id[1]=j1;
      for(int j2=0;j2<(int)c;++j2){
        id[2]=j2;
        for (int j3=0;j3<(int)d;++j3,++iall){
          id[3]=j3;
          long ib=id[l]+jd[l]*(id[k]+jd[k]*(id[j]+jd[j]*id[i]));
          sorted[ib]=unsorted[iall]*factor;
        }
      }
    } 
  } 
}


void CCR12_Info::sort_indices6(const double* unsorted,double* sorted,
                               const long a,const long b,const long c,
                               const long d,const long e,const long f,
                               const int i,const int j,const int k,
                               const int l,const int m,const int n,
                               const double factor)
{
  // prototype
  int id[6];
  int jd[6] = {a, b, c, d, e, f};

  int iall=0;
  for(int j0=0;j0<(int)a;++j0){
   id[0]=j0;
   for(int j1=0;j1<(int)b;++j1){
    id[1]=j1;
    for(int j2=0;j2<(int)c;++j2){
     id[2]=j2;
     for(int j3=0;j3<(int)d;++j3){
      id[3]=j3;
      for (int j4=0;j4<(int)e;++j4){
       id[4]=j4;
       for (int j5=0;j5<(int)f;++j5,++iall){
        id[5]=j5;
        long ib=id[n]+jd[n]*(id[m]+jd[m]*(id[l]+jd[l]*(id[k]+jd[k]*(id[j]+jd[j]*id[i]))));
        sorted[ib]=unsorted[iall]*factor;
       }
      }
     }
    } 
   }
  } 
}


void CCR12_Info::sort_indices8(const double* unsorted,double* sorted,
                               const long a,const long b,const long c,const long d,
                               const long e,const long f,const long g,const long h,
                               const int i,const int j,const int k,const int l,
                               const int m,const int n,const int o,const int p,
                               const double factor)
{
  // prototype
  int id[8];
  int jd[8] = {a, b, c, d, e, f, g, h};

  int iall=0;
  for(int j0=0;j0<(int)a;++j0){
   id[0]=j0;
   for(int j1=0;j1<(int)b;++j1){
    id[1]=j1;
    for(int j2=0;j2<(int)c;++j2){
     id[2]=j2;
     for(int j3=0;j3<(int)d;++j3){
      id[3]=j3;
      for(int j4=0;j4<(int)e;++j4){
       id[4]=j4;
       for(int j5=0;j5<(int)f;++j5){
        id[5]=j5;
        for(int j6=0;j6<(int)g;++j6){
         id[6]=j6;
         for(int j7=0;j7<(int)h;++j7,++iall){
          id[7]=j7;
          long ib=id[p]+jd[p]*(id[o]+jd[o]*(id[n]+jd[n]*(id[m]+jd[m]*(id[l]+jd[l]*(id[k]+jd[k]*(id[j]+jd[j]*id[i]))))));
          sorted[ib]=unsorted[iall]*factor;
         }
        }
       }
      }
     }
    } 
   }
  } 
}


void CCR12_Info::sort_indices_acc6(const double* unsorted,double* sorted,
                                   const long a,const long b,const long c,
                                   const long d,const long e,const long f,
                                   const int i,const int j,const int k,
                                   const int l,const int m,const int n,
                                   const double factor)
{
  // prototype
  int id[6];
  int jd[6] = {a, b, c, d, e, f};

  int iall=0;
  for(int j0=0;j0<(int)a;++j0){
   id[0]=j0;
   for(int j1=0;j1<(int)b;++j1){
    id[1]=j1;
    for(int j2=0;j2<(int)c;++j2){
     id[2]=j2;
     for(int j3=0;j3<(int)d;++j3){
      id[3]=j3;
      for (int j4=0;j4<(int)e;++j4){
       id[4]=j4;
       for (int j5=0;j5<(int)f;++j5,++iall){
        id[5]=j5;
        long ib=id[n]+jd[n]*(id[m]+jd[m]*(id[l]+jd[l]*(id[k]+jd[k]*(id[j]+jd[j]*id[i]))));
        sorted[ib]+=unsorted[iall]*factor;
       }
      }
     }
    } 
   }
  } 
}


void CCR12_Info::sort_indices_acc8(const double* unsorted,double* sorted,
                                   const long a,const long b,const long c,const long d,
                                   const long e,const long f,const long g,const long h,
                                   const int i,const int j,const int k,const int l,
                                   const int m,const int n,const int o,const int p,
                                   const double factor)
{
  // prototype
  int id[8];
  int jd[8] = {a, b, c, d, e, f, g, h};

  int iall=0;
  for(int j0=0;j0<(int)a;++j0){
   id[0]=j0;
   for(int j1=0;j1<(int)b;++j1){
    id[1]=j1;
    for(int j2=0;j2<(int)c;++j2){
     id[2]=j2;
     for(int j3=0;j3<(int)d;++j3){
      id[3]=j3;
      for(int j4=0;j4<(int)e;++j4){
       id[4]=j4;
       for(int j5=0;j5<(int)f;++j5){
        id[5]=j5;
        for(int j6=0;j6<(int)g;++j6){
         id[6]=j6;
         for(int j7=0;j7<(int)h;++j7,++iall){
          id[7]=j7;
          long ib=id[p]+jd[p]*(id[o]+jd[o]*(id[n]+jd[n]*(id[m]+jd[m]*(id[l]+jd[l]*(id[k]+jd[k]*(id[j]+jd[j]*id[i]))))));
          sorted[ib]+=unsorted[iall]*factor;
         }
        }
       }
      }
     }
    } 
   }
  } 
}


void CCR12_Info::transpose_1(double* d,const double* s,const int dim1,const int dim2){
  int dim1s8=(dim1/8)*8;
  int dim2s8=(dim2/8)*8;

  size_t dim1_2=dim1*2; 
  size_t dim1_3=dim1*3; 
  size_t dim1_4=dim1*4; 
  size_t dim1_5=dim1*5; 
  size_t dim1_6=dim1*6; 
  size_t dim1_7=dim1*7; 

  size_t dim2_2=dim2*2; 
  size_t dim2_3=dim2*3; 
  size_t dim2_4=dim2*4; 
  size_t dim2_5=dim2*5; 
  size_t dim2_6=dim2*6; 
  size_t dim2_7=dim2*7; 

//    d[i,j] <- s[j,i]
  for(int i=0;i<dim1s8;i+=8){
    size_t dim2_i=dim2*i;  

    for(int j=0;j<dim2s8;j+=8){
      size_t dim1_j=dim1*j;  

      d[i  +dim1_j       ]=s[j  +dim2_i       ];  
      d[i+1+dim1_j       ]=s[j  +dim2_i+dim2  ];  
      d[i+2+dim1_j       ]=s[j  +dim2_i+dim2_2];  
      d[i+3+dim1_j       ]=s[j  +dim2_i+dim2_3];  
      d[i+4+dim1_j       ]=s[j  +dim2_i+dim2_4];  
      d[i+5+dim1_j       ]=s[j  +dim2_i+dim2_5];  
      d[i+6+dim1_j       ]=s[j  +dim2_i+dim2_6];  
      d[i+7+dim1_j       ]=s[j  +dim2_i+dim2_7];  

      d[i  +dim1_j+dim1  ]=s[j+1+dim2_i       ];  
      d[i+1+dim1_j+dim1  ]=s[j+1+dim2_i+dim2  ];  
      d[i+2+dim1_j+dim1  ]=s[j+1+dim2_i+dim2_2];  
      d[i+3+dim1_j+dim1  ]=s[j+1+dim2_i+dim2_3];  
      d[i+4+dim1_j+dim1  ]=s[j+1+dim2_i+dim2_4];  
      d[i+5+dim1_j+dim1  ]=s[j+1+dim2_i+dim2_5];  
      d[i+6+dim1_j+dim1  ]=s[j+1+dim2_i+dim2_6];  
      d[i+7+dim1_j+dim1  ]=s[j+1+dim2_i+dim2_7];  

      d[i  +dim1_j+dim1_2]=s[j+2+dim2_i       ];  
      d[i+1+dim1_j+dim1_2]=s[j+2+dim2_i+dim2  ];  
      d[i+2+dim1_j+dim1_2]=s[j+2+dim2_i+dim2_2];  
      d[i+3+dim1_j+dim1_2]=s[j+2+dim2_i+dim2_3];  
      d[i+4+dim1_j+dim1_2]=s[j+2+dim2_i+dim2_4];  
      d[i+5+dim1_j+dim1_2]=s[j+2+dim2_i+dim2_5];  
      d[i+6+dim1_j+dim1_2]=s[j+2+dim2_i+dim2_6];  
      d[i+7+dim1_j+dim1_2]=s[j+2+dim2_i+dim2_7];  

      d[i  +dim1_j+dim1_3]=s[j+3+dim2_i       ];  
      d[i+1+dim1_j+dim1_3]=s[j+3+dim2_i+dim2  ];  
      d[i+2+dim1_j+dim1_3]=s[j+3+dim2_i+dim2_2];  
      d[i+3+dim1_j+dim1_3]=s[j+3+dim2_i+dim2_3];  
      d[i+4+dim1_j+dim1_3]=s[j+3+dim2_i+dim2_4];  
      d[i+5+dim1_j+dim1_3]=s[j+3+dim2_i+dim2_5];  
      d[i+6+dim1_j+dim1_3]=s[j+3+dim2_i+dim2_6];  
      d[i+7+dim1_j+dim1_3]=s[j+3+dim2_i+dim2_7];  

      d[i  +dim1_j+dim1_4]=s[j+4+dim2_i       ];  
      d[i+1+dim1_j+dim1_4]=s[j+4+dim2_i+dim2  ];  
      d[i+2+dim1_j+dim1_4]=s[j+4+dim2_i+dim2_2];  
      d[i+3+dim1_j+dim1_4]=s[j+4+dim2_i+dim2_3];  
      d[i+4+dim1_j+dim1_4]=s[j+4+dim2_i+dim2_4];  
      d[i+5+dim1_j+dim1_4]=s[j+4+dim2_i+dim2_5];  
      d[i+6+dim1_j+dim1_4]=s[j+4+dim2_i+dim2_6];  
      d[i+7+dim1_j+dim1_4]=s[j+4+dim2_i+dim2_7];  

      d[i  +dim1_j+dim1_5]=s[j+5+dim2_i       ];  
      d[i+1+dim1_j+dim1_5]=s[j+5+dim2_i+dim2  ];  
      d[i+2+dim1_j+dim1_5]=s[j+5+dim2_i+dim2_2];  
      d[i+3+dim1_j+dim1_5]=s[j+5+dim2_i+dim2_3];  
      d[i+4+dim1_j+dim1_5]=s[j+5+dim2_i+dim2_4];  
      d[i+5+dim1_j+dim1_5]=s[j+5+dim2_i+dim2_5];  
      d[i+6+dim1_j+dim1_5]=s[j+5+dim2_i+dim2_6];  
      d[i+7+dim1_j+dim1_5]=s[j+5+dim2_i+dim2_7];  

      d[i  +dim1_j+dim1_6]=s[j+6+dim2_i       ];  
      d[i+1+dim1_j+dim1_6]=s[j+6+dim2_i+dim2  ];  
      d[i+2+dim1_j+dim1_6]=s[j+6+dim2_i+dim2_2];  
      d[i+3+dim1_j+dim1_6]=s[j+6+dim2_i+dim2_3];  
      d[i+4+dim1_j+dim1_6]=s[j+6+dim2_i+dim2_4];  
      d[i+5+dim1_j+dim1_6]=s[j+6+dim2_i+dim2_5];  
      d[i+6+dim1_j+dim1_6]=s[j+6+dim2_i+dim2_6];  
      d[i+7+dim1_j+dim1_6]=s[j+6+dim2_i+dim2_7];  

      d[i  +dim1_j+dim1_7]=s[j+7+dim2_i       ];  
      d[i+1+dim1_j+dim1_7]=s[j+7+dim2_i+dim2  ];  
      d[i+2+dim1_j+dim1_7]=s[j+7+dim2_i+dim2_2];  
      d[i+3+dim1_j+dim1_7]=s[j+7+dim2_i+dim2_3];  
      d[i+4+dim1_j+dim1_7]=s[j+7+dim2_i+dim2_4];  
      d[i+5+dim1_j+dim1_7]=s[j+7+dim2_i+dim2_5];  
      d[i+6+dim1_j+dim1_7]=s[j+7+dim2_i+dim2_6];  
      d[i+7+dim1_j+dim1_7]=s[j+7+dim2_i+dim2_7];  
    } // end of loop j block 
    for(int j=dim2s8;j<dim2;++j) {
      size_t dim1_j=dim1*j;  

      d[i  +dim1_j]=s[j+dim2_i       ];
      d[i+1+dim1_j]=s[j+dim2_i+dim2  ];
      d[i+2+dim1_j]=s[j+dim2_i+dim2_2];
      d[i+3+dim1_j]=s[j+dim2_i+dim2_3];
      d[i+4+dim1_j]=s[j+dim2_i+dim2_4];
      d[i+5+dim1_j]=s[j+dim2_i+dim2_5];
      d[i+6+dim1_j]=s[j+dim2_i+dim2_6];
      d[i+7+dim1_j]=s[j+dim2_i+dim2_7];
    } // end of loop j
  } // end of loop i block

  for(int i=dim1s8;i<dim1;++i){
    size_t dim2_i=dim2*i;  

    for(int j=0;j<dim2s8;j+=8){
      size_t dim1_j=dim1*j;  

      d[i+dim1_j       ]=s[j  +dim2_i];
      d[i+dim1_j+dim1  ]=s[j+1+dim2_i];
      d[i+dim1_j+dim1_2]=s[j+2+dim2_i];
      d[i+dim1_j+dim1_3]=s[j+3+dim2_i];
      d[i+dim1_j+dim1_4]=s[j+4+dim2_i];
      d[i+dim1_j+dim1_5]=s[j+5+dim2_i];
      d[i+dim1_j+dim1_6]=s[j+6+dim2_i];
      d[i+dim1_j+dim1_7]=s[j+7+dim2_i];
    }
    for(int j=dim2s8;j<dim2;++j) d[i+dim1*j]=s[j+dim2_i];
  } // end of loop i
}

void CCR12_Info::transpose(double* d,const double* s,const int dim1,const int dim2,const double factor){
 if (::fabs(factor-1.0)<1.0e-12) { 
  transpose_1(d,s,dim1,dim2);
 } else {
  int dim1s8=(dim1/8)*8;
  int dim2s8=(dim2/8)*8;

  size_t dim1_2=dim1*2; 
  size_t dim1_3=dim1*3; 
  size_t dim1_4=dim1*4; 
  size_t dim1_5=dim1*5; 
  size_t dim1_6=dim1*6; 
  size_t dim1_7=dim1*7; 

  size_t dim2_2=dim2*2; 
  size_t dim2_3=dim2*3; 
  size_t dim2_4=dim2*4; 
  size_t dim2_5=dim2*5; 
  size_t dim2_6=dim2*6; 
  size_t dim2_7=dim2*7; 

//    d[i,j] <- s[j,i]
  for(int i=0;i<dim1s8;i+=8){
    size_t dim2_i=dim2*i;  

    for(int j=0;j<dim2s8;j+=8){
      size_t dim1_j=dim1*j;  

      d[i  +dim1_j       ]=s[j  +dim2_i       ]*factor;  
      d[i+1+dim1_j       ]=s[j  +dim2_i+dim2  ]*factor;  
      d[i+2+dim1_j       ]=s[j  +dim2_i+dim2_2]*factor;  
      d[i+3+dim1_j       ]=s[j  +dim2_i+dim2_3]*factor;  
      d[i+4+dim1_j       ]=s[j  +dim2_i+dim2_4]*factor;  
      d[i+5+dim1_j       ]=s[j  +dim2_i+dim2_5]*factor;  
      d[i+6+dim1_j       ]=s[j  +dim2_i+dim2_6]*factor;  
      d[i+7+dim1_j       ]=s[j  +dim2_i+dim2_7]*factor;  

      d[i  +dim1_j+dim1  ]=s[j+1+dim2_i       ]*factor;  
      d[i+1+dim1_j+dim1  ]=s[j+1+dim2_i+dim2  ]*factor;  
      d[i+2+dim1_j+dim1  ]=s[j+1+dim2_i+dim2_2]*factor;  
      d[i+3+dim1_j+dim1  ]=s[j+1+dim2_i+dim2_3]*factor;  
      d[i+4+dim1_j+dim1  ]=s[j+1+dim2_i+dim2_4]*factor;  
      d[i+5+dim1_j+dim1  ]=s[j+1+dim2_i+dim2_5]*factor;  
      d[i+6+dim1_j+dim1  ]=s[j+1+dim2_i+dim2_6]*factor;  
      d[i+7+dim1_j+dim1  ]=s[j+1+dim2_i+dim2_7]*factor;  

      d[i  +dim1_j+dim1_2]=s[j+2+dim2_i       ]*factor;  
      d[i+1+dim1_j+dim1_2]=s[j+2+dim2_i+dim2  ]*factor;  
      d[i+2+dim1_j+dim1_2]=s[j+2+dim2_i+dim2_2]*factor;  
      d[i+3+dim1_j+dim1_2]=s[j+2+dim2_i+dim2_3]*factor;  
      d[i+4+dim1_j+dim1_2]=s[j+2+dim2_i+dim2_4]*factor;  
      d[i+5+dim1_j+dim1_2]=s[j+2+dim2_i+dim2_5]*factor;  
      d[i+6+dim1_j+dim1_2]=s[j+2+dim2_i+dim2_6]*factor;  
      d[i+7+dim1_j+dim1_2]=s[j+2+dim2_i+dim2_7]*factor;  

      d[i  +dim1_j+dim1_3]=s[j+3+dim2_i       ]*factor;  
      d[i+1+dim1_j+dim1_3]=s[j+3+dim2_i+dim2  ]*factor;  
      d[i+2+dim1_j+dim1_3]=s[j+3+dim2_i+dim2_2]*factor;  
      d[i+3+dim1_j+dim1_3]=s[j+3+dim2_i+dim2_3]*factor;  
      d[i+4+dim1_j+dim1_3]=s[j+3+dim2_i+dim2_4]*factor;  
      d[i+5+dim1_j+dim1_3]=s[j+3+dim2_i+dim2_5]*factor;  
      d[i+6+dim1_j+dim1_3]=s[j+3+dim2_i+dim2_6]*factor;  
      d[i+7+dim1_j+dim1_3]=s[j+3+dim2_i+dim2_7]*factor;  

      d[i  +dim1_j+dim1_4]=s[j+4+dim2_i       ]*factor;  
      d[i+1+dim1_j+dim1_4]=s[j+4+dim2_i+dim2  ]*factor;  
      d[i+2+dim1_j+dim1_4]=s[j+4+dim2_i+dim2_2]*factor;  
      d[i+3+dim1_j+dim1_4]=s[j+4+dim2_i+dim2_3]*factor;  
      d[i+4+dim1_j+dim1_4]=s[j+4+dim2_i+dim2_4]*factor;  
      d[i+5+dim1_j+dim1_4]=s[j+4+dim2_i+dim2_5]*factor;  
      d[i+6+dim1_j+dim1_4]=s[j+4+dim2_i+dim2_6]*factor;  
      d[i+7+dim1_j+dim1_4]=s[j+4+dim2_i+dim2_7]*factor;  

      d[i  +dim1_j+dim1_5]=s[j+5+dim2_i       ]*factor;  
      d[i+1+dim1_j+dim1_5]=s[j+5+dim2_i+dim2  ]*factor;  
      d[i+2+dim1_j+dim1_5]=s[j+5+dim2_i+dim2_2]*factor;  
      d[i+3+dim1_j+dim1_5]=s[j+5+dim2_i+dim2_3]*factor;  
      d[i+4+dim1_j+dim1_5]=s[j+5+dim2_i+dim2_4]*factor;  
      d[i+5+dim1_j+dim1_5]=s[j+5+dim2_i+dim2_5]*factor;  
      d[i+6+dim1_j+dim1_5]=s[j+5+dim2_i+dim2_6]*factor;  
      d[i+7+dim1_j+dim1_5]=s[j+5+dim2_i+dim2_7]*factor;  

      d[i  +dim1_j+dim1_6]=s[j+6+dim2_i       ]*factor;  
      d[i+1+dim1_j+dim1_6]=s[j+6+dim2_i+dim2  ]*factor;  
      d[i+2+dim1_j+dim1_6]=s[j+6+dim2_i+dim2_2]*factor;  
      d[i+3+dim1_j+dim1_6]=s[j+6+dim2_i+dim2_3]*factor;  
      d[i+4+dim1_j+dim1_6]=s[j+6+dim2_i+dim2_4]*factor;  
      d[i+5+dim1_j+dim1_6]=s[j+6+dim2_i+dim2_5]*factor;  
      d[i+6+dim1_j+dim1_6]=s[j+6+dim2_i+dim2_6]*factor;  
      d[i+7+dim1_j+dim1_6]=s[j+6+dim2_i+dim2_7]*factor;  

      d[i  +dim1_j+dim1_7]=s[j+7+dim2_i       ]*factor;  
      d[i+1+dim1_j+dim1_7]=s[j+7+dim2_i+dim2  ]*factor;  
      d[i+2+dim1_j+dim1_7]=s[j+7+dim2_i+dim2_2]*factor;  
      d[i+3+dim1_j+dim1_7]=s[j+7+dim2_i+dim2_3]*factor;  
      d[i+4+dim1_j+dim1_7]=s[j+7+dim2_i+dim2_4]*factor;  
      d[i+5+dim1_j+dim1_7]=s[j+7+dim2_i+dim2_5]*factor;  
      d[i+6+dim1_j+dim1_7]=s[j+7+dim2_i+dim2_6]*factor;  
      d[i+7+dim1_j+dim1_7]=s[j+7+dim2_i+dim2_7]*factor;  
    } 
    for(int j=dim2s8;j<dim2;++j) {
      size_t dim1_j=dim2*j;  

      d[i  +dim1_j]=s[j+dim2_i       ]*factor;
      d[i+1+dim1_j]=s[j+dim2_i+dim2  ]*factor;
      d[i+2+dim1_j]=s[j+dim2_i+dim2_2]*factor;
      d[i+3+dim1_j]=s[j+dim2_i+dim2_3]*factor;
      d[i+4+dim1_j]=s[j+dim2_i+dim2_4]*factor;
      d[i+5+dim1_j]=s[j+dim2_i+dim2_5]*factor;
      d[i+6+dim1_j]=s[j+dim2_i+dim2_6]*factor;
      d[i+7+dim1_j]=s[j+dim2_i+dim2_7]*factor;
    }
  }

  for(int i=dim1s8;i<dim1;++i){
    size_t dim2_i=dim2*i;  

    for(int j=0;j<dim2s8;j+=8){
      size_t dim1_j=dim1*j;  

      d[i+dim1_j       ]=s[j  +dim2_i]*factor;
      d[i+dim1_j+dim1  ]=s[j+1+dim2_i]*factor;
      d[i+dim1_j+dim1_2]=s[j+2+dim2_i]*factor;
      d[i+dim1_j+dim1_3]=s[j+3+dim2_i]*factor;
      d[i+dim1_j+dim1_4]=s[j+4+dim2_i]*factor;
      d[i+dim1_j+dim1_5]=s[j+5+dim2_i]*factor;
      d[i+dim1_j+dim1_6]=s[j+6+dim2_i]*factor;
      d[i+dim1_j+dim1_7]=s[j+7+dim2_i]*factor;
    }
    for(int j=dim2s8;j<dim2;++j) {
      //std::cout << "*";
      d[i+dim1*j]=s[j+dim2_i]*factor;
    }
  }
 } 
}


