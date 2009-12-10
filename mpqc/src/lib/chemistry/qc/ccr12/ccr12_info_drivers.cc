//
// ccr12_info_drivers.cc
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@qtp.ufl.edu>
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

#include <algorithm>
#include <chemistry/qc/ccr12/ccr12_info.h>
#include <chemistry/qc/ccr12/tensor.h>

using namespace sc;

void CCR12_Info::jacobi_t1_(const Ref<Tensor>& d_r1_,Ref<Tensor>& d_t1_){
  for(long p1b=noab();p1b<noab()+nvab();++p1b){
   for(long h2b=0L;h2b<noab();++h2b){
    if(d_t1_->is_this_local(h2b+noab()*(p1b-noab()))){ 
     long size=get_range(p1b)*get_range(h2b);
     double* data=mem()->malloc_local_double(size); 
     d_r1_->get_block(h2b+noab()*(p1b-noab()),data); 
     long i=0;
     for(long p1=0;p1<get_range(p1b);++p1){ 
      double ep1=orbital_evl_sorted_[get_offset(p1b)+p1];
      for(long h2=0;h2<get_range(h2b);++h2,++i){ 
       double eh2=orbital_evl_sorted_[get_offset(h2b)+h2];
       data[i]/=(eh2-ep1); 
      }
     }
     d_t1_->add_block(h2b+noab()*(p1b-noab()),data); 
     mem()->free_local_double(data);
    }
   }
  }
}

void CCR12_Info::jacobi_t2_(const Ref<Tensor>& d_r2_,Ref<Tensor>& d_t2_){
  for(long p1b=noab();p1b<noab()+nvab();++p1b){
   for(long p2b=p1b;p2b<noab()+nvab();++p2b){
    for(long h3b=0L;h3b<noab();++h3b){
     for(long h4b=h3b;h4b<noab();++h4b){
      if(d_t2_->is_this_local(h4b+noab()*(h3b+noab()*(p2b-noab()+nvab()*(p1b-noab()))))){ 
       long size=get_range(p1b)*get_range(p2b)*get_range(h3b)*get_range(h4b);
       double* data=mem()->malloc_local_double(size); 
       d_r2_->get_block(h4b+noab()*(h3b+noab()*(p2b-noab()+nvab()*(p1b-noab()))),data); 
       long i=0;
       for(long p1=0;p1<get_range(p1b);++p1){ 
        double ep1=orbital_evl_sorted_[get_offset(p1b)+p1];
        for(long p2=0;p2<get_range(p2b);++p2){ 
         double ep2=orbital_evl_sorted_[get_offset(p2b)+p2];
         for(long h3=0;h3<get_range(h3b);++h3){ 
          double eh3=orbital_evl_sorted_[get_offset(h3b)+h3];
          for(long h4=0;h4<get_range(h4b);++h4,++i){ 
           double eh4=orbital_evl_sorted_[get_offset(h4b)+h4];
           data[i]/=(eh3+eh4-ep1-ep2); 
          }
         }
        }
       }
       d_t2_->add_block(h4b+noab()*(h3b+noab()*(p2b-noab()+nvab()*(p1b-noab()))),data); 
       mem()->free_local_double(data);
      }
     }
    }
   }
  }
}


void CCR12_Info::jacobi_t3_(const Ref<Tensor>& d_r3_,Ref<Tensor>& d_t3_){
  for(long p1b=noab();p1b<noab()+nvab();++p1b){
   for(long p2b=p1b;p2b<noab()+nvab();++p2b){
    for(long p3b=p2b;p3b<noab()+nvab();++p3b){
     for(long h4b=0L;h4b<noab();++h4b){
      for(long h5b=h4b;h5b<noab();++h5b){
       for(long h6b=h5b;h6b<noab();++h6b){
        long tag=h6b+noab()*(h5b+noab()*(h4b+noab()*(p3b-noab()+nvab()*(p2b-noab()+nvab()*(p1b-noab())))));
        if(d_t3_->is_this_local(tag)){ 
        long size=get_range(p1b)*get_range(p2b)*get_range(p3b)*get_range(h4b)*get_range(h5b)*get_range(h6b);
        double* data=mem()->malloc_local_double(size); 
        d_r3_->get_block(tag,data); 
        long i=0;
        for(long p1=0;p1<get_range(p1b);++p1){ 
         double ep1=orbital_evl_sorted_[get_offset(p1b)+p1];
         for(long p2=0;p2<get_range(p2b);++p2){ 
          double ep2=orbital_evl_sorted_[get_offset(p2b)+p2];
          for(long p3=0;p3<get_range(p3b);++p3){ 
           double ep3=orbital_evl_sorted_[get_offset(p3b)+p3];
           for(long h4=0;h4<get_range(h4b);++h4){ 
            double eh4=orbital_evl_sorted_[get_offset(h4b)+h4];
            for(long h5=0;h5<get_range(h5b);++h5){ 
             double eh5=orbital_evl_sorted_[get_offset(h5b)+h5];
             for(long h6=0;h6<get_range(h6b);++h6,++i){ 
              double eh6=orbital_evl_sorted_[get_offset(h6b)+h6];
              data[i]/=(eh4+eh5+eh6-ep1-ep2-ep3); 
             }
            }
           }
          }
         }
        }
        d_t3_->add_block(tag,data); 
        mem()->free_local_double(data);
       }
      }
     }
    }
   }
  }
 }
}


void CCR12_Info::jacobi_t4_(const Ref<Tensor>& d_r4_,Ref<Tensor>& d_t4_){
  for(long p1b=noab();p1b<noab()+nvab();++p1b){
   for(long p2b=p1b;p2b<noab()+nvab();++p2b){
    for(long p3b=p2b;p3b<noab()+nvab();++p3b){
     for(long p4b=p3b;p4b<noab()+nvab();++p4b){
      for(long h5b=0L;h5b<noab();++h5b){
       for(long h6b=h5b;h6b<noab();++h6b){
        for(long h7b=h6b;h7b<noab();++h7b){
         for(long h8b=h7b;h8b<noab();++h8b){
          long tag=h8b+noab()*(h7b+noab()*(h6b+noab()*(h5b+noab()*(p4b-noab()+nvab()*(p3b-noab()+nvab()*(p2b-noab()+nvab()*(p1b-noab())))))));
          if(d_t4_->is_this_local(tag)){ 
          long size=get_range(p1b)*get_range(p2b)*get_range(p3b)*get_range(p4b)*get_range(h5b)*get_range(h6b)*get_range(h7b)*get_range(h8b);
          double* data=mem()->malloc_local_double(size); 
          d_r4_->get_block(tag,data); 
          long i=0;
          for(long p1=0;p1<get_range(p1b);++p1){ 
           double ep1=orbital_evl_sorted_[get_offset(p1b)+p1];
           for(long p2=0;p2<get_range(p2b);++p2){ 
            double ep2=orbital_evl_sorted_[get_offset(p2b)+p2];
            for(long p3=0;p3<get_range(p3b);++p3){ 
             double ep3=orbital_evl_sorted_[get_offset(p3b)+p3];
             for(long p4=0;p4<get_range(p4b);++p4){ 
              double ep4=orbital_evl_sorted_[get_offset(p4b)+p4];
              for(long h5=0;h5<get_range(h5b);++h5){ 
               double eh5=orbital_evl_sorted_[get_offset(h5b)+h5];
               for(long h6=0;h6<get_range(h6b);++h6){ 
                double eh6=orbital_evl_sorted_[get_offset(h6b)+h6];
                for(long h7=0;h7<get_range(h7b);++h7){ 
                 double eh7=orbital_evl_sorted_[get_offset(h7b)+h7];
                 for(long h8=0;h8<get_range(h8b);++h8,++i){ 
                  double eh8=orbital_evl_sorted_[get_offset(h8b)+h8];
                  data[i]/=(eh5+eh6+eh7+eh8-ep1-ep2-ep3-ep4); 
                 }
                }
               }
              }
             }
            }
           }
          }
          d_t4_->add_block(tag,data); 
          mem()->free_local_double(data);
         }
        }
       }
      }
     }
    }
   }
  }
 }
}


void CCR12_Info::jacobi_l1_(const Ref<Tensor>& d_lr1_,Ref<Tensor>& d_l1_){
  for(long h2b=0L;h2b<noab();++h2b){
   for(long p1b=noab();p1b<noab()+nvab();++p1b){
    if(d_l1_->is_this_local(p1b-noab()+nvab()*h2b)){ 
     long size=get_range(p1b)*get_range(h2b);
     double* data=mem()->malloc_local_double(size); 
     d_lr1_->get_block(p1b-noab()+nvab()*h2b,data); 
     long i=0;
     for(long h2=0;h2<get_range(h2b);++h2){ 
      double eh2=orbital_evl_sorted_[get_offset(h2b)+h2];
      for(long p1=0;p1<get_range(p1b);++p1,++i){ 
       double ep1=orbital_evl_sorted_[get_offset(p1b)+p1];
       data[i]/=(eh2-ep1); 
      }
     }
     d_l1_->add_block(p1b-noab()+nvab()*h2b,data); 
     mem()->free_local_double(data);
    }
   }
  }
}

void CCR12_Info::jacobi_l2_(const Ref<Tensor>& d_lr2_,Ref<Tensor>& d_l2_){
  for(long h3b=0L;h3b<noab();++h3b){
   for(long h4b=h3b;h4b<noab();++h4b){
    for(long p1b=noab();p1b<noab()+nvab();++p1b){
     for(long p2b=p1b;p2b<noab()+nvab();++p2b){
      if(d_l2_->is_this_local(p2b-noab()+nvab()*(p1b-noab()+nvab()*(h4b+noab()*h3b)))){ 
       long size=get_range(p1b)*get_range(p2b)*get_range(h3b)*get_range(h4b);
       double* data=mem()->malloc_local_double(size); 
       d_lr2_->get_block(p2b-noab()+nvab()*(p1b-noab()+nvab()*(h4b+noab()*h3b)),data); 
       long i=0;
       for(long h3=0;h3<get_range(h3b);++h3){ 
        double eh3=orbital_evl_sorted_[get_offset(h3b)+h3];
        for(long h4=0;h4<get_range(h4b);++h4){ 
         double eh4=orbital_evl_sorted_[get_offset(h4b)+h4];
         for(long p1=0;p1<get_range(p1b);++p1){ 
          double ep1=orbital_evl_sorted_[get_offset(p1b)+p1];
          for(long p2=0;p2<get_range(p2b);++p2,++i){ 
           double ep2=orbital_evl_sorted_[get_offset(p2b)+p2];
           data[i]/=(eh3+eh4-ep1-ep2); 
          }
         }
        }
       }
       d_l2_->add_block(p2b-noab()+nvab()*(p1b-noab()+nvab()*(h4b+noab()*h3b)),data); 
       mem()->free_local_double(data);
      }
     }
    }
   }
  }
}

// read appropriate blocks from d_v2 and form the mp1 wavefunction
void CCR12_Info::guess_t2(Ref<Tensor>& d_t2_){
  for(long p1b=noab();p1b<noab()+nvab();++p1b){
   for(long p2b=p1b;p2b<noab()+nvab();++p2b){
    for(long h3b=0L;h3b<noab();++h3b){
     for(long h4b=h3b;h4b<noab();++h4b){
      if(d_t2_->is_this_local(h4b+noab()*(h3b+noab()*(p2b-noab()+nvab()*(p1b-noab()))))){ 
       long size=get_range(p1b)*get_range(p2b)*get_range(h3b)*get_range(h4b);
       double* data=mem()->malloc_local_double(size); 
       d_v2->get_block(h4b+nab()*(h3b+nab()*(p2b+nab()*p1b)),data); 
       long i=0;
       for(long p1=0;p1<get_range(p1b);++p1){ 
        double ep1=orbital_evl_sorted_[get_offset(p1b)+p1];
        for(long p2=0;p2<get_range(p2b);++p2){ 
         double ep2=orbital_evl_sorted_[get_offset(p2b)+p2];
         for(long h3=0;h3<get_range(h3b);++h3){ 
          double eh3=orbital_evl_sorted_[get_offset(h3b)+h3];
          for(long h4=0;h4<get_range(h4b);++h4,++i){ 
           double eh4=orbital_evl_sorted_[get_offset(h4b)+h4];
           data[i]/=(eh3+eh4-ep1-ep2); 
          }
         }
        }
       }
       d_t2_->add_block(h4b+noab()*(h3b+noab()*(p2b-noab()+nvab()*(p1b-noab()))),data); 
       mem()->free_local_double(data);
      }
     }
    }
   }
  }
}

// read appropriate blocks from d_v2, d_gt2 and form the mp1-r12 wavefunction
void CCR12_Info::guess_t2_r12(Ref<Tensor>& d_t2_, Ref<Tensor>& d_gt2_){
//if (r12world()->r12tech()->ansatz()->amplitudes() != LinearR12::GeminalAmplitudeAnsatz_fullopt) {
  // d_t2 += v / Denominator (MP2 amplitude)
  guess_t2(d_t2_);

  // ad = A^dagger
  Ref<Tensor> ad_local = new Tensor("ad_local", mem_);
  offset_l2(ad_local);
  form_ad(ad_local);

  // ca = A * gt2
  Ref<Tensor> ca_local = new Tensor("ca_local", mem_);
  offset_t2(ca_local, false);
  form_ca(d_gt2, ad_local, ca_local);

  // d_t2_ += A * gt2 / Denominator
  jacobi_t2_(ca_local, d_t2_);

}

// transpose t amplitudes to make a guess of lambda amplitudes
void CCR12_Info::guess_lambda1(Ref<Tensor>& d_lambda1_){
  for(long h2b=0L;h2b<noab();++h2b){
   for(long p1b=noab();p1b<noab()+nvab();++p1b){
    if(d_lambda1_->is_this_local(p1b-noab()+nvab()*h2b)){ 
     long size=get_range(p1b)*get_range(h2b);
     double* data1=mem()->malloc_local_double(size); 
     double* data2=mem()->malloc_local_double(size); 

     d_t1->get_block(h2b+noab()*(p1b-noab()),data1); 
     sort_indices2(data1,data2,get_range(p1b),get_range(h2b),1,0,1.0);
     d_lambda1_->put_block(p1b-noab()+nvab()*h2b,data2); 

     mem()->free_local_double(data2);
     mem()->free_local_double(data1);
    }
   }
  }
}

void CCR12_Info::guess_lambda2(Ref<Tensor>& d_lambda2_){
  for(long h3b=0L;h3b<noab();++h3b){
   for(long h4b=h3b;h4b<noab();++h4b){
    for(long p1b=noab();p1b<noab()+nvab();++p1b){
     for(long p2b=p1b;p2b<noab()+nvab();++p2b){
      if(d_lambda2_->is_this_local(p2b-noab()+nvab()*(p1b-noab()+nvab()*(h4b+noab()*h3b)))){ 
       long size=get_range(p1b)*get_range(p2b)*get_range(h3b)*get_range(h4b);
       double* data1=mem()->malloc_local_double(size); 
       double* data2=mem()->malloc_local_double(size); 

       d_t2->get_block(h4b+noab()*(h3b+noab()*(p2b-noab()+nvab()*(p1b-noab()))),data1); 
       sort_indices4(data1,data2,get_range(p1b),get_range(p2b),get_range(h3b),get_range(h4b),2,3,0,1,1.0);
       d_lambda2_->put_block(p2b-noab()+nvab()*(p1b-noab()+nvab()*(h4b+noab()*h3b)),data2); 

       mem()->free_local_double(data2);
       mem()->free_local_double(data1);
      }
     }
    }
   }
  }
}


void CCR12_Info::guess_glambda2(Ref<Tensor>& d_glambda2_){
  for(long h3b=0L;h3b<noab();++h3b){
   for(long h4b=h3b;h4b<noab();++h4b){
    for(long h1b=0L;h1b<noab();++h1b){
     for(long h2b=h1b;h2b<noab();++h2b){
      if(d_glambda2_->is_this_local(h2b+noab()*(h1b+noab()*(h4b+noab()*h3b)))){ 
       long size=get_range(h1b)*get_range(h2b)*get_range(h3b)*get_range(h4b);
       double* data1=mem()->malloc_local_double(size); 
       double* data2=mem()->malloc_local_double(size); 

       d_gt2->get_block(h4b+noab()*(h3b+noab()*(h2b+noab()*h1b)),data1); 
       sort_indices4(data1,data2,get_range(h1b),get_range(h2b),get_range(h3b),get_range(h4b),2,3,0,1,1.0);
       d_glambda2_->put_block(h2b+noab()*(h1b+noab()*(h4b+noab()*h3b)),data2); 

       mem()->free_local_double(data2);
       mem()->free_local_double(data1);
      }
     }
    }
   }
  }
}

