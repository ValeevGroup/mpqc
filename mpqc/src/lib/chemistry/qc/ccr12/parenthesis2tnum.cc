//
// parenthesis2tnum.cc
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

#include <util/misc/exenv.h>
#include <util/state/stateio.h>
#include <math/scmat/blocked.h>
#include <util/class/scexception.h>
#include <chemistry/qc/ccr12/parenthesis2tnum.h>

using namespace sc;

static ClassDesc Parenthesis2tNum_cd(
  typeid(Parenthesis2tNum),"Parenthesis2tNum",1,"virtual public RefCount"
  ,0,0,0);

Parenthesis2tNum::Parenthesis2tNum(CCR12_Info* info): z(info){
}
    
Parenthesis2tNum::~Parenthesis2tNum(){
}

void Parenthesis2tNum::compute_amp(double*,const long,const long,const long,
                                   const long,const long,const long,const long){
} 

void Parenthesis2tNum::compute_amp(double*,const long,const long,const long,const long,
                                   const long,const long,const long,const long,const long){
} 

