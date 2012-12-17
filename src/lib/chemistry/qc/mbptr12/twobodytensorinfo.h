//
// twobodytensorinfo.h
//
// Copyright (C) 2009 Martin Torheyden
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


#ifndef TWOBODYTENSORINFO_H_
#define TWOBODYTENSORINFO_H_

#include <util/ref/ref.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/basis/tbint.h>

namespace sc {
  /** Provides information about the type of a two body tensor. */
  class TwoBodyTensorInfo : public RefCount {
    public:
      enum tbtensor_type { tbint = 0, geminalcoeff = 1 };
    private:
      tbtensor_type twobodytensor_type_;
      TwoBodyOper::type twobodyint_type_;
    public:
      TwoBodyTensorInfo(tbtensor_type twobodytensor_type){
        if(twobodytensor_type==tbint) {
          throw ProgrammingError("Error in TwoBodyTensorInfo::TwoBodyTensorInfo -- if twobodytensor_type==tbint TwoBodyOper::type must be specified.",__FILE__,__LINE__);
        }
        twobodytensor_type_ = twobodytensor_type;
      }
      TwoBodyTensorInfo(TwoBodyOper::type twobodyint_type) { twobodytensor_type_ = tbint; twobodyint_type_ = twobodyint_type; }
      virtual ~TwoBodyTensorInfo(){}
      tbtensor_type twobodytensor_type() const { return(twobodytensor_type_); }
      TwoBodyOper::type twobodyint_type() const {
        if(twobodytensor_type_!=tbint)
          throw ProgrammingError("Error in TwoBodyTensorInfo::twobodyint_type() -- twobodyint_type_ can only be returned if twobodytensor_type_==tbint.",__FILE__,__LINE__);
        return(twobodyint_type_);
      }
  };

}

#endif /*TWOBODYTENSORINFO_H_*/
