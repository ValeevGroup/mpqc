//
// classda.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
// Maintainer: LPS
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

#ifndef CLASSNAME
#define CLASSNAME you_forgot_to_define_CLASSNAME
#endif

  private:
    static ClassDesc class_desc_;
    void* do_castdowns(void**,
                       const ClassDesc*cd);
  public:
    void* _castdown(const ClassDesc*);
    static CLASSNAME* require_castdown(DescribedClass*p,const char*,...);
    static const CLASSNAME* require_const_castdown(const DescribedClass*p,const char*,...);
    static CLASSNAME* castdown(DescribedClass*p);
    static const CLASSNAME* const_castdown(const DescribedClass*p);
    static CLASSNAME* castdown(const RefDescribedClass&p);
    static const ClassDesc* static_class_desc();
    const ClassDesc* class_desc() const;
#ifdef HAVE_CTOR
#undef HAVE_CTOR
#endif
#ifdef HAVE_KEYVAL_CTOR
#undef HAVE_KEYVAL_CTOR
#endif
#ifdef HAVE_STATEIN_CTOR
#undef HAVE_STATEIN_CTOR
#endif
  private:

#undef HAVE_KEYVAL_CTOR
#undef HAVE_STATEIN_CTOR
#undef HAVE_CTOR
#undef CLASSNAME
