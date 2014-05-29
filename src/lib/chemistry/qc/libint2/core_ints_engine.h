//
// core_ints_engine.h
//
// Copyright (C) 2014 Edward Valeev
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
#pragma interface
#endif

#ifndef _mpqc_src_lib_chemistry_qc_libint2_core_ints_engine_h
#define _mpqc_src_lib_chemistry_qc_libint2_core_ints_engine_h

#include <util/ref/ref.h>
#include <util/group/thread.h>

namespace sc {

  /// CoreIntsEngine manages Boys, and other core integral, engines.
  /// E.g. consider an engine that computes the Boys function, \f$ F_m(T) \f$.
  /// Since multiple users of a Boys function engine may exist, and they will require different ranges of params,
  /// it may be necessary to rebuild the engine when new users will require greater precision or larger values of parameter m.
  /// This class minimizes the number of rebuilds, and does the rebuilds in a thread-safe fashion.
  /// @tparam _Engine
  template <typename _Engine>
  class CoreIntsEngine {
    public:

      struct Engine :
          virtual public RefCount,
          public _Engine
      {
          Engine(int mmax): RefCount(), _Engine(mmax) {}
          Engine(int mmax, double prec): RefCount(), _Engine(mmax, prec) {}
      };

      template <typename Int>
      static Ref<Engine> instance(Int mmax) {
        if (default_engine_->max_m() >= mmax)
          return default_engine_;

        ThreadLockHolder lh(lock_);
        Ref<Engine> new_default_engine = new Engine(mmax);
        default_engine_ = new_default_engine;
        return default_engine_;
      }

      template <typename Int, typename Real>
      static Ref<Engine> instance(Int mmax, Real prec) {
        if (not (default_engine_->max_m() < mmax || default_engine_->precision() > prec))
          return default_engine_;

        ThreadLockHolder lh(lock_);
        Ref<Engine> new_default_engine = new Engine(mmax, prec);
        default_engine_ = new_default_engine;
        return default_engine_;
      }

    private:
      static Ref<Engine> default_engine_;
      static Ref<ThreadLock> lock_; //!< to serialize access to default_engine_
  };

  template <typename _Engine>
  Ref<ThreadLock> CoreIntsEngine<_Engine>::lock_ = ThreadGrp::get_default_threadgrp()->new_lock();

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
