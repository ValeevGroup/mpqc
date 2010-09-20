//
// tbint_batch.h
//
// Copyright (C) 2010 Edward Valeev
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

#ifndef _mpqc_src_lib_chemistry_qc_basis_tbint_batch_h
#define _mpqc_src_lib_chemistry_qc_basis_tbint_batch_h

#include <util/ref/ref.h>
#include <chemistry/qc/basis/gaussbas.h>

namespace sc {

  template <unsigned int NumCenters> class ShellRange;

  /** This is an abstract base type for classes that
      compute integrals involving two electrons and 2 functions per electron.
   */
  template <unsigned int NumCenters>
    class TwoBodyIntBatch : public RefCount {
    public:
      virtual ~TwoBodyIntBatch();

      /// prepare to iterate using seed s
      template <typename Seed> void init(const ShellRange<NumCenters>& i, Seed s = Seed());
      /// compute next batch, return true if have another
      bool next();

      template <typename T>
      struct tuple {
          T data[NumCenters];
          T operator[](size_t i) const { return data[i]; }
          T& operator[](size_t i) { return data[i]; }
      };

      /// returns the shell indices of the current batch
      const std::vector< tuple<unsigned int> >& current_batch() const;

      /** The computed shell integrals will be put in the buffer returned
          by this member.  Some TwoBodyInt specializations have more than
          one buffer:  The type arguments selects which buffer is returned.
          If the requested type is not supported, then 0 is returned. */
      virtual const double * buffer(TwoBodyOper::type type = TwoBodyOper::eri) const;

      /// Return the basis set on center c
      const Ref<GaussianBasisSet>& basis(unsigned int c = 0) const;

      /** Returns the type of the operator set that this object computes.
          this function is necessary to describe the computed integrals
          (their number, symmetries, etc.) and/or to implement cloning. */
      virtual TwoBodyOperSet::type type() const =0;
      /// return the operator set descriptor
      virtual const Ref<TwoBodyOperSetDescr>& descr() const =0;

      /// This storage is used to cache computed integrals.
      virtual void set_integral_storage(size_t storage);

      /** Return true if the clone member can be called.  The default
       * implementation returns false. */
      virtual bool cloneable();

      /** Returns a clone of this.  The default implementation throws an
       * exception. */
      virtual Ref<TwoBodyIntBatch> clone();

      /// Return the integral factory that was used to create this object.
      Integral *integral() const { return integral_; }

    private:
      double *log2_to_double_;

    protected:
      // this is who created me
      Integral *integral_;

      Ref<GaussianBasisSet> bs_[NumCenters];

      double *buffer_;

      TwoBodyIntBatch(Integral *integral,
                      const tuple<Ref<GaussianBasisSet> >& bs);

  };

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
