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

namespace detail {
template <int NumCenters, typename T>
struct tuple {
    T data[NumCenters];
    T operator[](size_t i) const {
      return data[i];
    }
    T& operator[](size_t i) {
      return data[i];
    }
    T operator= (const tuple<NumCenters, T> t )  {
        for(int i = 0; i < NumCenters; i++)
          data[i] = t[i];
    }
};
}

/// This is an abstract range of indices
template <int NDIM>
class IndexRangeIterator : public RefCount {
  public:
    typedef detail::tuple<NDIM, unsigned int> IntTuple;

    /// initialize the iterator
    virtual void init() =0;
    /// returns false if there are no more
    virtual bool in_range() const =0;
    /// update current
    virtual void next() =0;
    /// current shell set
    virtual const IntTuple& current() const =0;
};

/// TensorIndexRangeIterator is a direct product of shell ranges for each center
template <int NDIM>
class TensorIndexRangeIterator : public IndexRangeIterator<NDIM> {
  public:
    typedef typename IndexRangeIterator<NDIM>::IntTuple IntTuple;

    /// Constructs a [start, fence) IndexRangeIterator
    TensorIndexRangeIterator(const IntTuple& start, const IntTuple& fence) :
      start_(start), fence_(fence), in_range_(false) {
    }

    void init() {
      current_ = start_;
      next_ = current_;
      // if range is nonempty -- find next
      if (is_nonempty()) {
        in_range_ = true;
        next();
      }
    }

    bool in_range() const {
      return in_range_;
    }

    void next() {
      in_range_ = have_next_;
      current_ = next_;
      have_next_ = inc(next_);
    }

    const IntTuple& current() const {
      return current_;
    }

  private:
    // return true if can increment
    bool inc(IntTuple & s) {
      unsigned int d = NDIM - 1;
      ++s[d];
      while (s[d] >= fence_[d]) {
        if (d == 0)
          return false;
        s[d] = start_[d];
       --d;
        ++s[d];
      }
    }
    // return true if range is nonempty
    bool is_nonempty() {
      for (unsigned int d = 0; d < NDIM; ++d)
        if (fence_[d] <= start_[d])
          return false;
      return true;
    }

    IntTuple start_;
    IntTuple fence_;
    IntTuple current_;
    IntTuple next_;
    bool in_range_;
    bool have_next_;
};

/** This is an abstract base type for classes that
 compute integrals involving two electrons and 2 functions per electron.
 */
template <unsigned int NumCenters>
class TwoBodyIntBatch:public RefCount {


  public:
    typedef detail::tuple<NumCenters, unsigned int> IntTuple;

    //TwoBodyIntBatch(Integral *integral,
    //        const detail::tuple<NumCenters, Ref<GaussianBasisSet> >& bs); // breaks formatting in Eclipse
    //TwoBodyIntBatch(Integral *i, const detail::tuple<NumCenters, Ref<GaussianBasisSet> >& b);
    TwoBodyIntBatch(Ref<Integral> i) :
      integral_(i) {
      if (NumCenters > 0)
        bs_[0] = i->basis1();
      if (NumCenters > 1)
        bs_[1] = i->basis2();
      if (NumCenters > 2)
        bs_[2] = i->basis3();
      if (NumCenters > 3)
        bs_[3] - i->basis4();
    }


    virtual ~TwoBodyIntBatch(){}

    /// prepare to iterate using seed s
    /// TODO JTF implement
    template <typename Seed> void init(const Ref<IndexRangeIterator<NumCenters> >& n, Seed s = Seed());

    /// compute next batch, return true if have another
    /// may need to be split into have_next and next
    /// TODO JTF implement
    virtual bool next() = 0;

    /// returns the shell indices of the current batch
    virtual const std::vector<IntTuple>& current_batch() const = 0;

    /** The computed shell integrals will be put in the buffer returned
     by this member.  Some TwoBodyInt specializations have more than
     one buffer:  The type arguments selects which buffer is returned.
     If the requested type is not supported, then 0 is returned. */
    /// TODO JTF implement
    virtual const double * buffer(TwoBodyOper::type type = TwoBodyOper::eri) const = 0;

    /// Return the basis set on center c
    /// TODO JTF implement
    const Ref<GaussianBasisSet>& basis(unsigned int c = 0) const;

    /** Returns the type of the operator set that this object computes.
     this function is necessary to describe the computed integrals
     (their number, symmetries, etc.) and/or to implement cloning. */
    virtual TwoBodyOperSet::type type() const =0;

    /// return the operator set descriptor
    virtual const Ref<TwoBodyOperSetDescr>& descr() const =0;

    /// This storage is used to cache computed integrals.
    /// TODO: do we really want to allocate by size_t, instead of number of elements?
    virtual void set_integral_storage(size_t storage) = 0;

    /** Return true if the clone member can be called.  The default
     * implementation returns false. */
    /// TODO JTF implement
    virtual bool cloneable() = 0;

    /** Returns a clone of this.  The default implementation throws an
     * exception. */
    /// TODO JTF implement
    virtual Ref<TwoBodyIntBatch> clone() = 0;

    /// Return the integral factory that was used to create this object.
    Integral *integral() const {
      return integral_;
    }



  protected:
    // this is who created me
    Ref<Integral> integral_;
    Ref<GaussianBasisSet> bs_[NumCenters];


    std::vector<double> buffer_;
    std::vector<IntTuple> shells_in_buffer_;
    std::vector<IntTuple> start_;
    std::vector<IntTuple> fence_;
    size_t buf_cap_ ;
};

template <int NumCenters> struct TwoBodyIntEvalType;

/// This is a generic implementation of TwoBodyIntBatch in terms of a TwoBodyInt
template <int NumCenters>
class TwoBodyIntBatchGeneric:public TwoBodyIntBatch<NumCenters>  {

  typedef TensorIndexRangeIterator<NumCenters> SRange;
  typedef typename SRange::IntTuple IntTuple;

  public:
    typedef typename TwoBodyIntEvalType<NumCenters>::value TwoBodyIntEval;

    ~TwoBodyIntBatchGeneric() {
      this->buffer_.clear();
    }

    //    TwoBodyIntBatchGeneric(Integral *i, const detail::tuple<NumCenters, Ref<GaussianBasisSet> >& b) :
    TwoBodyIntBatchGeneric(Ref<TwoBodyIntEval> tbint) :
      TwoBodyIntBatch<NumCenters>(tbint->integral()), tbint_(tbint) {

      in_buf_ = tbint_->buffer();

    }

    bool cloneable(){
      return false;
    }

    Ref<TwoBodyIntBatch<NumCenters> > clone() {
      return new TwoBodyIntBatchGeneric<NumCenters>(this->tbint_->clone());
    }

    TwoBodyOperSet::type type() const {
      return tbint_->type();
    }
    const Ref<TwoBodyOperSetDescr>& descr() const {
      return tbint_->descr();
    }

    void set_integral_storage(size_t storage) {
      this->buf_cap_ = storage;
      this->buffer_.reserve(storage);
      this->shells_in_buffer_.reserve(storage); // this seems dumb; but no other way to guarantee enough room (all <ss|ss>?)
    }

    const std::vector<IntTuple>& current_batch() const {
      return this->shells_in_buffer_ ;
    }


    const std::vector<IntTuple>& start(){
      return this->start_ ;
    }

    const std::vector<IntTuple>& fence(){
      return this->fence_ ;
    }

    const double * buffer(TwoBodyOper::type type = TwoBodyOper::eri) const {
      return in_buf_;
    }

    template <unsigned int> void init(const Ref<IndexRangeIterator<NumCenters> >& i, int s = int()) {
      //
      int sp;
      sr_(i);
      sr_.start();
      if (s) {
        for (int sp = 0; sp < s, sr_.have_next(); sp++) {
          sr_.next();
        }
      }
    }

    bool next() {
      return 0;
    }


  private:

    const Ref<TwoBodyIntEval> tbint_;
    const Ref<SRange > sr_;
    const double *in_buf_;
    const Ref<GaussianBasisSet> basis_;

};




} // end of namespace sc

#endif // end of header guard
// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
