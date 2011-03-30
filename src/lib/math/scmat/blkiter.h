//
// blkiter.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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

#ifndef _math_scmat_blkiter_h
#define _math_scmat_blkiter_h

#include <math/scmat/block.h>

namespace sc {

class SCMatrixRectBlock;
class SCMatrixLTriBlock;
class SCMatrixDiagBlock;
class SCVectorSimpleBlock;

class SCElementOp;
class SCElementOp2;
class SCElementOp3;

/** The SCMatrixBlockIter class is used to described iterates that
    loop through the elements in a block. */
class SCMatrixBlockIter {
  public:
    SCMatrixBlockIter() {}
    virtual ~SCMatrixBlockIter();
    /// Returns the row index.
    virtual int i() = 0;
    /// Returns the column index.
    virtual int j() = 0;
    /// Set the current element to val.
    virtual void set(double val) = 0;
    /// Add val to the current element.
    virtual void accum(double val);
    /// Return the value of the current element.
    virtual double get() = 0;
    /// Return nonzero if there are more elements.
    virtual operator int() = 0;
    /// Move to the next element.
    virtual void operator++() = 0; // prefix ++
    void operator++(int) { operator++(); }
    /// Start the iteration over.
    virtual void reset() = 0;
};

class SCMatrixRectBlockIter: public SCMatrixBlockIter {
  private:
    SCMatrixRectBlock* block;
    int i_;
    int block_index;
    int j_;
  public:
    SCMatrixRectBlockIter(SCMatrixRectBlock*);
    virtual ~SCMatrixRectBlockIter();
    int i();
    int j();
    double get();
    void set(double);
    operator int();
    void operator++();
    void reset();
};

class SCMatrixRectSubBlockIter: public SCMatrixBlockIter {
  private:
    SCMatrixRectSubBlock* block;
    int i_;
    int block_index;
    int j_;
  public:
    SCMatrixRectSubBlockIter(SCMatrixRectSubBlock*);
    virtual ~SCMatrixRectSubBlockIter();
    int i();
    int j();
    double get();
    void set(double);
    operator int();
    void operator++();
    void reset();
};

class SCMatrixLTriBlockIter: public SCMatrixBlockIter {
  private:
    SCMatrixLTriBlock* block;
    int block_index;
    int i_;
    int j_;
  public:
    SCMatrixLTriBlockIter(SCMatrixLTriBlock*);
    virtual ~SCMatrixLTriBlockIter();
    int i();
    int j();
    double get();
    void set(double);
    operator int();
    void operator++();
    void reset();
};

class SCMatrixLTriSubBlockIter: public SCMatrixBlockIter {
  private:
    SCMatrixLTriSubBlock* block;
    int block_index;
    int i_;
    int j_;
  public:
    SCMatrixLTriSubBlockIter(SCMatrixLTriSubBlock*);
    virtual ~SCMatrixLTriSubBlockIter();
    int i();
    int j();
    double get();
    void set(double);
    operator int();
    void operator++();
    void reset();
};

class SCMatrixDiagBlockIter: public SCMatrixBlockIter {
  private:
    SCMatrixDiagBlock* block;
    int block_index;
    int i_;
  public:
    SCMatrixDiagBlockIter(SCMatrixDiagBlock*);
    virtual ~SCMatrixDiagBlockIter();
    int i();
    int j();
    double get();
    void set(double);
    operator int();
    void operator++();
    void reset();
};

class SCMatrixDiagSubBlockIter: public SCMatrixBlockIter {
  private:
    SCMatrixDiagSubBlock* block;
    int block_index;
    int i_;
  public:
    SCMatrixDiagSubBlockIter(SCMatrixDiagSubBlock*);
    virtual ~SCMatrixDiagSubBlockIter();
    int i();
    int j();
    double get();
    void set(double);
    operator int();
    void operator++();
    void reset();
};

class SCVectorSimpleBlockIter: public SCMatrixBlockIter {
  private:
    SCVectorSimpleBlock* block;
    int block_index;
    int i_;
  public:
    SCVectorSimpleBlockIter(SCVectorSimpleBlock*);
    virtual ~SCVectorSimpleBlockIter();
    int i();
    int j();
    double get();
    void set(double);
    operator int();
    void operator++();
    void reset();
};

class SCVectorSimpleSubBlockIter: public SCMatrixBlockIter {
  private:
    SCVectorSimpleSubBlock* block;
    int block_index;
    int i_;
  public:
    SCVectorSimpleSubBlockIter(SCVectorSimpleSubBlock*);
    virtual ~SCVectorSimpleSubBlockIter();
    int i();
    int j();
    double get();
    void set(double);
    operator int();
    void operator++();
    void reset();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
