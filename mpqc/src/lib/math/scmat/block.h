//
// block.h
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

#ifndef _math_scmat_block_h
#define _math_scmat_block_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>

namespace sc {

class SCElementOp;
class SCElementOp2;
class SCElementOp3;

/** SCMatrixBlock is the base clase for all types of blocks
    that comprise matrices and vectors. */
class SCMatrixBlock: public SavableState {
  public:
    int blocki, blockj;
  public:
    SCMatrixBlock();
    SCMatrixBlock(StateIn&s);
    virtual ~SCMatrixBlock();
    void save_data_state(StateOut&s);

    /** Return of copy of this.  A runtime error will be generated
        for blocks that cannot do a deepcopy.  These routines are only used
        internally in the matrix library. */
    virtual SCMatrixBlock *deepcopy() const;

    /** Return a pointer to the block's data and the number of elements
        in the block.  Some blocks cannot provide this information and
        a runtime error will be generated if these members are called.
        These routines are only used internally in the matrix library. */
    virtual double *dat();
    virtual int ndat() const;

    // These routines are obsolete.
    virtual void process(SCElementOp*) = 0;
    virtual void process(SCElementOp2*, SCMatrixBlock*) = 0;
    virtual void process(SCElementOp3*, SCMatrixBlock*, SCMatrixBlock*) = 0;
};


class SCMatrixBlockListLink {
  private:
    void operator = (const SCMatrixBlockListLink&) {}  // disallowed
    SCMatrixBlock* _block;
    SCMatrixBlockListLink* _next;
  public:
    SCMatrixBlockListLink(SCMatrixBlock*, SCMatrixBlockListLink* = 0);
    ~SCMatrixBlockListLink();
    void block(SCMatrixBlock*);
    void next(SCMatrixBlockListLink* link) { _next = link; }
    SCMatrixBlock* block() { return _block; }
    SCMatrixBlockListLink* next() { return _next; }
};

class SCMatrixBlockListIter {
  private:
    SCMatrixBlockListLink* link;
  public:
    SCMatrixBlockListIter(): link(0) {}
    SCMatrixBlockListIter(SCMatrixBlockListLink*l): link(l) {}
    int operator !=(const SCMatrixBlockListIter p) const {
        return link != p.link;
      }
    void operator ++() { link = link->next(); }
    void operator ++(int) { link = link->next(); }
    SCMatrixBlock* block() const { return link->block(); }
};

class SCMatrixBlockList: public SavableState {
  private:
    SCMatrixBlockListLink* _begin;
  public:
    SCMatrixBlockList();
    SCMatrixBlockList(StateIn&);
    ~SCMatrixBlockList();
    void save_data_state(StateOut&);
    void insert(SCMatrixBlock*);
    void append(SCMatrixBlock*);
    SCMatrixBlockListIter begin() { return _begin; }
    SCMatrixBlockListIter end() { return 0; }
    SCMatrixBlockList *deepcopy();
};


/** The SCVectorSimpleBlock describes a piece of a
vector.  The following bit of code illustrates the data layout:
fill(double *vector, SCVectorSimpleBlock &b)
{
  int i,offset=0;
  for (i=b.istart; i<b.iend; i++,offset++) {
      vector[i] = b.data[offset];
  }
}
*/
class SCVectorSimpleBlock: public SCMatrixBlock {
  public:
    SCVectorSimpleBlock(int istart,int iend);
    SCVectorSimpleBlock(StateIn&);
    virtual ~SCVectorSimpleBlock();
    void save_data_state(StateOut&);
    int istart;
    int iend;
    double* data;

    SCMatrixBlock *deepcopy() const;

    void process(SCElementOp*);
    void process(SCElementOp2*, SCMatrixBlock*);
    void process(SCElementOp3*, SCMatrixBlock*, SCMatrixBlock*);

    double *dat();
    int ndat() const;
};


/** The SCVectorSimpleSubBlock describes a subblock of a
vector.  The following bit of code illustrates the data layout:
fill(double *vector, SCVectorSimpleSubBlock &b)
{
  int i,offset=b.offset;
  for (i=b.istart; i<b.iend; i++,offset++) {
      vector[i] = b.data[offset];
  }
}
*/
class SCVectorSimpleSubBlock: public SCMatrixBlock {
  public:
    SCVectorSimpleSubBlock(int istart,int iend, int offset, double* data);
    SCVectorSimpleSubBlock(StateIn&);
    virtual ~SCVectorSimpleSubBlock();
    void save_data_state(StateOut&);
    int istart;
    int iend;
    int offset;
    double* data;

    void process(SCElementOp*);
    void process(SCElementOp2*, SCMatrixBlock*);
    void process(SCElementOp3*, SCMatrixBlock*, SCMatrixBlock*);
};


/** The SCMatrixRectBlock describes a rectangular piece of a
matrix.  The following bit of code illustrates the data layout:
fill(double **matrix, SCMatrixRectBlock &b)
{
  int offset=0;
  for (int i=b.istart; i<b.iend; i++) {
    for (int j=b.jstart; j<b.jend; j++,offset++) {
      matrix[i][j] = b.data[offset];
    }
  }
}
*/
class SCMatrixRectBlock: public SCMatrixBlock {
  public:
    SCMatrixRectBlock(int is, int ie, int js, int je);
    SCMatrixRectBlock(StateIn&);
    virtual ~SCMatrixRectBlock();
    void save_data_state(StateOut&);
    int istart;
    int jstart;
    int iend;
    int jend;
    double* data;

    SCMatrixBlock *deepcopy() const;

    void process(SCElementOp*);
    void process(SCElementOp2*, SCMatrixBlock*);
    void process(SCElementOp3*, SCMatrixBlock*, SCMatrixBlock*);

    double *dat();
    int ndat() const;
};


/** The SCMatrixRectSubBlock describes a rectangular piece of a
matrix.  The following bit of code illustrates the data layout:
fill(double **matrix, SCMatrixRectSubBlock &b)
{
  int offset=b.istart * b.istride + b.jstart;
  for (int i=b.istart; i<b.iend; i++) {
    for (int j=b.jstart; j<b.jend; j++,offset++) {
      matrix[i][j] = b.data[offset];
    }
  offset += b.istride - (b.jend - b.jstart);
  }
}
*/
class SCMatrixRectSubBlock: public SCMatrixBlock {
  public:
    SCMatrixRectSubBlock(int is, int ie, int istride, int js, int je,
                         double* data);
    SCMatrixRectSubBlock(StateIn&);
    // does not delete the data member
    virtual ~SCMatrixRectSubBlock();
    // does not save the data member
    void save_data_state(StateOut&);
    int istart;
    int jstart;
    int iend;
    int jend;
    int istride;
    double* data;

    void process(SCElementOp*);
    void process(SCElementOp2*, SCMatrixBlock*);
    void process(SCElementOp3*, SCMatrixBlock*, SCMatrixBlock*);
};


/** The SCMatrixLTriBlock describes a triangular piece of a
matrix.  The following bit of code illustrates the data layout:
fill(double **matrix, SCMatrixLTriBlock &b)
{
  int offset=0;
  for (int i=b.start; i<b.end; i++) {
    for (int j=b.start; j<=i; j++,offset++) {
      matrix[i][j] = b.data[offset];
    }
  }
}
*/
class SCMatrixLTriBlock: public SCMatrixBlock {
  public:
    SCMatrixLTriBlock(int s,int e);
    SCMatrixLTriBlock(StateIn&);
    virtual ~SCMatrixLTriBlock();
    void save_data_state(StateOut&);
    int start;
    int end;
    double* data;

    SCMatrixBlock *deepcopy() const;

    void process(SCElementOp*);
    void process(SCElementOp2*, SCMatrixBlock*);
    void process(SCElementOp3*, SCMatrixBlock*, SCMatrixBlock*);

    double *dat();
    int ndat() const;
};


/** The SCMatrixLTriSubBlock describes a triangular subblock of a
matrix.  The following bit of code illustrates the data layout:
fill(double **matrix, SCMatrixLTriSubBlock &b)
{
  int offset=(b.istart*(b.istart+1)>>1) + b.jstart;
  for (int i=b.start; i<b.end; i++) {
    for (int j=b.start; j<=i && j<b.jend; j++,offset++) {
      matrix[i][j] = b.data[offset];
    }
  if (j>i) offset += b.istart;
  else offset += i + b.jstart - b.jend;
  }
}
*/
class SCMatrixLTriSubBlock: public SCMatrixBlock {
  public:
    SCMatrixLTriSubBlock(int is,int ie,int js,int je,double*data);
    SCMatrixLTriSubBlock(StateIn&);
    // does not delete the data member
    virtual ~SCMatrixLTriSubBlock();
    // does not save the data member
    void save_data_state(StateOut&);
    int istart;
    int iend;
    int jstart;
    int jend;
    double* data;

    void process(SCElementOp*);
    void process(SCElementOp2*, SCMatrixBlock*);
    void process(SCElementOp3*, SCMatrixBlock*, SCMatrixBlock*);
};


/** The SCMatrixDiagBlock describes a diagonal piece of a
matrix.  The following bit of code illustrates the data layout:
fill(double **matrix, SCMatrixDiagBlock &b)
{
  int i,j,offset=0;
  for (i=b.istart,j=b.jstart; i<b.iend; i++,j++,offset++) {
      matrix[i][j] = b.data[offset];
  }
}
*/
class SCMatrixDiagBlock: public SCMatrixBlock {
  public:
    SCMatrixDiagBlock(int istart,int iend,int jstart);
    SCMatrixDiagBlock(int istart,int iend);
    SCMatrixDiagBlock(StateIn&);
    virtual ~SCMatrixDiagBlock();
    void save_data_state(StateOut&);
    int istart;
    int jstart;
    int iend;
    double* data;

    SCMatrixBlock *deepcopy() const;

    void process(SCElementOp*);
    void process(SCElementOp2*, SCMatrixBlock*);
    void process(SCElementOp3*, SCMatrixBlock*, SCMatrixBlock*);

    double *dat();
    int ndat() const;
};


/** The SCMatrixDiagSubBlock describes a diagonal subblock of a
matrix.  The following bit of code illustrates the data layout:
fill(double **matrix, SCMatrixDiagSubBlock &b)
{
  int i,j,offset=b.offset;
  for (i=b.istart,j=b.jstart; i<b.iend; i++,j++,offset++) {
      matrix[i][j] = b.data[offset];
  }
}
*/
class SCMatrixDiagSubBlock: public SCMatrixBlock {
  public:
    SCMatrixDiagSubBlock(int istart,int iend,int jstart, int offset,
                         double*data);
    SCMatrixDiagSubBlock(int istart,int iend, int offset, double*data);
    SCMatrixDiagSubBlock(StateIn&);
    // does not delete the data member
    virtual ~SCMatrixDiagSubBlock();
    // does not save the data member
    void save_data_state(StateOut&);
    int istart;
    int jstart;
    int iend;
    int offset;
    double* data;

    void process(SCElementOp*);
    void process(SCElementOp2*, SCMatrixBlock*);
    void process(SCElementOp3*, SCMatrixBlock*, SCMatrixBlock*);
};


// //////////////////////////////////////////////////////////////////
// Classes that iterate through the blocks of a matrix.

/** Objects of class SCMatrixSubblockIter are used to iterate through the
    blocks of a matrix.  The object must be deleted before using the matrix
    that owns the blocks that SCMatrixSubblockIter is iterating through. */
class SCMatrixSubblockIter: public RefCount {
  public:
    enum Access { Read, Write, Accum, None };
  protected:
    Access access_;
  public:
    /** The access variable should be one of Read, Write, Accum, and None,
        with the SCMatrixSubblockIter:: scope operator applied. */
    SCMatrixSubblockIter(Access access): access_(access) {}
    ~SCMatrixSubblockIter();
    /// Start at the beginning.
    virtual void begin() = 0;
    /// Returns nonzero if there is another block.
    virtual int ready() = 0;
    /// Proceed to the next block.
    virtual void next() = 0;
    /// Return the current block.
    virtual SCMatrixBlock *block() = 0;
    /// Return the type of Access allowed for these blocks.
    Access access() const { return access_; }
};


class SCMatrixSimpleSubblockIter: public SCMatrixSubblockIter {
  protected:
    Ref<SCMatrixBlock> block_;
    int ready_;
  public:
    SCMatrixSimpleSubblockIter(Access, const Ref<SCMatrixBlock> &b);
    void begin();
    int ready();
    void next();
    SCMatrixBlock *block();
};

class SCMatrixListSubblockIter: public SCMatrixSubblockIter {
  protected:
    Ref<SCMatrixBlockList> list_;
    SCMatrixBlockListIter iter_;
  public:
    SCMatrixListSubblockIter(Access, const Ref<SCMatrixBlockList> &list);
    void begin();
    int ready();
    void next();
    SCMatrixBlock *block();
};

class SCMatrixNullSubblockIter: public SCMatrixSubblockIter {
  public:
    SCMatrixNullSubblockIter();
    SCMatrixNullSubblockIter(Access);
    void begin();
    int ready();
    void next();
    SCMatrixBlock *block();
};

class SCMatrixCompositeSubblockIter: public SCMatrixSubblockIter {
  protected:
    int niters_;
    Ref<SCMatrixSubblockIter> *iters_;
    int iiter_;
  public:
    SCMatrixCompositeSubblockIter(Access, int niter);
    SCMatrixCompositeSubblockIter(Ref<SCMatrixSubblockIter>&,
                                  Ref<SCMatrixSubblockIter>&);
    ~SCMatrixCompositeSubblockIter();
    void set_iter(int i, const Ref<SCMatrixSubblockIter> &);
    void begin();
    int ready();
    void next();
    SCMatrixBlock *block();
    int current_block() const { return iiter_; }
};


class SCMatrixJointSubblockIter: public SCMatrixSubblockIter {
  protected:
    int niters_;
    Ref<SCMatrixSubblockIter> *iters_;
  public:
    SCMatrixJointSubblockIter(const Ref<SCMatrixSubblockIter>&,
                              const Ref<SCMatrixSubblockIter>&,
                              const Ref<SCMatrixSubblockIter>& = 0,
                              const Ref<SCMatrixSubblockIter>& = 0,
                              const Ref<SCMatrixSubblockIter>& = 0);
    ~SCMatrixJointSubblockIter();
    void begin();
    int ready();
    void next();
    SCMatrixBlock *block();
    SCMatrixBlock *block(int i);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
