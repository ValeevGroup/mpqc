
#ifndef _math_scmat_block_h
#define _math_scmat_block_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>

class SCElementOp;
class SCElementOp2;
class SCElementOp3;

class SCMatrixBlock: public SavableState {
#   define CLASSNAME SCMatrixBlock
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCMatrixBlock();
    SCMatrixBlock(StateIn&s): SavableState(s) {}
    virtual ~SCMatrixBlock();

    virtual void process(SCElementOp*) = 0;
    virtual void process(SCElementOp2*, SCMatrixBlock*) = 0;
    virtual void process(SCElementOp3*, SCMatrixBlock*, SCMatrixBlock*) = 0;
};
SavableState_REF_dec(SCMatrixBlock);

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
#   define CLASSNAME SCMatrixBlockList
#   define HAVE_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
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
};
SavableState_REF_dec(SCMatrixBlockList);

class SCVectorSimpleBlock: public SCMatrixBlock {
#   define CLASSNAME SCVectorSimpleBlock
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    SCVectorSimpleBlock(int istart,int iend);
    SCVectorSimpleBlock(StateIn&);
    virtual ~SCVectorSimpleBlock();
    void save_data_state(StateOut&);
    int istart;
    int iend;
    double* data;

    void process(SCElementOp*);
    void process(SCElementOp2*, SCMatrixBlock*);
    void process(SCElementOp3*, SCMatrixBlock*, SCMatrixBlock*);
};
SavableState_REF_dec(SCVectorSimpleBlock);

class SCVectorSimpleSubBlock: public SCMatrixBlock {
#   define CLASSNAME SCVectorSimpleSubBlock
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
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
SavableState_REF_dec(SCVectorSimpleSubBlock);

class SCMatrixRectBlock: public SCMatrixBlock {
#   define CLASSNAME SCMatrixRectBlock
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
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

    void process(SCElementOp*);
    void process(SCElementOp2*, SCMatrixBlock*);
    void process(SCElementOp3*, SCMatrixBlock*, SCMatrixBlock*);
};
SavableState_REF_dec(SCMatrixRectBlock);

class SCMatrixRectSubBlock: public SCMatrixBlock {
#   define CLASSNAME SCMatrixRectSubBlock
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
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
SavableState_REF_dec(SCMatrixRectSubBlock);

class SCMatrixLTriBlock: public SCMatrixBlock {
#   define CLASSNAME SCMatrixLTriBlock
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    SCMatrixLTriBlock(int s,int e);
    SCMatrixLTriBlock(StateIn&);
    virtual ~SCMatrixLTriBlock();
    void save_data_state(StateOut&);
    int start;
    int end;
    double* data;

    void process(SCElementOp*);
    void process(SCElementOp2*, SCMatrixBlock*);
    void process(SCElementOp3*, SCMatrixBlock*, SCMatrixBlock*);
};
SavableState_REF_dec(SCMatrixLTriBlock);

// off diagonal sub block of a lower triangular matrix
class SCMatrixLTriSubBlock: public SCMatrixBlock {
#   define CLASSNAME SCMatrixLTriSubBlock
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
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
SavableState_REF_dec(SCMatrixLTriSubBlock);

class SCMatrixDiagBlock: public SCMatrixBlock {
#   define CLASSNAME SCMatrixDiagBlock
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
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

    void process(SCElementOp*);
    void process(SCElementOp2*, SCMatrixBlock*);
    void process(SCElementOp3*, SCMatrixBlock*, SCMatrixBlock*);
};
SavableState_REF_dec(SCMatrixDiagBlock);

class SCMatrixDiagSubBlock: public SCMatrixBlock {
#   define CLASSNAME SCMatrixDiagSubBlock
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
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
SavableState_REF_dec(SCMatrixDiagSubBlock);

#endif
