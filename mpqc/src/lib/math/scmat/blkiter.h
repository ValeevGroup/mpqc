
#ifndef _math_scmat_blkiter_h
#define _math_scmat_blkiter_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/scmat/block.h>

class SCMatrixRectBlock;
class SCMatrixLTriBlock;
class SCMatrixDiagBlock;
class SCVectorSimpleBlock;

class SCElementOp;
class SCElementOp2;
class SCElementOp3;

class SCMatrixBlockIter {
  public:
    SCMatrixBlockIter() {}
    virtual ~SCMatrixBlockIter();
    virtual int i() = 0;
    virtual int j() = 0;
    virtual void set(double) = 0;
    virtual double get() = 0;
    virtual operator int() = 0;
    virtual void operator++() = 0; // prefix ++
    void operator++(int) { operator++(); }
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

#endif
