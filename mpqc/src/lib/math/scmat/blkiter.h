
#ifndef _math_scmat_blkiter_h
#define _math_scmat_blkiter_h

#include <math/scmat/matrix.h>

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

#endif
