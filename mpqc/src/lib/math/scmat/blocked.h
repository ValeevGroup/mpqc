
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _math_scmat_blocked_h
#define _math_scmat_blocked_h

#include <math/scmat/block.h>
#include <math/scmat/elemop.h>
#include <math/scmat/matrix.h>
#include <math/scmat/abstract.h>

class BlockedSCMatrixKit;
class BlockedSCVector;
class BlockedSCMatrix;
class BlockedSymmSCMatrix;
class BlockedDiagSCMatrix;

class BlockedSCMatrixKit: public SCMatrixKit {
#   define CLASSNAME BlockedSCMatrixKit
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  private:
    RefSCMatrixKit subkit_;
  public:
    BlockedSCMatrixKit(const RefSCMatrixKit& subkit);
    BlockedSCMatrixKit(const RefKeyVal&);
    ~BlockedSCMatrixKit();
    SCMatrix* matrix(const RefSCDimension&,const RefSCDimension&);
    SymmSCMatrix* symmmatrix(const RefSCDimension&);
    DiagSCMatrix* diagmatrix(const RefSCDimension&);
    SCVector* vector(const RefSCDimension&);

    RefSCMatrixKit subkit() { return subkit_; }
};
DescribedClass_REF_dec(BlockedSCMatrixKit);

class BlockedSCVector: public SCVector {
    friend class BlockedSCMatrix;
    friend class BlockedSymmSCMatrix;
    friend class BlockedDiagSCMatrix;
#   define CLASSNAME BlockedSCVector
#   include <util/class/classd.h>
  private:
    RefSCMatrixKit subkit;
    RefSCVector *vecs_;

    void resize(SCDimension*);

  public:
    BlockedSCVector(const RefSCDimension&, BlockedSCMatrixKit*);
    ~BlockedSCVector();

    // Save and restore this in an implementation independent way.
    void save(StateOut&);
    void restore(StateIn&);

    void assign(double);
    void assign(SCVector*);
    void assign(const double*);

    double get_element(int);
    void set_element(int,double);
    void accumulate_element(int,double);

    void accumulate_product(SCMatrix*,SCVector*);
    void accumulate_product(SymmSCMatrix*,SCVector*);

    void accumulate(SCVector*);
    void accumulate(SCMatrix*);
    double scalar_product(SCVector*);

    void element_op(const RefSCElementOp&);
    void element_op(const RefSCElementOp2&,
                    SCVector*);
    void element_op(const RefSCElementOp3&,
                    SCVector*,SCVector*);
    void print(const char* title=0,ostream& out=cout, int =10);

    // BlockedSCVector specific functions
    RefSCDimension dim() { return d; }
    RefSCDimension dim(int);
    int nblocks() const;
    RefSCVector block(int);

    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);
};

class BlockedSCMatrix: public SCMatrix {
    friend class BlockedSymmSCMatrix;
    friend class BlockedDiagSCMatrix;
    friend BlockedSCVector;
#   define CLASSNAME BlockedSCMatrix
#   include <util/class/classd.h>
  private:
    RefSCMatrixKit subkit;
    RefSCMatrix *mats_;
    int nblocks_;
    
    void resize(SCDimension*, SCDimension*);

  public:
    BlockedSCMatrix(const RefSCDimension&, const RefSCDimension&,
                    BlockedSCMatrixKit*);
    ~BlockedSCMatrix();

    // Save and restore this in an implementation independent way.
    void save(StateOut&);
    void restore(StateIn&);

    void assign(double);
    double get_element(int,int);
    void set_element(int,int,double);
    void accumulate_element(int,int,double);

    SCMatrix * get_subblock(int,int,int,int);
    void assign_subblock(SCMatrix*, int,int,int,int,int=0,int=0);
    void accumulate_subblock(SCMatrix*, int,int,int,int,int=0,int=0);

    SCVector * get_row(int i);
    SCVector * get_column(int i);
    void assign_row(SCVector *v, int i);
    void assign_column(SCVector *v, int i);
    void accumulate_row(SCVector *v, int i);
    void accumulate_column(SCVector *v, int i);

    void accumulate_outer_product(SCVector*,SCVector*);
    void accumulate_product(SCMatrix*,SCMatrix*);
    void accumulate_product(SCMatrix*,SymmSCMatrix*);
    void accumulate_product(SCMatrix*,DiagSCMatrix*);
    void accumulate(SCMatrix*);
    void accumulate(SymmSCMatrix*);
    void accumulate(DiagSCMatrix*);
    void accumulate(SCVector*);

    void transpose_this();
    double invert_this();
    void svd_this(SCMatrix *U, DiagSCMatrix *sigma, SCMatrix *V);
    double solve_this(SCVector*);
    double determ_this();
    double trace();
    void gen_invert_this();
    void schmidt_orthog(SymmSCMatrix*,int);

    void element_op(const RefSCElementOp&);
    void element_op(const RefSCElementOp2&,
                    SCMatrix*);
    void element_op(const RefSCElementOp3&,
                    SCMatrix*,SCMatrix*);

    void print(const char* title=0,ostream& out=cout, int =10);

    // BlockedSCMatrix specific functions
    RefSCDimension rowdim() { return d1; }
    RefSCDimension coldim() { return d2; }
    RefSCDimension rowdim(int);
    RefSCDimension coldim(int);
    int nblocks() const;
    RefSCMatrix block(int);

    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);
};

class BlockedSymmSCMatrix: public SymmSCMatrix {
    friend class BlockedSCMatrix;
    friend class BlockedDiagSCMatrix;
    friend BlockedSCVector;
#   define CLASSNAME BlockedSymmSCMatrix
#   include <util/class/classd.h>
  private:
    RefSCMatrixKit subkit;
    RefSymmSCMatrix *mats_;

    void resize(SCDimension*);

  public:
    BlockedSymmSCMatrix(const RefSCDimension&,BlockedSCMatrixKit*);
    ~BlockedSymmSCMatrix();

    // Save and restore this in an implementation independent way.
    void save(StateOut&);
    void restore(StateIn&);

    double get_element(int,int);
    void set_element(int,int,double);
    void accumulate_element(int,int,double);

    SCMatrix * get_subblock(int,int,int,int);
    SymmSCMatrix * get_subblock(int,int);
    void assign_subblock(SCMatrix*, int,int,int,int);
    void assign_subblock(SymmSCMatrix*, int,int);
    void accumulate_subblock(SCMatrix*, int,int,int,int);
    void accumulate_subblock(SymmSCMatrix*, int,int);
    SCVector * get_row(int i);
    void assign_row(SCVector *v, int i);
    void accumulate_row(SCVector *v, int i);

    double invert_this();
    double determ_this();
    double trace();
    double solve_this(SCVector*);
    void gen_invert_this();

    double scalar_product(SCVector*);
    void diagonalize(DiagSCMatrix*,SCMatrix*);

    void accumulate(SymmSCMatrix*);
    void accumulate_symmetric_outer_product(SCVector*);
    void accumulate_symmetric_product(SCMatrix*);
    void accumulate_symmetric_sum(SCMatrix*);
    void accumulate_transform(SCMatrix*,SymmSCMatrix*);
    void accumulate_transform(SCMatrix*,DiagSCMatrix*);
    //void accumulate_transform(SymmSCMatrix*,SymmSCMatrix*);

    void element_op(const RefSCElementOp&);
    void element_op(const RefSCElementOp2&,
                    SymmSCMatrix*);
    void element_op(const RefSCElementOp3&,
                    SymmSCMatrix*,SymmSCMatrix*);

    void print(const char* title=0,ostream& out=cout, int =10);

    // BlockedSymmSCMatrix specific functions
    RefSCDimension dim() { return d; }
    RefSCDimension dim(int);
    int nblocks() const;
    RefSymmSCMatrix block(int);

    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);
};

class BlockedDiagSCMatrix: public DiagSCMatrix {
    friend BlockedSCMatrix;
    friend BlockedSymmSCMatrix;
    friend BlockedSCVector;
#   define CLASSNAME BlockedDiagSCMatrix
#   include <util/class/classd.h>
  private:
    RefSCMatrixKit subkit;
    RefDiagSCMatrix *mats_;

    void resize(SCDimension*);

  public:
    BlockedDiagSCMatrix(const RefSCDimension&,BlockedSCMatrixKit*);
    ~BlockedDiagSCMatrix();

    // Save and restore this in an implementation independent way.
    void save(StateOut&);
    void restore(StateIn&);

    double get_element(int);
    void set_element(int,double);
    void accumulate_element(int,double);
    void accumulate(DiagSCMatrix*);

    double invert_this();
    double determ_this();
    double trace();
    void gen_invert_this();

    void element_op(const RefSCElementOp&);
    void element_op(const RefSCElementOp2&,
                    DiagSCMatrix*);
    void element_op(const RefSCElementOp3&,
                    DiagSCMatrix*,DiagSCMatrix*);
    void print(const char* title=0,ostream& out=cout, int =10);

    // BlockedDiagSCMatrix specific functions
    RefSCDimension dim() { return d; }
    RefSCDimension dim(int);
    int nblocks() const;
    RefDiagSCMatrix block(int);

    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);
};

class BlockedSCElementOp : public SCElementOp {
#   define CLASSNAME BlockedSCElementOp
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    int current_block_;
    
  public:
    void working_on(int);
    int current_block() const;
};

class BlockedSCElementOp2 : public SCElementOp2 {
#   define CLASSNAME BlockedSCElementOp2
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    int current_block_;
    
  public:
    void working_on(int);
    int current_block() const;
};

class BlockedSCElementOp3 : public SCElementOp3 {
#   define CLASSNAME BlockedSCElementOp3
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    int current_block_;
    
  public:
    void working_on(int);
    int current_block() const;
};

#endif
