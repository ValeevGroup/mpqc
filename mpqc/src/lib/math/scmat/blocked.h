
#ifndef _math_scmat_blocked_h
#define _math_scmat_blocked_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/scmat/block.h>
#include <math/scmat/elemop.h>
#include <math/scmat/matrix.h>
#include <math/scmat/abstract.h>

class BlockedSCMatrixKit: public SCMatrixKit {
#   define CLASSNAME BlockedSCMatrixKit
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    BlockedSCMatrixKit();
    BlockedSCMatrixKit(const RefKeyVal&);
    BlockedSCMatrixKit(StateIn&);
    ~BlockedSCMatrixKit();
    void save_data_state(StateOut&);

    SCDimension* dimension(int n, const char* name=0);
    SCDimension* dimension(int n, int *nelem, const char* name=0);

    SCDimension* dimension(const RefSCMatrixKit&, int n, const char* name=0);
    SCDimension* dimension(const RefSCMatrixKit&, int n, int *nelem,
                           const char* name=0);
};
SavableState_REF_dec(BlockedSCMatrixKit);

class BlockedSCDimension: public SCDimension {
#   define CLASSNAME BlockedSCDimension
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int n_;
    int nblocks_;
    RefSCDimension *dims_;
    
    int *first_;
    int *last_;

    void init(const RefSCMatrixKit&, int *, const char *name=0);
    
  public:
    BlockedSCDimension(const RefSCMatrixKit&,
                       int n, int *nelem, const char* name = 0);
    BlockedSCDimension(const RefSCMatrixKit&,
                       int n, const char* name = 0);
    BlockedSCDimension(const RefKeyVal&);
    BlockedSCDimension(StateIn&);
    ~BlockedSCDimension();

    void save_data_state(StateOut&);

    int equiv(SCDimension*) const;
    
    int n();
    int n(int) const;
    int nblocks() const;
    
    int first(int i) const;
    int last(int i) const;
    
    int block(int i) const;
    
    RefSCDimension dim(int i) const;
    
    SCMatrix* create_matrix(SCDimension*);
    SymmSCMatrix* create_symmmatrix();
    DiagSCMatrix* create_diagmatrix();
    SCVector* create_vector();
};
SavableState_REF_dec(BlockedSCDimension);

class BlockedSCVector: public SCVector {
    friend class BlockedSCMatrix;
    friend class BlockedSymmSCMatrix;
    friend class BlockedDiagSCMatrix;
#   define CLASSNAME BlockedSCVector
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefBlockedSCDimension d;
    RefSCVector *vecs_;

    void resize(BlockedSCDimension*);

  public:
    BlockedSCVector();
    BlockedSCVector(BlockedSCDimension*);
    BlockedSCVector(const RefKeyVal&);
    BlockedSCVector(StateIn&);
    ~BlockedSCVector();

    void save_data_state(StateOut&);

    RefSCDimension dim();

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
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefBlockedSCDimension d1;
    RefBlockedSCDimension d2;
    RefSCMatrix *mats_;
    int nblocks_;
    
    void resize(BlockedSCDimension*, BlockedSCDimension*);

  public:
    BlockedSCMatrix();
    BlockedSCMatrix(const RefKeyVal&);
    BlockedSCMatrix(StateIn&);
    BlockedSCMatrix(BlockedSCDimension*,BlockedSCDimension*);
    ~BlockedSCMatrix();

    void save_data_state(StateOut&);

    // implementations and overrides of virtual functions
    RefSCDimension rowdim();
    RefSCDimension coldim();

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
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefBlockedSCDimension d;
    RefSymmSCMatrix *mats_;

    void resize(BlockedSCDimension*);

  public:
    BlockedSymmSCMatrix();
    BlockedSymmSCMatrix(StateIn&);
    BlockedSymmSCMatrix(const RefKeyVal&);
    BlockedSymmSCMatrix(BlockedSCDimension*);
    ~BlockedSymmSCMatrix();

    // implementations and overrides of virtual functions
    void save_data_state(StateOut&);

    RefSCDimension dim();

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

    void element_op(const RefSCElementOp&);
    void element_op(const RefSCElementOp2&,
                    SymmSCMatrix*);
    void element_op(const RefSCElementOp3&,
                    SymmSCMatrix*,SymmSCMatrix*);

    void print(const char* title=0,ostream& out=cout, int =10);

    // BlockedSymmSCMatrix specific functions
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
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefBlockedSCDimension d;
    RefDiagSCMatrix *mats_;

    void resize(BlockedSCDimension*);

  public:
    BlockedDiagSCMatrix();
    BlockedDiagSCMatrix(const RefKeyVal&);
    BlockedDiagSCMatrix(StateIn&);
    BlockedDiagSCMatrix(BlockedSCDimension*);
    ~BlockedDiagSCMatrix();

    // implementations and overrides of virtual functions
    void save_data_state(StateOut&);

    RefSCDimension dim();

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
