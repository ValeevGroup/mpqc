
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _math_scmat_repl_h
#define _math_scmat_repl_h

#include <util/group/message.h>

#include <math/scmat/block.h>
#include <math/scmat/matrix.h>
#include <math/scmat/abstract.h>

class ReplSCMatrixKit: public SCMatrixKit {
#   define CLASSNAME ReplSCMatrixKit
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  public:
    ReplSCMatrixKit();
    ReplSCMatrixKit(const RefKeyVal&);
    ~ReplSCMatrixKit();
    SCMatrix* matrix(const RefSCDimension&,const RefSCDimension&);
    SymmSCMatrix* symmmatrix(const RefSCDimension&);
    DiagSCMatrix* diagmatrix(const RefSCDimension&);
    SCVector* vector(const RefSCDimension&);
};
DescribedClass_REF_dec(ReplSCMatrixKit);

class ReplSCMatrixListSubblockIter: public SCMatrixListSubblockIter {
  protected:
    RefMessageGrp grp_;
    double *data_;
    int ndata_;
  public:
    ReplSCMatrixListSubblockIter(Access,
                             const RefSCMatrixBlockList &list,
                             const RefMessageGrp &grp,
                             double *data, int ndata);
    ~ReplSCMatrixListSubblockIter();
};

class ReplSCVector: public SCVector {
    friend class ReplSCMatrix;
    friend class ReplSymmSCMatrix;
    friend class ReplDiagSCMatrix;
#   define CLASSNAME ReplSCVector
#   include <util/class/classd.h>
  protected:
    RefSCMatrixBlockList blocklist;
    double* vector;
    void init_blocklist();
    void before_elemop();
    void after_elemop();
  public:
    ReplSCVector(const RefSCDimension&,ReplSCMatrixKit*);
    ~ReplSCVector();
    void assign(double);
    void assign(SCVector*);
    void assign(const double*);

    void set_element(int,double);
    void accumulate_element(int,double);
    double get_element(int);
    void accumulate_product(SymmSCMatrix*,SCVector*);
    void accumulate_product(SCMatrix*,SCVector*);
    void accumulate(SCVector*);
    void accumulate(SCMatrix*);
    double scalar_product(SCVector*);
    void element_op(const RefSCElementOp&);
    void element_op(const RefSCElementOp2&,
                    SCVector*);
    void element_op(const RefSCElementOp3&,
                    SCVector*,SCVector*);
    void print(const char* title=0,ostream& out=cout, int =10);

    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);

    RefReplSCMatrixKit skit();
};

class ReplSCMatrix: public SCMatrix {
    friend class ReplSymmSCMatrix;
    friend class ReplDiagSCMatrix;
    friend ReplSCVector;
#   define CLASSNAME ReplSCMatrix
#   include <util/class/classd.h>
  protected:
    RefSCMatrixBlockList blocklist;
    double* matrix;
    double** rows;
  protected:
    // utility functions
    int compute_offset(int,int);
    void init_blocklist();

    void before_elemop();
    void after_elemop();
  public:
    ReplSCMatrix(const RefSCDimension&,const RefSCDimension&,
                 ReplSCMatrixKit*);
    ~ReplSCMatrix();

    // implementations and overrides of virtual functions
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

    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);

    RefReplSCMatrixKit skit();
};

class ReplSymmSCMatrix: public SymmSCMatrix {
    friend class ReplSCMatrix;
    friend class ReplDiagSCMatrix;
    friend ReplSCVector;
#   define CLASSNAME ReplSymmSCMatrix
#   include <util/class/classd.h>
  protected:
    RefSCMatrixBlockList blocklist;
    double* matrix;
    double** rows;
  protected:
    // utility functions
    int compute_offset(int,int);
    void init_blocklist();

    void before_elemop();
    void after_elemop();
  public:
    ReplSymmSCMatrix(const RefSCDimension&, ReplSCMatrixKit*);
    ~ReplSymmSCMatrix();

    // implementations and overrides of virtual functions
    void assign(double);
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

    void accumulate_product(SCMatrix*,SCMatrix*);
    void accumulate(SymmSCMatrix*);
    double invert_this();
    double solve_this(SCVector*);
    double trace();
    double determ_this();
    void gen_invert_this();

    double scalar_product(SCVector*);
    void diagonalize(DiagSCMatrix*,SCMatrix*);
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

    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);

    RefReplSCMatrixKit skit();
};

class ReplDiagSCMatrix: public DiagSCMatrix {
    friend ReplSCMatrix;
    friend ReplSymmSCMatrix;
    friend ReplSCVector;
#   define CLASSNAME ReplDiagSCMatrix
#   include <util/class/classd.h>
  protected:
    RefSCMatrixBlockList blocklist;
    void init_blocklist();
    double* matrix;

    void before_elemop();
    void after_elemop();
  public:
    ReplDiagSCMatrix(const RefSCDimension&, ReplSCMatrixKit*);
    ~ReplDiagSCMatrix();

    // implementations and overrides of virtual functions
    void assign(double);
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

    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);

    RefReplSCMatrixKit skit();
};

#endif
