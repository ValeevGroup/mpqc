
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _math_scmat_local_h
#define _math_scmat_local_h

#include <math/scmat/block.h>
#include <math/scmat/matrix.h>
#include <math/scmat/abstract.h>

class LocalSCMatrixKit;
class LocalSCVector;
class LocalSCMatrix;
class LocalSymmSCMatrix;
class LocalDiagSCMatrix;

class LocalSCMatrixKit: public SCMatrixKit {
#   define CLASSNAME LocalSCMatrixKit
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  public:
    LocalSCMatrixKit();
    LocalSCMatrixKit(const RefKeyVal&);
    ~LocalSCMatrixKit();
    SCMatrix* matrix(const RefSCDimension&,const RefSCDimension&);
    SymmSCMatrix* symmmatrix(const RefSCDimension&);
    DiagSCMatrix* diagmatrix(const RefSCDimension&);
    SCVector* vector(const RefSCDimension&);
};

class LocalSCVector: public SCVector {
    friend class LocalSCMatrix;
    friend class LocalSymmSCMatrix;
    friend class LocalDiagSCMatrix;
#   define CLASSNAME LocalSCVector
#   include <util/class/classd.h>
  private:
    RefSCVectorSimpleBlock block;

    void resize(int);
  public:
    LocalSCVector();
    LocalSCVector(const RefSCDimension&,LocalSCMatrixKit*);
    ~LocalSCVector();
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
    void vprint(const char* title=0,ostream& out=cout, int =10);

    // return a pointer to the data for fast access
    double *get_data();
    
    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);
};

class LocalSCMatrix: public SCMatrix {
    friend class LocalSymmSCMatrix;
    friend class LocalDiagSCMatrix;
    friend LocalSCVector;
#   define CLASSNAME LocalSCMatrix
#   include <util/class/classd.h>
  private:
    RefSCMatrixRectBlock block;
    double** rows;
  private:
    // utility functions
    int compute_offset(int,int);
    void resize(int,int);
  public:
    LocalSCMatrix(const RefSCDimension&,const RefSCDimension&,
                  LocalSCMatrixKit*);
    ~LocalSCMatrix();

    // implementations and overrides of virtual functions
    void assign(double);
    void assign(SCMatrix*);
    void assign(const double*);
    void assign(const double**);
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
    void accumulate_product(SymmSCMatrix*,SCMatrix*);
    void accumulate_product(DiagSCMatrix*,SCMatrix*);
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
    void vprint(const char* title=0,ostream& out=cout, int =10);

    // return a pointer to the data for fast access
    double *get_data();
    double **get_rows();
    
    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);
};

class LocalSymmSCMatrix: public SymmSCMatrix {
    friend class LocalSCMatrix;
    friend class LocalDiagSCMatrix;
    friend LocalSCVector;
#   define CLASSNAME LocalSymmSCMatrix
#   include <util/class/classd.h>
  private:
    RefSCMatrixLTriBlock block;
    double** rows;
  private:
    // utility functions
    int compute_offset(int,int);
    void resize(int n);
  public:
    LocalSymmSCMatrix(const RefSCDimension&, LocalSCMatrixKit*);
    ~LocalSymmSCMatrix();

    // implementations and overrides of virtual functions
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
    void accumulate_transform(SymmSCMatrix*,SymmSCMatrix*);
    void element_op(const RefSCElementOp&);
    void element_op(const RefSCElementOp2&,
                    SymmSCMatrix*);
    void element_op(const RefSCElementOp3&,
                    SymmSCMatrix*,SymmSCMatrix*);
    void vprint(const char* title=0,ostream& out=cout, int =10);

    // return a pointer to the data for fast access
    double *get_data();
    double **get_rows();
    
    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);
};

class LocalDiagSCMatrix: public DiagSCMatrix {
    friend LocalSCMatrix;
    friend LocalSymmSCMatrix;
    friend LocalSCVector;
#   define CLASSNAME LocalDiagSCMatrix
#   include <util/class/classd.h>
  private:
    RefSCMatrixDiagBlock block;
    void resize(int n);
  public:
    LocalDiagSCMatrix(const RefSCDimension&, LocalSCMatrixKit*);
    ~LocalDiagSCMatrix();

    // implementations and overrides of virtual functions
    void save_data_state(StateOut&);
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
    void vprint(const char* title=0,ostream& out=cout, int =10);

    // return a pointer to the data for fast access
    double *get_data();

    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);
};

#endif
