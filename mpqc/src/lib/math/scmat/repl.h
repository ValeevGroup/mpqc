
#ifndef _math_scmat_repl_h
#define _math_scmat_repl_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/group/message.h>

#include <math/scmat/block.h>
#include <math/scmat/matrix.h>
#include <math/scmat/abstract.h>

class ReplSCMatrixKit: public SCMatrixKit {
#   define CLASSNAME ReplSCMatrixKit
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefMessageGrp grp_;
  public:
    ReplSCMatrixKit();
    ReplSCMatrixKit(const RefKeyVal&);
    ~ReplSCMatrixKit();

    RefMessageGrp messagegrp() const { return grp_; }

    SCDimension* dimension(int n, const char* name=0);

    SCMatrix* restore_matrix(StateIn&,
                             const RefSCDimension&,const RefSCDimension&);
    SymmSCMatrix* restore_symmmatrix(StateIn&, const RefSCDimension&);
    DiagSCMatrix* restore_diagmatrix(StateIn&, const RefSCDimension&);
    SCVector* restore_vector(StateIn&, const RefSCDimension&);
};

class ReplSCDimension: public SCDimension {
#   define CLASSNAME ReplSCDimension
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int n_;
    int* blocks_;
    int nblocks_;
    RefMessageGrp grp_;
  public:
    ReplSCDimension(int n, const RefMessageGrp&, const char* name = 0);
    ~ReplSCDimension();
    int n();
    SCMatrix* create_matrix(SCDimension*);
    SymmSCMatrix* create_symmmatrix();
    DiagSCMatrix* create_diagmatrix();
    SCVector* create_vector();

    RefMessageGrp messagegrp() const { return grp_; }
    int blockstart(int b) { return blocks_[b]; }
    int blockfence(int b) { return blocks_[b+1]; }
    int nblock() { return nblocks_; }
};
SavableState_REF_dec(ReplSCDimension);

class ReplSCVector: public SCVector {
    friend class ReplSCMatrix;
    friend class ReplSymmSCMatrix;
    friend class ReplDiagSCMatrix;
#   define CLASSNAME ReplSCVector
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefReplSCDimension d;
    RefSCMatrixBlockList blocklist;
    double* vector;
    void init_blocklist();
    void before_elemop();
    void after_elemop();
  public:
    ReplSCVector(ReplSCDimension*);
    ~ReplSCVector();
    void assign(double);
    void assign(SCVector*);
    void assign(const double*);

    RefSCDimension dim();
    void set_element(int,double);
    void accumulate_element(int,double);
    double get_element(int);
    void accumulate_product(SymmSCMatrix*,SCVector*);
    void accumulate_product(SCMatrix*,SCVector*);
    void accumulate(SCVector*);
    double scalar_product(SCVector*);
    void element_op(const RefSCElementOp&);
    void element_op(const RefSCElementOp2&,
                    SCVector*);
    void element_op(const RefSCElementOp3&,
                    SCVector*,SCVector*);
    void print(const char* title=0,ostream& out=cout, int =10);

    RefMessageGrp messagegrp() { return d->messagegrp(); }
};

class ReplSCMatrix: public SCMatrix {
    friend class ReplSymmSCMatrix;
    friend class ReplDiagSCMatrix;
    friend ReplSCVector;
#   define CLASSNAME ReplSCMatrix
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefReplSCDimension d1;
    RefReplSCDimension d2;
    RefSCMatrixBlockList blocklist;
    double* matrix;
    double** rows;
  private:
    // utility functions
    int compute_offset(int,int);
    void init_blocklist();

    void before_elemop();
    void after_elemop();
  public:
    ReplSCMatrix(ReplSCDimension*,ReplSCDimension*);
    ~ReplSCMatrix();

    // implementations and overrides of virtual functions
    void assign(double);
    RefSCDimension rowdim();
    RefSCDimension coldim();
    double get_element(int,int);
    void set_element(int,int,double);
    void accumulate_element(int,int,double);
    SCMatrix * get_subblock(int,int,int,int);
    void assign_subblock(SCMatrix*, int,int,int,int);
    void accumulate_subblock(SCMatrix*, int,int,int,int);
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
    void transpose_this();
    double invert_this();
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

    RefMessageGrp messagegrp() { return d1->messagegrp(); }
};

class ReplSymmSCMatrix: public SymmSCMatrix {
    friend class ReplSCMatrix;
    friend class ReplDiagSCMatrix;
    friend ReplSCVector;
#   define CLASSNAME ReplSymmSCMatrix
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefReplSCDimension d;
    RefSCMatrixBlockList blocklist;
    double* matrix;
    double** rows;
  private:
    // utility functions
    int compute_offset(int,int);
    void init_blocklist();

    void before_elemop();
    void after_elemop();
  public:
    ReplSymmSCMatrix(ReplSCDimension*);
    ~ReplSymmSCMatrix();

    // implementations and overrides of virtual functions
    void assign(double);
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
    void element_op(const RefSCElementOp&);
    void element_op(const RefSCElementOp2&,
                    SymmSCMatrix*);
    void element_op(const RefSCElementOp3&,
                    SymmSCMatrix*,SymmSCMatrix*);
    void print(const char* title=0,ostream& out=cout, int =10);

    RefMessageGrp messagegrp() { return d->messagegrp(); }
};

class ReplDiagSCMatrix: public DiagSCMatrix {
    friend ReplSCMatrix;
    friend ReplSymmSCMatrix;
    friend ReplSCVector;
#   define CLASSNAME ReplDiagSCMatrix
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefReplSCDimension d;
    RefSCMatrixBlockList blocklist;
    void init_blocklist();
    double* matrix;

    void before_elemop();
    void after_elemop();
  public:
    ReplDiagSCMatrix(ReplSCDimension*);
    ~ReplDiagSCMatrix();

    // implementations and overrides of virtual functions
    void assign(double);
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

    RefMessageGrp messagegrp() { return d->messagegrp(); }
};

#endif
