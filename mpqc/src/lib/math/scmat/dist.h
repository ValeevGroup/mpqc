
#ifndef _math_scmat_dist_h
#define _math_scmat_dist_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/group/message.h>
#include <util/group/mstate.h>

#include <math/scmat/block.h>
#include <math/scmat/matrix.h>
#include <math/scmat/abstract.h>

class DistSCMatrixKit: public SCMatrixKit {
#   define CLASSNAME DistSCMatrixKit
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefMessageGrp grp_;
  public:
    DistSCMatrixKit(const RefMessageGrp &grp = 0);
    DistSCMatrixKit(const RefKeyVal&);
    ~DistSCMatrixKit();

    RefMessageGrp messagegrp() const { return grp_; }

    SCDimension* dimension(int n, const char* name=0);
    SCDimension* dimension(int n, int nblocks, const int *blocksizes,
                           const char* name=0);
};
SavableState_REF_dec(DistSCMatrixKit);

class SCBlockInfo: public VRefCount {
  private:
    int n_;
    int nblocks_;
    int *start_;
    int *size_;
  public:
    SCBlockInfo(int n_, int nblocks = 0, const int *blocksizes_ = 0);
    ~SCBlockInfo();
    int equiv(SCBlockInfo *);
    int nelem() const { return n_; }
    int nblock() const { return nblocks_; }
    int start(int i) const { return start_[i]; }
    int size(int i) const { return size_[i]; }
    int fence(int i) const { return start_[i] + size_[i]; }
    void elem_to_block(int i, int &block, int &offset);
};
REF_dec(SCBlockInfo);

class DistSCDimension: public SCDimension {
#   define CLASSNAME DistSCDimension
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int n_;
    RefSCBlockInfo blocks_;
    RefDistSCMatrixKit kit_;
  public:
    DistSCDimension(int n, const RefDistSCMatrixKit&,
                    const char* name = 0);
    DistSCDimension(int n, const RefDistSCMatrixKit&,
                    int nblocks, const int *blocksizes,
                    const char* name = 0);
    ~DistSCDimension();
    int equiv(SCDimension*) const;
    int n();
    SCMatrix* create_matrix(SCDimension*);
    SymmSCMatrix* create_symmmatrix();
    DiagSCMatrix* create_diagmatrix();
    SCVector* create_vector();

    RefMessageGrp messagegrp() { return kit_->messagegrp(); }

    RefSCBlockInfo blocks() { return blocks_; }
};
SavableState_REF_dec(DistSCDimension);

class DistSCVector: public SCVector {
    friend class DistSCMatrix;
    friend class DistSymmSCMatrix;
    friend class DistDiagSCMatrix;
#   define CLASSNAME DistSCVector
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefDistSCDimension d;
    RefSCMatrixBlockList blocklist;

    void init_blocklist();
    double *find_element(int i);
    int element_to_node(int i);
    int block_to_node(int);
    RefSCMatrixBlock block_to_block(int);
    void error(const char *);
  public:
    DistSCVector(DistSCDimension*);
    ~DistSCVector();
    void assign(SCVector*);
    void assign(const double*);
    void convert(double* v);

    RefSCDimension dim();
    void set_element(int,double);
    void accumulate_element(int,double);
    double get_element(int);
    void accumulate(SCVector*);
    void accumulate(SCMatrix*m);
    double scalar_product(SCVector*);
    void accumulate_product(SCMatrix *, SCVector *);
    void element_op(const RefSCElementOp&);
    void element_op(const RefSCElementOp2&,
                    SCVector*);
    void element_op(const RefSCElementOp3&,
                    SCVector*,SCVector*);
    void print(const char* title=0,ostream& out=cout, int =10);

    RefMessageGrp messagegrp() { return d->messagegrp(); }

    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);
};

class DistSCMatrix: public SCMatrix {
    friend class DistSymmSCMatrix;
    friend class DistDiagSCMatrix;
    friend DistSCVector;
#   define CLASSNAME DistSCMatrix
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefDistSCDimension d1;
    RefDistSCDimension d2;
    RefSCMatrixBlockList blocklist;

    int vecoff;
    int nvec;
    double **vec;
  protected:
    // utility functions
    void init_blocklist();
    void error(const char *);
    double *find_element(int i, int j);
    int element_to_node(int i, int j);
    int block_to_node(int,int);
    RefSCMatrixBlock block_to_block(int, int);
    RefSCBlockInfo rowblocks() { return d1->blocks(); }
    RefSCBlockInfo colblocks() { return d2->blocks(); }
    
    enum VecOp {CopyFromVec, CopyToVec, AccumFromVec, AccumToVec};
    enum Form { Row, Col } form;
    void create_vecform(Form, int nvec = -1);
    void delete_vecform();
    void vecform_op(VecOp op, int *ivec = 0);
    void vecform_zero();
  public:
    DistSCMatrix(DistSCDimension*,DistSCDimension*);
    ~DistSCMatrix();

    // implementations and overrides of virtual functions
    RefSCDimension rowdim();
    RefSCDimension coldim();
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
    void accumulate(SCMatrix*);
    void accumulate(SymmSCMatrix*);
    void accumulate(DiagSCMatrix*);
    void accumulate(SCVector*);
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

    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);
};

class DistSymmSCMatrix: public SymmSCMatrix {
    friend class DistSCMatrix;
    friend class DistDiagSCMatrix;
    friend DistSCVector;
#   define CLASSNAME DistSymmSCMatrix
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefDistSCDimension d;
    RefSCMatrixBlockList blocklist;
  protected:
    // utility functions
    void init_blocklist();
    double *find_element(int i, int j);
    int element_to_node(int i, int j);
    int block_to_node(int,int);
    RefSCMatrixBlock block_to_block(int, int);

    void error(const char *msg);
  public:
    DistSymmSCMatrix(DistSCDimension*);
    ~DistSymmSCMatrix();

    // implementations and overrides of virtual functions
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

    void diagonalize(DiagSCMatrix*,SCMatrix*);
    void accumulate_symmetric_sum(SCMatrix*);
    void element_op(const RefSCElementOp&);
    void element_op(const RefSCElementOp2&,
                    SymmSCMatrix*);
    void element_op(const RefSCElementOp3&,
                    SymmSCMatrix*,SymmSCMatrix*);

    RefMessageGrp messagegrp() { return d->messagegrp(); }

    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);
};

class DistDiagSCMatrix: public DiagSCMatrix {
    friend DistSCMatrix;
    friend DistSymmSCMatrix;
    friend DistSCVector;
#   define CLASSNAME DistDiagSCMatrix
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    RefDistSCDimension d;
    RefSCMatrixBlockList blocklist;

    void init_blocklist();
    double *find_element(int i);
    int element_to_node(int i);
    int block_to_node(int);
    RefSCMatrixBlock block_to_block(int);
    void error(const char *msg);
  public:
    DistDiagSCMatrix(DistSCDimension*);
    ~DistDiagSCMatrix();

    // implementations and overrides of virtual functions
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

    RefMessageGrp messagegrp() { return d->messagegrp(); }

    RefSCMatrixSubblockIter local_blocks(SCMatrixSubblockIter::Access);
    RefSCMatrixSubblockIter all_blocks(SCMatrixSubblockIter::Access);
};

class DistSCMatrixListSubblockIter: public SCMatrixListSubblockIter {
  protected:
    RefMessageGrp grp_;
    StateSend out_;
    StateRecv in_;
    int step_;
    RefSCMatrixBlockList locallist_;

    void maybe_advance_list();
    void advance_list();
  public:
    DistSCMatrixListSubblockIter(Access,
                                 const RefSCMatrixBlockList &locallist,
                                 const RefMessageGrp &grp);
    void begin();
    void next();
    ~DistSCMatrixListSubblockIter();
};

#endif
