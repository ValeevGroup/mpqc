
#ifndef _math_scmat_elemop_h
#define _math_scmat_elemop_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>
#include <util/group/message.h>

class SCMatrixBlock;
class SCMatrixBlockIter;
class SCMatrixRectBlock;
class SCMatrixLTriBlock;
class SCMatrixDiagBlock;
class SCVectorSimpleBlock;
class SCMatrixRectSubBlock;
class SCMatrixLTriSubBlock;
class SCMatrixDiagSubBlock;
class SCVectorSimpleSubBlock;

class SCMatrix;
class SymmSCMatrix;
class DiagSCMatrix;
class SCVector;

class SSRefSCElementOp;
typedef class SSRefSCElementOp RefSCElementOp;

class SSRefSCElementOp2;
typedef class SSRefSCElementOp2 RefSCElementOp2;

class SSRefSCElementOp3;
typedef class SSRefSCElementOp3 RefSCElementOp3;

//. Objects of class \clsnm{SCElementOp} are used to perform operations on
//the elements of matrices.  When the \clsnm{SCElementOp} object is given
//to the \srccd{element\_op} member of a matrix, each block the matrix is
//passed to one of the \srccd{process}, \srccd{process\_base}, or
//\srccd{process\_base} members.
class SCElementOp: public SavableState {
#   define CLASSNAME SCElementOp
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCElementOp();
    SCElementOp(StateIn&s): SavableState(s) {}
    virtual ~SCElementOp();
    //. If duplicates of the \clsnm{SCElementOp} exist (that is, there is
    //more than one node), then if \srccd{has\_collect} returns nonzero
    //then collect is called with a \clsnmref{MessageGrp} reference after
    //all of the blocks have been processed.  The default return value of
    //\srccd{has\_collect} is 0 and \srccd{collect}'s default action is do
    //nothing.  If \srccd{defer\_collect} member is called with nonzero,
    //\srccd{collect} will do nothing (this is only used by the blocked
    //matrices).
    virtual int has_collect();
    virtual void defer_collect(int);
    virtual void collect(const RefMessageGrp&);
    //. By default this returns nonzero.  If the \clsnm{ElementOp}
    //specialization will change any elements of the matrix, then
    //this must be overridden to return nonzero.
    virtual int has_side_effects();

    //. This is the fallback routine to process blocks and is called
    //by \srccd{process\_spec} members that are not overridden.
    virtual void process(SCMatrixBlockIter&) = 0;

    //. Lazy matrix implementors can call this member when the
    //type of block specialization is unknown.  However, this
    //will attempt to castdown \vrbl{block} to a block specialization
    //and will thus be less efficient.
    void process_base(SCMatrixBlock*block);

    //. Matrices should call these members when the type of block is known.
    //\clsnm{ElementOp} specializations should override these when
    //efficiency is important, since these give the most efficient access
    //to the elements of the block.
    virtual void process_spec(SCMatrixRectBlock*);
    virtual void process_spec(SCMatrixLTriBlock*);
    virtual void process_spec(SCMatrixDiagBlock*);
    virtual void process_spec(SCVectorSimpleBlock*);
    virtual void process_spec(SCMatrixRectSubBlock*);
    virtual void process_spec(SCMatrixLTriSubBlock*);
    virtual void process_spec(SCMatrixDiagSubBlock*);
    virtual void process_spec(SCVectorSimpleSubBlock*);
};
DCRef_declare(SCElementOp);
SSRef_declare(SCElementOp);

//. The \clsnm{SCElementOp2} class is very similar to the
//\clsnmref{SCElementOp} class except that pairs of blocks
//are treated simultaneously.  The two matrices involved must
//have identical storage layout, which will be the case if
//both matrices are of the same type and dimensions.
class SCElementOp2: public SavableState {
#   define CLASSNAME SCElementOp2
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCElementOp2();
    SCElementOp2(StateIn&s): SavableState(s) {}
    virtual ~SCElementOp2();
    virtual int has_collect();
    virtual void defer_collect(int);
    virtual int has_side_effects();
    virtual int has_side_effects_in_arg();
    virtual void collect(const RefMessageGrp&);
    virtual void process(SCMatrixBlockIter&,SCMatrixBlockIter&) = 0;
    void process_base(SCMatrixBlock*,SCMatrixBlock*);
    virtual void process_spec(SCMatrixRectBlock*,SCMatrixRectBlock*);
    virtual void process_spec(SCMatrixLTriBlock*,SCMatrixLTriBlock*);
    virtual void process_spec(SCMatrixDiagBlock*,SCMatrixDiagBlock*);
    virtual void process_spec(SCVectorSimpleBlock*,SCVectorSimpleBlock*);
};
DCRef_declare(SCElementOp2);
SSRef_declare(SCElementOp2);

//. The \clsnm{SCElementOp3} class is very similar to the
//\clsnmref{SCElementOp} class except that a triplet of blocks
//is treated simultaneously.  The three matrices involved must
//have identical storage layout, which will be the case if
//all matrices are of the same type and dimensions.
class SCElementOp3: public SavableState {
#   define CLASSNAME SCElementOp3
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCElementOp3();
    SCElementOp3(StateIn&s): SavableState(s) {}
    virtual ~SCElementOp3();
    virtual int has_collect();
    virtual void defer_collect(int);
    virtual int has_side_effects();
    virtual int has_side_effects_in_arg1();
    virtual int has_side_effects_in_arg2();
    virtual void collect(const RefMessageGrp&);
    virtual void process(SCMatrixBlockIter&,
                         SCMatrixBlockIter&,
                         SCMatrixBlockIter&) = 0;
    void process_base(SCMatrixBlock*,SCMatrixBlock*,SCMatrixBlock*);
    virtual void process_spec(SCMatrixRectBlock*,
                              SCMatrixRectBlock*,
                              SCMatrixRectBlock*);
    virtual void process_spec(SCMatrixLTriBlock*,
                              SCMatrixLTriBlock*,
                              SCMatrixLTriBlock*);
    virtual void process_spec(SCMatrixDiagBlock*,
                              SCMatrixDiagBlock*,
                              SCMatrixDiagBlock*);
    virtual void process_spec(SCVectorSimpleBlock*,
                              SCVectorSimpleBlock*,
                              SCVectorSimpleBlock*);
};
DCRef_declare(SCElementOp3);
SSRef_declare(SCElementOp3);

class SCElementScalarProduct: public SCElementOp2 {
#   define CLASSNAME SCElementScalarProduct
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int deferred_;
    double product;
  public:
    SCElementScalarProduct();
    SCElementScalarProduct(StateIn&);
    ~SCElementScalarProduct();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&,SCMatrixBlockIter&);
    int has_collect();
    void defer_collect(int);
    void collect(const RefMessageGrp&);
    double result();
    void init() { product = 0.0; }
};
SavableState_REF_dec(SCElementScalarProduct);

class SCDestructiveElementProduct: public SCElementOp2 {
#   define CLASSNAME SCDestructiveElementProduct
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    SCDestructiveElementProduct();
    SCDestructiveElementProduct(StateIn&);
    ~SCDestructiveElementProduct();
    int has_side_effects();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&,SCMatrixBlockIter&);
};

class SCElementScale: public SCElementOp {
#   define CLASSNAME SCElementScale
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    double scale;
  public:
    SCElementScale(double a);
    SCElementScale(StateIn&);
    ~SCElementScale();
    int has_side_effects();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementRandomize: public SCElementOp {
#   define CLASSNAME SCElementRandomize
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    double assign;
  public:
    SCElementRandomize();
    SCElementRandomize(StateIn&);
    ~SCElementRandomize();
    int has_side_effects();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementAssign: public SCElementOp {
#   define CLASSNAME SCElementAssign
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    double assign;
  public:
    SCElementAssign(double a);
    SCElementAssign(StateIn&);
    ~SCElementAssign();
    int has_side_effects();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementSquareRoot: public SCElementOp {
#   define CLASSNAME SCElementSquareRoot
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    SCElementSquareRoot();
    SCElementSquareRoot(double a);
    SCElementSquareRoot(StateIn&);
    ~SCElementSquareRoot();
    int has_side_effects();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementInvert: public SCElementOp {
#   define CLASSNAME SCElementInvert
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    double threshold_;
  public:
    SCElementInvert(double threshold = 0.0);
    SCElementInvert(StateIn&);
    ~SCElementInvert();
    int has_side_effects();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementScaleDiagonal: public SCElementOp {
#   define CLASSNAME SCElementScaleDiagonal
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    double scale_diagonal;
  public:
    SCElementScaleDiagonal(double a);
    SCElementScaleDiagonal(StateIn&);
    ~SCElementScaleDiagonal();
    int has_side_effects();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementShiftDiagonal: public SCElementOp {
#   define CLASSNAME SCElementShiftDiagonal
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    double shift_diagonal;
  public:
    SCElementShiftDiagonal(double a);
    SCElementShiftDiagonal(StateIn&);
    ~SCElementShiftDiagonal();
    int has_side_effects();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
};

class SCElementMaxAbs: public SCElementOp {
#   define CLASSNAME SCElementMaxAbs
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int deferred_;
    double r;
  public:
    SCElementMaxAbs();
    SCElementMaxAbs(StateIn&);
    ~SCElementMaxAbs();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
    int has_collect();
    void defer_collect(int);
    void collect(const RefMessageGrp&);
    double result();
};
SavableState_REF_dec(SCElementMaxAbs);

class SCElementDot: public SCElementOp {
#   define CLASSNAME SCElementDot
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    double** avects;
    double** bvects;
    int length;
  public:
    SCElementDot(StateIn&);
    void save_data_state(StateOut&);
    SCElementDot(double**a, double**b, int length);
    void process(SCMatrixBlockIter&);
    int has_side_effects();
};

class SCElementAccumulateSCMatrix: public SCElementOp {
#   define CLASSNAME SCElementAccumulateSCMatrix
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    SCMatrix *m;
  public:
    SCElementAccumulateSCMatrix(SCMatrix *);
    int has_side_effects();
    void process(SCMatrixBlockIter&);
};

class SCElementAccumulateSymmSCMatrix: public SCElementOp {
#   define CLASSNAME SCElementAccumulateSymmSCMatrix
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    SymmSCMatrix *m;
  public:
    SCElementAccumulateSymmSCMatrix(SymmSCMatrix *);
    int has_side_effects();
    void process(SCMatrixBlockIter&);
};

class SCElementAccumulateDiagSCMatrix: public SCElementOp {
#   define CLASSNAME SCElementAccumulateDiagSCMatrix
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    DiagSCMatrix *m;
  public:
    SCElementAccumulateDiagSCMatrix(DiagSCMatrix *);
    int has_side_effects();
    void process(SCMatrixBlockIter&);
};

class SCElementAccumulateSCVector: public SCElementOp {
#   define CLASSNAME SCElementAccumulateSCVector
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    SCVector *m;
  public:
    SCElementAccumulateSCVector(SCVector *);
    int has_side_effects();
    void process(SCMatrixBlockIter&);
};

#endif
