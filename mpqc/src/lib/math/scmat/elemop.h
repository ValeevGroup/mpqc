
#ifndef _math_scmat_elemop_h
#define _math_scmat_elemop_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>

class SCMatrixBlock;
class SCMatrixBlockIter;
class SCMatrixRectBlock;
class SCMatrixLTriBlock;
class SCMatrixDiagBlock;
class SCVectorSimpleBlock;

class SSRefSCElementOp;
typedef class SSRefSCElementOp RefSCElementOp;

class SSRefSCElementOp2;
typedef class SSRefSCElementOp2 RefSCElementOp2;

class SSRefSCElementOp3;
typedef class SSRefSCElementOp3 RefSCElementOp3;

class SCElementOp: public SavableState {
#   define CLASSNAME SCElementOp
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCElementOp();
    SCElementOp(StateIn&s): SavableState(s) {}
    virtual ~SCElementOp();
    // If duplicates of |SCElementOp| are made then if |has_collect| returns
    // true then collect is called in such a way that each duplicated
    // |SCElementOp| is the argument of |collect| once.  The default
    // return value of |has_collect| is 0 and |collect|'s default action
    // is do nothing.
    virtual int has_collect();
    virtual int has_side_effects();
    virtual void collect(RefSCElementOp&);
    virtual void process(SCMatrixBlockIter&) = 0;
    void process(SCMatrixBlock*);
    virtual void process(SCMatrixRectBlock*);
    virtual void process(SCMatrixLTriBlock*);
    virtual void process(SCMatrixDiagBlock*);
    virtual void process(SCVectorSimpleBlock*);
};
DCRef_declare(SCElementOp);
SSRef_declare(SCElementOp);

class SCElementOp2: public SavableState {
#   define CLASSNAME SCElementOp2
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCElementOp2();
    SCElementOp2(StateIn&s): SavableState(s) {}
    virtual ~SCElementOp2();
    // If duplicates of |SCElementOp2| are made then if |has_collect| returns
    // true then collect is called in such a way that each duplicated
    // |SCElementOp| is the argument of |collect| once.  The default
    // return value of |has_collect| is 0 and |collect|'s default action
    // is do nothing.
    virtual int has_collect();
    virtual int has_side_effects();
    virtual int has_side_effects_in_arg();
    virtual void collect(RefSCElementOp2&);
    virtual void process(SCMatrixBlockIter&,SCMatrixBlockIter&) = 0;
    void process(SCMatrixBlock*,SCMatrixBlock*);
    virtual void process(SCMatrixRectBlock*,SCMatrixRectBlock*);
    virtual void process(SCMatrixLTriBlock*,SCMatrixLTriBlock*);
    virtual void process(SCMatrixDiagBlock*,SCMatrixDiagBlock*);
    virtual void process(SCVectorSimpleBlock*,SCVectorSimpleBlock*);
};
DCRef_declare(SCElementOp2);
SSRef_declare(SCElementOp2);

class SCElementOp3: public SavableState {
#   define CLASSNAME SCElementOp3
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCElementOp3();
    SCElementOp3(StateIn&s): SavableState(s) {}
    virtual ~SCElementOp3();
    // If duplicates of |SCElementOp3| are made then if |has_collect| returns
    // true then collect is called in such a way that each duplicated
    // |SCElementOp| is the argument of |collect| once.  The default
    // return value of |has_collect| is 0 and |collect|'s default action
    // is do nothing.
    virtual int has_collect();
    virtual int has_side_effects();
    virtual int has_side_effects_in_arg1();
    virtual int has_side_effects_in_arg2();
    virtual void collect(RefSCElementOp3&);
    virtual void process(SCMatrixBlockIter&,
                         SCMatrixBlockIter&,
                         SCMatrixBlockIter&) = 0;
    void process(SCMatrixBlock*,SCMatrixBlock*,SCMatrixBlock*);
    virtual void process(SCMatrixRectBlock*,
                         SCMatrixRectBlock*,
                         SCMatrixRectBlock*);
    virtual void process(SCMatrixLTriBlock*,
                         SCMatrixLTriBlock*,
                         SCMatrixLTriBlock*);
    virtual void process(SCMatrixDiagBlock*,
                         SCMatrixDiagBlock*,
                         SCMatrixDiagBlock*);
    virtual void process(SCVectorSimpleBlock*,
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
    double product;
  public:
    SCElementScalarProduct();
    SCElementScalarProduct(StateIn&);
    ~SCElementScalarProduct();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&,SCMatrixBlockIter&);
    int has_collect();
    void collect(RefSCElementOp&);
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
  public:
    SCElementInvert();
    SCElementInvert(double a);
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
    double r;
  public:
    SCElementMaxAbs();
    SCElementMaxAbs(StateIn&);
    ~SCElementMaxAbs();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
    int has_collect();
    void collect(RefSCElementOp&);
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

#endif
