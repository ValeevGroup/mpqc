
#ifndef _math_scmat_block_h
#define _math_scmat_block_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>

class SCMatrixBlock: virtual public SavableState {
#   define CLASSNAME SCMatrixBlock
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  public:
    SCMatrixBlock();
    virtual ~SCMatrixBlock();
};

SavableState_REF_dec(SCMatrixBlock);

class SCVectorSimpleBlock: public SCMatrixBlock {
#   define CLASSNAME SCVectorSimpleBlock
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    SCVectorSimpleBlock(int istart,int iend);
    SCVectorSimpleBlock(StateIn&);
    virtual ~SCVectorSimpleBlock();
    void save_data_state(StateOut&);
    int istart;
    int iend;
    double* data;
};

SavableState_REF_dec(SCVectorSimpleBlock);

class SCMatrixRectBlock: public SCMatrixBlock {
#   define CLASSNAME SCMatrixRectBlock
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    SCMatrixRectBlock(int is, int ie, int js, int je);
    SCMatrixRectBlock(StateIn&);
    virtual ~SCMatrixRectBlock();
    void save_data_state(StateOut&);
    int istart;
    int jstart;
    int iend;
    int jend;
    double* data;
};
SavableState_REF_dec(SCMatrixRectBlock);

class SCMatrixLTriBlock: public SCMatrixBlock {
#   define CLASSNAME SCMatrixLTriBlock
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    SCMatrixLTriBlock(int s,int e);
    SCMatrixLTriBlock(StateIn&);
    virtual ~SCMatrixLTriBlock();
    void save_data_state(StateOut&);
    int start;
    int end;
    double* data;
};
SavableState_REF_dec(SCMatrixLTriBlock);

class SCMatrixDiagBlock: public SCMatrixBlock {
#   define CLASSNAME SCMatrixDiagBlock
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    SCMatrixDiagBlock(int istart,int iend,int jstart);
    SCMatrixDiagBlock(int istart,int iend);
    SCMatrixDiagBlock(StateIn&);
    virtual ~SCMatrixDiagBlock();
    void save_data_state(StateOut&);
    int istart;
    int jstart;
    int iend;
    double* data;
};
SavableState_REF_dec(SCMatrixDiagBlock);

#endif
