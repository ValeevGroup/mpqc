
#ifndef _chemistry_qc_scf_effh_h
#define _chemistry_qc_scf_effh_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/scmat/blocked.h>
#include <chemistry/qc/scf/scf.h>

class AccumEffectiveH: public BlockedSCElementOp2 {
  protected:
    SCF *scf_;
    double coef_[18];

    virtual void init() =0;
    
    // hindex is 0 for the closed and 1 for the open shell fock matrix
    // shelli and shellj are 0 for closed, 1 for open, and 2 for virtual
    int index(int hindex, int shelli, int shellj);

    // converts an occupation number to a shell number
    int shell(double);

    double& coef(int i, int j, int k) { return coef_[index(i,j,k)]; }

  public:
    AccumEffectiveH(SCF*);
    virtual ~AccumEffectiveH();

    virtual void process(SCMatrixBlockIter&,SCMatrixBlockIter&);
};

//  Guest & Saunders general form 
//        C        O         V
//    ----------
//    |        |
// C  |   fc   |
//    |        |
//    -------------------
//    |        |        |
// O  | 2fc-fo |   fc   |
//    |        |        |
//    ----------------------------
//    |        |        |        |
// V  |   fc   |   fo   |   fc   |
//    |        |        |        |
//    ----------------------------
class GSGeneralEffH: public AccumEffectiveH {
  protected:
    void init();
    
  public:
    GSGeneralEffH(SCF*);
    ~GSGeneralEffH();
};

//  Guest & Saunders' form for high spin
//        C        O         V
//    ----------
//    |        |
// C  | 2fc-fo |
//    |        |
//    -------------------
//    |        |        |
// O  | 2fc-fo | 2fc-fo |
//    |        |        |
//    ----------------------------
//    |        |        |        |
// V  |   fc   |   fo   | 2fc-fo |
//    |        |        |        |
//    ----------------------------
class GSHighSpinEffH: public AccumEffectiveH {
  protected:
    void init();

  public:
    GSHighSpinEffH(SCF*);
    ~GSHighSpinEffH();
};

//  test form
//        C        O         V
//    ----------
//    |        |
// C  |   fo   |
//    |        |
//    -------------------
//    |        |        |
// O  | 2fc-fo |   fo   |
//    |        |        |
//    ----------------------------
//    |        |        |        |
// V  |   fc   |   fo   |   fo   |
//    |        |        |        |
//    ----------------------------
class TestEffH: public AccumEffectiveH {
  protected:
    void init();

  public:
    TestEffH(SCF*);
    ~TestEffH();
};

//  form for converged wavefunction
//        C        O         V
//    ----------
//    |        |
// C  |   fc   |
//    |        |
//    -------------------
//    |        |        |
// O  | 2fc-fo |   fo   |
//    |        |        |
//    ----------------------------
//    |        |        |        |
// V  |   fc   |   fo   |   fo   |
//    |        |        |        |
//    ----------------------------
class PsiEffH: public AccumEffectiveH {
  protected:
    void init();

  public:
    PsiEffH(SCF*);
    ~PsiEffH();
};

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
