
#ifndef _chemistry_qc_wfn_effh_h
#define _chemistry_qc_wfn_effh_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/wfn/accum.h>

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
#   define CLASSNAME GSGeneralEffH
#   define HAVE_CTOR
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    void init();
    
  public:
    GSGeneralEffH();
    GSGeneralEffH(StateIn&);
    GSGeneralEffH(const RefKeyVal&);
    ~GSGeneralEffH();

    void save_data_state(StateOut&);
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
#   define CLASSNAME GSHighSpinEffH
#   define HAVE_CTOR
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    void init();

  public:
    GSHighSpinEffH();
    GSHighSpinEffH(StateIn&);
    GSHighSpinEffH(const RefKeyVal&);
    ~GSHighSpinEffH();

    void save_data_state(StateOut&);
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
#   define CLASSNAME TestEffH
#   define HAVE_CTOR
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    void init();

  public:
    TestEffH();
    TestEffH(StateIn&);
    TestEffH(const RefKeyVal&);
    ~TestEffH();

    void save_data_state(StateOut&);
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
#   define CLASSNAME PsiEffH
#   define HAVE_CTOR
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    void init();

  public:
    PsiEffH();
    PsiEffH(StateIn&);
    PsiEffH(const RefKeyVal&);
    ~PsiEffH();

    void save_data_state(StateOut&);
};

#endif
