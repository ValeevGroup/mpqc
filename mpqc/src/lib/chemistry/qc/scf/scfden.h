

#ifndef _chemistry_qc_scf_scfden_h
#define _chemistry_qc_scf_scfden_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/scmat/elemop.h>
#include <math/scmat/blocked.h>

#include <chemistry/qc/scf/scf.h>

class SCFDensity : public BlockedSCElementOp {
  private:
    SCF *scf_;
    RefSCMatrix vec;
    double occ;

  public:
    SCFDensity(SCF*, const RefSCMatrix&, double);
    ~SCFDensity();

    int has_side_effects();
    void set_occ(double);

    void process(SCMatrixBlockIter& bi);
};

class SCFEnergy : public SCElementOp2 {
  private:
    double eelec;
    int deferred_;
    
  public:
    SCFEnergy();
    ~SCFEnergy();

    int has_collect();
    void defer_collect(int h);
    void collect(const RefMessageGrp&grp);
    double result();
    void reset();

    void process(SCMatrixBlockIter&i, SCMatrixBlockIter&j);
};

class LevelShift : public BlockedSCElementOp {
  private:
    SCF *scf_;
    double shift;

  public:
    LevelShift(SCF*);
    ~LevelShift();

    int has_side_effects();
    void set_shift(double);
    
    void process(SCMatrixBlockIter&);
};

// MO lagrangian
//       c  o  v
//  c  |FC|FC| 0|
//     ----------
//  o  |FC|FO| 0|
//     ----------
//  v  | 0| 0| 0|
//
class MOLagrangian : public BlockedSCElementOp2 {
  private:
    SCF *scf_;

  public:
    MOLagrangian(SCF* s);
    ~MOLagrangian();

    int has_side_effects();

    void process(SCMatrixBlockIter& bi1, SCMatrixBlockIter& bi2);
};

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
