
#ifndef _chemistry_qc_wfn_density_h
#define _chemistry_qc_wfn_density_h

#include <math/isosurf/volume.h>
#include <chemistry/qc/wfn/wfn.h>

class ElectronDensity: public Volume {
  protected:
    Wavefunction& _wfn;
    virtual void compute();
  public:
    ElectronDensity(Wavefunction&);
    ~ElectronDensity();
    virtual void boundingbox(double valuemin,
                             double valuemax,
                             Point& p1, Point& p2);
};

#endif
