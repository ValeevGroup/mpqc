
#ifndef _chemistry_qc_wfn_density_h
#define _chemistry_qc_wfn_density_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/isosurf/volume.h>
#include <chemistry/qc/wfn/wfn.h>

class ElectronDensity: public Volume {
  protected:
    Wavefunction& wfn_;
    virtual void compute();
  public:
    ElectronDensity(Wavefunction&);
    ~ElectronDensity();
    virtual void boundingbox(double valuemin,
                             double valuemax,
                             RefSCVector& p1, RefSCVector& p2);
};

#endif
