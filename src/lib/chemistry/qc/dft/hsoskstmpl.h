
namespace sc {

class LocalHSOSKSContribution {
  private:
    double * const gmat;
    double * const gmato;
    double * const pmat;
    double * const pmato;
    double a0;

  public:
    LocalHSOSKSContribution(double *g, double *p, double *go, double *po,
                          double _a0) :
      gmat(g), gmato(go), pmat(p), pmato(po), a0(_a0) {}
    ~LocalHSOSKSContribution() {}

    void set_bound(double, double) {}

    inline void cont1(int ij, int kl, double val) {
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
    
    inline void cont2(int ij, int kl, double val) {
      val *= 0.25*a0;
      gmat[ij] -= val*pmat[kl];
      gmat[kl] -= val*pmat[ij];

      gmato[ij] += val*pmato[kl];
      gmato[kl] += val*pmato[ij];
    }
    
    inline void cont3(int ij, int kl, double val) {
      val *= 0.5*a0;
      gmat[ij] -= val*pmat[kl];
      gmat[kl] -= val*pmat[ij];

      gmato[ij] += val*pmato[kl];
      gmato[kl] += val*pmato[ij];
    }
    
    inline void cont4(int ij, int kl, double val) {
      gmat[ij] += (1.0 - 0.25*a0)*val*pmat[kl];
      gmat[kl] += (1.0 - 0.25*a0)*val*pmat[ij];

      gmato[ij] += a0*0.25*val*pmato[kl];
      gmato[kl] += a0*0.25*val*pmato[ij];
    }
    
    inline void cont5(int ij, int kl, double val) {
      gmat[ij] += (1.0 - 0.5*a0)*val*pmat[kl];
      gmat[kl] += (1.0 - 0.5*a0)*val*pmat[ij];

      gmato[ij] += 0.5*a0*val*pmato[kl];
      gmato[kl] += 0.5*a0*val*pmato[ij];
    }
};

class LocalHSOSKSEnergyContribution {
  private:
    double * const pmat;
    double * const pmato;
    double a0;

  public:
    double ec;
    double ex;
    
    LocalHSOSKSEnergyContribution(double *p, double *po,
                                double _a0) : pmat(p), pmato(po), a0(_a0) {
      ec=ex=0;
    }

    ~LocalHSOSKSEnergyContribution() {}

    void set_bound(double, double) {}

    inline void cont1(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
    }
    
    inline void cont2(int ij, int kl, double val) {
      ex -= a0*0.25*val*(pmat[ij]*pmat[kl] + pmato[ij]*pmato[kl]);
    }
    
    inline void cont3(int ij, int kl, double val) {
      ex -= a0*0.5*val*(pmat[ij]*pmat[kl] + pmato[ij]*pmato[kl]);
    }
    
    inline void cont4(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
      ex -= a0*0.25*val*(pmat[ij]*pmat[kl] + pmato[ij]*pmato[kl]);
    }
    
    inline void cont5(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
      ex -= a0*0.5*val*(pmat[ij]*pmat[kl] + pmato[ij]*pmato[kl]);
    }
};

}
