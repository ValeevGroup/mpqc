
namespace sc {

class LocalHSOSContribution {
  private:
    double * const gmat;
    double * const gmato;
    double * const pmat;
    double * const pmato;

  public:
    LocalHSOSContribution(double *g, double *p, double *go, double *po) :
      gmat(g), gmato(go), pmat(p), pmato(po) {}
    ~LocalHSOSContribution() {}

    void set_bound(double,double) {};

    inline void cont1(int ij, int kl, double val) {
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
    
    inline void cont2(int ij, int kl, double val) {
      val *= 0.25;
      gmat[ij] -= val*pmat[kl];
      gmat[kl] -= val*pmat[ij];

      gmato[ij] += val*pmato[kl];
      gmato[kl] += val*pmato[ij];
    }
    
    inline void cont3(int ij, int kl, double val) {
      val *= 0.5;
      gmat[ij] -= val*pmat[kl];
      gmat[kl] -= val*pmat[ij];

      gmato[ij] += val*pmato[kl];
      gmato[kl] += val*pmato[ij];
    }
    
    inline void cont4(int ij, int kl, double val) {
      gmat[ij] += 0.75*val*pmat[kl];
      gmat[kl] += 0.75*val*pmat[ij];

      gmato[ij] += 0.25*val*pmato[kl];
      gmato[kl] += 0.25*val*pmato[ij];
    }
    
    inline void cont5(int ij, int kl, double val) {
      val *= 0.5;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];

      gmato[ij] += val*pmato[kl];
      gmato[kl] += val*pmato[ij];
    }
};

class LocalHSOSEnergyContribution {
  private:
    double * const pmat;
    double * const pmato;

  public:
    double ec;
    double ex;

    void set_bound(double,double) {};
    
    LocalHSOSEnergyContribution(double *p, double *po) : pmat(p), pmato(po) {
      ec=ex=0;
    }

    ~LocalHSOSEnergyContribution() {}

    inline void cont1(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
    }
    
    inline void cont2(int ij, int kl, double val) {
      ex -= 0.25*val*(pmat[ij]*pmat[kl] + pmato[ij]*pmato[kl]);
    }
    
    inline void cont3(int ij, int kl, double val) {
      ex -= 0.5*val*(pmat[ij]*pmat[kl] + pmato[ij]*pmato[kl]);
    }
    
    inline void cont4(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
      ex -= 0.25*val*(pmat[ij]*pmat[kl] + pmato[ij]*pmato[kl]);
    }
    
    inline void cont5(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
      ex -= 0.5*val*(pmat[ij]*pmat[kl] + pmato[ij]*pmato[kl]);
    }
};
class LocalHSOSGradContribution {
  private:
    double * const pmat;
    double * const pmato;

  public:
    LocalHSOSGradContribution(double *p, double *po) : pmat(p), pmato(po) {}
    ~LocalHSOSGradContribution() {}

    inline double cont1(int ij, int kl) {
      return pmat[ij]*pmat[kl] +
        0.5*(pmato[ij]*pmat[kl] + pmat[ij]*pmato[kl]) +
        0.25*pmato[ij]*pmato[kl];
    }

    inline double cont2(int ij, int kl) {
      return pmat[ij]*pmat[kl] +
        0.5*(pmato[ij]*pmat[kl] + pmat[ij]*pmato[kl] + pmato[ij]*pmato[kl]);
    }
};

} // namespace sc
