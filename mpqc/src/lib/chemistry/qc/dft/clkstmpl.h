
namespace sc {

class LocalCLKSContribution {
  private:
    double * const gmat;
    double * const pmat;
    double a0;

  public:
    LocalCLKSContribution(double *g, double *p, double a) :
      gmat(g), pmat(p), a0(a) {}
    ~LocalCLKSContribution() {}

    void set_bound(double, double) {}

    inline void cont1(int ij, int kl, double val) {
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
    
    inline void cont2(int ij, int kl, double val) {
      val *= -0.25*a0;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
    
    inline void cont3(int ij, int kl, double val) {
      val *= -0.5*a0;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
    
    inline void cont4(int ij, int kl, double val) {
      val *= 1.0 - 0.25*a0;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
    
    inline void cont5(int ij, int kl, double val) {
      val *= 1.0 - 0.5*a0;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
};

class LocalCLKSEnergyContribution {
  private:
    double * const pmat;
    double a0;

  public:
    double ec;
    double ex;
    
    LocalCLKSEnergyContribution(double *p, double a) : pmat(p), a0(a) {
      ec=ex=0;
    }
    ~LocalCLKSEnergyContribution() {}

    void set_bound(double, double) {}

    inline void cont1(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
    }
    
    inline void cont2(int ij, int kl, double val) {
      ex -= a0*0.25*val*pmat[ij]*pmat[kl];
    }
    
    inline void cont3(int ij, int kl, double val) {
      ex -= a0*0.5*val*pmat[ij]*pmat[kl];
    }
    
    inline void cont4(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
      ex -= a0*0.25*val*pmat[ij]*pmat[kl];
    }
    
    inline void cont5(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
      ex -= a0*0.5*val*pmat[ij]*pmat[kl];
    }
};

}
