
class LocalCLHFContribution {
  private:
    double * const gmat;
    double * const pmat;

  public:
    LocalCLHFContribution(double *g, double *p) : gmat(g), pmat(p) {}
    ~LocalCLHFContribution() {}

    inline void cont1(int ij, int kl, double val) {
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
    
    inline void cont2(int ij, int kl, double val) {
      val *= -0.25;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
    
    inline void cont3(int ij, int kl, double val) {
      val *= -0.5;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
    
    inline void cont4(int ij, int kl, double val) {
      val *= 0.75;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
    
    inline void cont5(int ij, int kl, double val) {
      val *= 0.5;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
};

class LocalCLHFEnergyContribution {
  private:
    double * const pmat;

  public:
    double ec;
    double ex;
    
    LocalCLHFEnergyContribution(double *p) : pmat(p) {
      ec=ex=0;
    }
    ~LocalCLHFEnergyContribution() {}

    inline void cont1(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
    }
    
    inline void cont2(int ij, int kl, double val) {
      ex -= 0.25*val*pmat[ij]*pmat[kl];
    }
    
    inline void cont3(int ij, int kl, double val) {
      ex -= 0.5*val*pmat[ij]*pmat[kl];
    }
    
    inline void cont4(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
      ex -= 0.25*val*pmat[ij]*pmat[kl];
    }
    
    inline void cont5(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
      ex -= 0.5*val*pmat[ij]*pmat[kl];
    }
};

class LocalCLHFGradContribution {
  private:
    double * const pmat;

  public:
    LocalCLHFGradContribution(double *p) : pmat(p) {}
    ~LocalCLHFGradContribution() {}

    inline double cont1(int ij, int kl) {
      return pmat[ij]*pmat[kl];
    }

    inline double cont2(int ij, int kl) {
      return pmat[ij]*pmat[kl];
    }
};
