
#ifndef _chemistry_qc_scf_clcont_h
#define _chemistry_qc_scf_clcont_h

#ifdef __GNUC__
#pragma interface
#endif

///////////////////////////////////////////////////////////////////////////

class LocalCLContribution {
  private:
    double * const gmat;
    double * const pmat;

  public:
    LocalCLContribution(double *g, double *p) : gmat(g), pmat(p) {}
    ~LocalCLContribution() {}

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

class LocalCLGradContribution {
  private:
    double * const pmat;

  public:
    LocalCLGradContribution(double *p) : pmat(p) {}
    ~LocalCLGradContribution() {}

    inline double cont1(int ij, int kl) {
      return pmat[ij]*pmat[kl];
    }

    inline double cont2(int ij, int kl) {
      return pmat[ij]*pmat[kl];
    }
};

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
