
namespace sc {

class LocalOSSContribution {
  private:
    double * const gmat;
    double * const gmata;
    double * const gmatb;

    double * const pmat;
    double * const pmata;
    double * const pmatb;
  public:
    LocalOSSContribution(double *g, double *p, double *ga, double *pa,
                         double *gb, double *pb) :
      gmat(g), gmata(ga), gmatb(gb), pmat(p), pmata(pa), pmatb(pb) {}
    ~LocalOSSContribution() {}

    void set_bound(double,double) {}

    inline void cont1(int ij, int kl, double val) {
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
    }
    
    inline void cont2(int ij, int kl, double val) {
      val *= 0.25;
      gmat[ij] -= val*pmat[kl];
      gmat[kl] -= val*pmat[ij];

      gmata[ij] += val*pmata[kl];
      gmata[kl] += val*pmata[ij];

      gmatb[ij] += val*pmatb[kl];
      gmatb[kl] += val*pmatb[ij];

      val *= -3.0;
      gmatb[ij] += val*pmata[kl];
      gmatb[kl] += val*pmata[ij];

      gmata[ij] += val*pmatb[kl];
      gmata[kl] += val*pmatb[ij];
    }
    
    inline void cont3(int ij, int kl, double val) {
      val *= 0.5;
      gmat[ij] -= val*pmat[kl];
      gmat[kl] -= val*pmat[ij];

      gmata[ij] += val*pmata[kl];
      gmata[kl] += val*pmata[ij];

      gmatb[ij] += val*pmatb[kl];
      gmatb[kl] += val*pmatb[ij];

      val *= -3.0;
      gmata[ij] += val*pmatb[kl];
      gmata[kl] += val*pmatb[ij];

      gmatb[ij] += val*pmata[kl];
      gmatb[kl] += val*pmata[ij];
    }
    
    inline void cont4(int ij, int kl, double val) {
      gmat[ij] += 0.75*val*pmat[kl];
      gmat[kl] += 0.75*val*pmat[ij];

      gmata[ij] += 0.25*val*pmata[kl];
      gmata[kl] += 0.25*val*pmata[ij];

      gmatb[ij] += 0.25*val*pmatb[kl];
      gmatb[kl] += 0.25*val*pmatb[ij];

      gmata[ij] -= 0.75*val*pmatb[kl];
      gmata[kl] -= 0.75*val*pmatb[ij];

      gmatb[ij] -= 0.75*val*pmata[kl];
      gmatb[kl] -= 0.75*val*pmata[ij];
    }
    
    inline void cont5(int ij, int kl, double val) {
      val *= 0.5;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];

      gmata[ij] += val*pmata[kl];
      gmata[kl] += val*pmata[ij];

      gmatb[ij] += val*pmatb[kl];
      gmatb[kl] += val*pmatb[ij];

      val *= -3.0;
      gmata[ij] += val*pmatb[kl];
      gmata[kl] += val*pmatb[ij];

      gmatb[ij] += val*pmata[kl];
      gmatb[kl] += val*pmata[ij];
    }
};

class LocalOSSEnergyContribution {
  private:
    double * const pmat;
    double * const pmata;
    double * const pmatb;

  public:
    double ec;
    double ex;

    void set_bound(double,double) {};
    
    LocalOSSEnergyContribution(double *p, double *pa, double *pb) :
      pmat(p), pmata(pa), pmatb(pb) {
      ec=ex=0;
    }
    ~LocalOSSEnergyContribution() {}

    inline void cont1(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
    }
    
    inline void cont2(int ij, int kl, double val) {
      ex -= 0.25*val*(pmat[ij]*pmat[kl] + pmata[ij]*pmata[kl] +
                      pmatb[ij]*pmatb[kl]);
      ex += 0.75*val*(pmata[ij]*pmatb[kl] + pmatb[ij]*pmata[kl]);
    }
    
    inline void cont3(int ij, int kl, double val) {
      ex -= 0.5*val*(pmat[ij]*pmat[kl] + pmata[ij]*pmata[kl] +
                      pmatb[ij]*pmatb[kl]);
      ex += 1.5*val*(pmata[ij]*pmatb[kl] + pmatb[ij]*pmata[kl]);
    }
    
    inline void cont4(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
      ex -= 0.25*val*(pmat[ij]*pmat[kl] + pmata[ij]*pmata[kl] +
                      pmatb[ij]*pmatb[kl]);
      ex += 0.75*val*(pmata[ij]*pmatb[kl] + pmatb[ij]*pmata[kl]);
    }
    
    inline void cont5(int ij, int kl, double val) {
      ec += val*pmat[ij]*pmat[kl];
      ex -= 0.5*val*(pmat[ij]*pmat[kl] + pmata[ij]*pmata[kl] +
                      pmatb[ij]*pmatb[kl]);
      ex += 1.5*val*(pmata[ij]*pmatb[kl] + pmatb[ij]*pmata[kl]);
    }
};

class LocalOSSGradContribution {
  private:
    double * const pmat;
    double * const pmata;
    double * const pmatb;

  public:
    LocalOSSGradContribution(double *p, double *pa, double *pb) :
      pmat(p), pmata(pa), pmatb(pb) {}
    ~LocalOSSGradContribution() {}

    inline double cont1(int ij, int kl) {
      return pmat[ij]*pmat[kl] +
        0.5*(pmata[ij]*pmat[kl] + pmat[ij]*pmata[kl] +
             pmatb[ij]*pmat[kl] + pmat[ij]*pmatb[kl]) +
        0.25*(pmata[ij]*pmata[kl] + pmatb[ij]*pmatb[kl] +
              pmata[ij]*pmatb[kl] + pmatb[ij]*pmata[kl]);
    }

    inline double cont2(int ij, int kl) {
      return pmat[ij]*pmat[kl] +
        0.5*(pmata[ij]*pmat[kl] + pmat[ij]*pmata[kl] +
             pmatb[ij]*pmat[kl] + pmat[ij]*pmatb[kl] +
             pmata[ij]*pmata[kl] + pmatb[ij]*pmatb[kl] -
             pmata[ij]*pmatb[kl] - pmatb[ij]*pmata[kl]);
    }
};

}

