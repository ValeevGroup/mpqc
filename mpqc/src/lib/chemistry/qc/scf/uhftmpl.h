
namespace sc {

class LocalUHFContribution {
  private:
    double * const gmata;
    double * const gmatb;
    double * const pmata;
    double * const pmatb;

  public:
    LocalUHFContribution(double *ga, double *pa, double *gb, double *pb) :
      gmata(ga), gmatb(gb), pmata(pa), pmatb(pb) {}
    ~LocalUHFContribution() {}

    void set_bound(double,double) {};

    inline void cont1(int ij, int kl, double val) {
      gmata[ij] += val*(pmata[kl]+pmatb[kl]);
      gmata[kl] += val*(pmata[ij]+pmatb[ij]);

      gmatb[ij] += val*(pmata[kl]+pmatb[kl]);
      gmatb[kl] += val*(pmata[ij]+pmatb[ij]);
    }
    
    inline void cont2(int ij, int kl, double val) {
      val *= 0.5;
      gmata[ij] -= val*pmata[kl];
      gmata[kl] -= val*pmata[ij];

      gmatb[ij] -= val*pmatb[kl];
      gmatb[kl] -= val*pmatb[ij];
    }
    
    inline void cont3(int ij, int kl, double val) {
      gmata[ij] -= val*pmata[kl];
      gmata[kl] -= val*pmata[ij];

      gmatb[ij] -= val*pmatb[kl];
      gmatb[kl] -= val*pmatb[ij];
    }
    
    inline void cont4(int ij, int kl, double val) {
      cont1(ij,kl,val);
      cont2(ij,kl,val);
    }
    
    inline void cont5(int ij, int kl, double val) {
      cont1(ij,kl,val);
      cont3(ij,kl,val);
    }
};

class LocalUHFEnergyContribution {
  private:
    double * const pmata;
    double * const pmatb;

  public:
    double ec;
    double ex;
    
    LocalUHFEnergyContribution(double *a, double *b) : pmata(a), pmatb(b) {
      ec=ex=0;
    }

    ~LocalUHFEnergyContribution() {}

    void set_bound(double,double) {};

    inline void cont1(int ij, int kl, double val) {
      ec += val*(pmata[ij]+pmatb[ij])*(pmata[kl]+pmatb[kl]);
    }
    
    inline void cont2(int ij, int kl, double val) {
      ex -= 0.5*val*(pmata[ij]*pmata[kl]+pmatb[ij]*pmatb[kl]);
    }
    
    inline void cont3(int ij, int kl, double val) {
      ex -= val*(pmata[ij]*pmata[kl]+pmatb[ij]*pmatb[kl]);
    }
    
    inline void cont4(int ij, int kl, double val) {
      cont1(ij,kl,val);
      cont2(ij,kl,val);
    }
    
    inline void cont5(int ij, int kl, double val) {
      cont1(ij,kl,val);
      cont3(ij,kl,val);
    }
};

class LocalUHFGradContribution {
  private:
    double * const pmata;
    double * const pmatb;

  public:
    LocalUHFGradContribution(double *a, double *b) : pmata(a), pmatb(b) {}
    ~LocalUHFGradContribution() {}

    inline double cont1(int ij, int kl) {
      return (pmata[ij]*pmata[kl])+(pmatb[ij]*pmatb[kl]) +
             (pmata[ij]*pmatb[kl])+(pmatb[ij]*pmata[kl]);
    }

    inline double cont2(int ij, int kl) {
      return 2*((pmata[ij]*pmata[kl])+(pmatb[ij]*pmatb[kl]));
    }
};

}
