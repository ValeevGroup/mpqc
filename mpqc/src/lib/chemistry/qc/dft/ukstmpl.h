
class LocalUKSContribution {
  private:
    double * const gmata;
    double * const gmatb;
    double * const pmata;
    double * const pmatb;
    double a0;

  public:
    LocalUKSContribution(double *ga, double *pa, double *gb, double *pb,
                         double a) :
      gmata(ga), pmata(pa), gmatb(gb), pmatb(pb), a0(a) {}
    ~LocalUKSContribution() {}

    inline void cont1(int ij, int kl, double val) {
      gmata[ij] += val*(pmata[kl]+pmatb[kl]);
      gmata[kl] += val*(pmata[ij]+pmatb[ij]);

      gmatb[ij] += val*(pmata[kl]+pmatb[kl]);
      gmatb[kl] += val*(pmata[ij]+pmatb[ij]);
    }
    
    inline void cont2(int ij, int kl, double val) {
      val *= a0*0.5;
      gmata[ij] -= val*pmata[kl];
      gmata[kl] -= val*pmata[ij];

      gmatb[ij] -= val*pmatb[kl];
      gmatb[kl] -= val*pmatb[ij];
    }
    
    inline void cont3(int ij, int kl, double val) {
      val *= a0;
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

class LocalUKSEnergyContribution {
  private:
    double * const pmata;
    double * const pmatb;
    double a0;

  public:
    double ec;
    double ex;
    
    LocalUKSEnergyContribution(double *a, double *b, double an) :
      pmata(a), pmatb(b), a0(an) {
      ec=ex=0;
    }

    ~LocalUKSEnergyContribution() {}

    inline void cont1(int ij, int kl, double val) {
      ec += val*(pmata[ij]+pmatb[ij])*(pmata[kl]+pmatb[kl]);
    }
    
    inline void cont2(int ij, int kl, double val) {
      ex -= a0*0.5*val*(pmata[ij]*pmata[kl]+pmatb[ij]*pmatb[kl]);
    }
    
    inline void cont3(int ij, int kl, double val) {
      ex -= a0*val*(pmata[ij]*pmata[kl]+pmatb[ij]*pmatb[kl]);
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

class LocalUKSGradContribution {
  private:
    double * const pmata;
    double * const pmatb;

  public:
    LocalUKSGradContribution(double *a, double *b) : pmata(a), pmatb(b) {}
    ~LocalUKSGradContribution() {}

    inline double cont1(int ij, int kl) {
      return (pmata[ij]*pmata[kl])+(pmatb[ij]*pmatb[kl]) +
             (pmata[ij]*pmatb[kl])+(pmatb[ij]*pmata[kl]);
    }

    inline double cont2(int ij, int kl) {
      return 2*((pmata[ij]*pmata[kl])+(pmatb[ij]*pmatb[kl]));
    }
};
