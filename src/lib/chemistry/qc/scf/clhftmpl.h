
#include <mpqc_config.h>
#include <util/misc/exenv.h>

#undef SCF_CHECK_BOUNDS

#ifdef SCF_CHECK_BOUNDS
#define CHECK(ival,pval,ij,kl,con) check(ival,pval,ij,kl,con)
#else
#define CHECK(ival,pval,ij,kl,con)
#endif

namespace sc {

class LocalCLHFContribution {
  private:
    double * const RESTRICT gmat;
    double * const RESTRICT pmat;

    double ibound_;
    double pbound_;

  public:
    LocalCLHFContribution(double *g, double *p) : gmat(g), pmat(p) {}
    ~LocalCLHFContribution() {}

    void set_bound(double i, double p) { ibound_ = i; pbound_ = p; }
    void check(double ival, double pval, int ij, int kl, const char *contrib)
        {
          int bad = 0;
          if ( 1.000001 * ibound_ < (ival > 0 ? ival : -ival)) {
              ExEnv::errn() << "BAD INTEGRAL BOUND" << std::endl;
              ExEnv::errn() << " bound = " << ibound_ << std::endl;
              ExEnv::errn() << " value = " << ival << std::endl;
              bad = 1;
            }
          if ( 1.000001 * pbound_ < (pval > 0 ? pval : -pval)) {
              ExEnv::errn() << "BAD DENSITY BOUND" << std::endl;
              ExEnv::errn() << " bound = " << pbound_ << std::endl;
              ExEnv::errn() << " value = " << pval << std::endl;
              bad = 1;
            }
          if (bad) {
              ExEnv::errn() << " ij    = " << ij << std::endl;
              ExEnv::errn() << " kl    = " << kl << std::endl;
              ExEnv::errn() << " cont  = " << contrib << std::endl;
              abort();
            }
        }

    inline void cont1(int ij, int kl, double val) {
      gmat[ij] += val*pmat[kl]; CHECK(val,pmat[kl],ij,kl,"cont1a");
      gmat[kl] += val*pmat[ij]; CHECK(val,pmat[ij],ij,kl,"cont1b");
    }
    
    inline void cont2(int ij, int kl, double val) {
      val *= -0.25;
      gmat[ij] += val*pmat[kl]; CHECK(4*val,0.25*pmat[kl],ij,kl,"cont2a");
      gmat[kl] += val*pmat[ij]; CHECK(4*val,0.25*pmat[ij],ij,kl,"cont2b");
    }
    
    inline void cont3(int ij, int kl, double val) {
      val *= -0.5;
      gmat[ij] += val*pmat[kl]; CHECK(2*val,0.5*pmat[kl],ij,kl,"cont3a");
      gmat[kl] += val*pmat[ij]; CHECK(2*val,0.5*pmat[ij],ij,kl,"cont3b");
    }
    
    inline void cont4(int ij, int kl, double val) {
      val *= 0.75;
      gmat[ij] += val*pmat[kl]; CHECK(4./3.*val,0.75*pmat[kl],ij,kl,"cont4a");
      gmat[kl] += val*pmat[ij]; CHECK(4./3.*val,0.75*pmat[ij],ij,kl,"cont4b");
    }
    
    inline void cont5(int ij, int kl, double val) {
      val *= 0.5;
      gmat[ij] += val*pmat[kl]; CHECK(2*val,0.5*pmat[kl],ij,kl,"cont5a");
      gmat[kl] += val*pmat[ij]; CHECK(2*val,0.5*pmat[ij],ij,kl,"cont5b");
    }
};

class LocalCLHFEnergyContribution {
  private:
    double * const pmat;

  public:
    double ec;
    double ex;

    void set_bound(double,double) {}
    
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

}
