
// these provide integrals using Justin Fermann's cints routines

#ifndef _chemistry_qc_cints_int2jf_h
#define _chemistry_qc_cints_int2jf_h

#include <chemistry/qc/basis/basis.h>

struct iclass {
    int am[5];       // the 5th one is for m
    int m;           // so I screwed up...
    int operands[5];
    double *Val;
    int type;        // 0 for HRR build, n for number of operands in VRR build
    int done;
};

struct base_eri {
    int i;
    int j;
    int k;
    int l;
    double val;
};

struct shell {
    int shl;
    int atom;
    int nc;
    double Z;
    GaussianShell *gs;
    Point& p;

    shell(int a, int s, int n, GaussianBasisSet& gbs) :
      p(gbs.molecule()->atom(a).point())
    {
      atom=a;
      shl=s;
      nc=n;
      gs = &gbs(s);
      Z = gbs.molecule()->atom(a).element().charge();
    }

    int am() const { return gs->am(nc); }
    double exp(int prim) const { return gs->exponent(prim); }
    double coef(int prim) const { return gs->coefficient_unnorm(nc,prim); }
    double nprim() const { return gs->nprimitive(); }
};

class FJTable {
  private:
    int maxj;
    int itable_infinity;
    double wval_infinity;

    double **gtable;
    double *denomarray;
    double *int_fjttable;
    
    void done();
    
  public:
    FJTable();
    FJTable(int);
    ~FJTable();

    void reset(int);
    void fjt(int, double);
};

class TwoBodyIntJF
{
  protected:
    RefGaussianBasisSet gbs_;
    
    iclass *Classes;
    base_eri *tot_data;
    
    char *V_done;
    char *H_done;
    double *dp_use;
    
    FJTable fjt;
    
    void init();
    void List_HRR(iclass*, int[4], int&, double*&);
    void Top_VRR(int, iclass*, iclass*, int&, double*&);
    void Init_VRR(iclass*, iclass*, int);

  public:
    TwoBodyIntJF(const RefGaussianBasisSet&);
    ~TwoBodyIntJF();

    void compute_shell(int,int,int,int,double*);
};

#endif
