
// these provide integrals using Justin Fermann's cints routines

#ifndef _chemistry_qc_cints_int2jf_h
#define _chemistry_qc_cints_int2jf_h

#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/cints/cints.h>

/////////////////////////////////////////////////////////////////////////////

struct shell_stuff {
    int c;
    int am;
    int center;
    int func0;
    GaussianShell& gs;
    AtomicCenter& ac;

    shell_stuff(GaussianBasisSet& gbs, int shell, int nc) :
      c(nc),
      gs(gbs[shell]),
      ac(gbs.molecule()->atom(gbs.shell_to_center(shell)))
    {
      center = gbs.shell_to_center(shell);
      am = gs.am(nc);
      func0 = gbs.shell_to_function(shell);
      for (int i=0; i < nc; i++)
        func0 += gs.nfunction(i);
    }

    double coef(int prim) const { return gs.coefficient_unnorm(c,prim); }
};
      
struct iclass{
  int am[5]; /* the 5th one is for m */
  int m; /* so I screwed up... */
  int operands[5];
  double *Val;
  int type; /* 0 for HRR build, n for number of operands in VRR build */
  int done;
};

struct am_str{
   int index[4][3];
   int m;  /* auxiliary index */
};

struct shell_pair{
  int i, j;
  double P[3];
  double AB[3];
  double PA[3];
  double PB[3];
  double a1, a2, gamma;
  double inorm, jnorm;
  double Sovlp;
  double Smax;
};

/////////////////////////////////////////////////////////////////////////////

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

    double table_value(int i) const { return int_fjttable[i]; }
};

/////////////////////////////////////////////////////////////////////////////

class TwoBodyIntJF
{
  private:
    RefGaussianBasisSet gbs_;
    
    iclass *Classes;
    char *V_done;
    char *H_done;

    double *DP;
    double U[6][3];
    double F[MAXAM*4];

    FJTable fjt;
    
  private:
    double inorm;
    double jnorm;
    double knorm;
    double lnorm;
    
    double expi;
    double expj;
    double expk;
    double expl;

    double zeta;
    double eta;
    double rho;

    double oo2z;
    double oo2n;
    double oo2zn;
    double poz;
    double pon;

  private:
    int List_HRR(iclass*, int[5], int&, double*&);
    int List_VRR(int[5], iclass*, int&, double*&);

    void Init_VRR(iclass*, iclass*, int);
    void Top_VRR(int, iclass*, iclass*, int&, double*&);

    double * HRR_build(iclass *, int, double[3], double[3]);
    double * HRR_build_on1(iclass *, int, double[3], double[3]);
    double * HRR_build_on3(iclass *, int, double[3], double[3]);

    int Fill_data(iclass *, int);

    double *VRR_Build(iclass*, int);
    void Top_VRR_build(iclass*, iclass*, int);

  public:
    TwoBodyIntJF(const RefGaussianBasisSet&);
    ~TwoBodyIntJF();

    void compute_shell(int,int,int,int,double*);

  private:
    
    double *build_p000(iclass *, int);
    double *build_00p0(iclass *, int);
    double *build_d000(iclass *, int);
    double *build_00d0(iclass *, int);
    double *build_p0p0(iclass *, int);
    double *build_f000(iclass *, int);
    double *build_00f0(iclass *, int);
    double *build_d0p0(iclass *, int);
    double *build_p0d0(iclass *, int);
    double *build_f0p0(iclass *, int);
    double *build_p0f0(iclass *, int);
    double *build_d0d0(iclass *, int);
    double *build_f0d0(iclass *, int);
    double *build_d0f0(iclass *, int);
    double *build_f0f0(iclass *, int);
    double *build_g000(iclass *, int);
    double *build_00g0(iclass *, int);
    double *build_g0p0(iclass *, int);
    double *build_p0g0(iclass *, int);
    double *build_d0g0(iclass *, int);
    double *build_g0d0(iclass *, int);
    double *build_g0f0(iclass *, int);
    double *build_f0g0(iclass *, int);
    double *build_g0g0(iclass *, int);

    void top_build_p000(iclass*, iclass*);
    void top_build_00p0(iclass*, iclass*);
    void top_build_d000(iclass*, iclass*);
    void top_build_00d0(iclass*, iclass*);
    void top_build_p0p0(iclass*, iclass*);
    void top_build_f000(iclass*, iclass*);
    void top_build_00f0(iclass*, iclass*);
    void top_build_d0p0(iclass*, iclass*);
    void top_build_p0d0(iclass*, iclass*);
    void top_build_f0p0(iclass*, iclass*);
    void top_build_p0f0(iclass*, iclass*);
    void top_build_d0d0(iclass*, iclass*);
    void top_build_f0d0(iclass*, iclass*);
    void top_build_d0f0(iclass*, iclass*);
    void top_build_f0f0(iclass*, iclass*);
    void top_build_g000(iclass*, iclass*);
    void top_build_00g0(iclass*, iclass*);
    void top_build_g0p0(iclass*, iclass*);
    void top_build_p0g0(iclass*, iclass*);
    void top_build_d0g0(iclass*, iclass*);
    void top_build_g0d0(iclass*, iclass*);
    void top_build_g0f0(iclass*, iclass*);
    void top_build_f0g0(iclass*, iclass*);
    void top_build_g0g0(iclass*, iclass*);
};

#endif
