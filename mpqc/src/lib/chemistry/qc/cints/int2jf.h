
// these provide integrals using Justin Fermann's cints routines

#ifndef _chemistry_qc_cints_int2jf_h
#define _chemistry_qc_cints_int2jf_h

#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/cints.h>

/////////////////////////////////////////////////////////////////////////////

struct iclass{
  int am[5]; /* the 5th one is for m */
  int m; /* so I screwed up... */
  int operands[5];
  double *Val;
  int type; /* 0 for HRR build, n for number of operands in VRR build */
  int done;
};

struct coordinates{
   double x;  /*  what do you think these are? */
   double y;
   double z;
   double Z_nuc; /* nuclear charge */
};

struct am_str{
   int index[4][3];
   int m;  /* auxiliary index */
};

struct base_eri{
   int i, j, k, l;
   double val;
};

struct aux_eri{
   struct am_str L;
   double val;
   int flag;
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
    double table_value(int);
};

/////////////////////////////////////////////////////////////////////////////

class TwoBodyIntJF
{
  protected:
    RefGaussianBasisSet gbs_;
    
    iclass *Classes;
    char *V_done;
    char *H_done;

    base_eri *tot_data;

    double *DP;
    double U[6][3];
    double F[MAXAM*4];

    FJTable fjt;
    
  public:
    TwoBodyIntJF(const RefGaussianBasisSet&);
    ~TwoBodyIntJF();

    void compute_shell(int,int,int,int,double*);
};

#endif
