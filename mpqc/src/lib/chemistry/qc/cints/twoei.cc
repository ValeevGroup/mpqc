
#include <chemistry/qc/cints/cints.h>
#include <chemistry/qc/cints/int2jf.h>

void
TwoBodyIntJF::init()
{
  fjt.reset(20);

  V_done = new char[MAXCLASS];
  memset(V_done,0,sizeof(char)*MAXCLASS);
  
  H_done = new char[MAXCLASS];
  memset(H_done,0,sizeof(char)*MAXCLASS);

  Classes = new iclass[MAXCLASS];
  
  int sz=ioff(MAXAM);
  sz *= sz;
  sz *= sz;
  
  tot_data = new base_eri[sz];
  memset(tot_data,0,sizeof(base_eri)*sz);
  
  sz=0;
  for (int i=0; i < MAXAM*2; i++)
    for (int j=0; j < MAXAM*2; j++)
      sz += ioff(i+1)*ioff(j+1)*MAXAM*8;
  
  dp_use = new double[sz];
  
}

TwoBodyIntJF::TwoBodyIntJF(const RefGaussianBasisSet& gbs) :
  gbs_(gbs),
  V_done(0),
  H_done(0),
  dp_use(0),
  Classes(0),
  tot_data(0)
{
  init();
}

TwoBodyIntJF::~TwoBodyIntJF()
{
  if (V_done) {
    delete[] V_done;
    V_done = 0;
  }
  if (H_done) {
    delete[] H_done;
    H_done = 0;
  }
  if (dp_use) {
    delete[] dp_use;
    dp_use = 0;
  }
  if (Classes) {
    delete[] Classes;
    Classes = 0;
  }
  if (tot_data) {
    delete[] tot_data;
    tot_data = 0;
  }
}

void
TwoBodyIntJF::compute_shell(int sii, int sjj, int skk, int sll, double *buf)
{
  GaussianBasisSet& gbs = *gbs_.pointer();
  
  GaussianShell& gsi = gbs(sii);
  GaussianShell& gsj = gbs(sjj);
  GaussianShell& gsk = gbs(skk);
  GaussianShell& gsl = gbs(sll);

  int ati = gbs.shell_to_center(sii);
  int atj = gbs.shell_to_center(sjj);
  int atk = gbs.shell_to_center(skk);
  int atl = gbs.shell_to_center(sll);
  
  int ati_e_atj = (ati==atj);
  int atj_e_atk = (atj==atk);
  int atk_e_atl = (atk==atl);
  
  for (int ci=0; ci < gsi.ncontraction(); ci++) {
    int ami = gsi.am(ci);
    
    for (int cj=0; cj < gsj.ncontraction(); cj++) {
      int amj = gsj.am(cj);
      int amij = ami+amj;

      for (int ck=0; ck < gsk.ncontraction(); ck++) {
        int amk = gsk.am(ck);
        int amijk = amij + amk;

        for (int cl=0; cl < gsl.ncontraction(); cl++) {
          int aml = gsl.am(cl);
          int amijkl = amijk + aml;

          if (amijkl%2 && ati_e_atj && atj_e_atk && atk_e_atl)
            continue;
          
        }
      }
    }
  }
}


