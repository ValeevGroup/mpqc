
#include <chemistry/qc/cints/cints.h>
#include <chemistry/qc/cints/int2jf.h>

void
TwoBodyIntJF::init()
{
  fjt.reset(20);

  V_done = new char[MAXCLASS];
  H_done = new char[MAXCLASS];

  int sz=0;
  for (int i=0; i < MAXAM*2; i++)
    for (int j=0; j < MAXAM*2; j++)
      sz += ioff(i+1)*ioff(j+1)*MAXAM*8;
  
  dp_use = new double[sz];
  
}

TwoBodyIntJF::TwoBodyIntJF(const RefGaussianBasisSet& gbs) :
  gbs_(gbs), V_done(0), H_done(0), dp_use(0)
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
}

void
TwoBodyIntJF::compute_shell(int si, int sj, int sk, int sl, double *buf)
{
}


