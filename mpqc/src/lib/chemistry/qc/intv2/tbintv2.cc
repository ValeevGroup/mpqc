
#include <stdio.h>

#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/intv2/integralv2.h>

TwoBodyIntV2::TwoBodyIntV2(const RefGaussianBasisSet& b) :
  TwoBodyInt(b),
  same_center(1)
{
  c1 = c2 = c3 = c4 = int_centers_from_gbs(b);
  if (!c1) {
    fprintf(stderr,"TwoBodyIntV2::could not form centers\n");
    abort();
  }
  
  init();
}

TwoBodyIntV2::TwoBodyIntV2(const RefGaussianBasisSet& b1,
                           const RefGaussianBasisSet& b2,
                           const RefGaussianBasisSet& b3,
                           const RefGaussianBasisSet& b4) :
  TwoBodyInt(b1,b2,b3,b4),
  same_center(0)
{
  c1 = int_centers_from_gbs(b1);
  c2 = int_centers_from_gbs(b2);
  c3 = int_centers_from_gbs(b3);
  c4 = int_centers_from_gbs(b4);

  if (!c1 || !c2 || !c3 || !c4) {
    fprintf(stderr,"TwoBodyIntV2::could not form centers\n");
    abort();
  }

  init();
}

TwoBodyIntV2::~TwoBodyIntV2()
{
  int_done_bounds();
  int_done_erep();
  int_done_offsets2(c1,c2,c3,c4);
  int_done_storage();

  if (same_center) {
    free_centers(c1);
    free(c1);
    c1=c2=c3=c4=0;
  } else {
    free_centers(c1);
    free(c1);
    free_centers(c2);
    free(c2);
    free_centers(c3);
    free(c3);
    free_centers(c4);
    free(c4);
  }
}

void
TwoBodyIntV2::init()
{
  int_initialize_offsets2(c1,c2,c3,c4);

  //int flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2;
  int flags = INT_EREP|INT_NOSTRB;

  //if (!store1_)
    //flags |= INT_NOSTR1;

  //if (!store2_)
    //flags |= INT_NOSTR2;

  intbuf = int_initialize_erep(flags,0,c1,c2,c3,c4);

  int_storage(1000000);
  int_init_bounds();
}

void
TwoBodyIntV2::compute_shell(int is, int js, int ks, int ls)
{
  int_erep(INT_EREP|INT_NOBCHK|INT_NOPERM,&is,&js,&ks,&ls);

  int ni = INT_SH_NFUNC((c1),is);
  int nj = INT_SH_NFUNC((c2),js);
  int nk = INT_SH_NFUNC((c3),ks);
  int nl = INT_SH_NFUNC((c4),ls);

  memcpy(buffer_, intbuf, sizeof(double)*ni*nj*nk*nl);
}
