
#include <stdio.h>

extern "C" {
#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv2/atomsfree.h>
#include <chemistry/qc/intv2/int_flags.h>
#include <chemistry/qc/intv2/int_types.h>
#include <chemistry/qc/intv2/int_macros.h>

#include <chemistry/qc/intv2/bounds.gbl>
#include <chemistry/qc/intv2/init2e.gbl>
#include <chemistry/qc/intv2/offsets.gbl>
#include <chemistry/qc/intv2/comp_erep.gbl>
}

#include <chemistry/qc/intv2/int_cplus.h>
#include <chemistry/qc/intv2/tbintv2.h>
#include <chemistry/qc/intv2/storage.h>

TwoBodyIntV2::TwoBodyIntV2(const RefGaussianBasisSet& b1,
                           const RefGaussianBasisSet& b2,
                           const RefGaussianBasisSet& b3,
                           const RefGaussianBasisSet& b4) :
  TwoBodyInt(b1,b2,b3,b4)
{
  c1 = int_centers_from_gbs(b1);
  if (b1 != b2 || b1 != b3 || b1 != b4) {
      same_center = 0;
      c2 = int_centers_from_gbs(b2);
      c3 = int_centers_from_gbs(b3);
      c4 = int_centers_from_gbs(b4);
    }
  else {
      c2 = c3 = c4 = c1;
      same_center = 1;
    }

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

//////////////////////////////////////////////////////////////////////////

TwoBodyDerivIntV2::TwoBodyDerivIntV2(const RefGaussianBasisSet& b1,
                                     const RefGaussianBasisSet& b2,
                                     const RefGaussianBasisSet& b3,
                                     const RefGaussianBasisSet& b4) :
  TwoBodyDerivInt(b1,b2,b3,b4)
{
  c1 = int_centers_from_gbs(b1);
  if (b1 != b2 || b1 != b3 || b1 != b4) {
      same_center = 0;
      c2 = int_centers_from_gbs(b2);
      c3 = int_centers_from_gbs(b3);
      c4 = int_centers_from_gbs(b4);
    }
  else {
      c2 = c3 = c4 = c1;
      same_center = 1;
    }

  if (!c1 || !c2 || !c3 || !c4) {
    fprintf(stderr,"TwoBodyDerivIntV2::could not form centers\n");
    abort();
  }

  init();
}

TwoBodyDerivIntV2::~TwoBodyDerivIntV2()
{
  int_done_offsets2(c1,c2,c3,c4);
  int_done_erep();

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
TwoBodyDerivIntV2::init()
{
  int_initialize_offsets2(c1,c2,c3,c4);

  int flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2|INT_NODERB;

  intbuf = int_initialize_erep(flags,1,c1,c2,c3,c4);
}

void
TwoBodyDerivIntV2::compute_shell(int is, int js, int ks, int ls)
{
  der_centers_t dercenters;
  int sh[4], sz[4];
  
  sh[0]=is; sh[1]=js; sh[2]=ks; sh[3]=ls;

  int_erep_all1der_v(INT_EREP|INT_REDUND|INT_NOPERM,sh,sz,&dercenters);

  int ni = INT_SH_NFUNC((c1),is);
  int nj = INT_SH_NFUNC((c2),js);
  int nk = INT_SH_NFUNC((c3),ks);
  int nl = INT_SH_NFUNC((c4),ls);

  memcpy(buffer_, intbuf, sizeof(double)*ni*nj*nk*nl);
}
