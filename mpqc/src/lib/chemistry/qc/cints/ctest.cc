
#include <stdio.h>

#include <util/misc/timer.h>
#include <util/group/picl.h>

#include <chemistry/qc/basis/symgaussbas.h>
#include <chemistry/qc/cints/cints.h>
#include <chemistry/qc/cints/integraljf.h>
#include <chemistry/qc/cints/int2jf.h>
#include <chemistry/qc/integral/integralv2.h>
#include <chemistry/qc/intv2/int_libv2.h>

// force linkage...this MUST be fixed
const ClassDesc &fl0  = SymmGaussianBasisSet::class_desc_;

int
main(int argc, char *argv[])
{
  int i;
  
  char *filename = (argv[1]) ? argv[1] : SRCDIR "/ctest.kv";
  
  RefKeyVal keyval = new ParsedKeyVal(filename);

  int nproc, me, host, top, ord, dir;
  RefMessageGrp grp = MessageGrp::get_default_messagegrp();
  open0_messagegrp(&nproc, &me, &host, grp);
  setarc0(&nproc,&top,&ord,&dir);

  RefGaussianBasisSet gbs = keyval->describedclassvalue("basis");

  printf("\nMolecule:\n");
  gbs->molecule()->print();
  
  printf("\nnbasis = %d\n",gbs->nbasis());
  printf("nucrep = %20.15f\n",cints_nuclear_repulsion_energy(gbs));

  ////////////////////////////////////////////////////////////////////////////
  
  RefSymmSCMatrix s(gbs->basisdim());
  s.assign(0.0);

  tim_enter("intv2 overlap");
  RefSCElementOp op = new GaussianOverlapIntv2(gbs);
  s.element_op(op);
  op=0;
  tim_exit("intv2 overlap");

  if (argc > 1)
    s.print("overlap integrals v2");

  s.assign(0.0);

  tim_enter("jf overlap");
  op = new GaussianOverlapIntJF(gbs);
  s.element_op(op);
  op=0;
  tim_exit("jf overlap");

  if (argc > 1)
    s.print("overlap integrals jf");

  tim_print(0);

  ///////////////////////////////////////////////////////////////////////
  
#if 1
  s.assign(0.0);

  tim_enter("intv2 ke");
  op = new GaussianKineticIntv2(gbs);
  s.element_op(op);
  op=0;
  tim_exit("intv2 ke");

  if (argc > 1)
    s.print("ke integrals v2");

  tim_enter("intv2 pe");
  op = new GaussianNuclearIntv2(gbs);
  s.element_op(op);
  op=0;
  tim_exit("intv2 pe");

  if (argc > 1)
    s.print("hcore v2");

  RefSCMatrix vecs(gbs->basisdim(),gbs->basisdim());
  RefDiagSCMatrix vals(gbs->basisdim());
  s.diagonalize(vals,vecs);

  double sum=0;
  for (i=0; i < gbs->nbasis(); i++)
    sum += vals.get_element(i);
  
  printf("\n  sum of evals of hcore = %20.15f\n\n",sum);
  
  tim_print(0);
#endif
  
  ///////////////////////////////////////////////////////////////////////

  s.assign(0.0);

  tim_enter("jf ke");
  op = new GaussianKineticIntJF(gbs);
  s.element_op(op);
  op=0;
  tim_exit("jf ke");

  if (argc > 1)
    s.print("ke integrals jf");

  tim_enter("jf pe");
  op = new GaussianNuclearIntJF(gbs);
  s.element_op(op);
  op=0;
  tim_exit("jf pe");

  if (argc > 1)
    s.print("hcore jf");

#if 1
  s.diagonalize(vals,vecs);

  sum=0;
  for (i=0; i < gbs->nbasis(); i++)
    sum += vals.get_element(i);
  
  printf("\n  sum of evals of hcore = %20.15f\n\n",sum);
#endif

  tim_print(0);
  
  ///////////////////////////////////////////////////////////////////////

  tim_enter("2ei");
  tim_enter("init");
  TwoBodyIntJF twos(gbs);
  tim_exit("init");

  for (i=0; i < gbs->nshell(); i++) {
    for (int j=0; j <= i; j++) {
      for (int k=0; k <= i; k++) {
        for (int l=0; l <= ((i==k) ? j : k); l++) {
          twos.compute_shell(i,j,k,l,0);
        }
      }
    }
  }
  
  tim_exit("2ei");
  tim_print(0);

  tim_enter("2ei v2");
  tim_enter("init");

  centers_t *centers = gbs->convert_to_centers_t();
  
  if (!centers) {
    fprintf(stderr,"hoot man!  no centers\n");
    abort();
  }

  int_initialize_offsets2(centers,centers,centers,centers);

  tim_exit("init");

  int flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2;
  double *intbuf = 
    int_initialize_erep(flags,0,centers,centers,centers,centers);

  for (i=0; i < centers->nshell; i++) {
    for (int j=0; j <= i; j++) {
      for (int k=0; k <= i; k++) {
        for (int l=0; l <= ((k==i)?j:k); l++) {
          int s1=i, s2=j, s3=k, s4=l;
          int_erep(INT_REDUND|INT_EREP|INT_NOBCHK|INT_NOPERM,&s1,&s2,&s3,&s4);
        }
      }
    }
  }

  tim_exit("2ei v2");
  tim_print(0);
  
  return 0;
}
