
#include <stdio.h>

#include <util/misc/timer.h>
#include <util/group/picl.h>

#include <chemistry/qc/basis/symgaussbas.h>
#include <chemistry/qc/cints/cints.h>
#include <chemistry/qc/cints/integraljf.h>
#include <chemistry/qc/integral/integralv2.h>

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
  
  printf("\nnucrep = %20.15f\n",cints_nuclear_repulsion_energy(gbs));

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

  s.diagonalize(vals,vecs);

  sum=0;
  for (i=0; i < gbs->nbasis(); i++)
    sum += vals.get_element(i);
  
  printf("\n  sum of evals of hcore = %20.15f\n\n",sum);

  tim_print(0);
  
  return 0;
}
