
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
main(int, char *argv[])
{
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

  RefSymmSCMatrix s(gbs->basisdim());
  s.assign(0.0);

  tim_enter("intv2 overlap");
  RefSCElementOp op = new GaussianOverlapIntv2(gbs);
  s.element_op(op);
  op=0;
  tim_exit("intv2 overlap");

  //s.print("overlap integrals v2");

  s.assign(0.0);

  tim_enter("jf overlap");
  op = new GaussianOverlapIntJF(gbs);
  s.element_op(op);
  op=0;
  tim_exit("jf overlap");

  //s.print("overlap integrals jf");

  tim_print(0);

  ///////////////////////////////////////////////////////////////////////
  
  s.assign(0.0);

  tim_enter("intv2 ke");
  op = new GaussianKineticIntv2(gbs);
  s.element_op(op);
  op=0;
  tim_exit("intv2 ke");

  //s.print("ke integrals v2");

  s.assign(0.0);

  tim_enter("jf ke");
  op = new GaussianKineticIntJF(gbs);
  s.element_op(op);
  op=0;
  tim_exit("jf ke");

  //s.print("ke integrals jf");
  tim_print(0);
  
  return 0;
}
