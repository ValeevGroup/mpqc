
#include <stdio.h>

#include <chemistry/qc/wfn/obwfn.h>

// Force linkages:
const ClassDesc &fl0 = HCoreWfn::class_desc_;

main()
{
  char *input = "wfntest.in";

  RefKeyVal rpkv(new ParsedKeyVal(input));
  
  // the output stream is standard out
  SCostream o(stdout);

  RefOneBodyWavefunction wfn = rpkv->describedclassvalue("wavefunction");
  if (wfn.null()) {
    fprintf(stderr,"wfn is null\n");
    exit(1);
  }

  wfn->print(o);
  o << endl;

  RefOneBodyWavefunction oldwfn = rpkv->describedclassvalue("pwavefunction");
  
  RefSCMatrix evecs = wfn->projected_eigenvectors(oldwfn);

  evecs.print("projected wavefunction");

  StateOutText so("wfn.ckpt");
  wfn.save_state(so);
  so.close();

  RefMolecularEnergy me;
  StateInText si("wfn.ckpt");
  me.restore_state(si);
  
  me->print(o);
  o << me->value();
  
  return 0;
}
