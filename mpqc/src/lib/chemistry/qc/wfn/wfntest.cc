
#include <util/misc/formio.h>

#include <chemistry/qc/wfn/obwfn.h>

// Force linkages:
const ClassDesc &fl0 = HCoreWfn::class_desc_;

main(int argc, char *argv[])
{
  char *input = (argc > 1) ? argv[1] : SRCDIR "/wfntest.kv";

  RefKeyVal rpkv(new ParsedKeyVal(input));
  
  // the output stream is standard out
  ostream &o = cout;

  RefOneBodyWavefunction wfn = rpkv->describedclassvalue("wavefunction");
  if (wfn.null()) {
    cerr << node0 << "wfn is null\n";
    exit(1);
  }

  wfn->overlap()->print("overlap");
  wfn->core_hamiltonian()->print("Hcore");
  wfn->hcore_guess()->print("guess vector");

  //wfn->print(o);
  //o << endl;

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

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
