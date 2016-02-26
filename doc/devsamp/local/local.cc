//
// local.cc
//
// Copyright (C) Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
//

#include <mpqc_config.h>

#include <mpqcinit.h>
#include <chemistry/qc/scf/clhf.h>
#include <chemistry/qc/lmp2/pop_local.h>

using namespace std;
using namespace sc;

#ifdef HAVE_LIBINT2
#  include <chemistry/qc/libint2/linkage.h>
#endif
#include <util/state/linkage.h>

int main(int argc, char**argv)
{
  GetLongOpt opt;
  MPQCInit init(opt,argc,argv);
  int optind = opt.parse(argc, argv);

  Ref<ThreadGrp> thr = init.init_threadgrp();
  Ref<MessageGrp> msg = init.init_messagegrp();
  init.init_integrals();
  Ref<Integral> intf = Integral::get_default_integral();
  assert(std::string(intf->class_name()) == "IntegralLibint2");

  //
  // construct a Molecule object
  //
  Ref<Molecule> mol;
  {
    Ref<AssignedKeyVal> akv = new AssignedKeyVal;
//#define CONSTRUCT_MOLECULE_FROM_XYZ // undefine if want to construct programmatically
#ifdef CONSTRUCT_MOLECULE_FROM_XYZ
    akv->assign("xyz_file", "h2o.xyz");
    akv->assign("symmetry", "c1");
    mol = new Molecule(Ref<KeyVal>(akv));
#else
    mol = new Molecule;
    mol->add_atom(8, 0.0, 0.0, 0.0);
    mol->add_atom(1, 0.0, 1.0, 1.0);
    mol->add_atom(1, 0.0,-1.0, 1.0);
    Ref<PointGroup> c1_ptgrp = new PointGroup("C1");
    mol->set_point_group(c1_ptgrp);
#endif
    std::cout << "Molecule object:" << std::endl;
    mol->print(std::cout);
  }

  //
  // construct a GaussianBasisSet object
  //
  Ref<GaussianBasisSet> obs;
  {
    Ref<AssignedKeyVal> akv = new AssignedKeyVal;
    akv->assign("molecule", mol.pointer());
    akv->assign("name", "cc-pVDZ");
    obs = new GaussianBasisSet(Ref<KeyVal>(akv));
    std::cout << "GaussianBasisSet object:" << std::endl;
    obs->print(std::cout);
  }

  //
  // construct a CLHF object
  //
  Ref<CLHF> clhf;
  {
    Ref<AssignedKeyVal> akv = new AssignedKeyVal;
    akv->assign("molecule", mol.pointer());
    akv->assign("basis", obs.pointer());
    clhf = new CLHF(Ref<KeyVal>(akv));
    std::cout << "CLHF object:" << std::endl;
    clhf->print(std::cout);
  }

  const double e_hf = clhf->energy();
  const int nfzc = 1;
  RefSCMatrix Cloc = pop_local_mo(clhf, nfzc, clhf->overlap(), msg);
  Cloc.print("localized MOs", std::cout);

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
