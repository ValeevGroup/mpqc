
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/integral/integralv2.h>

#define CLASSNAME HCoreWfn
#define PARENTS public OneBodyWavefunction
//#include <util/state/statei.h>
#include <util/class/classia.h>
void *
HCoreWfn::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = OneBodyWavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

HCoreWfn::HCoreWfn(const RefKeyVal&keyval):
  OneBodyWavefunction(keyval)
{
}

HCoreWfn::~HCoreWfn()
{
}

RefSCMatrix
HCoreWfn::eigenvectors()
{
  fprintf(stderr,"HCoreWfn: not using smhalf yet\n");
  RefSymmSCMatrix h(basis_dimension());
  h.assign(0.0);
  RefSCElementOp op = new GaussianKineticIntv2(basis(), molecule());
  h.element_op(op);
  op = new GaussianNuclearIntv2(basis(), molecule());
  h.element_op(op);
  RefSCMatrix vec(basis_dimension(), basis_dimension());
  RefDiagSCMatrix val(basis_dimension());
  h.diagonalize(val,vec);
  return vec;
}

double
HCoreWfn::occupation(int)
{
  return 0.0;
}

void
HCoreWfn::compute()
{
  fprintf(stderr,"HCoreWfn::compute(): don't call me\n");
  abort();
}
