
#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
}

#include <iostream.h>

#include <util/keyval/keyval.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/symmint.h>
#include <chemistry/qc/intv3/intv3.h>

#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/wfn/hcore.h>

SavableState_REF_def(Wavefunction);

#define CLASSNAME Wavefunction
#define PARENTS public MolecularEnergy
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
Wavefunction::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = MolecularEnergy::_castdown(cd);
  return do_castdowns(casts,cd);
}

Wavefunction::Wavefunction(const RefKeyVal&keyval):
  // this will let molecule be retrieved from basis
  MolecularEnergy(new AggregateKeyVal(keyval,
                                      new PrefixKeyVal("basis", keyval))),
  overlap_(this),
  hcore_(this),
  natural_orbitals_(this),
  natural_density_(this),
  bs_values(0),
  bsg_values(0)
{
  overlap_.compute() = 0;
  hcore_.compute() = 0;
  natural_orbitals_.compute() = 0;
  natural_density_.compute() = 0;

  overlap_.computed() = 0;
  hcore_.computed() = 0;
  natural_orbitals_.computed() = 0;
  natural_density_.computed() = 0;

  gbs_ = GaussianBasisSet::require_castdown(
    keyval->describedclassvalue("basis").pointer(),
    "Wavefunction::Wavefunction\n"
    );

  integral_ = keyval->describedclassvalue("integrals");
  if (integral_.null())
    integral_ = new IntegralV3(gbs_);
  
  integral_->set_basis(gbs_);
  RefPetiteList pl = integral_->petite_list();

  basisdim_ = pl->SO_basisdim();
  basiskit_ = gbs_->so_matrixkit();
}

Wavefunction::Wavefunction(StateIn&s):
  MolecularEnergy(s),
  overlap_(s,this),
  hcore_(s,this),
  natural_orbitals_(s,this),
  natural_density_(s,this),
  bs_values(0),
  bsg_values(0)
  maybe_SavableState(s)
{
  overlap_.compute() = 0;
  hcore_.compute() = 0;
  natural_orbitals_.compute() = 0;
  natural_density_.compute() = 0;

  overlap_.computed() = 0;
  hcore_.computed() = 0;
  natural_orbitals_.computed() = 0;
  natural_density_.computed() = 0;

  gbs_.restore_state(s);
  integral_.restore_state(s);

  integral_->set_basis(gbs_);
  RefPetiteList pl = integral_->petite_list();

  basisdim_ = pl->SO_basisdim();
  basiskit_ = gbs_->so_matrixkit();
}

Wavefunction::~Wavefunction()
{
  if (bs_values) {
    delete[] bs_values;
    bs_values=0;
  }
  if (bsg_values) {
    delete[] bsg_values;
    bsg_values=0;
  }
}

void
Wavefunction::save_data_state(StateOut&s)
{
  MolecularEnergy::save_data_state(s);

  // overlap and hcore integrals are cheap so don't store them.
  // same goes for natural orbitals

  gbs_.save_state(s);
  integral_.save_state(s);
}

RefSCMatrix
Wavefunction::natural_orbitals()
{
  if (!natural_orbitals_.computed()) {
      RefSymmSCMatrix dens = density();

      // transform the density into an orthogonal basis
      RefSymmSCMatrix ortho = ao_to_orthog_ao();
      RefSymmSCMatrix orthoi = ortho.i();

      RefSymmSCMatrix densortho(basis_dimension(), basis_matrixkit());
      densortho.assign(0.0);
      densortho.accumulate_transform(orthoi,dens);

      RefSCMatrix natorb(basis_dimension(), basis_dimension(),
                         basis_matrixkit());
      RefDiagSCMatrix natden(basis_dimension(), basis_matrixkit());
      natural_orbitals_ = natorb;
      natural_density_ = natden;

      densortho.diagonalize(natural_density_,natural_orbitals_);

      // _natural_orbitals is the ortho to NO basis transform
      // make _natural_orbitals the AO to the NO basis transform
      natural_orbitals_ = ortho * natural_orbitals_;

      natural_orbitals_.computed() = 1;
      natural_density_.computed() = 1;
    }

  return natural_orbitals_.result_noupdate();
}

RefDiagSCMatrix
Wavefunction::natural_density()
{
  if (!natural_density_.computed()) {
      RefSymmSCMatrix dens = density();

      RefSCMatrix natorb(basis_dimension(), basis_dimension(),
                         basis_matrixkit());
      RefDiagSCMatrix natden(basis_dimension(), basis_matrixkit());
      natural_orbitals_ = natorb;
      natural_density_ = natden;

      dens.diagonalize(natural_density_,natural_orbitals_);

      natural_orbitals_.computed() = 1;
      natural_density_.computed() = 1;
    }

  return natural_density_.result_noupdate();
}

RefSymmSCMatrix
Wavefunction::overlap()
{
  if (!overlap_.computed()) {
    integral()->set_basis(gbs_);
    RefPetiteList pl = integral()->petite_list();

    // first form skeleton s matrix
    RefSymmSCMatrix s(basis()->basisdim(), basis()->matrixkit());
    RefSCElementOp ov =
      new OneBodyIntOp(new SymmOneBodyIntIter(integral()->overlap(), pl));

    s.assign(0.0);
    s.element_op(ov);
    ov=0;

    // then symmetrize it
    RefSymmSCMatrix sb(basis_dimension(), basis_matrixkit());
    pl->symmetrize(s,sb);

    overlap_ = sb;
    overlap_.computed() = 1;
  }

  return overlap_.result_noupdate();
}

RefSymmSCMatrix
Wavefunction::core_hamiltonian()
{
  if (!hcore_.computed()) {
    RefAccumDIH hc = new SymmAccumHCore();
    hc->init(basis(), integral());

    RefSymmSCMatrix h(basis_dimension(), basis_matrixkit());
    hc->accum(h);
    hc=0;

    hcore_ = h;
    hcore_.computed() = 1;
  }

  return hcore_.result_noupdate();
}

// at some point this will have to check for zero eigenvalues and not
// invert them
static void
form_m_half(RefSymmSCMatrix& M)
{
  // Diagonalize M to get m and U
  RefSCMatrix U(M.dim(), M.dim(), M.kit());
  RefDiagSCMatrix m(M.dim(), M.kit());
  M.diagonalize(m,U);

  // take square root of all elements of m
  RefSCElementOp op = new SCElementSquareRoot;
  m.element_op(op);

  // invert m
  op = new SCElementInvert;
  m.element_op(op);

  // back transform m^-1/2 to get M^-1/2 ( U*m^-1/2*U~ = M^-1/2)
  M.assign(0.0);
  M.accumulate_transform(U,m);
}

// returns a matrix containing S^-1/2
RefSymmSCMatrix
Wavefunction::ao_to_orthog_ao()
{
  // first calculate S
  RefSymmSCMatrix s = overlap().copy();
  
  // then form S^-1/2
  form_m_half(s);

  return s;
}

RefGaussianBasisSet
Wavefunction::basis()
{
  return gbs_;
}

RefIntegral
Wavefunction::integral()
{
  return integral_;
}

RefSCDimension
Wavefunction::basis_dimension()
{
  return basisdim_;
}

RefSCMatrixKit
Wavefunction::basis_matrixkit()
{
  return basiskit_;
}

void
Wavefunction::print(ostream&o)
{
  MolecularEnergy::print(o);
  // the other stuff is a wee bit too big to print
}
