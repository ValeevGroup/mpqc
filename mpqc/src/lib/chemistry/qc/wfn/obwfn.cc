
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/scmat/local.h>

#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/wfn/obwfn.h>

SavableState_REF_def(OneBodyWavefunction);

#define CLASSNAME OneBodyWavefunction
#define PARENTS public Wavefunction
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
OneBodyWavefunction::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Wavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

OneBodyWavefunction::OneBodyWavefunction(const OneBodyWavefunction& obwfn) :
  Wavefunction(obwfn),
  _eigenvectors(obwfn._eigenvectors,this),
  _density(obwfn._density,this)
{
}

OneBodyWavefunction::OneBodyWavefunction(const RefKeyVal&keyval):
  Wavefunction(keyval),
  _eigenvectors(this),
  _density(this)
{
  _eigenvectors.set_desired_accuracy(
    keyval->doublevalue("eigenvector_accuracy"));
  if (_eigenvectors.desired_accuracy() < 1.0e-7)
    _eigenvectors.set_desired_accuracy(1.0e-7);
}

OneBodyWavefunction::~OneBodyWavefunction()
{
}

OneBodyWavefunction::OneBodyWavefunction(StateIn&s):
  Wavefunction(s),
  _eigenvectors(s,this),
  _density(s,this)
  maybe_SavableState(s)
{
}

OneBodyWavefunction &
OneBodyWavefunction::operator=(const OneBodyWavefunction& obwfn)
{
  Wavefunction::operator=(obwfn);
  _eigenvectors = obwfn._eigenvectors;
  _density = obwfn._density;
  return *this;
}

void
OneBodyWavefunction::save_data_state(StateOut&s)
{
  Wavefunction::save_data_state(s);
  _eigenvectors.save_data_state(s);
  _density.save_data_state(s);
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

RefSCMatrix
OneBodyWavefunction::projected_eigenvectors(const RefOneBodyWavefunction& owfn)
{
  // find out how many occupied orbitals there should be
  int nocc = 0;
  while (occupation(nocc)) nocc++;
  
  // first we need the old eigenvectors
  RefSCMatrix ovec = owfn->eigenvectors().get_subblock(
                                         0,owfn->basis()->nbasis()-1,0,nocc-1);
  //ovec.print("old wavefunction");
  
  // now form the overlap between the old basis and the new one
  RefSCMatrix s2(basis_dimension(), ovec.rowdim(), matrixkit());
  RefSCElementOp op =
    new OneBodyIntOp(integral()->overlap_int(basis(), owfn->basis()));
  s2.assign(0.0);
  s2.element_op(op);
  op = 0;

  //s2.print("overlap between new basis and old");
  
  // form C' = S2 * Cold
  RefSCMatrix cprime(s2.rowdim(),ovec.coldim(),matrixkit());
  cprime.assign(0.0);
  cprime.accumulate_product(s2,ovec);
  //cprime.print("C' matrix");
  
  // we're done with s2 and ovec, free up some memory
  s2=0;
  ovec=0;
  
  // we'll need the inverse of S eventually
  RefSymmSCMatrix s = overlap().copy();
  s->invert_this();
  //s.print("s inverse");

  // now we need a matrix D = S^-1 * C' 
  RefSCMatrix D = cprime.clone();
  D = s * cprime;
  //D.print("D matrix");
  
  // we also need X = C'~ * S^-1 * C'
  RefSymmSCMatrix X(cprime.coldim(),matrixkit());
  X.assign(0.0);
  X.accumulate_transform(cprime.t(),s);
  //X.print("X matrix");
  
  // we're done with cprime, free up some memory
  cprime=0;

  // now form X^-1/2
  form_m_half(X);
  //X.print("X^-1/2 matrix");

  // and form C'' = D * X^-1/2
  RefSCMatrix Cpp = D * X;
  //Cpp.print("new vector (occupied bits)");
  
  // we're done with X, free up some memory
  X=0;
  D=0;
  
  // now form core guess in the present basis, and then stuff in the
  // projected bit
  HCoreWfn hcwfn(*this);
  RefSCMatrix vec = hcwfn.eigenvectors();
  vec.assign_subblock(Cpp,0,Cpp.nrow()-1,0,nocc-1);
  //vec.print("new vector");

  vec->schmidt_orthog(overlap().pointer(),vec.ncol());
  
  Cpp=0;
  
  return vec;
}

double
OneBodyWavefunction::density(const SCVector3 &c)
{
  return Wavefunction::density(c);
}

RefSymmSCMatrix
OneBodyWavefunction::density()
{
  if (!_density.computed()) {
    RefSCMatrix vec = eigenvectors();
    RefSymmSCMatrix newdensity(basis_dimension(), matrixkit());
    form_density(vec,newdensity,0,0,0);
    _density = newdensity;
    _density.computed() = 1;
  }

  return _density;
}

// Function for returning an orbital value at a point
double
OneBodyWavefunction::orbital(const SCVector3& r, int iorb)
{
  return Wavefunction::orbital(r,iorb,eigenvectors());
}

// Function for returning an orbital value at a point
double
OneBodyWavefunction::orbital_density(const SCVector3& r,
                                            int iorb,
                                            double* orbval)
{
  return Wavefunction::orbital_density(r,iorb,eigenvectors(),orbval);
}

void
OneBodyWavefunction::print(ostream&o)
{
  Wavefunction::print(o);
}

//////////////////////////////////////////////////////////////////////////

void
OneBodyWavefunction::form_density(const RefSCMatrix& vec,
                                  const RefSymmSCMatrix& density,
                                  const RefSymmSCMatrix& density_diff,
                                  const RefSymmSCMatrix& open_density,
                                  const RefSymmSCMatrix& open_density_diff)
{
  int nbasis = basis()->nbasis();

  // find out how many doubly occupied orbitals there should be
  int ndocc = 0, nsocc=0;
  for (int i=0; i < nbasis; i++) {
    if (occupation(i) > 1.9)
      ndocc++;
    else if (occupation(i) > 0.9)
      nsocc++;
  }
  
  // find out what type of matrices we're dealing with
  if (LocalSCMatrix::castdown(vec.pointer())) {
    LocalSCMatrix *lvec = LocalSCMatrix::require_castdown(
      vec.pointer(), "OneBodyWavefunction::form_density");
    LocalSymmSCMatrix *ldens = LocalSymmSCMatrix::require_castdown(
      density.pointer(), "OneBodyWavefunction::form_density");
    LocalSymmSCMatrix *ldensd = LocalSymmSCMatrix::require_castdown(
      density_diff.pointer(), "OneBodyWavefunction::form_density");
    LocalSymmSCMatrix *lodens = LocalSymmSCMatrix::require_castdown(
      open_density.pointer(), "OneBodyWavefunction::form_density");
    LocalSymmSCMatrix *lodensd = LocalSymmSCMatrix::require_castdown(
      open_density_diff.pointer(), "OneBodyWavefunction::form_density");

    for (int i=0; i < nbasis; i++) {
      for (int j=0; j <= i; j++) {
        double pt=0;
        int k;
        for (k=0; k < ndocc; k++)
          pt += 2.0*lvec->get_element(i,k)*lvec->get_element(j,k);

        double pto=0;
        for (k=ndocc; k < ndocc+nsocc; k++)
          pto += lvec->get_element(i,k)*lvec->get_element(j,k);

        if (lodens) {
          if (lodensd)
            lodensd->set_element(i,j,pto-lodens->get_element(i,j));
          lodens->set_element(i,j,pto);
        }
        if (ldensd)
          ldensd->set_element(i,j,pt+pto-ldens->get_element(i,j));
        ldens->set_element(i,j,pt+pto);
      }
    }
  } else {
    density.assign(0.0);
    open_density.assign(0.0);
    for (int i=0; i < nbasis; i++) {
      for (int j=0; j <= i; j++) {
        int k;
        for (k=0; k < ndocc; k++) {
          density.set_element(i,j, density.get_element(i,j)
                              + 2.0*vec.get_element(i,k)*vec.get_element(j,k));
        }
        for (k=ndocc; k < ndocc+nsocc; k++) {
          open_density.set_element(i,j, open_density.get_element(i,j)
                              + vec.get_element(i,k)*vec.get_element(j,k));
        }
      }
    }
  }
}
