
#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/integral/integralv2.h>

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
  RefSCMatrix U(M.dim(), M.dim());
  RefDiagSCMatrix m(M.dim());
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
  RefSCMatrix s2(basis_dimension(), ovec.rowdim());
  RefSCElementOp op = new GaussianOverlapIntv2(basis(), owfn->basis());
  s2.assign(0.0);
  s2.element_op(op);
  op = 0;

  //s2.print("overlap between new basis and old");
  
  // form C' = S2 * Cold
  RefSCMatrix cprime(s2.rowdim(),ovec.coldim());
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
  RefSymmSCMatrix X(cprime.coldim());
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

// returns a matrix containing S^-1/2
RefSymmSCMatrix
OneBodyWavefunction::ao_to_orthog_ao()
{
  // first calculate S
  RefSymmSCMatrix s = overlap().copy();
  
  // then form S^-1/2
  form_m_half(s);

  return s;
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
      RefSCMatrix ortho(basis_dimension(),basis_dimension());
      RefSCMatrix orthoi(basis_dimension(),basis_dimension());
      basis()->ortho(ortho,orthoi);
      int nbasis = basis()->nbasis();

      RefSymmSCMatrix newdensity(basis_dimension());
      _density = newdensity;
      newdensity.assign(0.0);
      for (int k=0; k<nbasis; k++) {
          double occ = occupation(k);
          if (occ == 0.0) continue;
          for (int i=0; i<nbasis; i++) {
              for (int j=0; j<=i; j++) {
                  newdensity.set_element(i,j,
                             newdensity.get_element(i,j)
                             + occ*vec.get_element(k,i)*vec.get_element(k,j));
                }
            }
        }

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
OneBodyWavefunction::print(SCostream&o)
{
  Wavefunction::print(o);
}
