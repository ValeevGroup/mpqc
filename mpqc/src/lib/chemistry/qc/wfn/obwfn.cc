
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/scmat/local.h>
#include <math/scmat/blocked.h>

#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/wfn/obwfn.h>


#ifndef DBL_EPSILON
#define DBL_EPSILON 1.0e-15
#endif

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

OneBodyWavefunction::OneBodyWavefunction(const RefKeyVal&keyval):
  Wavefunction(keyval),
  eigenvectors_(this),
  density_(this),
  nirrep_(0),
  nvecperirrep_(0)
{
  eigenvectors_.set_desired_accuracy(
    keyval->doublevalue("eigenvector_accuracy"));

  if (eigenvectors_.desired_accuracy() > 1.0e-7)
    eigenvectors_.set_desired_accuracy(1.0e-7);

  if (eigenvectors_.desired_accuracy() < DBL_EPSILON)
    eigenvectors_.set_desired_accuracy(DBL_EPSILON);
}

OneBodyWavefunction::OneBodyWavefunction(StateIn&s):
  Wavefunction(s),
  eigenvectors_(this),
  density_(this),
  nirrep_(0),
  nvecperirrep_(0)
  maybe_SavableState(s)
{
  eigenvectors_.result_noupdate() =
    basis_matrixkit()->matrix(basis_dimension(), basis_dimension());
  eigenvectors_.restore_state(s);
  eigenvectors_.result_noupdate().restore(s);

  density_.result_noupdate() =
    basis_matrixkit()->symmmatrix(basis_dimension());
  density_.restore_state(s);
  density_.result_noupdate().restore(s);
}

OneBodyWavefunction::~OneBodyWavefunction()
{
  delete[] nvecperirrep_;
}

void
OneBodyWavefunction::save_data_state(StateOut&s)
{
  Wavefunction::save_data_state(s);

  eigenvectors_.save_data_state(s);
  eigenvectors_.result_noupdate().save(s);

  density_.save_data_state(s);
  density_.result_noupdate().save(s);
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
  // first let's calculate the overlap between the new and old basis set
  integral()->set_basis(basis(), owfn->basis());
  RefSCMatrix s2(basis()->basisdim(), owfn->basis()->basisdim(),
                 basis()->matrixkit());

  RefSCElementOp op = new OneBodyIntOp(integral()->overlap());
  s2.assign(0.0);
  s2.element_op(op);
  op = 0;

  integral()->set_basis(basis());

  // now transform the AO s2 to the SO basis. this probably doesn't really
  // work.
  RefPetiteList pl = integral()->petite_list();
  RefPetiteList plo = integral()->petite_list(owfn->basis());
  RefSCMatrix s2b(pl->AO_basisdim(), plo->AO_basisdim(),
                  basis()->so_matrixkit());
  s2b->convert(s2);
  s2 = pl->aotoso().t() * s2b * plo->aotoso();

  //s2.print("overlap between new basis and old");
  s2b=0;

  BlockedSCMatrix *s2p = BlockedSCMatrix::require_castdown(s2.pointer(),
    "OneBodyWavefunction::projected_eigenvectors: s2p"
    );
  
  // now get the old eigenvectors and the new core hamiltonian guess
  RefSCMatrix ovec = owfn->eigenvectors();
  BlockedSCMatrix *ovecp = BlockedSCMatrix::require_castdown(ovec.pointer(),
    "OneBodyWavefunction::projected_eigenvectors: ovec"
    );
  //ovec.print("old wavefunction");
    
  RefSCMatrix vec = hcore_guess();
  BlockedSCMatrix *vecp = BlockedSCMatrix::require_castdown(vec.pointer(),
    "OneBodyWavefunction::projected_eigenvectors: vec"
    );
  //vec.print("hcore guess wavefunction");

  // we'll need the inverse of S eventually
  RefSymmSCMatrix s = overlap();
  BlockedSymmSCMatrix *sp = BlockedSymmSCMatrix::require_castdown(
    s.pointer(),
    "OneBodyWavefunction::projected_eigenvectors: s"
    );

  RefSymmSCMatrix sinv = s.i();
  BlockedSymmSCMatrix *sinvp = BlockedSymmSCMatrix::require_castdown(
    sinv.pointer(),
    "OneBodyWavefunction::projected_eigenvectors: sinv"
    );
  //sinv.print("Sinv");

  for (int irrep=0; irrep < vecp->nblocks(); irrep++) {
    // find out how many occupied orbitals there should be
    int nocc = 0;
    while (occupation(irrep,nocc) && nocc < plo->nfunction(irrep)) nocc++;
  
    if (!nocc)
      continue;
    
    // get subblock of old vector
    RefSCMatrix ovecsb = ovecp->block(irrep).get_subblock(
      0, plo->nfunction(irrep)-1, 0, nocc-1);
    //ovecsb.print("old vec subblock");
  
    // form C' = S2 * Cold
    RefSCMatrix cprime = s2p->block(irrep) * ovecsb;
    //cprime.print("C' matrix");
  
    // we're done with ovecsb, free up some memory
    ovecsb=0;
  
    // now we need a matrix D = S^-1 * C' 
    RefSCMatrix D = sinvp->block(irrep) * cprime;
    //D.print("D matrix");
  
    // we also need X = C'~ * S^-1 * C'
    RefSymmSCMatrix X(cprime.coldim(),basis()->matrixkit());
    X.assign(0.0);
    X.accumulate_transform(cprime.t(),sinvp->block(irrep));
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
  
    // stuff the projected bit into guess vector
    vecp->block(irrep).assign_subblock(Cpp,0,Cpp.nrow()-1,0,nocc-1);
    vecp->block(irrep)->schmidt_orthog(sp->block(irrep).pointer(),
                                       pl->nfunction(irrep));

    Cpp=0;
  }

  return vec;
}

RefSCMatrix
OneBodyWavefunction::hcore_guess()
{
  RefSymmSCMatrix hcore = core_hamiltonian();
  RefSCMatrix vec(basis_dimension(), basis_dimension(), basis_matrixkit());
  RefDiagSCMatrix val(basis_dimension(), basis_matrixkit());
  hcore.diagonalize(val,vec);

  vec = ao_to_orthog_ao() * vec;

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
#if 0
  if (!density_.computed()) {
    RefSCMatrix vec = eigenvectors();
    RefSymmSCMatrix newdensity(basis_dimension(), basis_matrixkit());

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
        newdensity.pointer(), "OneBodyWavefunction::form_density");

      for (int i=0; i < nbasis; i++) {
        for (int j=0; j <= i; j++) {
          double pt=0;
          int k;
          for (k=0; k < ndocc; k++)
            pt += 2.0*lvec->get_element(i,k)*lvec->get_element(j,k);

          double pto=0;
          for (k=ndocc; k < ndocc+nsocc; k++)
            pto += lvec->get_element(i,k)*lvec->get_element(j,k);

          ldens->set_element(i,j,pt+pto);
        }
      }
    } else {
      for (int i=0; i < nbasis; i++) {
        for (int j=0; j <= i; j++) {
          int k;
          double pt=0;
          for (k=0; k < ndocc; k++)
            pt += 2.0*vec.get_element(i,k)*vec.get_element(j,k);

          double pto=0;
          for (k=ndocc; k < ndocc+nsocc; k++)
            pto += vec.get_element(i,k)*vec.get_element(j,k);

          newdensity.set_element(i,j,pt+pto);
        }
      }
    }

    density_ = newdensity;
    density_.computed() = 1;
  }

  return density_.result_noupdate();
#else
  return 0;
#endif
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

void
OneBodyWavefunction::init_sym_info()
{
  RefPetiteList pl = integral()->petite_list();

  nirrep_ = pl->nirrep();
  nvecperirrep_ = new int[nirrep_];
  for (int i=0; i<nirrep_; i++) {
    nvecperirrep_[i] = pl->nfunction(i);
  }
}

double
OneBodyWavefunction::occupation(int vectornum)
{
  if (!nirrep_) init_sym_info();
  int svecnum = vectornum;
  int i;
  for (i=0; vectornum >= 0 && i<nirrep_; i++) {
    svecnum = vectornum;
    vectornum -= nvecperirrep_[i];
  }
  return occupation(i-1, svecnum);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
