
#ifndef _chemistry_qc_basis_gaussbas_h
#define _chemistry_qc_basis_gaussbas_h

#ifdef __GNUC__
#pragma interface
#endif

#include <stdio.h>
#include <util/container/ref.h>
//#include <util/container/pixmap.h>
#include <util/state/state.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>

// this is for centers_t*
#include <chemistry/qc/intv2/atoms.h>

class GaussianShell;
class RefKeyVal;
class cart_point;

class GaussianBasisSet: virtual public SavableState
{
#   define CLASSNAME GaussianBasisSet
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    char* name_;
    GaussianShell** shell;
    Arrayint shell_to_function_;
    Arrayint function_to_shell_;

    // Not using pixes anymore.  This is because they refer to molecule
    // pixes which are not really a part of the basis set class.  They
    // also make it difficult to restore state without molecule being a
    // member of GaussianBasisSet, which I don't really like.
    // Pix* shell_to_centerpix;
    // PixMap<int> centerpix_to_shellnum; // the number of the first shell is 0
    // PixMap<double*> centerpix_to_r;
    // PixMap<int> centerpix_to_nshell;

    // pix replacments (using center numbers instead of center pixes):
    int ncenter_;
    SSBArrayint shell_to_center_;
    SSBArrayint center_to_shell_;
    SSBArray2double center_to_r_;
    SSBArrayint center_to_nshell_;

    int nshell_;
    int nbasis_;
    int nprim_;

    void recursively_get_shell(int&,const RefKeyVal&,
                               const char*,const char*,int,int,int);

    void init(RefMolecule&,const RefKeyVal&,
              const char* basisname,
              int have_userkeyval,
              int pure);
    void init2();
  public:
    GaussianBasisSet(const RefKeyVal&);
    GaussianBasisSet(StateIn&);
    // pure is -1 if the pure from the basis set data is used.
    // Otherwise 0 is cartesian and 1 is pure am.
    GaussianBasisSet(RefMolecule&,const char*basisname,int pure = -1);
    virtual ~GaussianBasisSet();
    void save_data_state(StateOut&);
    const char* name() const;

    int ncenter() const;
    int nshell() const;
    int nshell_on_center(int icenter) const;
    int nbasis() const;
    int nprimitive() const;

    int max_nfunction_in_shell() const;

    int shell_to_function(int i) const;
    int function_to_shell(int i) const;

    // access to shells thru overall shell number
    const GaussianShell& operator()(int i) const;
    GaussianShell& operator()(int i);
    const GaussianShell& operator[](int i) const;
    GaussianShell& operator[](int i);

    // access to shells thru center number and relative shell number
    const GaussianShell& operator()(int icenter,int ishell) const;
    GaussianShell& operator()(int icenter,int ishell);

    // access to r thru center number
    const double& r(int icenter,int xyz) const;
    double& r(int icenter,int xyz);
    
    // converts the basis set to a centers_t for compatibility with libintv2
    // If the molecule is 0 then fake molecule information is included--
    // this is useful in computing the overlap, for example.
    centers_t* convert_to_centers_t(const Molecule*) const;
    centers_t* convert_to_centers_t(const RefMolecule&mol) const {
        convert_to_centers_t(mol.pointer());
      }
    //operator struct struct_centers*();

    // compute the value for this basis set at position r
    int values(cart_point& r, double* basis_values) const;
    int grad_values(cart_point& r,double*g_values,double* basis_values=0)const;

    // fill in matrix with a matrix that orthogonalizes the basis functions
    // note: this member is provided in the integrals library
    void ortho(const RefSCMatrix&ortho);
    void ortho(const RefSCMatrix&ortho,
               const RefSCMatrix&ortho_inverse);

    void print(FILE*fp=stdout) const;
};

SavableState_REF_dec(GaussianBasisSet);

#ifdef INLINE_FUNCTIONS
#include <chemistry/qc/basis/gaussbas_i.h>
#endif

#endif
