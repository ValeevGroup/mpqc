
#ifndef _chemistry_qc_basis_gaussbas_h
#define _chemistry_qc_basis_gaussbas_h

#ifdef __GNUC__
#pragma interface
#endif

#include <stdio.h>
#include <util/state/state.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>
#include <chemistry/molecule/molecule.h>

class GaussianShell;
class RefKeyVal;
class BasisFileSet;

SavableState_REF_fwddec(Integral)

class GaussianBasisSet: public SavableState
{
#   define CLASSNAME GaussianBasisSet
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    char* name_;
    GaussianShell** shell_;
    Arrayint shell_to_function_;
    Arrayint function_to_shell_;

    RefMolecule molecule_;

    RefSCMatrixKit matrixkit_;
    RefSCDimension basisdim_;
    
    int ncenter_;
    SSBArrayint shell_to_center_;
    SSBArrayint center_to_shell_;
    SSBArrayint center_to_nshell_;

    int nshell_;
    int nbasis_;
    int nprim_;

    void recursively_get_shell(int&,RefKeyVal&,
                               const char*,const char*,BasisFileSet&,
                               int,int,int);

    void init(RefMolecule&,RefKeyVal&,
              BasisFileSet&,
              int have_userkeyval,
              int pure);
    void init2();
    
  protected:
    GaussianBasisSet(const GaussianBasisSet&);
    virtual void set_matrixkit(const RefSCMatrixKit&);
    
  public:
    GaussianBasisSet(const RefKeyVal&);
    GaussianBasisSet(StateIn&);
    virtual ~GaussianBasisSet();

    void save_data_state(StateOut&);

    const char* name() const;

    RefMolecule molecule() const { return molecule_; }
    RefSCMatrixKit matrixkit() { return matrixkit_; }
    RefSCDimension basisdim() { return basisdim_; }

    int ncenter() const;
    int nshell() const;
    int nshell_on_center(int icenter) const;
    int shell_on_center(int icenter, int shell) const;
    int shell_to_center(int shell) const;
    int nbasis() const;
    int nprimitive() const;

    int max_nfunction_in_shell() const;
    int max_angular_momentum() const;

    int shell_to_function(int i) const;
    int function_to_shell(int i) const;

    // access to shells thru overall shell number
    const GaussianShell& operator()(int i) const;
    GaussianShell& operator()(int i);
    const GaussianShell& operator[](int i) const;
    GaussianShell& operator[](int i);
    const GaussianShell& shell(int i) const { return operator()(i); }
    GaussianShell& shell(int i) { return operator()(i); }

    // access to shells thru center number and relative shell number
    const GaussianShell& operator()(int icenter,int ishell) const;
    GaussianShell& operator()(int icenter,int ishell);

    // access to r thru center number
    double r(int icenter,int xyz) const;
    
    // compute the value for this basis set at position r
    int values(const RefIntegral&,
               const SCVector3& r, double* basis_values) const;
    int grad_values(const RefIntegral&, const SCVector3& r,
                    double*g_values,double* basis_values=0)const;

    // fill in matrix with a matrix that orthogonalizes the basis functions
    void ortho(const RefIntegral&, const RefSCMatrix&ortho);
    void ortho(const RefIntegral&,
               const RefSCMatrix&ortho, const RefSCMatrix&ortho_inverse);

    void print(FILE*fp=stdout) const;
};

SavableState_REF_dec(GaussianBasisSet);

#ifdef INLINE_FUNCTIONS
#include <chemistry/qc/basis/gaussbas_i.h>
#endif

#endif
