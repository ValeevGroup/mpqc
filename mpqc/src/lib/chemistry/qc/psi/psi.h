

#ifndef _chemistry_qc_psi_psi_h
#define _chemistry_qc_psi_psi_h

#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/qc/basis/basis.h>
#include "psiinput.h"

class PSISCF: public OneBodyWavefunction
{
#   define CLASSNAME PSISCF
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    PSI_Input_SCF psi_in;
  protected:
    void compute();
  public:
    PSISCF(const RefKeyVal&);
    PSISCF(StateIn&);
    virtual ~PSISCF();
    void save_data_state(StateOut&);
  
    void print(ostream& =cout);

    //int do_eigenvectors(int);
    //RefSCMatrix eigenvectors();
    double occupation(int vectornum);

    double occupation(int,int);
    RefDiagSCMatrix eigenvalues();
    RefSCMatrix eigenvectors();

    int gradient_implemented();
    int value_implemented();
};
  
#endif
