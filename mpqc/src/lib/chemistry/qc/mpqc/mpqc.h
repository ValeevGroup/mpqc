
#ifndef _chemistry_qc_mpqc_h
#define _chemistry_qc_mpqc_h

extern "C" {
#include <stdio.h>
#include <math/dmt/libdmt.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
};
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include <chemistry/qc/dmtscf/scf_dmt.h>
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/molecule/energy.h>
#include <math/scmat/matrix.h>
#include <math/topology/point.h>
#include <chemistry/qc/basis/basis.h>
class Molecule;
class MolecularCoor;

class MPSCF: public OneBodyWavefunction
{
#   define CLASSNAME MPSCF
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    char vecfile[512];
    int save_vector;
    static int active;
    sym_struct_t sym_info;
    scf_irreps_t irreps;
    scf_struct_t scf_info;
    int node_timings;
    centers_t centers;
    dmt_matrix SCF_VEC;
    Result _scf;
    dmt_matrix FOCK;
    dmt_matrix FOCKO;
    FILE* outfile;

    RefSCDimension _basisdim;
    Resultdouble _exchange_energy;
    ResultRefSCMatrix _eigenvectors;
    int _nocc;
    void compute();
    void init();
  public:
    MPSCF(KeyVal&);
    MPSCF(StateIn&);
    virtual ~MPSCF();
    void save_data_state(StateOut&);

    void print(SCostream& =SCostream::cout);
    int do_exchange_energy(int);
    int do_eigenvectors(int);
    RefSCMatrix eigenvectors();
    double exchange_energy();

    double occupation(int vectornum);
};

#endif
