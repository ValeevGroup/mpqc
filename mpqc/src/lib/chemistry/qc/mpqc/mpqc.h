
#ifndef _chemistry_qc_mpqc_h
#define _chemistry_qc_mpqc_h

#ifdef __GNUC__
#pragma interface
#endif

extern "C" {
#include <stdio.h>
#include <math/dmt/libdmt.h>
#include <math/array/math_lib.h>
};
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include <chemistry/qc/dmtscf/scf_dmt.h>
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/molecule/energy.h>
#include <math/scmat/matrix.h>
#include <math/topology/point.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/solvent/bem.h>
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
    int node_timings;
    sym_struct_t sym_info;
    scf_struct_t scf_info;
    centers_t centers;
    centers_t oldcenters;
    dmt_matrix Scf_Vec;
    Result _scf;
    dmt_matrix Fock;
    dmt_matrix FockO;
    FILE* outfile;

    //RefSCDimension _basisdim;
    Resultdouble _exchange_energy;
    ResultRefSCMatrix _eigenvectors;
    int _ndocc;
    int _nsocc;
    int nproc;
    int me;
    int throttle;
    int sync_loop;
    void compute();
    void init(int =0);

    void scfvec_to_eigenvectors();

    RefBEMSolvent solvent_;
  public:
    MPSCF(const RefKeyVal&);
    MPSCF(StateIn&);
    virtual ~MPSCF();
    void save_data_state(StateOut&);

    void print(ostream& =cout);
    int do_exchange_energy(int);
    int do_eigenvectors(int);
    RefSCMatrix eigenvectors();
    double exchange_energy();

    double occupation(int vectornum);

    int value_implemented();
    int gradient_implemented();
};

#endif
