
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _chemistry_qc_molfreq_h
#define _chemistry_qc_molfreq_h

#include <iostream.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/energy.h>

class MolecularFrequencies: public SavableState {
#   define CLASSNAME MolecularFrequencies
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefMolecularEnergy mole_;
    RefSCMatrixKit kit_;
    RefMolecule mol_;
    PointGroup original_point_group_;
    RefSCVector original_geometry_;
    // the cartesian displacement size in bohr
    double disp_;
    // displacements for each irrep
    int ndisp_;
    int nirrep_;
    RefSCMatrix *displacements_;
    RefSCVector *gradients_;

    // the number of external degrees of freedom
    int nexternal_;

    // the number of frequencies per irrep
    int *nfreq_;
    // the frequencies for each irrep
    double **freq_;

    RefSCDimension d3natom_;
    void get_disp(int disp, int &irrep, int &index, double &coef);
    void do_freq_for_irrep(int irrep,
                           const RefDiagSCMatrix &m,
                           const RefSymmSCMatrix &dhessian,
                           const RefSymmSCMatrix &xhessian);
    int debug_;
  public:
    MolecularFrequencies(const RefKeyVal &);
    MolecularFrequencies(StateIn &);
    ~MolecularFrequencies();
    void save_data_state(StateOut&);

    void compute_displacements();
    void compute_frequencies_from_gradients();
    int ndisplace() const;
    void displace(int disp);
    void set_gradient(int disp, const RefSCVector &grad);

    void thermochemistry(int degeneracy, double temp=298.15, double pres=1.0);

    RefSCMatrixKit matrixkit() { return kit_; }
    RefSCDimension d3natom() { return d3natom_; }
};

SavableState_REF_dec(MolecularFrequencies);

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
