
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _chemistry_qc_basis_dercent_h
#define _chemistry_qc_basis_dercent_h

#include <chemistry/qc/basis/basis.h>

//. \clsnm{DerivCenters} keeps track the centers that
//derivatives are taken with respect to.
class DerivCenters {
  private:
    int center_[4];
    int atom_[4];
    int ncenter_;
    int omitted_center_;
    int omitted_atom_;
  public:
    //. These are used by the \clsnm{Integral} specializations to
    //initializes the DerivCenters structure.
    DerivCenters();
    void clear();
    void add_center(int center, const RefGaussianBasisSet &, int shell);
    void add_omitted(int center, const RefGaussianBasisSet &, int shell);
    void add_center(int center, int atom);
    void add_omitted(int center, int atom);

    //. The number of unique centers minus one.
    int n() const { return ncenter_; }
    //. The center number.
    int center(int i) const { return center_[i]; }
    //. The atom number.
    int atom(int i) const { return atom_[i]; }
    //. The center that is omitted from the integral buffer.
    int omitted_center() const { return omitted_center_; }
    //. Returns 1 if there is an omitted center.
    int has_omitted_center() const { return omitted_center_ >= 0; }
    //. The atom that is omitted from the integral buffer.
    int omitted_atom() const { return omitted_atom_; }
};

#endif
