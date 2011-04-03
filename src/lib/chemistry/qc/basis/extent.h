//
// extent.h
//

#ifndef _chemistry_qc_basis_extent_h
#define _chemistry_qc_basis_extent_h

#include <vector>

#include <float.h>
#include <chemistry/qc/basis/basis.h>

namespace sc {

struct ExtentData {
    int shell;
    double bound;
    ExtentData() {}
    ExtentData(int s, double b): shell(s), bound(b) {}
};

class ShellExtent: public RefCount {
    double lower_[3];
    double resolution_;
    int n_[3];
    std::vector<ExtentData> *contributing_shells_;
    std::vector<ExtentData> null_;

    std::vector<ExtentData> &data(int *b);
    double distance(double loc, int axis, int origin, int point);
    std::vector<ExtentData> &data(int x, int y, int z);
  public:
    ShellExtent();
    ~ShellExtent();
    void init(const Ref<GaussianBasisSet>&,
              double resolution = 1.0, double tolerance = DBL_EPSILON);
    /** Returns the shells that are nonzero at coordinates x, y, z.
        The shells numbers are in ascending order. */
    const std::vector<ExtentData> &contributing_shells(int x, int y, int z)
        { return data(x,y,z); }
    const std::vector<ExtentData> &contributing_shells(double x, double y, double z);
    void print(std::ostream &o = ExEnv::out0());
    const int *n() const { return n_; }
    int n(int ixyz) const { return n_[ixyz]; }
    double lower(int ixyz) const { return lower_[ixyz]; }
    double upper(int ixyz) const { return resolution_*n_[ixyz] + lower_[ixyz]; }
    double resolution() const { return resolution_; }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
