
#ifndef _chemistry_qc_basis_h
#define _chemistry_qc_basis_h

#include <util/keyval/keyval.h>

class BasisFileSet {
  private:
    char *dir_[2];
    char **basissets_;
    int nbasissets_;
  public:
    BasisFileSet(const RefKeyVal&);
    ~BasisFileSet();
    RefKeyVal keyval(const RefKeyVal&, const char *name);
};

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
