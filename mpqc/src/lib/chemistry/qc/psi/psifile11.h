#ifdef __GNUG__
#pragma interface
#endif

#ifndef _chemistry_qc_psi_file11_h
#define _chemistry_qc_psi_file11_h

#include <fstream>
#include <string>
#include<util/ref/ref.h>
#include<chemistry/molecule/molecule.h>
#include<chemistry/qc/basis/basis.h>

namespace sc {

class PsiExEnv;

///////////////////////////////////////////////////
/// PsiFile11 is a Psi gradient file

class PsiFile11: public RefCount {

  std::string filename_;
  std::fstream file_;

  // No default constructor
  PsiFile11() {};

  void skip_lines(int n);
  void skip_entry();
  void rewind();

  public:
    PsiFile11(const std::string& name);
    ~PsiFile11();

    void open();
    void close();
    void remove();
    int get_natom(int entry);
    double get_energy(int entry);
    double get_coord(int entry, int atom, int xyz);
    double get_grad(int entry, int atom, int xyz);
};

}

#endif
