/*
** PSI Input Class
**
** This helper class will set up input decks for the PSI suite of
** ab initio quantum chemistry programs. 
**
** David Sherrill & Justin Fermann
** Center for Computational Quantum Chemistry, University of Georgia
**
*/

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _chemistry_qc_psi_input_h
#define _chemistry_qc_psi_input_h

using namespace std;

#include <fstream>
#include <string>
#include<util/ref/ref.h>
#include<chemistry/molecule/molecule.h>
#include<chemistry/qc/basis/basis.h>

namespace sc {

class PsiExEnv;
class CorrelationTable;

///////////////////////////////////////////////////
/// PsiInput is a Psi input file

class PsiInput: public RefCount {

  string filename_;
  std::ofstream file_;

  int indentation_;
  
  // No default constructor
  PsiInput() {};

  public:
    PsiInput(const string& name);
    ~PsiInput();
    void open();
    void close();
    void print(std::ostream&);

    void begin_section(const char * s);
    void end_section();
    void write_indent();
    void incindent(int);
    void decindent(int);
    void write_comment(const char *);
    void write_keyword(const char *, const char *);
    void write_keyword(const char *, bool);
    void write_keyword(const char *, int);
    void write_keyword(const char *, double);
    template <typename T> void write_keyword_array(const char *, int, const std::vector<T>&);
    void write_keyword_array(const char *, int, int *);
    void write_keyword_array(const char *, int, double *);
    void write_string(const char *);
    void write_key_wq(const char *, const char *);

    /// Construct the "basis" keyword for input. All functions with angular momentum >= 1 must be Cartesian or all must be sph. harm.
    void write_basis(const Ref<GaussianBasisSet>&);
    /// Write basis sets explicitly
    void write_basis_sets(const Ref<GaussianBasisSet>&);
    void write_geom(const Ref<Molecule>&);
    
    void write_defaults(const Ref<PsiExEnv>&, const char *dertype);
};

}

#endif
