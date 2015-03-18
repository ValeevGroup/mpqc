
#ifndef _mpqcin_h
#define _mpqcin_h

#include <vector>
#include <iostream>
#include <string.h>

#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>

class MPQCInFlexLexer;

namespace sc {

class IPV2;

template <class T>
class MPQCInDatum {
    int set_;
    T val_;
  public:
    MPQCInDatum(const T&v): val_(v), set_(0) {}
    const T &operator =(const T&v) { set_ = 1; val_ = v; return val_; }
    void reset(const T &val) { set_ = 0; val_ = val; }
    int set() const { return set_; }
    T val() const { return val_; }
};

/// @addtogroup Init
/// @{

/// Converts MPQC simple input to object-oriented input
class MPQCIn {

  public:
    struct Basis {
        Basis() : name(0), uc(0), split(0), puream(0) {}
        Basis(const char* n, bool u, bool s, bool p) : name(0), uc(u), split(s), puream(p) {
          if (n) name = strdup(n);
          if (uc.val()) split = false;  // uncontraction implies splitting
        }
        Basis(const Basis& other) : name(0), uc(0), split(0), puream(0) {
          if (other.name.set()) name = strdup(other.name.val());
          if (other.uc.set()) uc = other.uc;
          if (other.split.set()) split = other.split;
          if (other.puream.set()) puream = other.puream;
        }
        ~Basis() { if (name.val()) free(name.val()); }

        void set_name(char* c) { name = c; }
        void set_uc(bool b) { uc = b; }
        void set_split(bool b) { split = b; }
        void set_puream(bool b) { puream = b; }
        void write(std::ostream &ostrs,
                   const char *keyword) const;

        bool operator==(const Basis& other) {
          return strcmp(name.val(), other.name.val()) == 0 &&
                 uc.val() == other.uc.val() &&
                 split.val() == other.split.val() &&
                 puream.val() == other.puream.val();
        }
        bool operator!=(const Basis& other) {
          return ! (*this == other);
        }

        MPQCInDatum<char *> name; // name
        MPQCInDatum<int> uc;      // force uncontracted?
        MPQCInDatum<int> split;   // force split?
        MPQCInDatum<int> puream;  // force puream?
    };

  private:
    MPQCInFlexLexer *lexer_;
    Ref<Molecule> mol_;
    MPQCInDatum<int> gradient_;
    MPQCInDatum<int> frequencies_;
    MPQCInDatum<int> optimize_;
    MPQCInDatum<int> mult_;
    MPQCInDatum<int> restart_;
    MPQCInDatum<int> checkpoint_;
    MPQCInDatum<int> precise_findif_;
    MPQCInDatum<int> charge_;
    MPQCInDatum<int> atom_charge_;
    MPQCInDatum<int> molecule_bohr_;
    Basis basis_;
    Basis auxbasis_;
    Basis dfbasis_;
    MPQCInDatum<char *> method_;
    MPQCInDatum<char *> accuracy_;
    MPQCInDatum<char *> lindep_;
    // options for optimize
    MPQCInDatum<int> redund_coor_;
    MPQCInDatum<int> opt_type_;
    MPQCInDatum<char *> opt_convergence_;
    // options for frequencies
    MPQCInDatum<char *> freq_accuracy_;
    // options for SCF
    MPQCInDatum<char *> scf_maxiter_;
    // options for DFT methods
    MPQCInDatum<char *> dftmethod_xc_;
    MPQCInDatum<char *> dftmethod_grid_;
    // options for R12 methods
    MPQCInDatum<char *> r12method_f12_;
    MPQCInDatum<char *> r12method_app_;
    MPQCInDatum<char *> r12method_ri_;
    MPQCInDatum<char *> r12method_ansatz_;
    MPQCInDatum<char *> pccsd_alpha_;
    MPQCInDatum<char *> pccsd_beta_;
    MPQCInDatum<char *> pccsd_gamma_;
    MPQCInDatum<char *> symmetry_;
    MPQCInDatum<char *> memory_;
    MPQCInDatum<char *> tmpstore_;
    MPQCInDatum<char *> tmpdir_;
    MPQCInDatum<char *> debug_;
    MPQCInDatum<std::vector<int> *> alpha_;
    MPQCInDatum<std::vector<int> *> beta_;
    MPQCInDatum<std::vector<int> *> docc_;
    MPQCInDatum<std::vector<int> *> socc_;
    MPQCInDatum<std::vector<int> *> frozen_docc_;
    MPQCInDatum<std::vector<int> *> frozen_uocc_;

    int nirrep_;

    enum IntegralsFactoryType {
        IntV3,
        Libint2,
        Invalid
    };
    static std::string to_string(IntegralsFactoryType ifactory);
    static Basis guess_basis(IntegralsFactoryType ifactory);
    static bool psi_method(const char*);
    static bool r12_method(const char*);

    /// infer defaults for missing parameters
    void infer_defaults();

    void write_energy_object(std::ostream&, const char *keyword,
                             const char *method,
                             Basis const* basis, int coor,
                             IntegralsFactoryType& ifactory);
    void write_vector(std::ostream &ostrs,
                      const char *keyvalname,
                      const char *name,
                      MPQCInDatum<std::vector<int> *>&vec,
                      int require_nirrep);

    static int checking_;
  public:
    MPQCIn();
    ~MPQCIn();

    char *parse_string(const char *s);
    int check_string(const char *s);

    int ylex();
    int yparse();
    void error(const char* s);
    void error2(const char* s, const char* s2);
    void yerror(const char* s);
    void yerror2(const char* s, const char *);

    void begin_molecule();
    void end_molecule();
    void add_atom(char *, char *, char *, char *);
    void set_charge(char *);
    void set_method(char *);
    void set_multiplicity(char *);
    void set_memory(char *);
    void set_tmpstore(char *);
    void set_tmpdir(char *);
    void set_accuracy(char *);
    void set_lindep(char *);
    void set_optimize(int);
    void set_opt_type(int);
    void set_opt_convergence(char *);
    void set_atom_charge(char *);
    void set_molecule_unit(char *);
    void set_symmetry(char *);
    void set_redund_coor(int);
    void set_gradient(int);
    void set_frequencies(int);
    void set_freq_accuracy(char *);
    void set_restart(int);
    void set_checkpoint(int);
    void set_precise_findif(int);
    void set_molecule_bohr(int);
    void set_debug(char *);
    void set_pccsd(char *, char *, char *);
    void set_docc(std::vector<int> *);
    void set_socc(std::vector<int> *);
    void set_alpha(std::vector<int> *);
    void set_beta(std::vector<int> *);
    void set_frozen_docc(std::vector<int> *);
    void set_frozen_uocc(std::vector<int> *);
    std::vector<int> *make_nnivec(std::vector<int> *, char *);

    void set_scf_maxiter(char *);

    void set_dftmethod_xc(char *);
    void set_dftmethod_grid(char *);
    void set_r12method_f12(char *);
    void set_r12method_app(char *);
    void set_r12method_ri(char *);
    void set_r12method_ansatz(char *);

    static int checking() { return checking_; }
};

/// @}
// end of addtogroup Init

}

#endif
