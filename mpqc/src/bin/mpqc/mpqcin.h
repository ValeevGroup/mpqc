
#ifndef _mpqcin_h
#define _mpqcin_h

#include <vector>
#include <iostream>

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

class MPQCIn {
    MPQCInFlexLexer *lexer_;
    Ref<Molecule> mol_;
    MPQCInDatum<int> gradient_;
    MPQCInDatum<int> frequencies_;
    MPQCInDatum<int> optimize_;
    MPQCInDatum<int> mult_;
    MPQCInDatum<int> redund_coor_;
    MPQCInDatum<int> opt_type_;
    MPQCInDatum<int> restart_;
    MPQCInDatum<int> checkpoint_;
    MPQCInDatum<int> charge_;
    MPQCInDatum<int> atom_charge_;
    MPQCInDatum<int> molecule_bohr_;
    MPQCInDatum<char *> basis_;
    MPQCInDatum<char *> auxbasis_;
    MPQCInDatum<char *> dfbasis_;
    MPQCInDatum<char *> method_;
    MPQCInDatum<char *> accuracy_;
    MPQCInDatum<char *> lindep_;
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
    MPQCInDatum<char *> symmetry_;
    MPQCInDatum<char *> memory_;
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
        Cints,
        Libint2,
        Invalid
    };
    static std::string to_string(IntegralsFactoryType ifactory);
    static const char* guess_basis(IntegralsFactoryType ifactory);
    void write_energy_object(std::ostream&, const char *keyword,
                             const char *method,
                             const char *basis, int coor,
                             IntegralsFactoryType& ifactory);
    void write_basis_object(std::ostream&, const char *keyword,
                            const char *basis,
                            bool split = false,
                            bool uncontract = false);
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
    void set_basis(char *);
    void set_auxbasis(char *);
    void set_dfbasis(char *);
    void set_multiplicity(char *);
    void set_memory(char *);
    void set_accuracy(char *);
    void set_lindep(char *);
    void set_optimize(int);
    void set_opt_type(int);
    void set_atom_charge(char *);
    void set_molecule_unit(char *);
    void set_symmetry(char *);
    void set_redund_coor(int);
    void set_gradient(int);
    void set_frequencies(int);
    void set_restart(int);
    void set_checkpoint(int);
    void set_molecule_bohr(int);
    void set_debug(char *);
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

}

#endif
