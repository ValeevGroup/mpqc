
#ifndef _mpqcin_h
#define _mpqcin_h

#include <iostream.h>

#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>

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

class MPQCInFlexLexer;
class MPQCIn {
    MPQCInFlexLexer *lexer_;
    RefMolecule mol_;
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
    MPQCInDatum<char *> method_;
    MPQCInDatum<char *> method_xc_;
    MPQCInDatum<char *> method_grid_;
    MPQCInDatum<char *> symmetry_;
    MPQCInDatum<Arrayint *> alpha_;
    MPQCInDatum<Arrayint *> beta_;
    MPQCInDatum<Arrayint *> docc_;
    MPQCInDatum<Arrayint *> socc_;
    MPQCInDatum<Arrayint *> frozen_docc_;
    MPQCInDatum<Arrayint *> frozen_uocc_;

    int nirrep_;

    void write_energy_object(ostream&, const char *keyword, const char *method,
                             const char *basis, int coor);
    void write_basis_object(ostream&, const char *keyword, const char *basis);
    void write_vector(ostream &ostrs,
                      const char *keyvalname,
                      const char *name,
                      MPQCInDatum<Arrayint *>&vec,
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
    void set_multiplicity(char *);
    void set_optimize(int);
    void set_opt_type(int);
    void set_atom_charge(char *);
    void set_molecule_unit(char *);
    void set_method_xc(char *);
    void set_method_grid(char *);
    void set_symmetry(char *);
    void set_redund_coor(int);
    void set_gradient(int);
    void set_frequencies(int);
    void set_restart(int);
    void set_checkpoint(int);
    void set_molecule_bohr(int);
    void set_docc(Arrayint *);
    void set_socc(Arrayint *);
    void set_alpha(Arrayint *);
    void set_beta(Arrayint *);
    void set_frozen_docc(Arrayint *);
    void set_frozen_uocc(Arrayint *);
    Arrayint *make_nnivec(Arrayint *, char *);

    static int checking() { return checking_; }
};

#endif
