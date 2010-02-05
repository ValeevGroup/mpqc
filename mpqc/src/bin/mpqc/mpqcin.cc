
#include <scconfig.h>
#include <sstream>
#include <stdlib.h>
#include <stdexcept>

#include <util/misc/formio.h>
#include <util/class/scexception.h>

using namespace std;
using namespace sc;

#undef yyFlexLexer
#define yyFlexLexer MPQCInFlexLexer
#include <FlexLexer.h>

#include "mpqcin.h"
#include "parse.h"

int MPQCIn::checking_ = 0;

MPQCIn::MPQCIn():
  nirrep_(0),
  mult_(1),
  charge_(0),
  basis_(0),
  auxbasis_(0),
  dfbasis_(0),
  method_(0),
  optimize_(0),
  gradient_(0),
  frequencies_(0),
  accuracy_(0),
  lindep_(0),
  opt_type_(T_INTERNAL),
  redund_coor_(0),
  restart_(0),
  checkpoint_(1),
  atom_charge_(0),
  symmetry_(0),
  memory_(0),
  molecule_bohr_(0),
  alpha_(0),
  beta_(0),  
  docc_(0),
  socc_(0),  
  frozen_docc_(0),
  frozen_uocc_(0),
  debug_(0),
  scf_maxiter_(0),
  dftmethod_xc_(0),
  dftmethod_grid_(0),
  r12method_f12_(0),
  r12method_app_(0),
  r12method_ri_(0),
  r12method_ansatz_(0)
{
  lexer_ = new MPQCInFlexLexer;
}

MPQCIn::~MPQCIn()
{
  delete lexer_;
  if (basis_.val()) free(basis_.val());
  if (auxbasis_.val()) free(auxbasis_.val());
  if (method_.val()) free(method_.val());
  if (symmetry_.val()) free(symmetry_.val());
  if (memory_.val()) free(memory_.val());
  if (alpha_.val()) free(alpha_.val());
  if (beta_.val()) free(beta_.val());
  if (docc_.val()) free(docc_.val());
  if (socc_.val()) free(socc_.val());
  if (accuracy_.val()) free(accuracy_.val());
  if (lindep_.val()) free(lindep_.val());
  if (debug_.val()) free(debug_.val());
  if (frozen_docc_.val()) free(frozen_docc_.val());
  if (frozen_uocc_.val()) free(frozen_uocc_.val());
  if (scf_maxiter_.val()) free(scf_maxiter_.val());
  if (dftmethod_xc_.val()) free(dftmethod_xc_.val());
  if (dftmethod_grid_.val()) free(dftmethod_grid_.val());
  if (r12method_f12_.val()) free(r12method_f12_.val());
  if (r12method_app_.val()) free(r12method_app_.val());
  if (r12method_ansatz_.val()) free(r12method_ansatz_.val());
  if (r12method_ri_.val()) free(r12method_ri_.val());
}

void
MPQCIn::error(const char* s)
{
  ExEnv::outn() << ExEnv::program_name()
               << ": error: " << s
               << endl;
  abort();
}

void
MPQCIn::error2(const char* s, const char *s2)
{
  ExEnv::outn() << ExEnv::program_name()
               << ": error: " << s << "\"" << s2 << "\""
               << endl;
  abort();
}

void
MPQCIn::yerror(const char* s)
{
  ExEnv::outn() << ExEnv::program_name()
               << ": " << s
               << " at line " << lexer_->lineno()+1
               << endl;
  abort();
}

void
MPQCIn::yerror2(const char* s, const char *s2)
{
  ExEnv::outn() << ExEnv::program_name()
               << ": " << s
               << " \"" << s2 << "\" at line " << lexer_->lineno()+1
               << endl;
  abort();
}

int
MPQCIn::ylex()
{
  return lexer_->yylex();
}

void
MPQCIn::begin_molecule()
{
  if (mol_.nonnull()) {
      ExEnv::outn() << ExEnv::program_name()
                   << ": error: second molecule given at line "
                   << lexer_->lineno()+1
                   << endl;
      abort();
    }
  mol_ = new Molecule;
}

void
MPQCIn::end_molecule()
{
  //double symtol = 1e-4;
  //mol_->set_point_group(mol_->highest_point_group(symtol), symtol*10.0);
}

void
MPQCIn::add_atom(char *sym, char *xs, char *ys, char *zs)
{
  int Z = mol_->atominfo()->string_to_Z(sym, 0);
  if (Z == 0) yerror2("bad element", sym);
  free(sym);

  char *xse;
  double x = strtod(xs,&xse);
  if (xse == xs) yerror2("bad x coordinate", xs);
  free(xs);

  char *yse;
  double y = strtod(ys,&yse);
  if (yse == ys) yerror2("bad y coordinate", ys);
  free(ys);

  char *zse;
  double z = strtod(zs,&zse);
  if (zse == zs) yerror2("bad z coordinate", zs);
  free(zs);

  mol_->add_atom(Z, x, y, z, "", 0, atom_charge_.set(), atom_charge_.val());
  atom_charge_.reset(0);
}

void
MPQCIn::set_charge(char *cs)
{
  char *cse;
  int c = strtol(cs,&cse,10);
  if (cse == cs) yerror2("bad charge", cs);
  charge_ = c;
  free(cs);
}

void
MPQCIn::set_atom_charge(char *cs)
{
  char *cse;
  int c = strtol(cs,&cse,10);
  if (cse == cs) yerror2("bad atom charge", cs);
  atom_charge_ = c;
  free(cs);
}

void
MPQCIn::set_method(char *m)
{
  method_ = m;
}

void
MPQCIn::set_molecule_bohr(int i)
{
  molecule_bohr_ = i;
}

void
MPQCIn::set_basis(char *b)
{
  basis_ = b;
}

void
MPQCIn::set_auxbasis(char *b)
{
  auxbasis_ = b;
}

void
MPQCIn::set_dfbasis(char *b)
{
  dfbasis_ = b;
}

void
MPQCIn::set_symmetry(char *s)
{
  symmetry_ = s;
  if (strcmp(s,"auto")) {
      Ref<PointGroup> p = new PointGroup(s);
      nirrep_ = p->char_table().nirrep();
    }
}

void
MPQCIn::set_memory(char *s)
{
  memory_ = s;
}

void
MPQCIn::set_multiplicity(char *ms)
{
  char *mse;
  int m = strtol(ms,&mse,10);
  if (mse == ms || m <= 0) yerror2("bad multiplicity", ms);
  mult_ = m;
  free(ms);
}

void
MPQCIn::set_optimize(int i)
{
  optimize_ = i;
}

void
MPQCIn::set_gradient(int i)
{
  gradient_ = i;
}

void
MPQCIn::set_frequencies(int i)
{
  frequencies_ = i;
}

void
MPQCIn::set_restart(int i)
{
  restart_ = i;
}

void
MPQCIn::set_checkpoint(int i)
{
  checkpoint_ = i;
}

void
MPQCIn::set_accuracy(char* c)
{
  accuracy_ = c;
}

void
MPQCIn::set_lindep(char* c)
{
  lindep_ = c;
}

void
MPQCIn::set_redund_coor(int i)
{
  redund_coor_ = i;
}

void
MPQCIn::set_opt_type(int i)
{
  opt_type_ = i;
}

void
MPQCIn::set_debug(char* c)
{
  debug_ = c;
}

void
MPQCIn::set_docc(std::vector<int> *a)
{
  docc_ = a;
}

void
MPQCIn::set_socc(std::vector<int> *a)
{
  socc_ = a;
}

void
MPQCIn::set_alpha(std::vector<int> *a)
{
  alpha_ = a;
}

void
MPQCIn::set_beta(std::vector<int> *a)
{
  beta_ = a;
}

void
MPQCIn::set_frozen_docc(std::vector<int> *a)
{
  frozen_docc_ = a;
}

void
MPQCIn::set_frozen_uocc(std::vector<int> *a)
{
  frozen_uocc_ = a;
}

void
MPQCIn::set_scf_maxiter(char* c)
{
  scf_maxiter_ = c;
}

void
MPQCIn::set_dftmethod_xc(char *m)
{
  dftmethod_xc_ = m;
}

void
MPQCIn::set_dftmethod_grid(char *m)
{
  dftmethod_grid_ = m;
}

void
MPQCIn::set_r12method_f12(char *m)
{
  r12method_f12_ = m;
}

void
MPQCIn::set_r12method_app(char *m)
{
  r12method_app_ = m;
}

void
MPQCIn::set_r12method_ri(char *m)
{
  r12method_ri_ = m;
}

void
MPQCIn::set_r12method_ansatz(char* m)
{
  r12method_ansatz_ = m;
}

std::vector<int> *
MPQCIn::make_nnivec(std::vector<int> *a, char *ms)
{
  if (ms == 0) return new std::vector<int>;

  char *mse;
  int m = strtol(ms,&mse,10);
  if (mse == ms || m < 0) yerror2("bad positive integer", ms);
  free(ms);

  if (a == 0) a = new std::vector<int>;
  a->push_back(m);
  return a;
}

int
MPQCIn::check_string(const char *s)
{
  checking_ = 1;
  istringstream in(s);
  lexer_->switch_streams(&in, &ExEnv::outn());
  int token;
  while ((token = ylex())) {
      if (token == T_OO_INPUT_KEYWORD) return 0;
    }
  checking_ = 0;
  return 1;
}

char *
MPQCIn::parse_string(const char *s)
{
  // read in easy input
  istringstream in(s);
  lexer_->switch_streams(&in, &ExEnv::outn());
  yparse();

  // form the oo input
  ostringstream ostrs;
  SCFormIO::init_ostream(ostrs);
  ostrs << decindent;
  if (mol_.null()) error("no molecule given");
  if (symmetry_.set() && strcmp(symmetry_.val(),"auto") != 0) {
      mol_->symmetrize(new PointGroup(symmetry_.val()));
    }
  ostrs << indent << "molecule<Molecule>: (" << endl;
  ostrs << incindent;
  ostrs << indent << "symmetry = "
        << (symmetry_.set()?symmetry_.val():"auto") << endl;
  ostrs << indent << "unit = \""
        << (molecule_bohr_.val()?"bohr":"angstrom")
        << "\"" << endl;
  mol_->print_parsedkeyval(ostrs, 0, 0, 0);
  ostrs << decindent;
  ostrs << indent << ")" << endl;
  write_basis_object(ostrs, "basis", basis_.val());
  ostrs << indent << "mpqc: (" << endl;
  ostrs << incindent;
  ostrs << indent << "do_gradient = " << gradient_.val() << endl;
  ostrs << indent << "optimize = " << optimize_.val() << endl;
  ostrs << indent << "restart = " << restart_.val() << endl;
  ostrs << indent << "checkpoint = " << checkpoint_.val() << endl;
  ostrs << indent << "savestate = " << checkpoint_.val() << endl;
  IntegralsFactoryType ifactory = IntV3;    // this is set by write_energy_object()
  write_energy_object(ostrs, "mole", method_.val(), 0, optimize_.val(),
                      ifactory);
  ostrs << indent << "integrals<" << to_string(ifactory) << ">: ()" << std::endl;
  if (optimize_.val()) {
      const char *coortype = "SymmMolecularCoor";
      if (opt_type_.val() == T_CARTESIAN) coortype = "CartMolecularCoor";
      else if (redund_coor_.val()) coortype = "RedundMolecularCoor";
      ostrs << indent << "coor<" << coortype << ">: (" << endl;
      ostrs << indent << "  molecule = $:molecule" << endl;
      if (opt_type_.val() == T_INTERNAL) {
          ostrs << indent << "  generator<IntCoorGen>: (" << endl;
          ostrs << indent << "    molecule = $:molecule" << endl;
          ostrs << indent << "  )" << endl;
        }
      ostrs << indent << ")" << endl;
      ostrs << indent << "opt<QNewtonOpt>: (" << endl;
      ostrs << indent << "  function = $:mpqc:mole" << endl;
      ostrs << indent << "  update<BFGSUpdate>: ()" << endl;
      ostrs << indent << "  convergence<MolEnergyConvergence>: (" << endl;
      ostrs << indent << "    cartesian = yes" << endl;
      ostrs << indent << "    energy = $:mpqc:mole" << endl;
      ostrs << indent << "  )" << endl;
      ostrs << indent << ")" << endl;
    }

  if (frequencies_.val()) {
      ostrs << indent << "freq<MolecularFrequencies>: (" << endl;
      ostrs << indent << "  molecule = $:molecule" << endl;
      ostrs << indent << ")" << endl;
    }

  ostrs << decindent;
  ostrs << indent << ")" << endl;
  ostrs << ends;

  int n = 1 + strlen(ostrs.str().c_str());
  char *in_char_array = strcpy(new char[n],ostrs.str().c_str());
  return in_char_array;
}

void
MPQCIn::write_vector(ostream &ostrs,
                     const char *keyvalname,
                     const char *name, MPQCInDatum<std::vector<int> *>&vec,
                     int require_nirrep)
{
  if (vec.set()) {
      ostrs << indent << keyvalname << " = ";
      if (!require_nirrep && vec.val()->size() == 1) {
          ostrs << (*vec.val())[0] << endl;
        }
      else if (nirrep_ && vec.val()->size() == nirrep_) {
          ostrs << "[";
              for (int i=0; i<nirrep_; i++) {
                  ostrs << " " << (*vec.val())[i];
                }
          ostrs << "]" << endl;
        }
      else {
          if (!require_nirrep) error2("need 1 or n_irrep elements in ", name);
          else {
              error2("need n_irrep (must give symmetry) elements in ", name);
            }
        }
    }
}

namespace {

  struct R12TechDescr {

      enum CorrFactorType {
          NullCorrFactor,
          R12CorrFactor,
          GTGCorrFactor
      };

      R12TechDescr(const std::string& app,  // approximation
                   const std::string& an,   // ansatz
                   const std::string& cf,   // correlation factor
                   const std::string& r     // RI approach
                  ) :
                    approx(app),
                    ansatz(an),
                    corrfactor(cf),
                    ri(r) { }
      static R12TechDescr* default_instance() {
        return new R12TechDescr("C", "diag", "stg-6g[1.3]", "cabs+");
      }

      void validate() {
        if (approx != "C" &&
            approx != "A'" &&
            approx != "A''" &&
            approx != "B") {
          throw sc::InputError("invalid R12 approximation",__FILE__,__LINE__);
        }
      }
      void write(ostream& os) {
        validate();
        os << indent << "stdapprox = " << approx << endl;
        os << indent << "abs_method = " << ri << endl;
        os << indent << "ansatz<R12Ansatz>: (" << endl << incindent;
        if (ansatz == "diag" || ansatz == "ijij") {
          os << indent << "diag = true" << endl;
        }
        else if (ansatz == "ijpq") {
          os << indent << "orbital_product_GG = pq" << endl;
        }
        os << decindent << indent << ")" << endl;
        // parse correlation factor
        std::string cf_label, cf_param;
        parse_corrfactor(corrfactor, cf_label, cf_param);
        os << indent << "corr_factor = " << cf_label << endl;
        if (!cf_param.empty()) os << indent << "corr_param = " << cf_param << endl;
      }

      CorrFactorType corrfactor_type() {
        std::string cf_label, cf_param;
        parse_corrfactor(corrfactor, cf_label, cf_param);
        CorrFactorType result;
        if (cf_label == "none" || cf_label == "NONE")
          result = NullCorrFactor;
        else if (cf_label == "r12" || cf_label == "R12")
          result = R12CorrFactor;
        else if (cf_label.find_first_of("stg-") == 0 ||
                 cf_label.find_first_of("STG-") == 0)
          result = GTGCorrFactor;
        else
          throw InputError("invalid correlation factor given to an R12 method");
        return result;
      }

      static void parse_corrfactor(const std::string& cf,
                                   std::string& cf_label,
                                   std::string& cf_param) {
        const std::string::size_type lpar_pos = cf.find_first_of('[');
        if (lpar_pos != std::string::npos) { // parameter given in brakets? must be an STG
          cf_label = cf.substr(0, lpar_pos);
          const std::string::size_type rpar_pos = cf.find_first_of(']');
          if (rpar_pos == std::string::npos) { // oops
            throw InputError("parentheses mismatch", __FILE__, __LINE__);
          }
          cf_param = cf.substr(lpar_pos+1, rpar_pos-lpar_pos-1);
        }
        else {
          cf_label = cf;
          cf_param.clear();
        }
      }

      std::string approx;
      std::string ansatz;
      std::string corrfactor;
      std::string ri;
  };
}

void
MPQCIn::write_energy_object(ostream &ostrs,
                            const char *keyword,
                            const char *method,
                            const char *basis,
                            int coor, IntegralsFactoryType& ifactory)
{
  int nelectron = int(mol_->nuclear_charge()+1e-6) - charge_.val();
  if (nelectron < 0) {
      error("charge is impossibly large");
    }
  if (nelectron%2 == 0 && mult_.val()%2 == 0
      ||nelectron%2 == 1 && mult_.val()%2 == 1) {
      error("given multiplicity is not possible");
    }

  const char *method_object = 0;
  const char *reference_method = 0;
  const char *guess_method = method;
  bool ask_auxbasis = false;
  bool ask_dfbasis = false;
  R12TechDescr* r12descr = 0;
  int dft = 0;
  int uscf = 0;
  bool scf = false;
  ostringstream o_extra;
  SCFormIO::init_ostream(o_extra);
  o_extra << incindent;
  if (method) {
      // Hartree-Fock methods
      if (!strcmp(method, "HF")) {
          if (mult_.val() == 1) method_object = "CLHF";
          else { uscf = 1; method_object = "UHF"; }
          scf = true;
        }
      else if (!strcmp(method, "RHF")) {
          if (mult_.val() == 1) method_object = "CLHF";
          else method_object = "HSOSHF";
          scf = true;
        }
      else if (!strcmp(method, "UHF")) {
          method_object = "UHF";
          uscf = 1;
          scf = true;
        }
      // Density Functional Methods
      else if (!strcmp(method, "KS")) {
          guess_method = "HF";
          if (mult_.val() == 1) method_object = "CLKS";
          else { uscf = 1; method_object = "UKS"; }
          dft = 1;
          scf = true;
        }
      else if (!strcmp(method, "RKS")) {
          guess_method = "RHF";
          if (mult_.val() == 1) method_object = "CLKS";
          else method_object = "HSOSKS";
          dft = 1;
          scf = true;
        }
      else if (!strcmp(method, "UKS")) {
          guess_method = "UHF";
          method_object = "UKS";
          dft = 1;
          uscf = 1;
          scf = true;
        }
      // Perturbation Theory
      else if (!strcmp(method, "MP2") ||
               !strcmp(method, "UMP2") ||
               !strcmp(method, "RMP2") ) {
          guess_method = 0;
          if (mult_.val() == 1) { // closed-shell
            reference_method = "RHF";
            method_object = "MBPT2";
          }
          else { // open-shell will use MP2-R12 code
            method_object = "MBPT2_R12";
            r12descr = R12TechDescr::default_instance();
            r12descr->corrfactor = "none";
            if (!strcmp(method, "RMP2")) {
              reference_method = "HSOSHF";
            }
            else {
              reference_method = "UHF";
            }
          }
        }
      else if (!strcmp(method, "ZAPT2")) {
          guess_method = 0;
          method_object = "MBPT2";
          reference_method = "RHF";
          if (mult_.val() == 1) {
              error("ZAPT2 can only be used with multiplicity != 1: try MP2");
            }
          if (optimize_.val() || gradient_.val() || frequencies_.val()) {
              error("cannot do a gradient or optimization with ZAPT2");
            }
        }
      // MP2-R12
      else if (strncmp(method, "MP2-R12", 7) == 0 ||
               strncmp(method, "MP2-F12", 7) == 0) {
        guess_method = 0;
        ask_auxbasis = true;
        ask_dfbasis = true;
        method_object = "MBPT2_R12";
        reference_method = "HF";
        r12descr = R12TechDescr::default_instance();
      }
      else error2("invalid method: ", method);
    }
  else error("no method given");
  ostrs << indent << keyword << "<" << method_object << ">: (" << endl;
  ostrs << incindent;
  ostrs << o_extra.str();

  if (ask_auxbasis
      && auxbasis_.val() != 0
      && strcmp(auxbasis_.val(),basis_.val()) != 0)
    write_basis_object(ostrs, "aux_basis", auxbasis_.val());
  if (ask_dfbasis
      && dfbasis_.val() != 0)
    write_basis_object(ostrs, "df_basis", dfbasis_.val());

  if (r12descr) {
    if (r12method_f12_.set()) r12descr->corrfactor = r12method_f12_.val();
    if (r12method_app_.set()) r12descr->approx = r12method_app_.val();
    if (r12method_ansatz_.set()) r12descr->ansatz = r12method_ansatz_.val();
    if (r12method_ri_.set()) r12descr->ri = r12method_ri_.val();
    r12descr->write(ostrs);
    switch (r12descr->corrfactor_type()) {
      case R12TechDescr::R12CorrFactor:  ifactory = Cints; break;
      case R12TechDescr::GTGCorrFactor:  ifactory = Libint2; break;
    }
  }

  ostrs << indent << "integrals<" << to_string(ifactory) << ">: ()" << endl;
  ostrs << indent << "total_charge = " << charge_.val() << endl;
  ostrs << indent << "multiplicity = " << mult_.val() << endl;
  if (accuracy_.set()) ostrs << indent << "value_accuracy = " << accuracy_.val() << endl;
  ostrs << indent << "molecule = $:molecule" << endl;
  if (lindep_.set()) ostrs << indent << "lindep_tol = " << lindep_.val() << endl;
  if (memory_.set()) ostrs << indent << "memory = " << memory_.val() << endl;
  if (debug_.set()) ostrs << indent << "debug = " << debug_.val() << endl;
  if (scf && scf_maxiter_.set())
    ostrs << indent << "maxiter = " << scf_maxiter_.val() << endl;
  if (!strcmp(keyword, "mole") && !reference_method) {
      ostrs << indent << "print_npa = 1" << endl;
    }
  if (reference_method) {
      write_vector(ostrs, "nfzc", "frozen_docc", frozen_docc_, 0);
      write_vector(ostrs, "nfzv", "frozen_uocc", frozen_uocc_, 0);
    }
  else {
      if (uscf && (docc_.set() || socc_.set())) {
          error("cannot set docc or socc for unrestricted methods"
                " (use alpha and beta)");
        }
      else if (uscf) {
          write_vector(ostrs, "alpha", "alpha", alpha_, 1);
          write_vector(ostrs, "beta", "beta", beta_, 1);
        }
      else if (alpha_.set() || beta_.set()) {
          error("cannot set alpha or beta for restricted methods"
                " (use docc and socc)");
        }
      else {
          write_vector(ostrs, "docc", "docc", docc_, 1);
          write_vector(ostrs, "socc", "socc", socc_, 1);
        }
    }
  if (coor) ostrs << indent << "coor = $:mpqc:coor" << endl;
  if (basis) {
      write_basis_object(ostrs, "basis", basis);
    }
  else {
      ostrs << indent << "basis = $:basis" << endl;
    }
  if (dft) {
      if (dftmethod_xc_.set()) {
          ostrs << indent << "functional<StdDenFunctional>: ( name = \""
                << dftmethod_xc_.val()
                << "\" )" << endl;
        }
      else error("no exchange-correlation functional given");
      if (dftmethod_grid_.set()) {
          ostrs << indent << "integrator<RadialAngularIntegrator>: ( grid = \""
                << dftmethod_grid_.val()
                << "\" )" << endl;
        }
    }
  if (dft || (!(basis
                && !strcmp(guess_basis(ifactory),basis))
              && strcmp(guess_basis(ifactory),basis_.val())
              && guess_method)) {
      if (frequencies_.val()) {
          ostrs << indent << "keep_guess_wavefunction = 1" << endl;;
        }
      write_energy_object(ostrs, "guess_wavefunction",
                          guess_method, guess_basis(ifactory), 0, ifactory);
    }
  if (reference_method) {
      ostrs << indent << "nfzc = auto" << endl;;
      write_energy_object(ostrs, "reference",
                          reference_method, 0, 0, ifactory);
    }
  ostrs << decindent;
  ostrs << indent << ")" << endl;
}

void
MPQCIn::write_basis_object(ostream &ostrs,
                           const char *keyword,
                           const char *basis,
                           bool split,
                           bool uncontract)
{
  // validate input
  if (!basis) error("no basis given");
  if (uncontract) split = false;  // uncontraction implies splitting

  if (!split && !uncontract) {
    ostrs << indent << keyword << "<GaussianBasisSet>: (" << endl;
    ostrs << incindent;
    ostrs << indent << "molecule = $:molecule" << endl;
    ostrs << indent << "name = \"" << basis << "\"" << endl;
    ostrs << decindent;
    ostrs << indent << ")" << endl;
  }
  else {
    std::ostringstream oss;
    oss << "m" << keyword;
    const char* mkeyword = oss.str().c_str();
    this->write_basis_object(ostrs, mkeyword, basis, false, false);
    if (uncontract)
      ostrs << indent << keyword << "<UncontractedBasisSet>: (" << endl;
    if (split)
      ostrs << indent << keyword << "<SplitBasisSet>: (" << endl;
    ostrs << incindent;
    ostrs << indent << "basis = $..:" << mkeyword << endl;
    ostrs << decindent;
    ostrs << indent << ")" << endl;
  }
}

std::string
MPQCIn::to_string(IntegralsFactoryType ifactory) {
  std::string result;
  switch (ifactory) {
    case IntV3:
      result = "IntegralV3"; break;
    case Cints:
      result = "IntegralCints"; break;
    case Libint2:
      result = "IntegralLibint2"; break;
    default:
      throw std::logic_error("Invalid integral factory");
  }
  return result;
}

const char*
MPQCIn::guess_basis(IntegralsFactoryType ifactory) {
  if (ifactory == IntV3)
    return "STO-3G";
  else
    return "DZ (Dunning)";
}
