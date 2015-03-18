
#include <mpqc_config.h>
#include <sstream>
#include <stdlib.h>
#include <stdexcept>

#include <util/misc/formio.h>
#include <util/misc/scexception.h>

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
  method_(0),
  optimize_(0),
  gradient_(0),
  frequencies_(0),
  precise_findif_(0),
  accuracy_(0),
  lindep_(0),
  opt_type_(T_INTERNAL),
  redund_coor_(0),
  opt_convergence_(0),
  freq_accuracy_(0),
  restart_(0),
  checkpoint_(1),
  atom_charge_(0),
  symmetry_(0),
  memory_(0),
  tmpstore_(0),
  tmpdir_(0),
  molecule_bohr_(0),
  alpha_(0),
  beta_(0),
  docc_(0),
  socc_(0),
  frozen_docc_(0),
  frozen_uocc_(0),
  debug_(0),
  scf_maxiter_(0),
  pccsd_alpha_(0),
  pccsd_beta_(0),
  pccsd_gamma_(0),
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
  if (method_.val()) free(method_.val());
  if (symmetry_.val()) free(symmetry_.val());
  if (memory_.val()) free(memory_.val());
  if (alpha_.val()) free(alpha_.val());
  if (beta_.val()) free(beta_.val());
  if (docc_.val()) free(docc_.val());
  if (socc_.val()) free(socc_.val());
  if (accuracy_.val()) free(accuracy_.val());
  if (lindep_.val()) free(lindep_.val());
  if (opt_convergence_.val()) free(opt_convergence_.val());
  if (freq_accuracy_.val()) free(freq_accuracy_.val());
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
  if (mol_) {
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
MPQCIn::set_precise_findif(int i)
{
  precise_findif_ = i;
}

void
MPQCIn::set_accuracy(char* c)
{
  accuracy_ = c;
}

void
MPQCIn::set_tmpstore(char* c)
{
  if (strcmp(c,"mem") &&
      strcmp(c,"disk"))
    yerror2("bad tmpstore", c);
  tmpstore_ = c;
}

void
MPQCIn::set_tmpdir(char* c)
{
  tmpdir_ = c;
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
MPQCIn::set_opt_convergence(char* c)
{
  opt_convergence_ = c;
}

void
MPQCIn::set_freq_accuracy(char* c)
{
  freq_accuracy_ = c;
}

void
MPQCIn::set_debug(char* c)
{
  debug_ = c;
}

void
MPQCIn::set_pccsd(char *a, char *b, char *c)
{
  pccsd_alpha_ = a;
  pccsd_beta_ = b;
  pccsd_gamma_ = c;
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

namespace {

  struct Status {
      static bool need_psi_exenv;
  };
  bool Status::need_psi_exenv = false;

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

  infer_defaults();

  basis_.write(ostrs, "basis");

  ostrs << indent << "mpqc: (" << endl;
  ostrs << incindent;
  ostrs << indent << "do_gradient = " << gradient_.val() << endl;
  ostrs << indent << "optimize = " << optimize_.val() << endl;
  ostrs << indent << "do_freq = " << frequencies_.val() << endl;
  ostrs << indent << "restart = " << restart_.val() << endl;
  ostrs << indent << "checkpoint = " << checkpoint_.val() << endl;
  ostrs << indent << "savestate = " << checkpoint_.val() << endl;
  ostrs << indent << "precise_findif = " << precise_findif_.val() << endl;
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
      ostrs << indent << "  convergence<MolEnergyConvergence>: (" << endl;
      ostrs << indent << "    cartesian = yes" << endl;
      ostrs << indent << "    energy = $:mpqc:mole" << endl;
      if (opt_convergence_.val()) {
        ostrs << indent << "    max_disp = " << opt_convergence_.val() << endl;
        ostrs << indent << "    max_grad = " << opt_convergence_.val() << endl;
        ostrs << indent << "    graddisp = " << opt_convergence_.val() << endl;
      }
      ostrs << indent << "  )" << endl;
      ostrs << indent << ")" << endl;
    }

  if (frequencies_.val()) {
      if (freq_accuracy_.val()) {
        ostrs << indent << "freq_accuracy = " << freq_accuracy_.val() << endl;
      }
      ostrs << indent << "freq<MolecularFrequencies>: (" << endl;
      ostrs << indent << "  molecule = $:molecule" << endl;
      ostrs << indent << ")" << endl;
    }

  ostrs << decindent;
  ostrs << indent << ")" << endl; // end of mpqc section

  if (Status::need_psi_exenv) {
    if (tmpdir_.set()) {
      ostrs << indent << "psienv<PsiExEnv>: ( scratch = [\""
          << tmpdir_.val()
          << "\"] )" << endl;
    }
  }

  ostrs << ends;

  int n = 1 + strlen(ostrs.str().c_str());
  char *in_char_array = strcpy(new char[n],ostrs.str().c_str());
  return in_char_array;
}

namespace {

  int pVXXF12_to_X(const char *name) {
    int X = 0;
    if (strncmp("cc-pV", name, 5) == 0 &&
        strncmp("Z-F12", name+6, 5) == 0) {
      switch (name[5]) {
        case 'D': X = 2; break;
        case 'T': X = 3; break;
        case 'Q': X = 4; break;
        default: throw InputError("unknown basis set", __FILE__, __LINE__);
      }
    }

    return X;
  }

  std::string X_to_XZ(int X) {
    if (X < 2 && X > 9) throw ProgrammingError("X in cc-pVXZ out of range",__FILE__,__LINE__);
    const char* letters = "DTQ56789";
    const char letter = letters[X-2];
    std::ostringstream oss;
    oss << letter << 'Z';
    return oss.str();
  }

  char* X_to_pVXZF12CABS(int X) {
    std::ostringstream oss;
    oss << "cc-pV" << X_to_XZ(X) << "-F12-CABS";
    return strdup(oss.str().c_str());
  }

  char* X_to_pVXZRI(int X) {
    std::ostringstream oss;
    oss << "cc-pV" << X_to_XZ(X) << "-RI";
    return strdup(oss.str().c_str());
  }

  char* X_pVXZF12_STG(int X) {
    switch (X) {
      case 2: return strdup("0.9");
      case 3: return strdup("1.0");
      case 4: return strdup("1.1");
      default:
        throw ProgrammingError("invalid basis",__FILE__,__LINE__);
    }
  }

}

void
MPQCIn::infer_defaults() {

  // Psi can only use all harmonics or all cartesians -> set to all harmonics by default
  const bool force_puream = psi_method(method_.val()) ? true : false;
  basis_.set_puream(force_puream);

  // R12 caculations using cc-pVXZ-F12 basis sets should use the corresponding RI and DF basis sets
  // as well as geminal exponent
  if (r12_method(method_.val())) {
    const int X = pVXXF12_to_X(basis_.name.val());
    if (X) {
      if (!auxbasis_.name.set()) {
        auxbasis_.name = X_to_pVXZF12CABS(X);
      }
      if (!dfbasis_.name.set()) {
        dfbasis_.name = X_to_pVXZRI(X+1);
      }
      if (!r12method_f12_.set()) {
        std::ostringstream oss;
        oss << "stg-6g[" << X_pVXZF12_STG(X) << "]";
        set_r12method_f12( strdup(oss.str().c_str()) );
      }
    }
  }

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
                   const std::string& r,    // RI approach
                   bool c = false,          // include coupling in H0?
                   bool e = false           // assume EBC?
                  ) :
                    approx(app),
                    ansatz(an),
                    corrfactor(cf),
                    ri(r),
                    coupling(c),
                    ebc(e)
      {
      }
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
        os << indent << "coupling = " << (coupling ? "true" : "false") << endl;
        os << indent << "ebc = " << (ebc ? "true" : "false") << endl;
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
      bool coupling;
      bool ebc;
  };
}

void
MPQCIn::write_energy_object(ostream &ostrs,
                            const char *keyword,
                            const char *method,
                            Basis const* basis,
                            int coor, IntegralsFactoryType& ifactory)
{
  const int nelectron = mol_->total_Z() - charge_.val();
  if (nelectron < 0) {
      error("charge is impossibly large");
    }
  if ((nelectron%2 == 0 && mult_.val()%2 == 0) ||
      (nelectron%2 == 1 && mult_.val()%2 == 1) ) {
      error("given multiplicity is not possible");
    }

  const char *method_object = 0;
  const char *reference_method = 0;
  const char *guess_method = method;
  bool ask_auxbasis = false;
  R12TechDescr* r12descr = 0;
  int dft = 0;
  int uscf = 0;
  bool scf = false;
  bool psi = false;
  bool psi_ccr12 = false;
  bool need_wfnworld = false;
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
            // MBPT2 can do regular MP2 only
            // if density fitting is requested use MP2-R12
            if (dfbasis_.name.set() == false) {
              method_object = "MBPT2";
            }
            else {
              method_object = "MBPT2_R12";
              r12descr = R12TechDescr::default_instance();
              r12descr->corrfactor = "none";
              need_wfnworld = true;
            }
          }
          else { // open-shell will use MP2-R12 code
            method_object = "MBPT2_R12";
            r12descr = R12TechDescr::default_instance();
            r12descr->corrfactor = "none";
            need_wfnworld = true;
            if (!strcmp(method, "RMP2")) {
              reference_method = "RHF";
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
        }
      // Local Perturbation Theory
      else if (!strcmp(method, "LMP2")) {
        guess_method = 0;
        if (mult_.val() != 1) // only closed-shell allowed
          throw InputError("LMP2 calculations are only allowed on closed-shell molecules",
                           __FILE__, __LINE__);
        reference_method = "RHF";
        method_object = "LMP2";
      }
      // MP2-R12
      else if (strncmp(method,   "MP2-R12", 7) == 0 ||
               strncmp(method,   "MP2-F12", 7) == 0 ||
               strncmp(method+1, "MP2-R12", 7) == 0 || // R/U
               strncmp(method+1, "MP2-F12", 7) == 0) { // R/U
        guess_method = 0;
        ask_auxbasis = true;
        need_wfnworld = true;
        method_object = "MBPT2_R12";
        if (method[0] == 'U')
          reference_method = "UHF";
        else if (method[0] == 'M' || method[0] == 'R')
          reference_method = "RHF";
        else
          error2("invalid method: ", method);  // XMP2-R12, X!=U && X!=R
        r12descr = R12TechDescr::default_instance();
        r12descr->coupling = true;  // include coupling in MP2-R12 by default
        r12descr->ebc = false;      // do not assume EBC in MP2-R12 by default
      }
      // CCSD(2)_R12 / CCSD(T)_R12
      else if (strncmp(method,   "CCSD(2)_R12", 11) == 0 ||
               strncmp(method,   "CCSD(2)_F12", 11) == 0 ||
               strncmp(method+1, "CCSD(2)_R12", 11) == 0 || // R/U
               strncmp(method+1, "CCSD(2)_F12", 11) == 0 || // R/U
               strncmp(method,   "CCSD(T)_R12", 11) == 0 ||
               strncmp(method,   "CCSD(T)_F12", 11) == 0 ||
               strncmp(method+1, "CCSD(T)_R12", 11) == 0 || // R/U
               strncmp(method+1, "CCSD(T)_F12", 11) == 0 || // R/U
               strncmp(method,   "CC3(2)_R12", 11)  == 0 ||
               strncmp(method,   "CC3(2)_F12", 11)  == 0 ||
               strncmp(method+1, "CC3(2)_R12", 11)  == 0 || // R/U
               strncmp(method+1, "CC3(2)_F12", 11)  == 0) { //R/U
        guess_method = 0;
        ask_auxbasis = true;
        need_wfnworld = true;
        psi = true;
        if (method[0] == 'U')
          reference_method = "PsiUHF";
        else if (method[0] == 'C' || method[0] == 'R')
          reference_method = "PsiRHF";
        else
          error2("invalid method: ", method);

        if ((method[2] == 'S' || method[3] == 'S')) { // CCSD
          if (method[5] == '2' || method[6] == '2') // CCSD(2)
            method_object = "PsiCCSD_PT2R12";
          else  // CCSD(T)
            method_object = "PsiCCSD_PT2R12T";
        }
        else { // CC3
          method_object = "PsiCC3_PT2R12";
        }

        r12descr = R12TechDescr::default_instance();
        r12descr->coupling = false;  // inclusion of coupling in CC-R12 not implemented
        r12descr->ebc = false;       // do not assume EBC in CC-R12 by default
        psi_ccr12 = true;

      }
      else if (!strcmp(method, "PsiHF")) {
        if (mult_.val() == 1) method_object = "PsiCLHF";
        else { uscf = 1; method_object = "PsiUHF"; }
        scf = true;
        psi = true;
        guess_method = "HF";
      }
      else if (!strcmp(method, "PsiRHF")) {
        guess_method = "RHF";
        if (mult_.val() == 1) method_object = "PsiCLHF";
        else method_object = "PsiHSOSHF";
        scf = true;
        psi = true;
        guess_method = "RHF";
      }
      else if (!strcmp(method, "PsiUHF")) {
        guess_method = "UHF";
        method_object = "PsiUHF";
        uscf = 1;
        scf = true;
        psi = true;
        guess_method = "UHF";
      }
      else error2("invalid method: ", method);
    }
  else error("no method given");

  ostrs << indent << keyword << "<" << method_object << ">: (" << endl;
  ostrs << incindent;
  ostrs << o_extra.str();

  if (ask_auxbasis
      && auxbasis_.name.val() != 0
      && auxbasis_ != basis_)
    auxbasis_.write(ostrs, "aux_basis");

  if (r12descr) {
    if (r12method_f12_.set()) r12descr->corrfactor = r12method_f12_.val();
    if (r12method_app_.set()) r12descr->approx = r12method_app_.val();
    if (r12method_ansatz_.set()) r12descr->ansatz = r12method_ansatz_.val();
    if (r12method_ri_.set()) r12descr->ri = r12method_ri_.val();
    r12descr->write(ostrs);
    switch (r12descr->corrfactor_type()) {
      case R12TechDescr::R12CorrFactor:  ifactory = Libint2; break;
      case R12TechDescr::GTGCorrFactor:  ifactory = Libint2; break;
      default: break;
    }

    if (need_wfnworld) {
      if (tmpstore_.set()) {
        ostrs << indent << "store_ints = " << (!strcmp(tmpstore_.val(),"mem") ? "mem" : "posix") << endl;
      }
      if (tmpdir_.set())
        ostrs << indent << "ints_file = " << tmpdir_.val() << endl;
      if (dfbasis_.name.val() != 0)
        dfbasis_.write(ostrs, "df_basis");
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

  // basis
  if (basis) { // if basis is given explicitly, use it
    if (psi) { // Psi only allows all puream or all cartesians. Default to all puream
      Basis puream_basis(*basis);
      puream_basis.write(ostrs, "basis");
    }
    else
      basis->write(ostrs, "basis");
    }
  else { // if basis is not given, refer to the top keyword
    ostrs << indent << "basis = $:basis" << endl;
  }

  // dft
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

  // guess basis
  {
    Basis gbasis(guess_basis(ifactory));
    bool gbasis_eq_basis;
    if (basis)
      gbasis_eq_basis = (gbasis == *basis);
    else
      gbasis_eq_basis = (gbasis == basis_);
    if (dft || (guess_method && !gbasis_eq_basis) ) {
      if (frequencies_.val()) {
          ostrs << indent << "keep_guess_wavefunction = 1" << endl;;
        }
      write_energy_object(ostrs, "guess_wavefunction",
                          guess_method, &gbasis, 0, ifactory);
    }
  }

  // reference wfn
  if (reference_method) {
    ostrs << indent << "nfzc = auto" << endl;;
    write_energy_object(ostrs, "reference",
                        reference_method, 0, 0, ifactory);
  }

  // Psi wfn? need psi environment
  if (psi) {
    Status::need_psi_exenv = true;
    ostrs << indent << "psienv = $:psienv" << endl;
  }

  // a Psi-based CC object
  if (psi_ccr12) {
    // TODO make sure this is a closed-shell
    if (pccsd_alpha_.set())
      ostrs << indent << "pccsd_alpha = " << pccsd_alpha_.val() << endl;;
    if (pccsd_beta_.set())
      ostrs << indent << "pccsd_beta = " << pccsd_beta_.val() << endl;;
    if (pccsd_gamma_.set())
      ostrs << indent << "pccsd_gamma = " << pccsd_gamma_.val() << endl;;
  }

  ostrs << decindent;
  ostrs << indent << ")" << endl;
}

void
MPQCIn::Basis::write(ostream &ostrs,
                     const char *keyword) const
{
  // validate input
  if (!name.val()) throw InputError("no basis given", __FILE__, __LINE__);

  if (!split.val() && !uc.val()) {
    ostrs << indent << keyword << "<GaussianBasisSet>: (" << endl;
    ostrs << incindent;
    ostrs << indent << "molecule = $:molecule" << endl;
    ostrs << indent << "name = \"" << name.val() << "\"" << endl;
    if (puream.val()) ostrs << indent << "puream = true" << endl;
    ostrs << decindent;
    ostrs << indent << ")" << endl;
  }
  else {
    std::ostringstream oss;
    oss << "m" << keyword;
    const char* mkeyword = oss.str().c_str();
    Basis mother(*this);
    mother.set_uc(false);
    mother.set_split(false);
    mother.write(ostrs, mkeyword);
    if (uc.val())
      ostrs << indent << keyword << "<UncontractedBasisSet>: (" << endl;
    if (split.val())
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
    case Libint2:
      result = "IntegralLibint2"; break;
    default:
      throw std::logic_error("Invalid integral factory");
  }
  return result;
}

MPQCIn::Basis
MPQCIn::guess_basis(IntegralsFactoryType ifactory) {
  // split STO-3G basis if factory is not IntV3
  return Basis("STO-3G", false, (ifactory != IntV3), false);
}

bool
MPQCIn::psi_method(const char* method) {
  return (strncmp(method,   "Psi",         3) == 0 ||
          strncmp(method,   "CCSD(2)_R12", 11) == 0 ||
          strncmp(method,   "CCSD(2)_F12", 11) == 0 ||
          strncmp(method+1, "CCSD(2)_R12", 11) == 0 || // R/U
          strncmp(method+1, "CCSD(2)_F12", 11) == 0 || // R/U
          strncmp(method,   "CC3(2)_R12", 11)  == 0 ||
          strncmp(method,   "CC3(2)_F12", 11)  == 0 ||
          strncmp(method+1, "CC3(2)_R12", 11)  == 0 || // R/U
          strncmp(method+1, "CC3(2)_F12", 11)  == 0 || // R/U
          strncmp(method,   "CCSD(T)_R12", 11) == 0 ||
          strncmp(method,   "CCSD(T)_F12", 11) == 0 ||
          strncmp(method+1, "CCSD(T)_R12", 11) == 0 || // R/U
          strncmp(method+1, "CCSD(T)_F12", 11) == 0 );
}

bool
MPQCIn::r12_method(const char* method) {
  return (strstr(method, "R12") != 0 ||
          strstr(method, "F12") != 0);
}
