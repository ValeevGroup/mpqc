
#include <scconfig.h>
#ifdef HAVE_SSTREAM
#  include <sstream>
#else
#  include <strstream.h>
#endif
#include <stdlib.h>

#include <util/misc/formio.h>
#include <util/container/array.h>

using namespace std;

#undef yyFlexLexer
#define yyFlexLexer MPQCInFlexLexer
#include <FlexLexer.h>

#include "parse.h"
#include "mpqcin.h"

int MPQCIn::checking_ = 0;

MPQCIn::MPQCIn():
  nirrep_(0),
  mult_(1),
  charge_(0),
  basis_(0),
  method_(0),
  optimize_(0),
  gradient_(0),
  frequencies_(0),
  opt_type_(T_INTERNAL),
  redund_coor_(0),
  restart_(0),
  checkpoint_(1),
  atom_charge_(0),
  method_xc_(0),
  method_grid_(0),
  symmetry_(0),
  molecule_bohr_(0),
  alpha_(0),
  beta_(0),  
  docc_(0),
  socc_(0),  
  frozen_docc_(0),
  frozen_uocc_(0)
{
  lexer_ = new MPQCInFlexLexer;
}

MPQCIn::~MPQCIn()
{
  delete lexer_;
  free(basis_.val());
  free(method_.val());
  free(method_xc_.val());
  free(method_grid_.val());
  free(symmetry_.val());
  free(alpha_.val());
  free(beta_.val());
  free(docc_.val());
  free(socc_.val());
  free(frozen_docc_.val());
  free(frozen_uocc_.val());
}

void
MPQCIn::error(const char* s)
{
  ExEnv::out() << ExEnv::program_name()
               << ": error: " << s
               << endl;
  abort();
}

void
MPQCIn::error2(const char* s, const char *s2)
{
  ExEnv::out() << ExEnv::program_name()
               << ": error: " << s << "\"" << s2 << "\""
               << endl;
  abort();
}

void
MPQCIn::yerror(const char* s)
{
  ExEnv::out() << ExEnv::program_name()
               << ": " << s
               << " at line " << lexer_->lineno()+1
               << endl;
  abort();
}

void
MPQCIn::yerror2(const char* s, const char *s2)
{
  ExEnv::out() << ExEnv::program_name()
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
      ExEnv::out() << ExEnv::program_name()
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
  int Z = AtomInfo::string_to_Z(sym, 0);
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

  mol_->add_atom(Z, x, y, z, 0, 0, atom_charge_.set(), atom_charge_.val());
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
MPQCIn::set_method_xc(char *m)
{
  method_xc_ = m;
}

void
MPQCIn::set_method_grid(char *m)
{
  method_grid_ = m;
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
MPQCIn::set_symmetry(char *s)
{
  symmetry_ = s;
  if (strcmp(s,"auto")) {
      RefPointGroup p = new PointGroup(s);
      nirrep_ = p->char_table().nirrep();
    }
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
MPQCIn::set_docc(Arrayint *a)
{
  docc_ = a;
}

void
MPQCIn::set_socc(Arrayint *a)
{
  socc_ = a;
}

void
MPQCIn::set_alpha(Arrayint *a)
{
  alpha_ = a;
}

void
MPQCIn::set_beta(Arrayint *a)
{
  beta_ = a;
}

void
MPQCIn::set_frozen_docc(Arrayint *a)
{
  frozen_docc_ = a;
}

void
MPQCIn::set_frozen_uocc(Arrayint *a)
{
  frozen_uocc_ = a;
}

Arrayint *
MPQCIn::make_nnivec(Arrayint *a, char *ms)
{
  if (ms == 0) return new Arrayint;

  char *mse;
  int m = strtol(ms,&mse,10);
  if (mse == ms || m < 0) yerror2("bad positive integer", ms);
  free(ms);

  if (a == 0) a = new Arrayint;
  a->push_back(m);
  return a;
}

int
MPQCIn::check_string(const char *s)
{
  checking_ = 1;
#ifdef HAVE_SSTREAM
  istringstream in(s);
#else
  istrstream in(s);
#endif
  lexer_->switch_streams(&in, &ExEnv::out());
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
#ifdef HAVE_SSTREAM
  istringstream in(s);
#else
  istrstream in(s);
#endif
  lexer_->switch_streams(&in, &ExEnv::out());
  yparse();

  // form the oo input
#ifdef HAVE_SSTREAM
  ostringstream ostrs;
#else
  ostrstream ostrs;
#endif
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
  write_energy_object(ostrs, "mole", method_.val(), 0, optimize_.val());
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
      ostrs << indent << "freq<MolecularFrequencies>: (";
      ostrs << indent << "  molecule = $:molecule";
      ostrs << indent << ")";
    }

  ostrs << decindent;
  ostrs << indent << ")" << endl;
  ostrs << ends;

#ifdef HAVE_SSTREAM
  int n = 1 + strlen(ostrs.str().c_str());
  char *in_char_array = strcpy(new char[n],ostrs.str().c_str());
#else
  char *in_char_array = ostrs.str();
#endif
  return in_char_array;
}

void
MPQCIn::write_vector(ostream &ostrs,
                     const char *keyvalname,
                     const char *name, MPQCInDatum<Arrayint *>&vec,
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

void
MPQCIn::write_energy_object(ostream &ostrs,
                            const char *keyword,
                            const char *method,
                            const char *basis,
                            int coor)
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
  int dft = 0;
  int uscf = 0;
  if (method) {
      // Hartree-Fock methods
      if (!strcmp(method, "HF")) {
          if (mult_.val() == 1) method_object = "CLHF";
          else { uscf = 1; method_object = "UHF"; }
        }
      else if (!strcmp(method, "RHF")) {
          if (mult_.val() == 1) method_object = "CLHF";
          else method_object = "HSOSHF";
        }
      else if (!strcmp(method, "UHF")) {
          method_object = "UHF";
          uscf = 1;
        }
      // Density Functional Methods
      else if (!strcmp(method, "KS")) {
          guess_method = "HF";
          if (mult_.val() == 1) method_object = "CLKS";
          else { uscf = 1; method_object = "UKS"; }
          dft = 1;
        }
      else if (!strcmp(method, "RKS")) {
          guess_method = "RHF";
          if (mult_.val() == 1) method_object = "CLKS";
          else method_object = "HSOSKS";
          dft = 1;
        }
      else if (!strcmp(method, "UKS")) {
          guess_method = "UHF";
          method_object = "UKS";
          dft = 1;
          uscf = 1;
        }
      // Perturbation Theory
      else if (!strcmp(method, "MP2")) {
          guess_method = 0;
          method_object = "MBPT2";
          reference_method = "HF";
          if (mult_.val() != 1) {
              error("MP2 can only be used with multiplicity 1: try ZAPT2");
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
      else error2("invalid method: ", method);
    }
  else error("no method given");
  ostrs << indent << keyword << "<" << method_object << ">: (" << endl;
  ostrs << incindent;
  ostrs << indent << "total_charge = " << charge_.val() << endl;
  ostrs << indent << "molecule = $:molecule" << endl;
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
      if (method_xc_.set()) {
          ostrs << indent << "functional<StdDenFunctional>: ( name = \""
                << method_xc_.val()
                << "\" )" << endl;
        }
      else error("no exchange-correlation functional given");
      if (method_grid_.set()) {
          ostrs << indent << "integrator<RadialAngularIntegrator>: ( grid = \""
                << method_grid_.val()
                << "\" )" << endl;
        }
    }
  if (dft || (!(basis
                && !strncmp("STO",basis,3))
              && strncmp("STO",basis_.val(),3)
              && guess_method)) {
      write_energy_object(ostrs, "guess_wavefunction",
                          guess_method, "STO-3G", 0);
    }
  if (reference_method) {
      ostrs << indent << "nfzc = auto" << endl;;
      write_energy_object(ostrs, "reference",
                          reference_method, 0, 0);
    }
  ostrs << decindent;
  ostrs << indent << ")" << endl;
}

void
MPQCIn::write_basis_object(ostream &ostrs,
                           const char *keyword,
                           const char *basis)
{
  if (!basis) error("no basis given");
  ostrs << indent << keyword << "<GaussianBasisSet>: (" << endl;
  ostrs << incindent;
  ostrs << indent << "molecule = $:molecule" << endl;
  ostrs << indent << "name = \"" << basis << "\"" << endl;
  ostrs << decindent;
  ostrs << indent << ")" << endl;
}
