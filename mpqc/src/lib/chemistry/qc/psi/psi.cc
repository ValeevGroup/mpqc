
#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
#include <stdlib.h>
#include <math.h>
}

#include <util/keyval/keyval.h>
#include <util/misc/libmisc.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <math/scmat/matrix.h>
#include <math/symmetry/pointgrp.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/psi/psi.h>
#include <chemistry/qc/psi/file11.h>

//////////////////////////////////////////////////////////////////////////

#define CLASSNAME PsiWfn
#define PARENTS public Wavefunction
#include <util/state/statei.h>
#include <util/class/classia.h>

void *
PsiWfn::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Wavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

PsiWfn::PsiWfn(const RefKeyVal&keyval):
  Wavefunction(keyval), psi_in(keyval)
{
}

PsiWfn::~PsiWfn()
{
}

PsiWfn::PsiWfn(StateIn&s):
  SavableState(s),
  Wavefunction(s)
{
  abort();
}

void
PsiWfn::save_data_state(StateOut&s)
{
  Wavefunction::save_data_state(s);
  abort();
}

void
PsiWfn::print(ostream&o) const
{
  Wavefunction::print(o);
}

void
PsiWfn::compute()
{  
  double energy_acc = desired_value_accuracy();
  double grad_acc = desired_gradient_accuracy();
  if (energy_acc > 1.0e-6) energy_acc = 1.0e-6;
  if (grad_acc > 1.0e-7) grad_acc = 1.0e-7;
  if (gradient_needed() && energy_acc > grad_acc/10.0)
      energy_acc = grad_acc/10.0;

  if (!psi_in.test()){
      write_input((int)-log10(energy_acc));
      if (system("inputth")) {
          cout << "PsiWfn: inputth failed" << endl;
          abort();
        }
      if (system("psi")) {
          cout << "PsiWfn: psi failed" << endl;
          abort();
        }
    }

  // read output
  if (gradient_needed()) {
      int i, j, ii;
      int *reorder = new int[molecule()->natom()];
      double tol = 1e-6;
      for (i=0; i<molecule()->natom(); i++) reorder[i] = -1;
      FILE11 file11(0);
      for(i=0; i<molecule()->natom(); i++) {
          int gotj = 0;
          for(j=0; j<molecule()->natom(); j++) {
              if( fabs(file11.coordinate(0,i) - molecule()->r(j,0)) < tol &&
                  fabs(file11.coordinate(1,i) - molecule()->r(j,1)) < tol && 
                  fabs(file11.coordinate(2,i) - molecule()->r(j,2)) < tol){
                  reorder[j] = i;
                  gotj = 1;
                  break;
                }
            }
          if (!gotj) {
              cerr << "ERROR: Psi: atom "<<i<<" not found in file11" << endl;
              abort();
            }
        }
      cout << node0 << indent << " Psi<->MPQC atom reordering:";
      for(i=0; i<molecule()->natom(); i++) cout << node0 << " " << reorder[i];
      cout << node0 << endl;

      RefSCVector g(moldim(),matrixkit());

      for (ii=0, i=0; i<molecule()->natom(); i++) {
          for (j=0; j<3; j++, ii++) {
              g(ii) = file11.coordinate(j,reorder[i]);
            }
        }
      print_natom_3(g, "Reordered FILE11 Geometry");

      for (ii=0, i=0; i<molecule()->natom(); i++) {
          for (j=0; j<3; j++, ii++) {
              g(ii) = file11.gradient(j,reorder[i]);
            }
        }
      delete[] reorder;
      print_natom_3(g, "Reordered FILE11 Gradient");
      set_gradient(g);
      set_actual_gradient_accuracy(grad_acc);
      set_energy(file11.energy());
      set_actual_value_accuracy(energy_acc);
    }
  else {
      double r = read_energy();
      set_energy(r);
      set_actual_value_accuracy(energy_acc);
    }
}

double
PsiWfn::density(const SCVector3&d)
{
  abort();
  return 0.0;
}

RefSymmSCMatrix
PsiWfn::density()
{
  abort();
  return 0;
}

int
PsiWfn::spin_polarized()
{
  return 1;
}

int
PsiWfn::nelectron()
{
  abort();
  return 0;
}

int
PsiWfn::value_implemented() const
{
  return 1;
}

int
PsiWfn::gradient_implemented() const
{
  return 1;
}

void
PsiWfn::write_basic_input(int conv, const char *wfn)
{
  const char *dertype = gradient_needed() ? "FIRST" : "NONE";
  psi_in.write_defaults(dertype, wfn);
  psi_in.write_input();
  psi_in.write_basis();
}

//////////////////////////////////////////////////////////////////////////

#define CLASSNAME PsiCCSD
#define PARENTS public PsiWfn
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
PsiCCSD::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = PsiWfn::_castdown(cd);
  return do_castdowns(casts,cd);
}

PsiCCSD::PsiCCSD(const RefKeyVal&keyval):
  PsiWfn(keyval)
{
}

PsiCCSD::~PsiCCSD()
{
}

PsiCCSD::PsiCCSD(StateIn&s):
  SavableState(s),
  PsiWfn(s)
{
  abort();
}

void
PsiCCSD::save_data_state(StateOut&s)
{
  PsiWfn::save_data_state(s);
  abort();
}

void
PsiCCSD::write_input(int convergence)
{
  psi_in.open("input.dat");
  PsiWfn::write_basic_input(convergence, "CCSD");
  psi_in.begin_section("cceg");
  if(convergence != 0){
      psi_in.write_keyword("convergence", convergence);
    }
  psi_in.write_keyword("restart", "no");
  psi_in.write_keyword("maxiter", 50);
  psi_in.end_section();
  psi_in.begin_section("cczv");
  if(convergence != 0){
      psi_in.write_keyword("convergence", convergence);
    }
  psi_in.write_keyword("restart", "no");
  psi_in.write_keyword("maxiter", 50);
  psi_in.end_section();

  psi_in.close();
}

double
PsiCCSD::read_energy()
{
  FILE *in;
  remove("psitmp.energy");
  if(!psi_in.test()){
      if (system("grep \"CCSD \" energy.dat > psitmp.energy")) {
          cout << "PsiWfn: could not find CCSD energy in output file" << endl;
          abort();
        }
    }
  in = fopen("psitmp.energy","r");
  if (!in) {
      cerr << "PsiCCSD::read_energy(): cannot open psitmp.energy" << endl;
      abort();
    }
  double r;
  while(fscanf(in,"%*s %lf", &r)>0){ }
  fclose(in);

  return r;
}

//////////////////////////////////////////////////////////////////////////

#define CLASSNAME PsiCCSD_T
#define PARENTS public PsiWfn
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
PsiCCSD_T::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = PsiWfn::_castdown(cd);
  return do_castdowns(casts,cd);
}

PsiCCSD_T::PsiCCSD_T(const RefKeyVal&keyval):
  PsiWfn(keyval)
{
}

PsiCCSD_T::~PsiCCSD_T()
{
}

PsiCCSD_T::PsiCCSD_T(StateIn&s):
  SavableState(s),
  PsiWfn(s)
{
  abort();
}

void
PsiCCSD_T::save_data_state(StateOut&s)
{
  PsiWfn::save_data_state(s);
  abort();
}

void
PsiCCSD_T::write_input(int convergence)
{
  psi_in.open("input.dat");
  PsiWfn::write_basic_input(convergence, "CCSD_T");
  psi_in.begin_section("cceg");
  if(convergence != 0){
      psi_in.write_keyword("convergence", convergence);
    }
  psi_in.write_keyword("restart", "no");
  psi_in.write_keyword("maxiter", 50);
  psi_in.end_section();
  psi_in.begin_section("cczv");
  if(convergence != 0){
      psi_in.write_keyword("convergence", convergence);
    }
  psi_in.write_keyword("restart", "no");
  psi_in.write_keyword("maxiter", 50);
  psi_in.end_section();

  psi_in.close();
}

double
PsiCCSD_T::read_energy()
{
  FILE *in;
  remove("psitmp.energy");
  if(!psi_in.test()){
      if (system("grep \"CCT \" energy.dat > psitmp.energy")) {
          cout << "PsiWfn: could not find CCSD(T) energy in output file"
               << endl;
          abort();
        }
    }
  in = fopen("psitmp.energy","r");
  if (!in) {
      cerr << "PsiCCSD_T::read_energy(): cannot open psitmp.energy" << endl;
      abort();
    }
  double r;
  while(fscanf(in,"%*s %lf", &r)>0){ }
  fclose(in);

  return r;
}

//////////////////////////////////////////////////////////////////////////

#define CLASSNAME PsiCCSDT
#define PARENTS public PsiWfn
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
PsiCCSDT::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = PsiWfn::_castdown(cd);
  return do_castdowns(casts,cd);
}

PsiCCSDT::PsiCCSDT(const RefKeyVal&keyval):
  PsiWfn(keyval)
{
}

PsiCCSDT::~PsiCCSDT()
{
}

PsiCCSDT::PsiCCSDT(StateIn&s):
  SavableState(s),
  PsiWfn(s)
{
  abort();
}

void
PsiCCSDT::save_data_state(StateOut&s)
{
  PsiWfn::save_data_state(s);
  abort();
}

void
PsiCCSDT::write_input(int convergence)
{
  psi_in.open("input.dat");
  PsiWfn::write_basic_input(convergence, "CCSDT");
  psi_in.begin_section("cceg");
  if(convergence != 0){
      psi_in.write_keyword("convergence", convergence);
    }
  psi_in.write_keyword("restart", "no");
  psi_in.write_keyword("maxiter", 50);
  psi_in.end_section();
  psi_in.begin_section("cczv");
  if(convergence != 0){
      psi_in.write_keyword("convergence", convergence);
    }
  psi_in.write_keyword("restart", "no");
  psi_in.write_keyword("maxiter", 50);
  psi_in.end_section();

  psi_in.close();
}

double
PsiCCSDT::read_energy()
{
  FILE *in;
  remove("psitmp.energy");
  if(!psi_in.test()){
      if (system("grep \"FSDT \" energy.dat > psitmp.energy")) {
          cout << "PsiWfn: could not find CCSDT energy in output file"
               << endl;
          abort();
        }
    }
  in = fopen("psitmp.energy","r");
  if (!in) {
      cerr << "PsiCCSDT::read_energy(): cannot open psitmp.energy" << endl;
      abort();
    }
  double r;
  while(fscanf(in,"%*s %lf", &r)>0){ }
  fclose(in);

  return r;
}

int
PsiCCSDT::gradient_implemented() const
{
  return 0;
}

//////////////////////////////////////////////////////////////////////////

#define CLASSNAME PsiCI
#define PARENTS public PsiWfn
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
PsiCI::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = PsiWfn::_castdown(cd);
  return do_castdowns(casts,cd);
}

PsiCI::PsiCI(const RefKeyVal&keyval):
  PsiWfn(keyval)
{
}

PsiCI::~PsiCI()
{
}

PsiCI::PsiCI(StateIn&s):
  SavableState(s),
  PsiWfn(s)
{
  abort();
}

void
PsiCI::save_data_state(StateOut&s)
{
  PsiWfn::save_data_state(s);
  abort();
}

void
PsiCI::write_input(int convergence)
{
  psi_in.open("input.dat");
  PsiWfn::write_basic_input(convergence, "CI");
  if(convergence != 0){
      psi_in.begin_section("gugaci");
      psi_in.write_keyword("convergence", convergence);
      psi_in.end_section();
    }
  psi_in.close();
}

double
PsiCI::read_energy()
{
  FILE *in;
  remove("psitmp.energy");
  if(!psi_in.test()){
      if (system("grep \"1 ECI \" output.dat > psitmp.energy")) {
          cout << "PsiWfn: could not find CI energy in output file"
               << endl;
          abort();
        }
    }
  in = fopen("psitmp.energy","r");
  if (!in) {
      cerr << "PsiCI::compute(): cannot open psitmp.energy" << endl;
      abort();
    }
  double r;
  fscanf(in,"%*s %*s %*s %*s %lf", &r);
  fclose(in);
  return r;
}

//////////////////////////////////////////////////////////////////////////

#define CLASSNAME PsiHF
#define PARENTS public PsiWfn
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
PsiHF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = PsiWfn::_castdown(cd);
  return do_castdowns(casts,cd);
}

PsiHF::PsiHF(const RefKeyVal&keyval):
  PsiWfn(keyval)
{
}

PsiHF::~PsiHF()
{
}

PsiHF::PsiHF(StateIn&s):
  SavableState(s),
  PsiWfn(s)
{
  abort();
}

void
PsiHF::save_data_state(StateOut&s)
{
  PsiWfn::save_data_state(s);
  abort();
}

void
PsiHF::write_input(int convergence)
{
  psi_in.open("input.dat");
  PsiWfn::write_basic_input(convergence, "SCF");
  if(convergence != 0){
      psi_in.begin_section("scf");
      psi_in.write_keyword("convergence", convergence);
      psi_in.end_section();
    }
  psi_in.close();
}

double
PsiHF::read_energy()
{
  FILE *in;
  remove("psitmp.energy");
  if(!psi_in.test()){
      if (system("grep \"total energy *=\" output.dat > psitmp.energy")) {
          cout << "PsiWfn: could not find Hartree-Fock energy in output file"
               << endl;
          abort();
        }
    }
  in = fopen("psitmp.energy","r");
  if (!in) {
      fprintf(stderr,"PsiHF::compute(): cannot open psitmp.energy\n");
      abort();
    }
  double r;
  fscanf(in,"%*s %*s %*s %lf", &r);
  fclose(in);
  return r;
}

//////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
