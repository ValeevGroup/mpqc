
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdlib.h>
#include <math.h>

#include <util/keyval/keyval.h>
#include <util/misc/libmisc.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <math/scmat/matrix.h>
#include <math/symmetry/pointgrp.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/psi/psi.h>
#include <chemistry/qc/psi/file11.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiWfn_cd(
  typeid(PsiWfn),"PsiWfn",1,"public Wavefunction",
  0, 0, 0);

PsiWfn::PsiWfn(const Ref<KeyVal>&keyval):
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
          ExEnv::out() << "PsiWfn: inputth failed" << endl;
          abort();
        }
      if (system("psi")) {
          ExEnv::out() << "PsiWfn: psi failed" << endl;
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
              ExEnv::err() << "ERROR: Psi: atom "<<i<<" not found in file11" << endl;
              abort();
            }
        }
      ExEnv::out() << node0 << indent << " Psi<->MPQC atom reordering:";
      for(i=0; i<molecule()->natom(); i++) ExEnv::out() << node0 << " " << reorder[i];
      ExEnv::out() << node0 << endl;

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

static ClassDesc PsiCCSD_cd(
  typeid(PsiCCSD),"PsiCCSD",1,"public PsiWfn",
  0, create<PsiCCSD>, create<PsiCCSD>);

PsiCCSD::PsiCCSD(const Ref<KeyVal>&keyval):
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
          ExEnv::out() << "PsiWfn: could not find CCSD energy in output file" << endl;
          abort();
        }
    }
  in = fopen("psitmp.energy","r");
  if (!in) {
      ExEnv::err() << "PsiCCSD::read_energy(): cannot open psitmp.energy" << endl;
      abort();
    }
  double r;
  while(fscanf(in,"%*s %lf", &r)>0){ }
  fclose(in);

  return r;
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCCSD_T_cd(
  typeid(PsiCCSD_T),"PsiCCSD_T",1,"public PsiWfn",
  0, create<PsiCCSD_T>, create<PsiCCSD_T>);

PsiCCSD_T::PsiCCSD_T(const Ref<KeyVal>&keyval):
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
          ExEnv::out() << "PsiWfn: could not find CCSD(T) energy in output file"
               << endl;
          abort();
        }
    }
  in = fopen("psitmp.energy","r");
  if (!in) {
      ExEnv::err() << "PsiCCSD_T::read_energy(): cannot open psitmp.energy" << endl;
      abort();
    }
  double r;
  while(fscanf(in,"%*s %lf", &r)>0){ }
  fclose(in);

  return r;
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCCSDT_cd(
  typeid(PsiCCSDT),"PsiCCSDT",1,"public PsiWfn",
  0, create<PsiCCSDT>, create<PsiCCSDT>);

PsiCCSDT::PsiCCSDT(const Ref<KeyVal>&keyval):
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
          ExEnv::out() << "PsiWfn: could not find CCSDT energy in output file"
               << endl;
          abort();
        }
    }
  in = fopen("psitmp.energy","r");
  if (!in) {
      ExEnv::err() << "PsiCCSDT::read_energy(): cannot open psitmp.energy" << endl;
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

static ClassDesc PsiCI_cd(
  typeid(PsiCI),"PsiCI",1,"public PsiWfn",
  0, create<PsiCI>, create<PsiCI>);

PsiCI::PsiCI(const Ref<KeyVal>&keyval):
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
          ExEnv::out() << "PsiWfn: could not find CI energy in output file"
               << endl;
          abort();
        }
    }
  in = fopen("psitmp.energy","r");
  if (!in) {
      ExEnv::err() << "PsiCI::compute(): cannot open psitmp.energy" << endl;
      abort();
    }
  double r;
  fscanf(in,"%*s %*s %*s %*s %lf", &r);
  fclose(in);
  return r;
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiHF_cd(
  typeid(PsiHF),"PsiHF",1,"public PsiWfn",
  0, create<PsiHF>, create<PsiHF>);

PsiHF::PsiHF(const Ref<KeyVal>&keyval):
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
          ExEnv::out() << "PsiWfn: could not find Hartree-Fock energy in output file"
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
