//
// psiexenv.cc
//
// Copyright (C) 2002 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUG__
#pragma implementation
#endif

#include <string>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <errno.h>
#include <scconfig.h>
#include <util/class/scexception.h>
#include <util/ref/ref.h>
#include <util/keyval/keyval.h>
#include <util/misc/formio.h>
#include <chemistry/qc/psi/psiexenv.h>
#include <psifiles.h>
#include <chemistry/qc/mbptr12/print.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/shellrot.h>

using namespace std;
using namespace sc;

static ClassDesc PsiExEnv_cd(
  typeid(PsiExEnv),"PsiExEnv",1,"public DescribedClass",
  create<PsiExEnv>, create<PsiExEnv>, 0);

string PsiExEnv::defaultinputname_("input.dat");
string PsiExEnv::defaultoutputname_("output.dat");
string PsiExEnv::file11name_("file11.dat");
int PsiExEnv::ckptfile_(PSIF_CHKPT);
string PsiExEnv::defaultcwd_(".");
string PsiExEnv::defaultfileprefix_("psi");
string PsiExEnv::defaultpsiprefix_(PSIROOTDIR "/bin");
string PsiExEnv::defaultstdout_("stdout");
string PsiExEnv::defaultstderr_("stderr");

PsiExEnv::PsiExEnv(const Ref<KeyVal>& keyval) :
	psio_(), chkpt_(0), me_(MessageGrp::get_default_messagegrp()->me())
{
  const std::string prefix(SCFormIO::fileext_to_filename("."));

  // Find Psi
  char *psibin = getenv("PSIBIN");
  if (psibin)
    psiprefix_ = string(psibin);
  else
    psiprefix_ = string(defaultpsiprefix_);
  add_to_path(psiprefix_);

  cwd_ = keyval->stringvalue("cwd", KeyValValuestring(defaultcwd_));
  fileprefix_ = keyval->stringvalue("fileprefix", KeyValValuestring(prefix + defaultfileprefix_));

  inputname_ = fileprefix_ + "." + defaultinputname_;
  outputname_ = fileprefix_ + "." + defaultoutputname_;

  stdout_ = keyval->stringvalue("stdout", KeyValValuestring(fileprefix_ + "." + defaultstdout_));
  stderr_ = keyval->stringvalue("stderr", KeyValValuestring(fileprefix_ + "." + defaultstderr_));

  if (keyval->exists("scratch")) {
    nscratch_ = keyval->count("scratch");
    scratch_ = new string[nscratch_];
    for (int i=0; i<nscratch_; i++)
      scratch_[i] = keyval->stringvalue("scratch", i);
  }
  else {
    nscratch_ = 1;
    scratch_ = new string[nscratch_];
    scratch_[0] = cwd_ + "/";
  }

  char *s = new char[cwd_.size() + inputname_.size() + 2];
  sprintf(s,"%s/%s",cwd_.c_str(),inputname_.c_str());
  psiinput_ = new PsiInput(s);
  delete[] s;

  s = new char[cwd_.size() + fileprefix_.size() + file11name_.size() + 3];
  sprintf(s,"%s/%s.%s",cwd_.c_str(),fileprefix_.c_str(),file11name_.c_str());
  psifile11_ = new PsiFile11(s);
  delete[] s;

  config_psio();
}

PsiExEnv::PsiExEnv() :
    psio_(), chkpt_(0), me_(MessageGrp::get_default_messagegrp()->me()),
    cwd_(defaultcwd_), nscratch_(1)
{
  // Find Psi
  char *psibin = getenv("PSIBIN");
  if (psibin)
    psiprefix_ = string(psibin);
  else
    psiprefix_ = string(defaultpsiprefix_);
  add_to_path(psiprefix_);

  scratch_ = new string[nscratch_];
  for (int i=0; i<nscratch_; i++)
    scratch_[i] = cwd_ + '/';

  fileprefix_ = SCFormIO::fileext_to_filename(".") + defaultfileprefix_;
  inputname_ = fileprefix_ + "." + defaultinputname_;
  outputname_ = fileprefix_ + "." + defaultoutputname_;
  stdout_ = fileprefix_ + "." + defaultstdout_;
  stderr_ = fileprefix_ + "." + defaultstderr_;

  config_psio();
}

PsiExEnv::~PsiExEnv()
{
  run_psiclean();
  delete[] scratch_;
  if (chkpt_) delete chkpt_;
}

void PsiExEnv::config_psio()
{
	// configure libpsio object
	{
		psio_.filecfg_kwd("DEFAULT", "NAME", -1, fileprefix_.c_str());
		std::ostringstream oss;
		oss << nscratch_;
		psio_.filecfg_kwd("DEFAULT", "NVOLUME", -1, oss.str().c_str());
		for (int i=0; i<nscratch_; i++) {
			std::ostringstream oss;
			oss << "VOLUME"<< (i+1);
			psio_.filecfg_kwd("DEFAULT", oss.str().c_str(), -1, scratch_[i].c_str());
		}
		psio_.filecfg_kwd("DEFAULT", "NVOLUME", PSIF_CHKPT, "1");
		psio_.filecfg_kwd("DEFAULT", "VOLUME1", PSIF_CHKPT, "./");
	}
	// this libpsio object is the default object
	psi::_default_psio_lib_ = &psio_;
}

void PsiExEnv::add_to_path(const string& dir)
{
  if (dir.size()) {
    char *path = getenv("PATH");
    int newpath_len = strlen(path) + dir.size() + 2;
    char *newpath = new char[newpath_len];
    sprintf(newpath,"%s:%s",dir.c_str(),path);
#ifdef HAVE_SETENV
    setenv("PATH",newpath,1);
#else
    string putenvstr("PATH=");
    putenvstr += newpath;
    char *putenvcstr = strcpy(new char[putenvstr.size()+1], putenvstr.c_str());
    putenv(putenvcstr);
#endif
    delete[] newpath;
  }
}

Ref<PsiInput>
PsiExEnv::get_psi_input() {
  if (psiinput_.null()) {
    char *s = new char[cwd_.size() + inputname_.size() + 2];
    sprintf(s,"%s/%s",cwd_.c_str(),inputname_.c_str());
    psiinput_ = new PsiInput(s);
    delete[] s;
  }
  return psiinput_;
}

Ref<PsiFile11>
PsiExEnv::get_psi_file11() {
  if (psifile11_.null()) {
    char* s = new char[cwd_.size() + fileprefix_.size() + file11name_.size() + 3];
    sprintf(s,"%s/%s.%s",cwd_.c_str(),fileprefix_.c_str(),file11name_.c_str());
    psifile11_ = new PsiFile11(s);
    delete[] s;
  }
  return psifile11_;
}

void PsiExEnv::run_psi()
{
  std::ostringstream oss;
  oss << "psi3 --messy";
  run_psi_module(oss.str().c_str());
}

void PsiExEnv::run_psi_module(const char *module)
{
  // can only run on node 0
  if (me_ != 0) return;

  // can't run unless input file has been created
  if (psiinput_.null()) throw ProgrammingError("input file has not been created",
                                               __FILE__, __LINE__,
                                               class_desc());

  // delete chkpt file in case it gets overwritten
  if (chkpt_) { delete chkpt_;  chkpt_ = 0; }

  std::ostringstream oss;
  oss << "cd " << cwd_ << "; pwd >> " << stdout_ << "; env >> " << stdout_ << "; " << psiprefix_ << "/" << module << " -f " << inputname_ << " -o " << outputname_
      << " -p " << fileprefix_ << " 1>> " << stdout_ << " 2>> " << stderr_;
  const int errcod = system(oss.str().c_str());
  if (errcod) {
      // errcod == -1 means fork() failed. check errno
      if (errcod == -1) {
        std::string errmsg;
        switch(errno) {
          case EAGAIN:
            errmsg = "limit on the number of processors reached";
            break;
          case ENOMEM:
            errmsg = "insufficient memory";
            break;
          default:
            errmsg = "unknown error";
        }
        perror("system(psi3)");
        throw SystemException((std::string("fork() failed: ")+errmsg).c_str(),
                              __FILE__,__LINE__);
      }
      std::ostringstream oss; oss << "PsiExEnv::run_psi_module -- module " << module << " failed";
      // clean up if wasn't a cleanup attempt already
      if (strcmp(module,"psiclean"))
        run_psiclean();
      throw SystemException(oss.str().c_str(),__FILE__,__LINE__);
  }
}

void PsiExEnv::run_psiclean()
{
  // can only run on node 0
  if (me_ != 0) return;

  // can't run unless input file has been created
  if (psiinput_.nonnull())
    psio_.purge(true);
}

void PsiExEnv::print(std::ostream&o) const
{
  o << endl;
  o << indent << "PsiExEnv:" << endl << incindent;
  o << indent << "Location of Psi: " << psiprefix_ << endl;
  o << indent << "Current Psi Working Directory: " << cwd_ << endl;
  o << indent << "Current Psi File Prefix: " << fileprefix_ << endl;
  o << indent << "Number of Scratch Groups: " << nscratch_ << endl;
  for(int i=0; i<nscratch_; i++)
    o << indent << "Scratch Group " << i << ": " << scratch_[i] << endl;
  o << endl << decindent;
}

psi::Chkpt&
PsiExEnv::chkpt() {
  if (chkpt_ == 0)
    chkpt_ = new psi::Chkpt(&psio_,PSIO_OPEN_OLD);
  return *chkpt_;
}

extern "C" char* gprgid() { return "MPQC"; }

//////////////////////////////////////////////////

PsiChkpt::PsiChkpt(const Ref<PsiExEnv>& exenv,
                   const Ref<Integral>& integral,
                   int debug) :
                   exenv_(exenv),
                   integral_(integral),
                   debug_(debug) {
}

PsiChkpt::~PsiChkpt() {
}

RefDiagSCMatrix
PsiChkpt::evals(SpinCase1 spin,
                bool spin_restricted) const {

  // grab orbital info
  const int num_mo = exenv()->chkpt().rd_nmo();
  const int nirrep = exenv()->chkpt().rd_nirreps();
  int* mopi = exenv()->chkpt().rd_orbspi();
  // get the eigenvalues
  double* E;
  if (spin_restricted)
    E = exenv()->chkpt().rd_evals();
  else {
    E = (spin == Alpha) ? exenv()->chkpt().rd_alpha_evals() : exenv()->chkpt().rd_beta_evals();
    if (E == 0)
      E = exenv()->chkpt().rd_evals();
  }

  // convert raw matrices to SCMatrices
  RefSCDimension modim = new SCDimension(num_mo,nirrep,mopi);
  for (unsigned int h=0; h<nirrep; ++h)
    modim->blocks()->set_subdim(h, new SCDimension(mopi[h]));
  RefDiagSCMatrix result = integral()->basis1()->so_matrixkit()->diagmatrix(modim);
  result.assign(E);
  if (debug() >= DefaultPrintThresholds::mostN)
    result.print(prepend_spincase(spin,"Psi3 SCF eigenvalues").c_str());

  psi::Chkpt::free(E);
  psi::Chkpt::free(mopi);

  return result;
}

RefSCMatrix
PsiChkpt::coefs(SpinCase1 spin,
                bool spin_restricted) const {

  const Ref<GaussianBasisSet>& bs = integral()->basis1();
  Ref<SCMatrixKit> sokit = bs->so_matrixkit();

  psi::PSIO& psio = exenv()->psio();
  // grab orbital info
  const int num_so = exenv()->chkpt().rd_nso();
  const int num_mo = exenv()->chkpt().rd_nmo();
  const int nirrep = exenv()->chkpt().rd_nirreps();
  int* mopi = exenv()->chkpt().rd_orbspi();
  int* sopi = exenv()->chkpt().rd_sopi();
  // get MO coefficients in SO basis
  double** C;
  if (spin_restricted)
    C = exenv()->chkpt().rd_scf();
  else {
    C = (spin == Alpha) ? exenv()->chkpt().rd_alpha_scf() : exenv()->chkpt().rd_beta_scf();
    if (C == 0)
      C = exenv()->chkpt().rd_scf();
  }

  // get AO->SO matrix (MPQC AO equiv PSI3 BF)
  double** ao2so = exenv()->chkpt().rd_usotbf();

  // convert raw matrices to SCMatrices
  RefSCDimension sodim_nb = new SCDimension(num_so,1);
  sodim_nb->blocks()->set_subdim(0, new SCDimension(num_so));
  RefSCDimension sodim = new SCDimension(num_so,nirrep,sopi);
  for (unsigned int h=0; h<nirrep; ++h)
    sodim->blocks()->set_subdim(h, new SCDimension(sopi[h]));
  RefSCDimension modim = new SCDimension(num_mo,nirrep,mopi);
  for (unsigned int h=0; h<nirrep; ++h)
    modim->blocks()->set_subdim(h, new SCDimension(mopi[h]));
  RefSCMatrix C_so = sokit->matrix(sodim, modim);
  C_so.assign(C[0]);
  if (debug() >= DefaultPrintThresholds::allN2)
    C_so.print(prepend_spincase(spin,"Psi3 eigenvector in SO basis").c_str());
  RefSCMatrix aotoso = sokit->matrix(sodim, sodim_nb);
  aotoso.assign(ao2so[0]);
  Ref<PetiteList> plist = integral()->petite_list();
  if (debug() >= DefaultPrintThresholds::allN2) {
    aotoso.print("Psi3 SO->AO matrix");
    plist->sotoao().print("MPQC SO->AO matrix");
    plist->aotoso().print("MPQC AO->SO matrix");
  }
  RefSCMatrix result = aotoso.t() * C_so;
  if (debug() >= DefaultPrintThresholds::allN2)
    result.print(prepend_spincase(spin,"Psi3 eigenvector in AO basis (Psi-ordered)").c_str());

  // shells in Psi3 do not have to follow same order as atoms, as they do in MPQC
  // resort AOs from Psi to MPQC order using Psi3->MPQC shell map
  {
    std::vector<unsigned int> aomap = ao_map();
    const int nao = aomap.size();
    const char* name = "PsiSCF::evecs";
    RefSCMatrix coefs_mpqc = result.clone();
    BlockedSCMatrix* coefs_psi_blkd = require_dynamic_cast<BlockedSCMatrix*>(result.pointer(),name);
    BlockedSCMatrix* coefs_mpqc_blkd = require_dynamic_cast<BlockedSCMatrix*>(coefs_mpqc.pointer(),name);
    for (unsigned int h=0; h<nirrep; ++h) {
      RefSCMatrix coefs_mpqc_blk = coefs_mpqc_blkd->block(h);
      if (coefs_mpqc_blk.null()) continue;
      RefSCMatrix coefs_psi_blk = coefs_psi_blkd->block(h);

      for (unsigned int aopsi=0; aopsi<nao; ++aopsi) {
        RefSCVector row = coefs_psi_blk.get_row(aopsi);
        const unsigned int aompqc = aomap[aopsi];
        coefs_mpqc_blk.assign_row(row, aompqc);
      }
    }
    result = coefs_mpqc;
  }
  if (debug() >= DefaultPrintThresholds::allN2)
    result.print(prepend_spincase(spin,"Psi3 eigenvector in AO basis (MPQC-ordered)").c_str());

  // Psi3 also uses the symmetry frame for the molecule, whereas MPQC may use a different frame
  // rotate the Psi3 eigenvector from the symmetry frame to the MPQC frame
  {
    SymmetryOperation rr = bs->molecule()->point_group()->symm_frame();  rr.transpose();
    const int nshell = bs->nshell();
    const char* name = "PsiSCF::evecs";
    BlockedSCMatrix* coefs_blkd = require_dynamic_cast<BlockedSCMatrix*>(result.pointer(),name);
    double* tmpvec_orig = new double[bs->nbasis()];
    double* tmpvec_tformed = new double[bs->nbasis()];
    for(int s=0; s<nshell; ++s) {
      const GaussianShell& shell = bs->shell(s);
      const int ncontr = shell.ncontraction();
      // transform each contraction separately
      for(int c=0; c<ncontr; ++c) {
        const int am = shell.am(c);
        // am=0 functions are invariant to rotations
        if (am == 0) continue;
        ShellRotation sr = integral()->shell_rotation(am,rr,shell.is_pure(c));
        const int nf = sr.dim();
        const int foff = bs->shell_to_function(s) + shell.contraction_to_function(c);
        // in each block
        for (unsigned int h=0; h<nirrep; ++h) {
          RefSCMatrix coefs_blk = coefs_blkd->block(h);
          if (coefs_blk.null()) continue;
          const int ncol = coefs_blk.coldim().n();
          // transform each vector
          for(int col=0; col<ncol; ++col) {
            // initialize original vector
            for(int f=0; f<nf; ++f)
              tmpvec_orig[f] = coefs_blk.get_element(f+foff,col);
            // transform
            for(int f=0; f<nf; ++f) {
              double tmp = 0.0;
              for(int g=0; g<nf; ++g) {
                tmp += sr(f,g) * tmpvec_orig[g];
              }
              tmpvec_tformed[f] = tmp;
            }
            // copy to the original location
            for(int f=0; f<nf; ++f)
              coefs_blk.set_element(f+foff,col,tmpvec_tformed[f]);
          }
        }
      }
    }
    delete[] tmpvec_orig; delete[] tmpvec_tformed;
  }

  // lastly, change the dimensions to match those used by SCF classes (AO dimension much have subdimension blocked by shells)
  {
    RefSCMatrix coefs_redim = result.kit()->matrix(plist->AO_basisdim(), result.coldim());
    RefSCMatrix coefs = result;
    const int nrow = coefs.nrow();
    const int ncol = coefs.ncol();
    for(int r=0; r<nrow; ++r) {
      for(int c=0; c<ncol; ++c) {
        coefs_redim.set_element(r, c, coefs.get_element(r,c) );
      }
    }
    result = coefs_redim;
  }

  if (debug() >= DefaultPrintThresholds::allN2)
    result.print(prepend_spincase(spin,"Psi3 eigenvector in AO basis (MPQC-ordered, in MPQC frame)").c_str());

  using psi::Chkpt;
  Chkpt::free(mopi);
  Chkpt::free(sopi);
  Chkpt::free(C);
  Chkpt::free(ao2so);

  return result;
}

std::vector< std::pair<unsigned int, unsigned int> > PsiChkpt::shell_map() const {
  typedef std::vector< std::pair<unsigned int, unsigned int> > ShellMap;
  const Ref<GaussianBasisSet>& bs = integral()->basis1();
  // # of MPQC contractions = # of Psi3 shells
  const unsigned int ncontr = exenv()->chkpt().rd_nshell();
  // # of MPQC shells
  const unsigned int nshells = bs->nshell();
  ShellMap map(ncontr);
  int* snuc = exenv()->chkpt().rd_snuc();
  // ordering of shells on an atom is the same in Psi3 and MPQC
  // but shells si and sj from different atoms (i<j) may be ordered differently (si > sj)
  int first_shell_on_curr_atom = 0;
  int atom_curr = snuc[0] - 1;
  int contr = 0;
  for(unsigned int s=0; s<nshells; ++s) {
    int atom = snuc[contr] - 1;
    if (atom != atom_curr) {
      atom_curr = atom;
      first_shell_on_curr_atom = s;
    }
    const unsigned int shell_mpqc = bs->shell_on_center(atom,s-first_shell_on_curr_atom);
    const unsigned int ncontr_in_shell = bs->shell(shell_mpqc).ncontraction();
    for(unsigned int c=0; c<ncontr_in_shell; ++c, ++contr) {
      map[contr] = make_pair(shell_mpqc,c);
    }
  }
  psi::Chkpt::free(snuc);
  return map;
}

std::vector<unsigned int> PsiChkpt::ao_map() const {
  typedef std::pair<unsigned int, unsigned int> ShellContrPair;
  typedef std::vector<ShellContrPair> ShellMap;
  const Ref<GaussianBasisSet>& bs = integral()->basis1();
  const unsigned int nao = bs->nbasis();
  psi::Chkpt& chkpt = exenv()->chkpt();
  ShellMap smap = shell_map();

  std::vector<unsigned int> map(nao);
  const unsigned int nshells_psi = chkpt.rd_nshell();
  int* first_ao_from_shell_psi = chkpt.rd_puream() ? chkpt.rd_sloc_new() : chkpt.rd_sloc();
  for(unsigned int spsi=0; spsi<nshells_psi; ++spsi) {
    const ShellContrPair& shellcontr = smap[spsi];
    unsigned int smpqc = shellcontr.first;
    unsigned int cmpqc = shellcontr.second;
    const GaussianShell& shell = bs->shell(smpqc);
    unsigned int nbf = shell.nfunction(cmpqc);
    int first_bf_psi = first_ao_from_shell_psi[spsi] - 1;
    int first_bf_mpqc = bs->shell_to_function(smpqc) + shell.contraction_to_function(cmpqc);
    for(unsigned int bf=0; bf<nbf; ++bf) {
      map[first_bf_psi + bf] = first_bf_mpqc + bf;
    }
  }
  psi::Chkpt::free(first_ao_from_shell_psi);
  return map;
}

