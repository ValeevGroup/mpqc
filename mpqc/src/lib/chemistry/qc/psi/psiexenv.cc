
#ifdef __GNUG__
#pragma implementation
#endif

#include <string>
#include <string.h>
#include <stdlib.h>
#include <scconfig.h>
#include <util/ref/ref.h>
#include <util/keyval/keyval.h>
#include <util/misc/formio.h>
#include <chemistry/qc/psi/psiexenv.h>
#include <psifiles.h>

using namespace std;

namespace sc {

static ClassDesc PsiExEnv_cd(
  typeid(PsiExEnv),"PsiExEnv",1,"public DescribedClass",
  0, create<PsiExEnv>, 0);

string PsiExEnv::defaultinputname_("input.dat");
string PsiExEnv::defaultoutputname_("output.dat");
string PsiExEnv::file11name_("file11.dat");
int PsiExEnv::ckptfile_(PSIF_CHKPT);
string PsiExEnv::defaultcwd_("/tmp");
string PsiExEnv::defaultfileprefix_("psi");
string PsiExEnv::defaultpsiprefix_("/usr/local/psi/bin");
string PsiExEnv::defaultstdout_("stdout");
string PsiExEnv::defaultstderr_("stderr");

PsiExEnv::PsiExEnv(const Ref<KeyVal>& keyval) :
	psio_()
{
  const std::string prefix(SCFormIO::fileext_to_filename("."));

  // Find Psi
  char *psibin = getenv("PSIBIN");
  if (psibin)
    psiprefix_ = string(psibin);
  else
    psiprefix_ = string(defaultpsiprefix_);
  add_to_path(psiprefix_);

  char *cwdchar = keyval->pcharvalue("cwd");
  if (cwdchar)
    cwd_ = string(cwdchar);
  else
    cwd_ = string(defaultcwd_);
  char *fileprefixchar = keyval->pcharvalue("fileprefix");
  if (fileprefixchar)
    fileprefix_ = string(fileprefixchar);
  else
    fileprefix_ = prefix + defaultfileprefix_;

  inputname_ = fileprefix_ + "." + defaultinputname_;
  outputname_ = fileprefix_ + "." + defaultoutputname_;

  char *stdout_char = keyval->pcharvalue("stdout");
  if (stdout_char)
    stdout_ = string(stdout_char);
  else
    stdout_ = fileprefix_ + "." + defaultstdout_;
  delete[] stdout_char;
  char *stderr_char = keyval->pcharvalue("stderr");
  if (stderr_char)
    stderr_ = string(stderr_char);
  else
    stderr_ = fileprefix_ + "." + defaultstderr_;
  delete[] stderr_char;

  nscratch_ = keyval->intvalue("nscratch");
  if (nscratch_ != keyval->count("scratch")) {
    ExEnv::err0() << indent
		  << "PsiExEnv: number of scratch directories != nscratch\n";
    abort();
  }
  scratch_ = new string[nscratch_];
  for (int i=0; i<nscratch_; i++)
    scratch_[i] = string(keyval->pcharvalue("scratch",i));

  char *s = new char[cwd_.size() + inputname_.size() + 2];
  sprintf(s,"%s/%s",cwd_.c_str(),inputname_.c_str());
  psiinput_ = new PsiInput(s);
  delete[] s;

  s = new char[cwd_.size() + fileprefix_.size() + file11name_.size() + 3];
  sprintf(s,"%s/%s.%s",cwd_.c_str(),fileprefix_.c_str(),file11name_.c_str());
  psifile11_ = new PsiFile11(s);
  delete[] s;

  config_psio();
  chkpt_ = new psi::Chkpt(&psio_,PSIO_OPEN_OLD);
}

PsiExEnv::PsiExEnv(char *cwd, char *fileprefix, int nscratch, char **scratch):
    cwd_(cwd), fileprefix_(fileprefix), nscratch_(nscratch)
{
  const std::string prefix(SCFormIO::fileext_to_filename("."));

  // Find Psi
  char *psibin = 0;
  psibin = getenv("PSIBIN");
  if (psibin)
    psiprefix_ = string(psibin);
  else
    psiprefix_ = prefix + defaultpsiprefix_;
  add_to_path(psiprefix_);

  scratch_ = new string[nscratch_];
  for(int i=0; i<nscratch_; i++)
    scratch_[i] = string(scratch[i]);
  
  char *s = new char[cwd_.size() + inputname_.size() + 2];
  sprintf(s,"%s/%s",cwd_.c_str(),inputname_.c_str());
  psiinput_ = new PsiInput(s);
  delete[] s;

  s = new char[cwd_.size() + fileprefix_.size() + file11name_.size() + 3];
  sprintf(s,"%s/%s.%s",cwd_.c_str(),fileprefix_.c_str(),file11name_.c_str());
  psifile11_ = new PsiFile11(s);
  delete[] s;

  config_psio();
  chkpt_ = new psi::Chkpt(&psio_,PSIO_OPEN_OLD);
}

PsiExEnv::~PsiExEnv()
{
  delete[] scratch_;
  delete chkpt_;
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

int PsiExEnv::run_psi()
{
  std::ostringstream oss;
  oss << "psi3 --messy";
  int errcod;
  if (errcod = run_psi_module(oss.str().c_str())) {
    return errcod;
  }
  return 0;
}

int PsiExEnv::run_psi_module(const char *module)
{
  std::ostringstream oss;
  oss << "cd " << cwd_ << "; " << psiprefix_ << "/" << module << " -f " << inputname_ << " -o " << outputname_
      << " -p " << fileprefix_ << " 1>> " << stdout_ << " 2>> " << stderr_;
  int errcod;
  if (errcod = system(oss.str().c_str())) {
      ExEnv::outn() << "PsiExEnv::run_psi_module -- module " << module << " failed" << endl;
      abort();
  }
  return errcod;
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

}

extern "C" char* gprgid() { return "MPQC"; }
