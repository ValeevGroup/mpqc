/*
**
** PSI Input Class
**
** This helper class will set up input decks for the PSI suite of
** ab initio quantum chemistry programs.
**
** David Sherrill & Justin Fermann
** Center for Computational Quantum Chemistry, University of Georgia
**
*/

#ifdef __GNUG__
#pragma implementation
#endif

#include <iostream>

#include <util/misc/formio.h>
#include <math/symmetry/corrtab.h>
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/atominfo.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/psi/psiexenv.h>
#include <chemistry/qc/psi/psiinput.h>

using namespace std;

PsiInput::PsiInput(const string& name) : file_()
{
  filename_ = string(name);
  indentation_ = 0;
}

PsiInput::~PsiInput()
{
}

void
PsiInput::open()
{
  file_.open(filename_.c_str(),ios::out);
  indentation_ = 0;
}

void
PsiInput::close()
{
  file_.close();
  indentation_ = 0;
}


void 
PsiInput::write_indent()
{
  for (int i=0; i<indentation_; i++)
    file_ << " ";
}

void
PsiInput::begin_section(const char * s)
{
   write_indent();
   indentation_ += 2;
   file_ << s << ":(" << endl;
}

void
PsiInput::end_section(void)
{
   indentation_ -= 2;
   write_indent();
   file_ << ")" << endl;
}

void
PsiInput::write_keyword(const char *keyword, const char *value)
{
   write_indent();
   file_ << scprintf("%s = %s",keyword,value) << endl;
}

void
PsiInput::write_keyword(const char *keyword, int value)
{
   write_indent();
   file_ << scprintf("%s = %d",keyword,value) << endl;
}

void
PsiInput::write_keyword(const char *keyword, double value)
{
   write_indent();
   file_ << scprintf("%s = %20.15lf",keyword,value) << endl;
}

void
PsiInput::write_keyword_array(const char *keyword, int num, int *values)
{
  write_indent();
  file_ << scprintf("%s = (", keyword);
  for (int i=0; i<num; i++) {
    file_ << scprintf(" %d", values[i]);
  }
  file_ << ")" << endl;
}

void
PsiInput::write_keyword_array(const char *keyword, int num, double *values)
{
   write_indent();
   file_ << scprintf("%s = (", keyword);
  for (int i=0; i<num; i++) {
    file_ << scprintf(" %20.15lf", values[i]);
  }
  file_ << ")" << endl;
}

void 
PsiInput::write_string(const char *s)
{
   write_indent();
   file_ << s;
}

void
PsiInput::write_key_wq(const char *keyword, const char *value)
{
   write_indent();
   file_ << scprintf("%s = \"%s\"", keyword, value) << endl;
}


void
PsiInput::write_geom(const Ref<Molecule>& mol)
{
  write_string("geometry = (\n");
  for (int i=0; i < mol->natom(); i++) {
    write_string("  (");
    char *s;
    file_ << AtomInfo::symbol(mol->Z(i)) <<
	scprintf(" %14.12lf %14.12lf %14.12lf",mol->r(i,0),mol->r(i,1),mol->r(i,2))
	  << ")" << endl;
  } 
  write_string(")\n");
}


void
PsiInput::write_basis(const Ref<GaussianBasisSet>& basis)
{
  begin_section("basis");
  Ref<AtomInfo> atominfo = basis->molecule()->atominfo(); 
  end_section();
}

void
PsiInput::write_defaults(const Ref<PsiExEnv>& exenv, const char *wfn,
			 const char *dertype)
{
  begin_section("default");
  
  write_key_wq("label"," ");
  write_keyword("wfn",wfn);
  write_keyword("reference","rhf");
  write_keyword("dertype",dertype);
  begin_section("files");
  begin_section("default");
  write_key_wq("name",(exenv->get_fileprefix()).c_str());
  int nscratch = exenv->get_nscratch();
  write_keyword("nvolume",nscratch);
  char *scrname; scrname = new char[10];
  for(int i=0; i<nscratch; i++) {
    sprintf(scrname,"volume%d",i+1);
    write_key_wq(scrname,(exenv->get_scratch(i)).c_str());
  }
  delete[] scrname;
  end_section();
  write_string("file30: ( nvolume = 1 volume1 = \"./\" )");
  end_section();

  end_section();
}


void
PsiInput::print(ostream& o)
{
}

void
PsiInput::write_input_file(const char *dertype, const char *wavefn, 
    const int convergence, const char *fname )
{
}

void
PsiInput::write_input(void)
{
}

