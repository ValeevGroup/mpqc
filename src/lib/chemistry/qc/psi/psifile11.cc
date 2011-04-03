#include <iostream>

#include <util/misc/formio.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/atominfo.h>
#include <chemistry/qc/psi/psiexenv.h>
#include <chemistry/qc/psi/psifile11.h>

using namespace std;

namespace sc {

PsiFile11::PsiFile11(const string& name) : file_()
{
  filename_ = string(name);
}

PsiFile11::~PsiFile11()
{
}

void
PsiFile11::rewind()
{
  file_.seekg(0,ios::beg);
}

void
PsiFile11::skip_lines(int n)
{
  // Lines in File11 are guaranteed to be 80 characters
  char line[100];
  for(int i=0; i<n; i++)
    file_.getline(line,100);
}

void
PsiFile11::skip_entry()
{
  skip_lines(1);
  int natom;
  file_ >> natom;
  double energy;
  file_ >> energy;
  skip_lines(2*natom);
}

void
PsiFile11::open()
{
  file_.open(filename_.c_str(),ios::in);
}

void
PsiFile11::close()
{
  if (!file_.is_open())
    file_.close();
}

void
PsiFile11::remove()
{
  if (file_.is_open())
    file_.close();
  file_.open(filename_.c_str(),ios::out | ios::trunc);
  file_.close();
}

int
PsiFile11::get_natom(int entry)
{
  skip_lines(1);

  int natom;
  file_ >> natom;
  rewind();
  return natom;
}

double
PsiFile11::get_energy(int entry)
{
  skip_lines(1);

  int natom;
  file_ >> natom;
  double energy;
  file_ >> energy;
  rewind();
  return energy;
}

double
PsiFile11::get_coord(int entry, int atom, int xyz)
{
  skip_lines(1);
  int natom;
  file_ >> natom;
  if (natom <= atom)
    abort();
  double energy;
  file_ >> energy;

  skip_lines(atom+1);
  double charge;
  file_ >> charge;
  double trash;
  for(int i=0; i<xyz; i++)
    file_ >> trash;
  double coord;
  file_ >> coord;
  
  rewind();
  return coord;
}

double
PsiFile11::get_grad(int entry, int atom, int xyz)
{
  skip_lines(1);
  int natom;
  file_ >> natom;
  if (natom <= atom)
    abort();
  double energy;
  file_ >> energy;

  skip_lines(natom+atom+1);
  double trash;
  for(int i=0; i<xyz; i++)
    file_ >> trash;
  double grad;
  file_ >> grad;

  rewind();
  return grad;
}

}
