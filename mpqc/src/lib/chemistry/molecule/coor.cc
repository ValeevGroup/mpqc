
#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
#include <math.h>
};

#include <math/scmat/matrix.h>
#include <math/scmat/local.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>

///////////////////////////////////////////////////////////////////////////
// members of IntCoor

double IntCoor::bohr_conv = 0.52917706;
double IntCoor::radian_conv = 180.0/3.14159265358979323846;

SavableState_REF_def(IntCoor);
ARRAY_def(RefIntCoor);
SET_def(RefIntCoor);
ARRAYSET_def(RefIntCoor);

#define CLASSNAME IntCoor
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
IntCoor::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

IntCoor::IntCoor(const char *re):
  value_(0.0), label_(0)
{
  if (!re) re = "noname";
  label_=new char[strlen(re)+1]; strcpy(label_,re);
}

IntCoor::IntCoor(const IntCoor& c):
  label_(0)
{
  value_ = c.value_;
  if (c.label_) label_ = strcpy(new char[strlen(c.label_)+1],c.label_);
}

IntCoor::IntCoor(KeyVal&keyval)
{
  label_ = keyval.pcharvalue("label");
  value_ = keyval.doublevalue("value");
}

IntCoor::IntCoor(StateIn& si):
  SavableState(si,IntCoor::class_desc_)
{
  si.get(value_);
  si.getstring(label_);
}

IntCoor::~IntCoor()
{
  if (label_) delete[] label_;
}

void
IntCoor::save_data_state(StateOut& so)
{
  so.put(value_);
  so.putstring(label_);
}

const char*
IntCoor::label() const
{
  return label_;
}

double
IntCoor::value() const
{
  return value_;
}

#ifndef __GNUC__
void
IntCoor::print()
{
  print(0);
}
#endif

void
IntCoor::print(RefMolecule mol, SCostream& os)
{
  os.setf(ios::fixed,ios::floatfield);
  os.precision(10);
  os.setf(ios::left,ios::adjustfield);
  os.width(10);

  os.indent() << ctype()
              << " \""
              << label()
              << "\" "
              << preferred_value()
              << endl;
}

double
IntCoor::preferred_value() const
{
  return value_;
}

///////////////////////////////////////////////////////////////////////////
// members of SetIntCoor

SavableState_REF_def(SetIntCoor);

#define CLASSNAME SetIntCoor
#define PARENTS virtual public SavableState
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
SetIntCoor::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

SetIntCoor::SetIntCoor()
{
}

SetIntCoor::SetIntCoor(KeyVal& keyval)
{
  int n = keyval.count();
  if (!n) {
      fprintf(stderr,"SetIntCoor::SetIntCoor: bad input\n");
      abort();
    }

  for (int i=0; i<n; i++) {
      coor_.add(keyval.describedclassvalue(i));
    }
}

SetIntCoor::SetIntCoor(StateIn& s):
  SavableState(s,SetIntCoor::class_desc_)
{
  int n;
  s.get(n);

  RefIntCoor tmp;
  for (int i=0; i<n; i++) {
      tmp.restore_state(s);
      coor_.add(tmp);
    }
}

SetIntCoor::~SetIntCoor()
{
}

void
SetIntCoor::save_data_state(StateOut& s)
{
  int n = coor_.length();
  s.put(n);

  for (int i=0; i<n; i++) {
      coor_[i].save_state(s);
    }
}

void
SetIntCoor::add(const RefIntCoor& coor)
{
  coor_.add(coor);
}

void
SetIntCoor::add(const RefSetIntCoor& coor)
{
  for (int i=0; i<coor->n(); i++) {
      coor_.add(coor->coor(i));
    }
}

int
SetIntCoor::n() const
{
  return coor_.length();
}

RefIntCoor
SetIntCoor::coor(int i) const
{
  return coor_[i];
}

// compute the bmatrix by finite displacements
void
SetIntCoor::fd_bmat(RefMolecule& mol,RefSCMatrix& fd_bmatrix)
{
  fd_bmatrix.assign(0.0);
  
  int i;
  Molecule& m = * mol.pointer();

  const double cart_disp = 0.01;

  RefSCDimension dn3(fd_bmatrix.coldim());
  RefSCDimension dnc(fd_bmatrix.rowdim());
  int n3 = dn3.n();
  int nc = dnc.n();
  RefSCVector internal(dnc);
  RefSCVector internal_p(dnc);
  RefSCVector internal_m(dnc);

  // the internal coordinates
  update_values(mol);
  for (i=0; i<nc; i++) {
      internal(i) = coor_[i]->value();
    }

  // the finite displacement bmat
  for (i=0; i<n3; i++) {
      // the plus displacement
      m[i/3][i%3] += cart_disp;
      update_values(mol);
      for (int j=0; j<nc; j++) {
          internal_p(j) = coor_[j]->value();
        }
      // the minus displacement
      m[i/3][i%3] -= 2.0*cart_disp;
      update_values(mol);
      for (j=0; j<nc; j++) {
          internal_m(j) = coor_[j]->value();
        }
      // reset the cartesian coordinate to its original value
      m[i/3][i%3] += cart_disp;

      // construct the entries in the finite displacement bmat
      for (j=0; j<nc; j++) {
          fd_bmatrix(j,i) = (internal_p(j)-internal_m(j))/(2.0*cart_disp);
        }
    }
}

void
SetIntCoor::bmat(RefMolecule& mol, RefSCMatrix& bmat)
{
  bmat.assign(0.0);

  int i, ncoor = n(), ncart = bmat.coldim().n();

  RefSCVector bmatrow(bmat.coldim());
  // send the rows of the b matrix to each of the coordinates
  for (i=0; i<ncoor; i++) {
      bmatrow.assign(0.0);
      coor_[i]->bmat(mol,bmatrow);
      for (int j=0; j<ncart; j++) bmat(i,j) = bmatrow(j);
    }
}

void
SetIntCoor::guess_hessian(RefMolecule& mol,RefSymmSCMatrix& hessian)
{
  int ncoor = hessian.n();

  hessian.assign(0.0);
  for (int i=0; i<ncoor; i++) {
      hessian(i,i) = coor_[i]->force_constant(mol);
    }
}

#ifndef __GNUC__
void
SetIntCoor::print()
{
  print(0);
}
#endif

void
SetIntCoor::print(RefMolecule mol, SCostream& os)
{
  int i;

  for(i=0; i<coor_.length(); i++) {
      coor_[i]->print(mol,os);
    }
}

void
SetIntCoor::update_values(RefMolecule&mol)
{
  for (int i=0; i<coor_.length(); i++) {
      coor_[i]->update_value(mol);
    }
}

void
SetIntCoor::values_to_vector(RefSCVector&v)
{
  for (int i=0; i<coor_.length(); i++) {
      v(i) = coor_[i]->value();
    }
}  

void
SetIntCoor::clear()
{
  coor_.clear();
}

///////////////////////////////////////////////////////////////////////////
// members of SumIntCoor

#define CLASSNAME SumIntCoor
#define PARENTS public IntCoor
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
SumIntCoor::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = IntCoor::_castdown(cd);
  return do_castdowns(casts,cd);
}

SumIntCoor::SumIntCoor(const char* label):
  IntCoor(label)
{
}

SumIntCoor::SumIntCoor(KeyVal&keyval):
  IntCoor(keyval)
{
  static const char* coor = "coor";
  static const char* coef = "coef";
  int n = keyval.count(coor);
  int ncoef = keyval.count(coef);
  if (n != ncoef || !n) {
      fprintf(stderr,"SumIntCoor::SumIntCoor: bad input\n");
      abort();
    }

  for (int i=0; i<n; i++) {
      double coe = keyval.doublevalue(coef,i);
      RefIntCoor coo = keyval.describedclassvalue(coor,i);
      add(coo,coe);
    }
}

SumIntCoor::SumIntCoor(StateIn&s):
  SavableState(s,SumIntCoor::class_desc_),
  IntCoor(s)
{
  int n;
  s.get(n);

  coef_.set_length(n);
  coor_.set_length(n);
  for (int i=0; i<n; i++) {
      s.get(coef_[i]);
      coor_[i].restore_state(s);
    }
}

SumIntCoor::~SumIntCoor()
{
}

void
SumIntCoor::save_data_state(StateOut&s)
{
  int n = coef_.length();
  IntCoor::save_data_state(s);
  s.put(coef_.length());

  for (int i=0; i<n; i++) {
      s.put(coef_[i]);
      coor_[i].save_state(s);
    }
}

int
SumIntCoor::n()
{
  return coef_.length();
}

void
SumIntCoor::add(RefIntCoor&coor,double coef)
{
  // if a sum is added to a sum, unfold the nested sum
  SumIntCoor* scoor = SumIntCoor::castdown(coor.pointer());
  if (scoor) {
      int l = scoor->coor_.length();
      for (int i=0; i<l; i++) {
          add(scoor->coor_[i],coef * scoor->coef_[i]);
        }
    }
  else {
      int l = coef_.length();
      for (int i=0; i<l; i++) {
          if (coor_[i]->equivalent(coor)) {
              coef_[i] += coef;
              return;
            }
        }
      coef_.reset_length(l+1);
      coor_.reset_length(l+1);
      coef_[l] = coef;
      coor_[l] = coor;
    }
}

int
SumIntCoor::equivalent(RefIntCoor&c)
{
  return 0;
}

// this normalizes and makes the biggest coordinate positive
void
SumIntCoor::normalize()
{
  int i;
  int n = coef_.length();
  double norm = 0.0;

  double biggest = 0.0;
  for (i=0; i<n; i++) {
      norm += coef_[i] * coef_[i];
      if (fabs(biggest) < fabs(coef_[i])) biggest = coef_[i];
    }
  norm = (biggest < 0.0? -1.0:1.0)/sqrt(norm);

  for (i=0; i<n; i++) {
      coef_[i] = coef_[i]*norm;
    }
}

double
SumIntCoor::preferred_value() const
{
  return value_;
}

const char*
SumIntCoor::ctype() const
{
  return "SUM";
}

#ifndef __GNUC__
void
SumIntCoor::print()
{
  print(0);
}
#endif

void
SumIntCoor::print(RefMolecule mol, SCostream& os)
{
  os.setf(ios::fixed,ios::floatfield);
  os.precision(10);
  os.setf(ios::left,ios::adjustfield);

  int initial_indent = os.get_indent();
  int i;

  os.width(5);
  os.indent() << ctype()
              << " ";
  os.width(10);
  os          << label()
              << " ";
  os.width(14);
  os.setf(ios::right,ios::adjustfield);
  os          << preferred_value()
              << endl;

  for(i=0; i<coor_.length(); i++) {
      os++;
      os.width(14);
      os.indent() << coef_[i] << " ";
      os.set_indent_to_column(); os.skip_next_indent();
      coor_[i]->print(mol,os);
      os.set_indent(initial_indent);
    }
}

// the SumIntCoor should be normalized before this is called.
double
SumIntCoor::force_constant(RefMolecule&molecule)
{
  double fc = 0.0;
  
  for (int i=0; i<n(); i++) {
      fc += coef_[i] * coef_[i] * coor_[i]->force_constant(molecule);
    }

  return fc;
}

void
SumIntCoor::update_value(RefMolecule&molecule)
{
  int i, l = n();

  value_ = 0.0;
  for (i=0; i<l; i++) {
      coor_[i]->update_value(molecule);
      value_ += coef_[i] * coor_[i]->value();
    }

}

void
SumIntCoor::bmat(RefMolecule&molecule,RefSCVector&bmat,double coef)
{
  int i, l = n();
  
  for (i=0; i<l; i++) {
      coor_[i]->bmat(molecule,bmat,coef*coef_[i]);
    }
}

///////////////////////////////////////////////////////////////////////////
// members of MolecularCoor

SavableState_REF_def(MolecularCoor);

#define CLASSNAME MolecularCoor
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
MolecularCoor::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

MolecularCoor::MolecularCoor(RefMolecule&mol):
  molecule_(mol)
{
}

MolecularCoor::MolecularCoor(KeyVal&keyval)
{
  molecule_ = keyval.describedclassvalue("molecule");
}

MolecularCoor::MolecularCoor(StateIn&s):
  SavableState(s,MolecularCoor::class_desc_)
{
  molecule_.restore_state(s);
}

MolecularCoor::~MolecularCoor()
{
}

void
MolecularCoor::save_data_state(StateOut&s)
{
  molecule_.save_state(s);
}
