
#ifndef _chemistry_qc_basis_gaussbas_h
#define _chemistry_qc_basis_gaussbas_h

#include <stdio.h>
#include <math/newmat7/newmat.h>
#include <util/container/pixmap.h>

// this is for centers_t*
#include <chemistry/qc/intv2/atoms.h>

class GaussianShell;
class KeyVal;
class Molecule;

// This class is complicated by the use of several different sorts of Pixes:
//   Centerpixes are the pixes for the Molecule class and are used
//     as keys into the shell PixMap.  Some centerpixes may be owned
//     by the molecule class (the basis set may have basis functions that
//     are not on a shell).
//   Shellpixes select a particular shell.  These are used by the r and
//     operator[] members.
//
//   The r and operator[] members also take integer arguments for those
//     not interested in using pixes.  The pix_to_ix and ix_to_pix
//     members map between the representations.
class GaussianBasisSet
{
 private:
  char* name_;
  GaussianShell** shell;
  Pix* shell_to_centerpix;
  int* shell_to_function_;
  PixMap<int> centerpix_to_shellnum; // the number of the first shell is 0
  PixMap<double*> centerpix_to_r;
  PixMap<int> centerpix_to_nshell;
  int nshell_;
  int nbasis_;
  void recursively_get_shell(int&,KeyVal&,const char*,const char*,int,int,int);
 public:
  GaussianBasisSet(KeyVal&,Molecule&);
  virtual ~GaussianBasisSet();
  inline const char* name() const { return name_; }
  inline int nshell() const { return nshell_; }
  inline int nbasis() const { return nbasis_; }
  int max_nfunction_in_shell() const;
  int shell_to_function(Pix shellpix) const;
  inline int shell_to_function(int i) const
  { shell_to_function(ix_to_pix(i)); }

  // iteration over centers with centerpixes
  int ncenter() const;
  Pix first_center() const;
  void next_center(Pix& centerpix) const;
  int nshell_on_center(Pix centerpix) const;
  // the centershell pixes are equivalent to the shell pixes, below
  Pix first_shell_on_center(Pix centerpix) const;
  void next_shell_on_center(Pix centerpix,Pix& shellpix) const;

  // iteration over shells with shellpixes
  Pix first() const;
  void next(Pix& shellpix) const;

  // access to shells and r thru pixes
  const GaussianShell& operator[](Pix shellpix) const;
  const double& r_center(Pix centerpix,int) const;
  const double& r_shell(Pix shellpix,int) const;
  GaussianShell& operator[](Pix shellpix);
  double& r_center(Pix centerpix,int);
  double& r_shell(Pix shellpix,int);

  // map between shellpixes and integer indices for shells
  inline Pix ix_to_pix(int i) const { return (Pix) (i+1); }
  inline int pix_to_ix(Pix i) const { return ((int)i)-1; }

  // access to shells and r thru ints
  inline const GaussianShell& operator[](int i) const
  { return (*this)[ix_to_pix(i)]; }
  inline const double& r(int i,int xyz) const
  { return this->r_shell(ix_to_pix(i),xyz); }
  inline GaussianShell& operator[](int i)
  { return (*this)[ix_to_pix(i)]; }
  inline double& r(int i,int xyz)
  { return this->r_shell(ix_to_pix(i),xyz); }

  // this operator converts a shell pix to a center pix
  Pix shell_to_center(Pix) const;

  // converts the basis set to a centers_t for compatibility with libintv2
  // If the molecule is 0 then fake molecule information is included--
  // this is useful in computing the overlap, for example.
  centers_t* convert_to_centers_t(const Molecule*) const;
  //operator struct struct_centers*();

  // compute the value for this basis set at position r
  int values(cart_point& r, double* basis_values) const;
  int grad_values(cart_point& r, double*g_values, double* basis_values=0)const;

  // fill in matrix with a matrix that orthogonalizes the basis functions
  // note: this member is provided in the integrals library
  void ortho(Matrix&ortho) const;
  void ortho(Matrix&ortho,Matrix&ortho_inverse) const;

  void print(FILE*fp=stdout) const;
};

#endif
