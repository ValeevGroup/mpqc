
#ifndef _chemistry_qc_basis_petite_h
#define _chemistry_qc_basis_petite_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/gaussbas.h>

class PetiteList {
  private:
    int _natom;
    int _nshell;
    int _ng;

    RefGaussianBasisSet _gbs;
    
    char *_p1;        // p1[n] is 1 if shell n is in the group P1
    int **_atom_map;  // atom_map[n][g] is the atom that symop g maps atom n
                     // into
    int **_shell_map; // shell_map[n][g] is the shell that symop g maps shell n
                     // into
    char *_lamij;     // see Dupuis & King, IJQC 11,613,(1977)

    inline int ioff(int i) const { return i*(i+1)>>1; }
    inline int ioff(int i, int j) const
      { return (i>=j) ? ioff(i)+j : ioff(j)+i; }
    
  public:
    PetiteList();
    PetiteList(const RefGaussianBasisSet&);
    ~PetiteList();

    void init(const RefGaussianBasisSet&);

    int in_p1(int n) const { return (int) _p1[n]; }
    int atom_map(int n, int g) const { return _atom_map[n][g]; }
    int shell_map(int n, int g) const { return _shell_map[n][g]; }
    int lambda(int ij) const { return (int) _lamij[ij]; }
    int lambda(int i, int j) const { return (int) _lamij[ioff(i,j)]; }

    int in_p4(int ij, int kl, int i, int j, int k, int l) const;
    
    void print(FILE* =stdout);

    RefSCMatrix r(int g);
};

inline int
PetiteList::in_p4(int ij, int kl, int i, int j, int k, int l) const
{
  int ijkl = ioff(ij)+kl;
  int nijkl=0;

  for (int g=0; g < _ng; g++) {
    int gij = ioff(_shell_map[i][g],_shell_map[j][g]);
    int gkl = ioff(_shell_map[k][g],_shell_map[l][g]);
    int gijkl = ioff(gij,gkl);

    if (gijkl > ijkl)
      return 0;
    else if (gijkl == ijkl)
      nijkl++;
  }

  return _ng/nijkl;
}

#endif
    
