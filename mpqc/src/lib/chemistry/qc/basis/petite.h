
#ifndef _chemistry_qc_basis_petite_h
#define _chemistry_qc_basis_petite_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/ref/ref.h>
#include <util/container/array.h>
#include <math/scmat/blocked.h>
#include <math/scmat/offset.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/integral.h>

////////////////////////////////////////////////////////////////////////////

struct contribution {
    int bfn;
    double coef;

    contribution();
    contribution(int b, double c);
    ~contribution();
};

struct SO {
    int len;
    contribution *cont;

    SO();
    SO(int);
    ~SO();

    SO& operator=(const SO&);
    
    void set_length(int);
    void reset_length(int);
    
    // is this equal to so to within a sign
    int equiv(const SO& so);
};

struct SO_block {
    int len;
    SO *so;

    SO_block();
    SO_block(int);
    ~SO_block();

    void set_length(int);
    void reset_length(int);

    int add(SO& s, int i);
    void print(const char *title);
};

////////////////////////////////////////////////////////////////////////////
// this should only be used from within a SymmGaussianBasisSet

class PetiteList : public VRefCount {
  private:
    int natom_;
    int nshell_;
    int ng_;
    int nirrep_;
    int nblocks_;

    RefGaussianBasisSet gbs_;
    RefIntegral ints_;
    
    char *p1_;        // p1[n] is 1 if shell n is in the group P1
    int **atom_map_;  // atom_map[n][g] is the atom that symop g maps atom n
                      // into
    int **shell_map_; // shell_map[n][g] is the shell that symop g maps shell n
                      // into
    char *lamij_;     // see Dupuis & King, IJQC 11,613,(1977)

    int *nbf_in_ir_;

    void init();

  public:
    PetiteList(const RefGaussianBasisSet&, const RefIntegral&);
    ~PetiteList();

    int order() const { return ng_; }
    int atom_map(int n, int g) const { return atom_map_[n][g]; }
    int shell_map(int n, int g) const { return shell_map_[n][g]; }
    int lambda(int ij) const { return (int) lamij_[ij]; }
    int lambda(int i, int j) const { return (int) lamij_[ij_offset(i,j)]; }

    int in_p1(int n) const { return (int) p1_[n]; }
    int in_p2(int ij) const { return (int) lamij_[ij]; }
    int in_p4(int ij, int kl, int i, int j, int k, int l) const;
    
    int nfunction(int i) const { return nbf_in_ir_[i]; }

    int nblocks() const { return nblocks_; }

    void print(FILE* =stdout, int verbose=1);

    // these return blocked dimensions
    RefSCDimension AO_basisdim();
    RefSCDimension SO_basisdim();

    // return the basis function rotation matrix R(g)
    RefSCMatrix r(int g);

    // return information about the transformation from AOs to SOs
    SO_block * aotoso_info();
    RefSCMatrix aotoso();
    RefSCMatrix sotoao();

    // given a skeleton matrix, form the symmetrized matrix in the SO basis
    void symmetrize(const RefSymmSCMatrix& skel, const RefSymmSCMatrix& sym);

    // transform a matrix from AO->SO or SO->AO.
    // this can take either a blocked or non-blocked AO basis matrix.
    RefSymmSCMatrix to_SO_basis(const RefSymmSCMatrix&);

    // this returns a non-blocked AO basis matrix.
    RefSymmSCMatrix to_AO_basis(const RefSymmSCMatrix&);
};

inline int
PetiteList::in_p4(int ij, int kl, int i, int j, int k, int l) const
{
  int ijkl = i_offset(ij)+kl;
  int nijkl=0;

  for (int g=0; g < ng_; g++) {
    int gij = ij_offset(shell_map_[i][g],shell_map_[j][g]);
    int gkl = ij_offset(shell_map_[k][g],shell_map_[l][g]);
    int gijkl = ij_offset(gij,gkl);

    if (gijkl > ijkl)
      return 0;
    else if (gijkl == ijkl)
      nijkl++;
  }

  return ng_/nijkl;
}

REF_dec(PetiteList);

#endif
    
