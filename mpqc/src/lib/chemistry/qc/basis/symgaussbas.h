
#ifndef _chemistry_qc_basis_symgaussbas_h
#define _chemistry_qc_basis_symgaussbas_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/petite.h>

class SymmGaussianBasisSet: public GaussianBasisSet {
#   define CLASSNAME SymmGaussianBasisSet
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    PetiteList pl;
    
  public:
    SymmGaussianBasisSet(const GaussianBasisSet&);
    SymmGaussianBasisSet(const RefKeyVal&);
    SymmGaussianBasisSet(StateIn&);
    ~SymmGaussianBasisSet();
    void save_data_state(StateOut&);

    PetiteList& petite_list();
};

SavableState_REF_dec(SymmGaussianBasisSet);

////////////////////////////////////////////////////////////////////////////

class SymmetryOrbitals {
    friend class AOSO_Transformation;
    friend PetiteList;
    
  private:
    SO_block *sos_;
    SO_block **som_;
    RefSymmGaussianBasisSet gbs_;

  public:
    SymmetryOrbitals();
    SymmetryOrbitals(const RefSymmGaussianBasisSet&);
    ~SymmetryOrbitals();

    void print(FILE* =stdout);

    RefBlockedSCDimension AO_basisdim();
    RefBlockedSCDimension SO_basisdim();
};
    
////////////////////////////////////////////////////////////////////////////

class AOSO_Unit : public BlockedSCElementOp {
  private:
    RefBlockedSCDimension d1;
    RefBlockedSCDimension d2;

  public:
    AOSO_Unit(const RefBlockedSCDimension&,const RefBlockedSCDimension&);

    void process(SCMatrixBlockIter&);
    void process(SCMatrixRectBlock*);
};

class AOSO_Transformation : public BlockedSCElementOp {
  private:
    SymmetryOrbitals sos;
    CharacterTable ct;

  public:
    AOSO_Transformation(const RefSymmGaussianBasisSet&);

    void process(SCMatrixBlockIter&);
    void process(SCMatrixRectBlock*);
};
#endif
