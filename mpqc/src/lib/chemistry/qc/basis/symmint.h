
#ifndef _chemistry_qc_integral_symmint_h
#define _chemistry_qc_integral_symmint_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/petite.h>

////////////////////////////////////////////////////////////////////////////

class SymmOneBodyIntIter : public OneBodyIntIter {
  protected:
    RefPetiteList pl;
    
  public:
    SymmOneBodyIntIter(const RefOneBodyInt&, const RefPetiteList&);
    ~SymmOneBodyIntIter();

    void start(int ist=0, int jst=0, int ien=0, int jen=0);
    void next();

    double scale() const;
};

////////////////////////////////////////////////////////////////////////////

class SymmetryOrbitals {
    friend class AOSO_Transformation;
    friend PetiteList;
    
  private:
    SO_block *sos_;
    SO_block **som_;
    RefGaussianBasisSet gbs_;
    RefPetiteList pl_;

  public:
    SymmetryOrbitals();
    SymmetryOrbitals(const RefGaussianBasisSet&, const RefPetiteList&);
    ~SymmetryOrbitals();

    void print(FILE* =stdout);

    RefSCDimension AO_basisdim();
    RefSCDimension SO_basisdim();
};
    
////////////////////////////////////////////////////////////////////////////

class AOSO_Unit : public BlockedSCElementOp {
  private:
    RefSCDimension d1;
    RefSCDimension d2;

  public:
    AOSO_Unit(const RefSCDimension&,const RefSCDimension&);

    void process(SCMatrixBlockIter&);
    void process_spec(SCMatrixRectBlock*);
};

class AOSO_Transformation : public BlockedSCElementOp {
  private:
    SymmetryOrbitals sos;
    CharacterTable ct;

  public:
    AOSO_Transformation(const RefGaussianBasisSet&, const RefPetiteList&);

    void process(SCMatrixBlockIter&);
    void process_spec(SCMatrixRectBlock*);
};
#endif
