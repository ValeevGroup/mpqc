
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
    SymmOneBodyIntIter(const RefPetiteList&);
    ~SymmOneBodyIntIter();

    void start();
    void next();

    void start_ltri();
    void next_ltri();

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
    AOSO_Transformation(const RefGaussianBasisSet&, const RefPetiteList&);

    void process(SCMatrixBlockIter&);
    void process(SCMatrixRectBlock*);
};
#endif
