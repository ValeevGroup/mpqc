
#include <stream.h>
#include "moleculei.h"
#include "keyvali.h"
#include "energyi.h"

class ChemistryLoader: public CORBA(LoaderClass) {
  public:
    ChemistryLoader();
    ~ChemistryLoader() {}
    CORBA(ObjectRef) load (const char *interface,
                           const char *marker,
                           unsigned char isBind,
                           CORBA(Environment)&); 
};

ChemistryLoader::ChemistryLoader():
  CORBA(LoaderClass)(1)
{
}

CORBA(ObjectRef)
ChemistryLoader::load(const char *interface,
                 const char *marker,   
                 unsigned char isBind, 
                 CORBA(Environment)&)
{
  if (!strcmp(interface, "C_MolecularEnergy")) {
      TIE_C_MolecularEnergy(C_MolecularEnergyImpl) *g
          = new TIE_C_MolecularEnergy(C_MolecularEnergyImpl)
                (new C_MolecularEnergyImpl(), marker, this);
      return g;
    }
  else if (!strcmp(interface, "C_KeyValCreatable")) {
      TIE_C_KeyValCreatable(C_KeyValCreatableImpl) *g
          = new TIE_C_KeyValCreatable(C_KeyValCreatableImpl)
                (new C_KeyValCreatableImpl(), marker, this);
      return g;
    }
  else if (!strcmp(interface, "C_Molecule")) {
      TIE_C_Molecule(C_MoleculeImpl) *g
          = new TIE_C_Molecule(C_MoleculeImpl)
                (new C_MoleculeImpl(), marker, this);
      return g;
    }
  else return 0;
}


main() 
{      
  ChemistryLoader *mel = new ChemistryLoader();

  CORBA(Orbix).impl_is_ready("Chemistry");

  cout << "server terminating" << endl;

  return 0;
}
