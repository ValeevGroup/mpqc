
#include <stdio.h>
#include <stream.h>
#include "moleculei.h"
#include "optimizei.h"
#include "functioni.h"
#include "keyvali.h"
#include "energyi.h"

#include <util/class/proxy.h>
#include <util/keyval/keyval.h>

///////////////////////////////////////////////////////////////////
// ChemistryLoader creates implementation objects.

class ChemistryLoader: public CORBA(LoaderClass) {
  public:
    ChemistryLoader();
    ~ChemistryLoader() {}
    CORBA(ObjectRef) load (const char *interface,
                           const char *marker,
                           unsigned char isBind,
                           CORBA(Environment)&); 

  private:
    enum {MaxMarkers = 100};
    const char *markers_[MaxMarkers];
    C_KeyValCreatableImpl *objects_[MaxMarkers];
    int n_;
  public:
    C_KeyValCreatableImpl *lookup(const char *marker);
    void store(const char *marker, C_KeyValCreatableImpl *);
};

C_KeyValCreatableImpl *
ChemistryLoader::lookup(const char *marker)
{
  for (int i=0; i<n_; i++) {
      if (!strcmp(markers_[i], marker)) return objects_[i];
    }
  return 0;
}

void
ChemistryLoader::store(const char *marker, C_KeyValCreatableImpl *g)
{
  if (n_ >= MaxMarkers) {
      cerr << "out of storagage space for marker->describedclass map" << endl;
    }
  else {
      markers_[n_] = strcpy(new char[strlen(marker)+1], marker);
      objects_[n_] = g;
      n_++;
    }
}

ChemistryLoader::ChemistryLoader():
  CORBA(LoaderClass)(1)
{
  n_ = 0;
}

CORBA(ObjectRef)
ChemistryLoader::load(const char *interface,
                 const char *marker,   
                 unsigned char isBind, 
                 CORBA(Environment)&)
{
  fflush(stdout); cout.flush();
  cout << "server loading object " << interface << endl;
  fflush(stdout); cout.flush();
  C_KeyValCreatableImpl *o;
  CORBA(ObjectRef) g;
  if (!strcmp(interface, "C_MolecularEnergy")) {
      C_MolecularEnergyImpl *s = new C_MolecularEnergyImpl();
      g = new TIE_C_MolecularEnergy(C_MolecularEnergyImpl) (s , marker, this);
      o = s;
    }
  else if (!strcmp(interface, "C_KeyValCreatable")) {
      C_KeyValCreatableImpl *s = new C_KeyValCreatableImpl();
      g = new TIE_C_KeyValCreatable(C_KeyValCreatableImpl) (s , marker, this);
      o = s;
    }
  else if (!strcmp(interface, "C_Molecule")) {
      C_MoleculeImpl *s = new C_MoleculeImpl();
      g = new TIE_C_Molecule(C_MoleculeImpl) (s , marker, this);
      o = s;
    }
  else if (!strcmp(interface, "C_Optimize")) {
      C_OptimizeImpl *s = new C_OptimizeImpl();
      g = new TIE_C_Optimize(C_OptimizeImpl) (s , marker, this);
      o = s;
    }
  else if (!strcmp(interface, "C_Function")) {
      C_FunctionImpl *s = new C_FunctionImpl();
      g = new TIE_C_Function(C_FunctionImpl) (s , marker, this);
      o = s;
    }
  else {
      return 0;
    }
  store(marker, o);
  return g;
}

///////////////////////////////////////////////////////////////////////
// ORBProxy objects can be given in keyval input to grab an object
// that is already in in the server by its marker

class ORBProxy: public DescribedClassProxy {
  private:
    char *marker_;
  public:
    static ChemistryLoader *loader;

    ORBProxy(const Ref<KeyVal>&);
    ~ORBProxy();
    Ref<DescribedClass> object();
};

ChemistryLoader *ORBProxy::loader = 0;

static ClassDesc ORBProxy_cd(
  typeid(ORBProxy),"ORBProxy",1,"public DescribedClassProxy",
  0, create<ORBProxy>, 0);

ORBProxy::ORBProxy(const Ref<KeyVal> &keyval)
{
  marker_ = keyval->pcharvalue("marker");
}

ORBProxy::~ORBProxy()
{
  delete[] marker_;
}

Ref<DescribedClass>
ORBProxy::object()
{
  if (!loader) return 0;

  C_KeyValCreatableImpl *kvc = loader->lookup(marker_);

  if (kvc) return kvc->object();
  else return 0;
}

///////////////////////////////////////////////////////////////////////
// main creates a loader and starts serving objects.

main() 
{      
  ChemistryLoader *mel = new ChemistryLoader();

  ORBProxy::loader = mel;

  CORBA(Orbix).impl_is_ready("Chemistry");

  cout << "server terminating" << endl;

  return 0;
}
