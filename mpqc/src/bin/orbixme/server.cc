
#include <stream.h>
#include "molei.h"

class MolELoader: public CORBA(LoaderClass) {
  public:
    MolELoader();
    ~MolELoader() {}
    CORBA(ObjectRef) load (const char *interface,
                           const char *marker,
                           unsigned char isBind,
                           CORBA(Environment)&); 
};

MolELoader::MolELoader():
  CORBA(LoaderClass)(1)
{
}

CORBA(ObjectRef)
MolELoader::load(const char *interface,
                 const char *marker,   
                 unsigned char isBind, 
                 CORBA(Environment)&)
{
  if (!strcmp(interface, "MolE")) {
      MolEImpl *g = new MolEImpl(marker,this);
      return g;
    }
  else return 0;
}


main() 
{      
  MolELoader *mel = new MolELoader();

  CORBA(Orbix).impl_is_ready("MolE");

  cout << "server terminating" << endl;

  return 0;
}
