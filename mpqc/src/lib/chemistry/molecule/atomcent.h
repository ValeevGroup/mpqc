
#ifndef _chemistry_molecule_atomcent_h
#define _chemistry_molecule_atomcent_h

#include <stdio.h>

#include <chemistry/molecule/chemelem.h>
#include <math/topology/point.h>

class AtomicCenter :
  virtual public SavableState
{
#   define CLASSNAME AtomicCenter
#   define HAVE_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    Point p;
    ChemicalElement element_;
    char *label_;
  public:
    AtomicCenter();
    AtomicCenter(const AtomicCenter&);
    AtomicCenter(const char*symbol,double x,double y,double z, const char* =0);
    AtomicCenter(StateIn&);
    ~AtomicCenter();

    AtomicCenter& operator=(const AtomicCenter&ac);

    inline double& operator[](int i) { return p[i]; };
    inline operator int() { return element_.number(); };
    inline ChemicalElement& element() { return element_; };
    inline const ChemicalElement& element() const { return element_; };
    inline Point& point() { return p; };
    inline const Point& point() const { return p; };
    inline const char * label() const { return label_; }

    void save_data_state(StateOut&);

    void print(FILE*fp=stdout);
};
DescribedClass_REF_dec(AtomicCenter);

#endif
