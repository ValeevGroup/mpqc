
#ifndef _libQC_atomcent_h
#define _libQC_atomcent_h

#include <stdio.h>

#include <chemistry/molecule/chemelem.h>
#include <math/topology/point.h>

class AtomicCenter :
  virtual public DescribedClass, virtual public SavableState
{
  DescribedClass_DECLARE(AtomicCenter)
  SavableState_DECLARE(AtomicCenter)
  private:
    Point p;
    ChemicalElement element_;
  public:
    AtomicCenter();
    AtomicCenter(const AtomicCenter&);
    AtomicCenter(const char*symbol,double x,double y,double z);
    ~AtomicCenter();

    AtomicCenter& operator=(const AtomicCenter&ac);

    inline double& operator[](int i) { return p[i]; };
    inline operator int() { return element_.number(); };
    inline ChemicalElement& element() { return element_; };
    inline Point& point() { return p; };

    void save_data_state(StateOut&);
    void restore_data_state(int,StateIn&);

    void print(FILE*fp=stdout);
};

#endif
