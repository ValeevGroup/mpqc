
#ifndef _chemistry_molecule_molinfo_h
#define _chemistry_molecule_molinfo_h

#include <util/class/class.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/chemelem.h>

class MolInfo: public DescribedClass {
#   define CLASSNAME MolInfo
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  protected:
    RefKeyVal keyval;
  public:
    MolInfo(const RefKeyVal&);
    virtual ~MolInfo();
    double doublevalue(const char * sym, const char * property);
    double doublevalue(const char * sym, const char * property, int i);
    KeyVal::KeyValError error();
};
DescribedClass_REF_dec(MolInfo);

#define ATOMINFO_MAXZ 100
class AtomInfo: public MolInfo {
#   define CLASSNAME AtomInfo
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  private:
    // certain values are cached here for fast access
    double radius_vals[ATOMINFO_MAXZ];
    int have_radius[ATOMINFO_MAXZ];
    double rgb_vals[ATOMINFO_MAXZ][3];
    int have_rgb[ATOMINFO_MAXZ];
    int get_zindex(const ChemicalElement&);

    // these items are cached for quick lookup
    double radius_scale_factor_;
    double radius_to_bohr_;

    void preload_values();
  public:
    AtomInfo(const RefKeyVal&);
    ~AtomInfo();
    double doublevalue(const ChemicalElement& atom,
                       const char * property);
    double doublevalue(const ChemicalElement& atom,
                       const char * property, int i);

    // routines for common properties
    double radius(const ChemicalElement&);
    double rgb(const ChemicalElement&, int color);
    double red(const ChemicalElement&);
    double green(const ChemicalElement&);
    double blue(const ChemicalElement&);
};
DescribedClass_REF_dec(AtomInfo);

#endif

