
#ifndef _util_render_material_h
#define _util_render_material_h

#include <stdio.h>

#include <util/class/class.h>
#include <util/keyval/keyval.h>
#include <util/render/parameter.h>

class Color {
  private:
    double red_;
    double green_;
    double blue_;
  public:
    Color(): red_(0.0), green_(0.0), blue_(0.0) {}
    Color(double r, double g, double b): red_(r), green_(g), blue_(b) {}
    Color(const RefKeyVal&);
    double red() const { return red_; }
    double green() const { return green_; }
    double blue() const { return blue_; }
    void red(double r) { red_ = r; }
    void green(double g) { green_ = g; }
    void blue(double b) { blue_ = b; }
};
    
    

class Material: public DescribedClass {
#   define CLASSNAME Material
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  private:
    Parameter<Color> diffuse_;
    Parameter<Color> ambient_;
  public:
    Material();
    Material(const RefKeyVal&);
    ~Material();
    Parameter<Color>& diffuse() { return diffuse_; }
    Parameter<Color>& ambient() { return ambient_; }
    void print(FILE*fp = stdout);
};
DescribedClass_REF_dec(Material);

#endif
