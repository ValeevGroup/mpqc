
#ifndef _util_render_color_h
#define _util_render_color_h

#include <util/keyval/keyval.h>

class Color {
  private:
    double red_;
    double green_;
    double blue_;
  public:
    Color() {}
    Color(double r, double g, double b): red_(r), green_(g), blue_(b) {}
    Color(const RefKeyVal&);
    double red() const { return red_; }
    double green() const { return green_; }
    double blue() const { return blue_; }
    void set_rgb(double r, double g, double b) {
        red_ = r;
        green_ = g;
        blue_ = b;
      }
};

#endif
