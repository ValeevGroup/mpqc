
#ifndef _util_render_parameter_h
#define _util_render_parameter_h

#include <util/keyval/keyval.h>

template <class T>
class Parameter {
    T parameter_;
    int is_set_;
    int overrides_;
  public:
    Parameter(): is_set_(0), overrides_(0) {}
    void set(const T& a) { parameter_ = a; is_set_ = 1; }
    void override(int overrides = 1) { overrides_ = overrides; }
    const T& value() const { return parameter_; }
    int overrides() const { return overrides_; }
    int is_set() const { return is_set_; }
};

#endif
