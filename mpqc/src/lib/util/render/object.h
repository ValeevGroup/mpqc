
#ifndef _util_render_object_h
#define _util_render_object_h

#include <stdio.h>

#include <util/keyval/keyval.h>
#include <util/render/material.h>
#include <util/render/appearance.h>
#include <util/render/transform.h>

DescribedClass_REF_fwddec(Render);

class RenderedObject: public DescribedClass {
#   define CLASSNAME RenderedObject
#   include <util/class/classda.h>
  protected:
    char* name_;
    RefMaterial material_;
    RefAppearance appearance_;
    RefTransform transform_;

    friend class Render;
  public:
    RenderedObject(const RefMaterial& = 0);
    RenderedObject(const RefKeyVal&);
    ~RenderedObject();
    const char* name() const { return name_; }
    RefMaterial material() const { return material_; }
    RefAppearance appearance() const { return appearance_; }
    RefTransform transform() const { return transform_; }
    void material(const RefMaterial&m) { material_ = m; }
    void appearance(const RefAppearance&a) { appearance_ = a; }
    void transform(const RefTransform&t) { transform_ = t; }

    virtual void print(FILE* fp = stdout);

    // to be called only by derivatives of Render
    virtual void render(const RefRender&) = 0;
};
DescribedClass_REF_dec(RenderedObject);

class RenderedObjectSet: public RenderedObject {
#   define CLASSNAME RenderedObjectSet
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  private:
    int capacity_;
    int n_;
    RefRenderedObject *array_;
  protected:
    void render(const RefRender&);
  public:
    RenderedObjectSet(int capacity = 10);
    RenderedObjectSet(const RefKeyVal&);
    ~RenderedObjectSet();
    int n() const { return n_; }
    const RefRenderedObject& element(int i) const { return array_[i]; }
    virtual void add(const RefRenderedObject&);
};
DescribedClass_REF_dec(RenderedObjectSet);

#endif
