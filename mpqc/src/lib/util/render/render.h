
#ifndef _util_render_render_h
#define _util_render_render_h

#include <util/class/class.h>
#include <util/container/ref.h>
#include <util/render/appearance.h>
#include <util/render/material.h>
#include <util/render/transform.h>
#include <util/render/stack.h>

DescribedClass_REF_fwddec(RenderedObject);
DescribedClass_REF_fwddec(RenderedObjectSet);
DescribedClass_REF_fwddec(RenderedSphere);
DescribedClass_REF_fwddec(RenderedPolygons);

class Render: public DescribedClass {
#   define CLASSNAME Render
#   include <util/class/classda.h>
  protected:
    RefMaterial default_material_;
    RefAppearance default_appearance_;
    RefTransform default_transform_;

    Stack<RefMaterial> material_stack_;
    Stack<RefAppearance> appearance_stack_;
    Stack<RefTransform> transform_stack_;

    virtual void push_material(const RefMaterial& m);
    virtual void push_appearance(const RefAppearance& a);
    virtual void push_transform(const RefTransform& t);
    virtual RefMaterial pop_material();
    virtual RefAppearance pop_appearance();
    virtual RefTransform pop_transform();

  public:
    Render();
    Render(const RefKeyVal&);
    virtual ~Render();

    RefMaterial default_material() { return default_material_; }
    RefAppearance default_appearance() { return default_appearance_; }
    RefTransform default_transform() { return default_transform_; }
    void default_material(const RefMaterial& m) { default_material_ = m; }
    void default_appearance(const RefAppearance& a) {default_appearance_ = a;}
    void default_transform(const RefTransform& t) {default_transform_ = t;}

    virtual void clear() = 0;

    virtual void render(const RefRenderedObject&);

    virtual void set(const RefRenderedObjectSet&);
    virtual void sphere(const RefRenderedSphere&);
    virtual void polygons(const RefRenderedPolygons&) = 0;
};
DescribedClass_REF_dec(Render);

#endif
