
#include "mkclasses.h"

void
Class::write_members(Out &p)
{
  // required function
  p("static %sClassDesc %s_class_desc_;", name().c_str(), name().c_str());
  p("const ClassDesc* %s::class_desc() const", name().c_str());
  p("{");
  p("  return &%s_class_desc_;", name().c_str());
  p("}");

  // required function
  p("void *");
  p("%s::_castdown(const ClassDesc*cd)", name().c_str());
  p("{");
  p++;
  p("if (!cd) return 0;");
  p("if (cd == &%s_class_desc_) {", name().c_str());
  p("    return this;");
  p("  }");
  p("const ParentClasses& parents = %s_class_desc_.parents();",
    name().c_str());
  p("void *tmp,*p=0;");
  for (vector<Parent>::const_iterator i=parents().begin();
       i != parents().end();
       i++) {
      const Parent &parent = *i;
      if (parent.access() == Parent::Private) continue;
      p("tmp = %s::_castdown(cd);", parent.name().c_str());
      p("if (tmp) {");
      p++;
      if (i != parents().begin()) {
          p("if (p && tmp != p) {");
          p++;
          p("fprintf(stderr,\"%s: castdown to %%s ambiguous\\n\",",
            name().c_str());
          p("        cd->name());");
          p("fprintf(stderr,\" tmp = 0x%%lx p = 0x%%lx\\n\",");
          p("        (long)tmp,(long)p);");
          p("}");
          p--;
        }
      p("p = tmp;");
      p("}");
      p--;
    }
  p("return p;");
  p--;
  p("}");

  if (static_class_desc()) {
      p("const ClassDesc* %s::static_class_desc()", name().c_str());
      p("{");
      p("  return &%s_class_desc_;", name().c_str());
      p("}");
    }

  if (castdown()) {
      p("%s*", name().c_str());
      p("%s::castdown(const RefDescribedClass&p)", name().c_str());
      p("{");
      p("  if (p.null()) return 0;");
      p("  return (%s*) p->_castdown(%s::static_class_desc());",
        name().c_str(), name().c_str());
      p("}");
    }

  if (require_castdown()) {
      p("%s*", name().c_str());
      p("%s::require_castdown(const RefDescribedClass&p,",
        name().c_str());
      p("                     const char * errmsg,...)");
      p("{");
      p++;
      p("if (p.null()) return 0;");
      p("%s* t = (%s*) p->_castdown(%s::static_class_desc());",
        name().c_str(), name().c_str(), name().c_str());
      p("if (!t) {");
      p("  va_list args;");
      p("  va_start(args,errmsg);");
      p("  fprintf(stderr,\"A required castdown failed in: \");");
      p("  vfprintf(stderr,errmsg,args);");
      p("  fprintf(stderr,\"\\n\");");
      p("  fprintf(stderr,");
      p("          \"wanted type \\\"%%s\\\" but got \\\"%%s\\\"\\n\",");
      p("          \"%s\",p.nonnull()?p->class_name():\"(null)\");",
        name().c_str());
      p("  va_end(args);");
      p("  abort();");
      p("  }");
      p("return t;");
      p--;
      p("}");
    }
  
  if (void_ctor()) {
      p("DescribedClass*");
      p("%sClassDesc::create()", name().c_str());
      p("{");
      p("  return (DescribedClass*) new %s();", name().c_str());
      p("}");
    }

  if (keyval_ctor()) {
      p("DescribedClass*");
      p("%sClassDesc::create(const RefKeyVal& keyval)", name().c_str());
      p("{");
      p("  return (DescribedClass*) new %s(keyval);", name().c_str());
      p("}");
    }

  if (statein_ctor()) {
      p("DescribedClass*");
      p("%sClassDesc::create(StateIn& statein)", name().c_str());
      p("{");
      p("  return (DescribedClass*) new %s(statein);", name().c_str());
      p("}");
    }
}
