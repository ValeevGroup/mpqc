
#include "mkclasses.h"

void
Class::write_ctor(Out &p)
{
  p("%sClassDesc::%sClassDesc():", name().c_str(), name().c_str());
  p++;
  for (vector<DataMember>::const_iterator mi = member_data().begin();
       mi != member_data().end();
       mi++) {
      const DataMember &m = *mi;
      p("%s(&%s::%s),", m.name().c_str(), name().c_str(), m.name().c_str());
    }
  p("ClassDesc(");
  p++;
  p("\"%s\",", name().c_str());
  p("%d,", version());
  p("\"%s\",", parentstring().c_str());
  p("0,0,0)");
  p--;
  p--;
  p("{");
  p("}");
}

