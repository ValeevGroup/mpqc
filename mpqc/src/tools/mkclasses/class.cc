
#include "mkclasses.h"

Class::Class()
{
}

Class::~Class()
{
}

void
Class::write_declaration(Out &p)
{
  p("class %sClassDesc: public ClassDesc {", name().c_str());

  p++;

  p("private:");
  p++;

  for (vector<DataMember>::const_iterator mi = member_data().begin();
       mi != member_data().end();
       mi++) {
      write_member_declaration(p, *mi);
    }

  p--;
  p("public:");
  p++;

  p("%sClassDesc();", name().c_str());
  if (ctor_void()) {
      p("DescribedClass* create();");
    }
  if (ctor_keyval()) {
      p("DescribedClass* create(const RefKeyVal& keyval);");
    }
  if (ctor_statein()) {
      p("DescribedClass* create(StateIn& statein);");
    }

  p--;
  p--;

  p("};");
}

void
Class::write_member_declaration(Out& p, const DataMember &m)
{
  if (m.type().empty() || name().empty() || m.name().empty()) return;
  p("DescribedMemberDatum<%s,%s> %s;",
    m.type().c_str(), name().c_str(), m.name().c_str());
}

string
Class::parentstring() const
{
  string result;
  for (vector<Parent>::const_iterator p = parents().begin();
       p != parents().end();
       p++) {
      if (!result.empty()) result += ",";
      result += (*p).stringrep();
    }
  return result;
}

void
Class::set_parents(const vector<string> &p)
{
  for (vector<string>::const_iterator i = p.begin();
       i != p.end();
       i++) {
      Parent p(*i);
      parents_.insert(parents_.end(), p);
    }
}

void
Class::start_member_datum()
{
  DataMember c;
  current_member_datum_ = member_data_.insert(member_data_.end(), c);
}

void
Class::end_member_datum()
{
}
