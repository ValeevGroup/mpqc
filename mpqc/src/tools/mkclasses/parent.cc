
#include "mkclasses.h"

Parent::Parent(const string &dec):
  access_(Private),
  virtual_(false)
{
  // if parents is empty then we are done
  if (dec.empty()) return;

  char *tokens = ::strcpy(new char[dec.length()+1], dec.c_str());
  const char* whitesp = "\t\n,() ";
  char* token;
  for (token = ::strtok(tokens,whitesp);
       token;
       token = ::strtok(0,whitesp)) {
      if (!strcmp(token,"virtual")) {
          virtual_ = true;
        }
      else if (!strcmp(token,"public")) {
          access_ = Parent::Public;
        }
      else if (!strcmp(token,"protected")) {
          access_ = Parent::Protected;
        }
      else if (!strcmp(token,"private")) {
          access_ = Parent::Private;
        }
      else {
          name_ = token;
          break;
        }
    }
  delete[] tokens;
}

Parent::~Parent()
{
}

string
Parent::stringrep() const
{
  string result;
  if (virtual_) {
      result = "virtual ";
    }
  switch (access_) {
  case Private:
      result += "private ";
      break;
  case Protected:
      result += "protected ";
      break;
  case Public:
      result += "public ";
      break;
    };
  result += name_;
  return result;
}
