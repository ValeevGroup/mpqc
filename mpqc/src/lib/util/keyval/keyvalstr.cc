
extern "C" {
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
}
#include <iostream.h>
#include <util/keyval/keyval.h>
#include <util/keyval/keyvalImplMap.h>
#include <util/class/proxy.h>


////////////////////////////////////////////////////////////////////////
// StringKeyVal



StringKeyVal::StringKeyVal()
{
  _map = new MAPCTOR;

}

StringKeyVal::~StringKeyVal()
{
  delete _map;
}

const char*
StringKeyVal::classname(const char * key)
{
  return 0;
}

const char*
StringKeyVal::truekeyword(const char * key)
{
  return key;
}

// This does not cause objects to be constructed.
int
StringKeyVal::key_exists(const char* key)
{
  stringvalue(key);
  if (error() == OK || error() == HasNoValue) {
      return 1;
    }
  return 0;
}

RefKeyValValue
StringKeyVal::key_value(const char* key)
{
  // convert the key to the true key so variable assignments in the
  // input will effectively be done by reference
  // check to see if the datum is a described class
  const char* tkw = truekeyword(key);
  //cout << "truekeyword = "<< tkw << '\n';
  if (!tkw) return 0;

  // if a classname exists then read in the datum as an object
  const char* classn = classname(tkw);
  //cout << "classname = " << classn << '\n';
  if (classn) {
      KeyValKeyword truekey(tkw);
      //cout << "truekey = " << tkw << '\n';
      // see if a reference to this datum already exists
      //cout << "(*_map).contains(truekey) = "
      //     << (*_map).contains(truekey)
      //     << '\n';
      if ((*_map).contains(truekey)) {
          return (*_map)[truekey];
        }
      else {
          // create a new instance of this datum
          RefKeyVal pkv = new PrefixKeyVal(this, tkw);
          const ClassDesc* cd = ClassDesc::name_to_class_desc(classn);
          if (!cd) {
              ClassDesc::load_class(classn);
              cd = ClassDesc::name_to_class_desc(classn);
            }
          // the original error status must be saved
          KeyValError original_error = error();
          RefDescribedClass newdc(cd->create(pkv));
          if (newdc.null()) {
              fprintf(stderr,"StringKeyVal::value: create failed for:\n");
              fprintf(stderr," keyword = \"%s\" class = \"%s\"\n",
                      tkw,classn);
              fprintf(stderr," either the KeyVal create operator doesn't\n");
              fprintf(stderr," exist or memory was exhausted\n");
            }
          DescribedClassProxy *proxy
              = DescribedClassProxy::castdown(newdc.pointer());
          if (proxy) {
              newdc = proxy->object();
            }
          seterror(original_error);
          //pkv.dump(stderr);
          KeyValValueRefDescribedClass* keyvalvalue =
            new KeyValValueRefDescribedClass(newdc);
          (*_map)[truekey] = keyvalvalue;
          return keyvalvalue;
        }
    }
  else {
      const char* string = stringvalue(tkw);
      if (string) return new KeyValValueString(string);
      return 0;
    }
}

void
StringKeyVal::errortrace(ostream&fp,int n)
{
  offset(fp,n); fp << "StringKeyVal: error: \"" << errormsg() << "\"" << endl;
}

void
StringKeyVal::dump(ostream&fp,int n)
{
  offset(fp,n); fp << "StringKeyVal: error: \"" << errormsg() << "\"" << endl;
}
