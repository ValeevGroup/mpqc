
#include <stdio.h>
#include "ipv2.h"
#include "keyval.h"

ParsedKeyVal::ParsedKeyVal()
{
  ipv2 = new IPV2;
}

ParsedKeyVal::ParsedKeyVal(const char* name):
nfp(0),
nfile(0),
file(0)
{
  ipv2 = new IPV2;
  read(name);
}

ParsedKeyVal::ParsedKeyVal(FILE* fp):
nfp(0),
nfile(0),
file(0)
{
  ipv2 = new IPV2;
  read(fp);
}

void
ParsedKeyVal::read(const char* name)
{
  FILE *infp;

  infp = fopen(name,"r");
  if (!infp) {
    fprintf(stderr,"ParsedKeyVal couldn't open %s\n",name);
    exit(1);
    }

  int i;
  char**newfile = new char*[nfile+1];
  for (i=0; i<nfile; i++) newfile[i] = file[i];
  if (file) delete[] file;
  file = newfile;
  newfile[nfile] = strdup(name);
  nfile++;

  read(infp);
  nfp--; // read(infp) will incr nfp, but this isn't desired so undo
  fclose(infp);
}

void ParsedKeyVal::read(FILE*infp)
{
  nfp++;
  ipv2->read(infp,stderr);
}

ParsedKeyVal::~ParsedKeyVal()
{
  delete ipv2;
  for (int i=0; i<nfile; i++) free(file[i]);
  delete[] file;
}

static KeyVal::KeyValError maperr(IPV2::Status err)
  {
  if (err == IPV2::OK            ) return KeyVal::OK;
  if (err == IPV2::KeyNotFound ) return KeyVal::UnknownKeyword;
  if (err == IPV2::OutOfBounds ) return KeyVal::UnknownKeyword;
  if (err == IPV2::Malloc        ) return KeyVal::OperationFailed;
  if (err == IPV2::NotAnArray  ) return KeyVal::UnknownKeyword;
  if (err == IPV2::NotAScalar  ) return KeyVal::HasNoValue;
  if (err == IPV2::Type          ) return KeyVal::WrongType;
  if (err == IPV2::HasNoValue  ) return KeyVal::HasNoValue;
  if (err == IPV2::ValNotExpd  ) return KeyVal::OperationFailed;
  return KeyVal::OperationFailed;
  }

const char* ParsedKeyVal::stringvalue(const char* key)
{
  char* result;
  seterror(maperr(ipv2->value_v((char *)key,&result,0,0)));
  if (error() != OK) {
      result = 0;
    }
  return result;
}

const char*
ParsedKeyVal::classname(const char* key)
{
  char* result;
  seterror(maperr(ipv2->classname_v((char *)key,&result,0,0)));
  return result;
}

const char*
ParsedKeyVal::truekeyword(const char*key)
{
  char* result;
  seterror(maperr(ipv2->truekeyword_v((char *)key,&result,0,0)));
  if (!result && error() == OK) return key;
  else return result;
}

void ParsedKeyVal::errortrace(FILE*fp,int n)
{
  offset(fp,n); fprintf(fp,"ParsedKeyVal: error: \"%s\"\n",errormsg());
  if (nfp) {
      offset(fp,n);
      fprintf(fp,"    reading from %d files with unknown names\n", nfp);
    }
  for (int i=0; i<nfile; i++) {
      offset(fp,n); fprintf(fp,"    reading from \"%s\"\n",file[i]);
    }
}

void ParsedKeyVal::dump(FILE*fp,int n)
{
  offset(fp,n); fprintf(fp,"ParsedKeyVal: error: \"%s\"\n",errormsg());
  if (nfp) {
      offset(fp,n);
      fprintf(fp,"    reading from %d files with unknown names\n", nfp);
    }
  for (int i=0; i<nfile; i++) {
      offset(fp,n); fprintf(fp,"    reading from \"%s\"\n",file[i]);
    }
  offset(fp,n); fprintf(fp,"The IPV2 tree:\n");
  ipv2->print_tree(fp);

}
