
#include <util/keyval/keyval.h>

KeyValKeyword::KeyValKeyword() :
  keyword_(0)
{
}

KeyValKeyword::KeyValKeyword(const char* name)
{
  // get rid of the leading ':' character
  if (name[0] == ':') {
      keyword_ = ::strcpy(new char[strlen(name)],&name[1]);
    }
  else {
      keyword_ = ::strcpy(new char[strlen(name)+1],name);
    }
}

KeyValKeyword::KeyValKeyword(KeyValKeyword& key):
  keyword_(::strcpy(new char[strlen(key.keyword_)+1],key.keyword_))
{
}

KeyValKeyword::~KeyValKeyword()
{
  if (keyword_) delete[] keyword_;
}

KeyValKeyword& KeyValKeyword::operator=(const KeyValKeyword& key)
{
  if (keyword_ && keyword_ != key.keyword_) {
      delete[] keyword_;
      keyword_ = ::strcpy(new char[strlen(key.keyword_)+1],key.keyword_);
    }
  return *this;
}

int KeyValKeyword::operator==(KeyValKeyword& ck)
{
  if (!keyword_) {
      if (!ck.keyword_) return 1;
      else return 0;
    }
  else if (!ck.keyword_) return 0;
  
  return !::strcmp(keyword_,ck.keyword_);
}

int
KeyValKeyword::hash() const
{
  int r=0;
  int i;

  // Even numbered bytes make up the lower part of the hash index
  for (i=0; i<::strlen(keyword_); i+=2) {
      r ^= keyword_[i];
    }

  // Odd numbered bytes make up the upper part of the hash index
  for (i=1; i<::strlen(keyword_); i+=2) {
      r ^= keyword_[i]<<8;
    }

  return r;
}
