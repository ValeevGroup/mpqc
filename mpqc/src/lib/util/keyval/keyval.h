
#ifndef _KeyVal_h
#define _KeyVal_h
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include <util/class/class.h>

class KeyValKeyword {
  private:
    char* keyword_;
  public:
    KeyValKeyword();
    KeyValKeyword(const char* name);
    KeyValKeyword(KeyValKeyword&);
    ~KeyValKeyword();
    KeyValKeyword& operator=(const KeyValKeyword&);
    int operator==(KeyValKeyword& ck);
    int hash() const;
    inline int cmp(KeyValKeyword&ck) const
    {
      return strcmp(keyword_,ck.keyword_);
    }
    inline const char* name() const {return keyword_;}
  };

class RefKeyValValue;

class KeyVal {
  public:
    enum {MaxKeywordLength = 256};
    enum KeyValError { OK, HasNoValue, WrongType,
                       UnknownKeyword, OperationFailed };
  private:
    KeyValError errcod;
  protected:
    KeyVal();
    inline void seterror(KeyValError err) { errcod = err; }

    void offset(FILE* fp,int n); // Put n ' ' into fp.
    enum {OffsetDelta=4};
  public:
    virtual ~KeyVal();

    // These access data given a keyword.
    virtual int    key_exists(const char*) = 0;
    virtual int    key_count(const char* =0);
    virtual RefKeyValValue key_value(const char*) = 0;
    virtual int    key_booleanvalue(const char*);
    virtual double key_doublevalue(const char* key);
    virtual float  key_floatvalue(const char* key);
    virtual char   key_charvalue(const char* key);
    virtual int    key_intvalue(const char* key);
    virtual char*  key_pcharvalue(const char* key);
    virtual RefDescribedClass key_describedclassvalue(const char* key);

    // For nonindexed things.   If a subclass defines one of these,
    // then the overloaded functions will be hidden.  The key_... functionss
    // should be overridden instead.
    int    exists(const char*);
    int    count(const char* =0);
    RefKeyValValue value(const char*);
    int    booleanvalue(const char*);
    double doublevalue(const char* key);
    float  floatvalue(const char* key);
    char   charvalue(const char* key);
    int    intvalue(const char* key);
    char*  pcharvalue(const char* key);
    RefDescribedClass describedclassvalue(const char* key);

    // For vectors:
    int    exists(const char*,int);
    int    count(const char*,int);
    int    booleanvalue(const char*,int);
    double doublevalue(const char* key,int);
    float  floatvalue(const char* key,int);
    char   charvalue(const char* key,int);
    int    intvalue(const char* key,int);
    char*  pcharvalue(const char* key,int);
    RefDescribedClass describedclassvalue(const char* key,int);

    inline int    exists(int i) { return exists((const char*)0,i); };
    inline int    count(int i) { return count((const char*)0,i); };
    inline int    booleanvalue(int i) {
        return booleanvalue((const char*)0,i);
      };
    inline double doublevalue(int i) { return doublevalue((const char*)0,i); };
    inline float  floatvalue(int i) { return floatvalue((const char*)0,i); };
    inline char   charvalue(int i) { return charvalue((const char*)0,i); };
    inline int    intvalue(int i) { return intvalue((const char*)0,i); };
    inline char*  pcharvalue(int i) { return pcharvalue((const char*)0,i); };
    inline RefDescribedClass describedclassvalue(int i) {
        return describedclassvalue((const char*)0,i);
      };

    // For arrays:
    int    exists(const char*,int,int);
    int    count(const char*,int,int);
    int    booleanvalue(const char*,int,int);
    double doublevalue(const char* key,int,int);
    float  floatvalue(const char* key,int,int);
    char   charvalue(const char* key,int,int);
    int    intvalue(const char* key,int,int);
    char*  pcharvalue(const char* key,int,int);
    RefDescribedClass describedclassvalue(const char* key,int,int);

    inline int    exists(int i,int j) { return exists((const char*)0,i,j); };
    inline int    count(int i,int j) { return count((const char*)0,i,j); };
    inline int    booleanvalue(int i,int j) {
        return booleanvalue((const char*)0,i,j);
      };
    inline double doublevalue(int i,int j) {
        return doublevalue((const char*)0,i,j);
      };
    inline float  floatvalue(int i,int j) {
        return floatvalue((const char*)0,i,j);
      };
    inline char   charvalue(int i,int j) {
        return charvalue((const char*)0,i,j);
      };
    inline int    intvalue(int i,int j) {
        return intvalue((const char*)0,i,j);
      };
    inline char*  pcharvalue(int i,int j) {
        return pcharvalue((const char*)0,i,j);
      };
    inline RefDescribedClass describedclassvalue(int i,int j) {
        return describedclassvalue((const char*)0,i,j);
      };

    // For all else:
    int    Va_exists(const char*,int,...);
    int    Va_count(const char*,int,...);
    int    Va_booleanvalue(const char*,int,...);
    double Va_doublevalue(const char* key,int,...);
    float  Va_floatvalue(const char* key,int,...);
    char   Va_charvalue(const char* key,int,...);
    int    Va_intvalue(const char* key,int,...);
    char*  Va_pcharvalue(const char* key,int,...);
    RefDescribedClass Va_describedclassvalue(const char* key,int,...);

    // default values
    inline static double Defaultdouble() { return 0.0; };
    inline static int    Defaultint() { return 0; };
    inline static float  Defaultfloat() { return 0.0; };
    inline static char   Defaultchar() { return 0; };
    inline static char*  Defaultpchar() { return 0; };
    inline static int    Defaultboolean() { return 0; };
    inline static RefDescribedClass DefaultRefDescribedClass() { return 0; };

    inline KeyValError error() { return errcod; }
    char*  errormsg(KeyValError);
    inline char*  errormsg() { return errormsg(errcod); }
    virtual void errortrace(FILE*fp=stderr,int offset = 0);
    virtual void dump(FILE*fp=stderr,int offset = 0);
};

class KeyValValue: public VRefCount {
  public:
    virtual ~KeyValValue();
    virtual KeyVal::KeyValError value(double&);
    virtual KeyVal::KeyValError booleanvalue(int&);
    virtual KeyVal::KeyValError value(float&);
    virtual KeyVal::KeyValError value(char&);
    virtual KeyVal::KeyValError value(int&);
    virtual KeyVal::KeyValError value(char*&);
    virtual KeyVal::KeyValError value(RefDescribedClass&);
};

REF_dec(KeyValValue);

class KeyValValuedouble: public KeyValValue {
  private:
    double _val;
  public:
    KeyValValuedouble(double);
    KeyVal::KeyValError value(double&);
};

class KeyValValueboolean: public KeyValValue {
  private:
    int _val;
  public:
    KeyValValueboolean(int);
    KeyVal::KeyValError booleanvalue(int&);
};

class KeyValValuefloat: public KeyValValue {
  private:
    float _val;
  public:
    KeyValValuefloat(float);
    KeyVal::KeyValError value(float&);
};

class KeyValValuechar: public KeyValValue {
  private:
    char _val;
  public:
    KeyValValuechar(char);
    KeyVal::KeyValError value(char&);
};

class KeyValValueint: public KeyValValue {
  private:
    int _val;
  public:
    KeyValValueint(int);
    KeyVal::KeyValError value(int&);
};

class KeyValValuepchar: public KeyValValue {
  private:
    char* _val;
  public:
    KeyValValuepchar(const char*);
    ~KeyValValuepchar();
    KeyVal::KeyValError value(char*&);
};

class KeyValValueRefDescribedClass: public KeyValValue {
  private:
    RefDescribedClass _val;
  public:
    KeyValValueRefDescribedClass(RefDescribedClass&);
    KeyVal::KeyValError value(RefDescribedClass&);
};

class KeyValValueString: public KeyValValue {
  private:
    const char* _val;
  public:
    KeyValValueString(const char*);
    KeyVal::KeyValError value(double&);
    KeyVal::KeyValError booleanvalue(int&);
    KeyVal::KeyValError value(float&);
    KeyVal::KeyValError value(char&);
    KeyVal::KeyValError value(int&);
    KeyVal::KeyValError value(char*&);
};

// this class allows keyval associations to be set up by the program,
// rather than determined by an external file
class KeyValKeywordRefKeyValValueMap;
class AssignedKeyVal: public KeyVal {
  private:
    KeyValKeywordRefKeyValValueMap* _map;
  public:
    AssignedKeyVal();
    ~AssignedKeyVal();
    int    key_exists(const char*);
    RefKeyValValue key_value(const char*);

    void assign(const char*, RefKeyValValue&);
    void assign(const char*, double);
    void assignboolean(const char*, int);
    void assign(const char*, float);
    void assign(const char*, char);
    void assign(const char*, int);
    void assign(const char*, const char*);
    void assign(const char*, RefDescribedClass&);
};

class StringKeyVal: public KeyVal {
  private:
    // once a described class is found it is kept here so
    // multiple references to it return the same instance
    KeyValKeywordRefKeyValValueMap* _map;
protected:
    StringKeyVal();
  public:
    virtual ~StringKeyVal();
    RefKeyValValue key_value(const char*);
    virtual const char* stringvalue(const char *) = 0;
    // returns the name of the exact class the object at the keyword
    virtual const char* classname(const char*);
    // returns a string which is the actual keyword if some sort
    // of variable substitution takes place (needed to make multiple
    // references to the same object work in input files)
    virtual const char* truekeyword(const char*);
    int    key_exists(const char*);

    virtual void errortrace(FILE*fp=stderr,int offset = 0);
    virtual void dump(FILE*fp=stderr,int offset = 0);
};

class AggregateKeyVal : public KeyVal {
  private:
    enum { MaxKeyVal = 4 };
    KeyVal* kv[MaxKeyVal];
    KeyVal* getkeyval(const char*key);
  public:
    AggregateKeyVal(KeyVal&);
    AggregateKeyVal(KeyVal&,KeyVal&);
    AggregateKeyVal(KeyVal&,KeyVal&,KeyVal&);
    AggregateKeyVal(KeyVal&,KeyVal&,KeyVal&,KeyVal&);
    ~AggregateKeyVal();
    int    key_exists(const char*);
    RefKeyValValue key_value(const char*);
    void errortrace(FILE*fp=stderr, int offset = 0);
    void dump(FILE*fp=stderr,int offset = 0);
};

class PrefixKeyVal : public KeyVal {
  private:
    int nprefix;
    char** prefices;
    KeyVal* keyval;
    void setup(const char*,int,int,int,int,int);
    int getnewprefixkey(const char*key,char*newkey);
  public:
    PrefixKeyVal(const char*,KeyVal&);
    PrefixKeyVal(const char*,KeyVal&,int);
    PrefixKeyVal(const char*,KeyVal&,int,int);
    PrefixKeyVal(const char*,KeyVal&,int,int,int);
    PrefixKeyVal(const char*,KeyVal&,int,int,int,int);
    ~PrefixKeyVal();
    int    key_exists(const char*);
    RefKeyValValue key_value(const char*);
    void errortrace(FILE*fp=stderr, int offset = 0);
    void dump(FILE*fp=stderr,int offset = 0);
};

class IPV2;
class ParsedKeyVal : public StringKeyVal {
  private:
    int nfile;
    char**file;
    int nfp;
    IPV2* ipv2;
  public:
    ParsedKeyVal();
    ParsedKeyVal(const char*);
    ParsedKeyVal(FILE*);
    ~ParsedKeyVal();
    const char* stringvalue(const char*);
    virtual const char* classname(const char*);
    virtual const char* truekeyword(const char*);
    void read(const char*);
    void read(FILE*);
    void errortrace(FILE*fp=stderr, int offset = 0);
    void dump(FILE*fp=stderr,int offset = 0);
};

#endif /* _KeyVal_h */
