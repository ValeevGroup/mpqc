
#ifndef _KeyVal_h
#define _KeyVal_h
#ifdef __GNUG__
#pragma interface
#endif
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include <util/class/class.h>
#include <util/container/ref.h>

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

class KeyVal: public VRefCount {
  public:
    enum {MaxKeywordLength = 256};
    enum KeyValError { OK, HasNoValue, WrongType,
                       UnknownKeyword, OperationFailed };
  private:
    KeyValError errcod;
    // do not allow a copy constructor or assignment
    KeyVal(const KeyVal&);
    operator=(const KeyVal&);
  protected:
    KeyVal();
    void seterror(KeyValError err);

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
    // then the overloaded functions will be hidden.  The key_... functions
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

    int    exists(int i);
    int    count(int i);
    int    booleanvalue(int i);
    double doublevalue(int i);
    float  floatvalue(int i);
    char   charvalue(int i);
    int    intvalue(int i);
    char*  pcharvalue(int i);
    RefDescribedClass describedclassvalue(int i);

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

    int    exists(int i,int j);
    int    count(int i,int j);
    int    booleanvalue(int i,int j);
    double doublevalue(int i,int j);
    float  floatvalue(int i,int j);
    char   charvalue(int i,int j);
    int    intvalue(int i,int j);
    char*  pcharvalue(int i,int j);
    RefDescribedClass describedclassvalue(int i,int j);

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
    static double Defaultdouble();
    static int    Defaultint();
    static float  Defaultfloat();
    static char   Defaultchar();
    static char*  Defaultpchar();
    static int    Defaultboolean();
    static RefDescribedClass DefaultRefDescribedClass();

    KeyValError error();
    char*  errormsg(KeyValError);
    char*  errormsg();
    virtual void errortrace(FILE*fp=stderr,int offset = 0);
    virtual void dump(FILE*fp=stderr,int offset = 0);
};

REF_dec(KeyVal);

class KeyValValue: public VRefCount {
  protected:
    KeyValValue();
    KeyValValue(KeyValValue&);
  public:
    virtual ~KeyValValue();
    virtual KeyVal::KeyValError doublevalue(double&);
    virtual KeyVal::KeyValError booleanvalue(int&);
    virtual KeyVal::KeyValError floatvalue(float&);
    virtual KeyVal::KeyValError charvalue(char&);
    virtual KeyVal::KeyValError intvalue(int&);
    virtual KeyVal::KeyValError pcharvalue(char*&);
    virtual KeyVal::KeyValError describedclassvalue(RefDescribedClass&);
};

REF_dec(KeyValValue);

class KeyValValuedouble: public KeyValValue {
  private:
    double _val;
  public:
    KeyValValuedouble(double);
    KeyValValuedouble(const KeyValValuedouble&);
    ~KeyValValuedouble();
    KeyVal::KeyValError doublevalue(double&);
};

class KeyValValueboolean: public KeyValValue {
  private:
    int _val;
  public:
    KeyValValueboolean(int);
    KeyValValueboolean(const KeyValValueboolean&);
    ~KeyValValueboolean();
    KeyVal::KeyValError booleanvalue(int&);
};

class KeyValValuefloat: public KeyValValue {
  private:
    float _val;
  public:
    KeyValValuefloat(float);
    KeyValValuefloat(const KeyValValuefloat&);
    ~KeyValValuefloat();
    KeyVal::KeyValError floatvalue(float&);
};

class KeyValValuechar: public KeyValValue {
  private:
    char _val;
  public:
    KeyValValuechar(char);
    KeyValValuechar(const KeyValValuechar&);
    ~KeyValValuechar();
    KeyVal::KeyValError charvalue(char&);
};

class KeyValValueint: public KeyValValue {
  private:
    int _val;
  public:
    KeyValValueint(int);
    KeyValValueint(const KeyValValueint&);
    ~KeyValValueint();
    KeyVal::KeyValError intvalue(int&);
};

class KeyValValuepchar: public KeyValValue {
  private:
    char* _val;
  public:
    KeyValValuepchar(const char*);
    KeyValValuepchar(const KeyValValuepchar&);
    ~KeyValValuepchar();
    KeyVal::KeyValError pcharvalue(char*&);
};

class KeyValValueRefDescribedClass: public KeyValValue {
  private:
    RefDescribedClass _val;
  public:
    KeyValValueRefDescribedClass(RefDescribedClass&);
    KeyValValueRefDescribedClass(KeyValValueRefDescribedClass&);
    ~KeyValValueRefDescribedClass();
    KeyVal::KeyValError describedclassvalue(RefDescribedClass&);
};

class KeyValValueString: public KeyValValue {
  private:
    const char* _val;
  public:
    KeyValValueString(const char*);
    KeyValValueString(const KeyValValueString&);
    ~KeyValValueString();
    KeyVal::KeyValError doublevalue(double&);
    KeyVal::KeyValError booleanvalue(int&);
    KeyVal::KeyValError floatvalue(float&);
    KeyVal::KeyValError charvalue(char&);
    KeyVal::KeyValError intvalue(int&);
    KeyVal::KeyValError pcharvalue(char*&);
};

// this class allows keyval associations to be set up by the program,
// rather than determined by an external file
class KeyValKeywordRefKeyValValueMap;
class AssignedKeyVal: public KeyVal {
  private:
    KeyValKeywordRefKeyValValueMap* _map;
    // do not allow a copy constructor or assignment
    AssignedKeyVal(const AssignedKeyVal&);
    operator=(const AssignedKeyVal&);
  public:
    AssignedKeyVal();
    ~AssignedKeyVal();
    int    key_exists(const char*);
    RefKeyValValue key_value(const char*);

    void assign(const char*, const RefKeyValValue&);
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
    // do not allow a copy constructor or assignment
    StringKeyVal(const StringKeyVal&);
    operator=(const StringKeyVal&);
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
    RefKeyVal kv[MaxKeyVal];
    RefKeyVal getkeyval(const char*key);
    // do not allow a copy constructor or assignment
    AggregateKeyVal(const AggregateKeyVal&);
    operator=(const AggregateKeyVal&);
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
    RefKeyVal keyval;
    void setup(const char*,int,int,int,int,int);
    int getnewprefixkey(const char*key,char*newkey);
    // do not allow a copy constructor or assignment
    PrefixKeyVal(const PrefixKeyVal&);
    operator=(const PrefixKeyVal&);
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
    // do not allow a copy constructor or assignment
    ParsedKeyVal(const ParsedKeyVal&);
    operator=(const ParsedKeyVal&);
  public:
    ParsedKeyVal();
    ParsedKeyVal(const char*);
    ParsedKeyVal(FILE*);
    ParsedKeyVal(IPV2*);
    // This ctor is given a string which is used to form keywords
    // that are sought in the keyval argument.  The associated values
    // are used to construct file names that are used to initialize
    // the parsedkeyval.  The keywords sought are string'dir' for the
    // directory prefix and string'files' for an array of file names.
    ParsedKeyVal(const char*,KeyVal&);
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
