
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

//texi
// The @code{KeyVal} class is designed to simplify the process of allowing
// a user to specify keyword/value associations to a C++ program.  A
// flexible input style and ease of use for the programmer is achieved with
// this method.  Keywords are represented by null terminated character arrays.
// The keywords are organized hierarchially, in a manner similar to the way
// that many file systems are organized.  One character is special,
// '@code{:}', which is used to separate the various hierarchial labels,
// which are referred to as ``segments'', in the keyword.
//
// A convention for specifying arrays is provided by @code{KeyVal}.  Each
// index of the array is given by appending a segment containing the
// character representation of the index.  Thus, @code{array:3:4} would be
// a the keyword corresponding to third row and fourth column of
// @code{array}.
//
// To allow the @code{KeyVal} class to have associations that can represent
// data for classes, the keyword can be associated with a class as well as
// a value.  This permits polymorphic data to be unambiguously represented
// by keyword/value associations.  Most use of @code{KeyVal} need not be
// concerned with this.
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

    //texi Set the current error condition.
    void seterror(KeyValError err);

    void offset(FILE* fp,int n); // Put n ' ' into fp.
    enum {OffsetDelta=4};

    //texi Ultimately called by @code{exists}.
    virtual int    key_exists(const char*) = 0;
    //texi Ultimately called by @code{count}.
    virtual int    key_count(const char* =0);
    //texi Ultimately called by @code{value}.
    virtual RefKeyValValue key_value(const char*) = 0;
    //texi Ultimately called by @code{booleavalue}.
    virtual int    key_booleanvalue(const char*);
    //texi Ultimately called by @code{doublevalue}.
    virtual double key_doublevalue(const char* key);
    //texi Ultimately called by @code{floatvalue}.
    virtual float  key_floatvalue(const char* key);
    //texi Ultimately called by @code{charvalue}.
    virtual char   key_charvalue(const char* key);
    //texi Ultimately called by @code{intvalue}.
    virtual int    key_intvalue(const char* key);
    //texi Ultimately called by @code{pcharvalue}.
    virtual char*  key_pcharvalue(const char* key);
    //texi Ultimately called by @code{describedclassvalue}.
    virtual RefDescribedClass key_describedclassvalue(const char* key);

  public:
    virtual ~KeyVal();

    // For nonindexed things.   If a subclass defines one of these,
    // then the overloaded functions will be hidden.  The key_... functions
    // should be overridden instead.

    //texi This takes as its only argument a keyword.
    // Returns 1 if the keyword has a value and 0 otherwise.
    int    exists(const char*);
    //texi If the value of a keyword is an array, then return its length.
    // If no arguments are given then the top level will be checked to
    // see if it is an array and, if so, the number of elements will be
    // counted.
    int    count(const char* =0);
    //texi Return the value associated with the keyword.
    RefKeyValValue value(const char*);
    //texi Returns the boolean value (0 = false, 1 = true) of @var{key}.
    int    booleanvalue(const char* key);
    //texi Returns the double value of @var{key}.
    double doublevalue(const char* key);
    //texi Returns the float value of @var{key}.
    float  floatvalue(const char* key);
    //texi Returns the char value of @var{key}.
    char   charvalue(const char* key);
    //texi Returns the int value of @var{key}.
    int    intvalue(const char* key);
    //texi Returns a copy of the string representation of the @var{key}'s
    // value. Storage for the copy is obtained with new.
    char*  pcharvalue(const char* key);
    //texi Returns a reference to an object of type DescribedClass
    // (@pxref{The DescribedClass Class}).
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

    //texi Return the current error condition.
    KeyValError error();
    //texi Return a textual representation of @var{err}.  The current error
    // condition will be used if the argument is omitted.
    char*  errormsg(KeyValError err);
    //texi Return a textual representation of the current error.
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
    KeyValValueRefDescribedClass(const RefDescribedClass&);
    KeyValValueRefDescribedClass(const KeyValValueRefDescribedClass&);
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
    void assign(const char*, const RefDescribedClass&);

    void clear();
};

REF_dec(AssignedKeyVal);

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
    AggregateKeyVal(const RefKeyVal&);
    AggregateKeyVal(const RefKeyVal&,const RefKeyVal&);
    AggregateKeyVal(const RefKeyVal&,const RefKeyVal&,const RefKeyVal&);
    AggregateKeyVal(const RefKeyVal&,const RefKeyVal&,const RefKeyVal&,
                    const RefKeyVal&);
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
    PrefixKeyVal(const char*,const RefKeyVal&);
    PrefixKeyVal(const char*,const RefKeyVal&,int);
    PrefixKeyVal(const char*,const RefKeyVal&,int,int);
    PrefixKeyVal(const char*,const RefKeyVal&,int,int,int);
    PrefixKeyVal(const char*,const RefKeyVal&,int,int,int,int);
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
    ParsedKeyVal(const char*,const RefKeyVal&);
    ~ParsedKeyVal();
    const char* stringvalue(const char*);
    virtual const char* classname(const char*);
    virtual const char* truekeyword(const char*);
    void read(const char*);
    void read(FILE*);
    void errortrace(FILE*fp=stderr, int offset = 0);
    void dump(FILE*fp=stderr,int offset = 0);
};

REF_dec(ParsedKeyVal);

#endif /* _KeyVal_h */
