//
// keyval.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef _util_keyval_keyval_h
#define _util_keyval_keyval_h
#ifdef __GNUG__
#pragma interface
#endif

#include <iostream.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include <util/container/avlmap.h>
#include <util/class/class.h>
#include <util/keyval/keyvalval.h>

class KeyValKeyword {
  private:
    char* keyword_;
  public:
    KeyValKeyword();
    KeyValKeyword(const char* name);
    KeyValKeyword(const KeyValKeyword&);
    ~KeyValKeyword();
    KeyValKeyword& operator=(const KeyValKeyword&);
    int operator==(const KeyValKeyword& ck) const;
    int operator<(const KeyValKeyword& ck) const;
    int hash() const;
    inline int cmp(const KeyValKeyword&ck) const
    {
      if (!keyword_) {
          if (!ck.keyword_) return 0;
          return -1;
        }
      if (!ck.keyword_) return 1;
      return strcmp(keyword_,ck.keyword_);
    }
    inline const char* name() const {return keyword_;}
  };

REF_fwddec(KeyValValue);

//.
// The \clsnm{KeyVal} class is designed to simplify the process of allowing
// a user to specify keyword/value associations to a C++ program.  A
// flexible input style and ease of use for the programmer is achieved with
// this method.  Keywords are represented by null terminated character arrays.
// The keywords are organized hierarchially, in a manner similar to the way
// that many file systems are organized.  One character is special,
// '\srccd{:}', which is used to separate the various hierarchial labels,
// which are referred to as ``segments'', in the keyword.
//
// A convention for specifying arrays is provided by \clsnm{KeyVal}.  Each
// index of the array is given by appending a segment containing the
// character representation of the index.  Thus, \srccd{array:3:4} would be
// a the keyword corresponding to third row and fourth column of
// \srccd{array}.
//
// To allow the \clsnm{KeyVal} class to have associations that can represent
// data for classes, the keyword can be associated with a class as well as
// a value.  This permits polymorphic data to be unambiguously represented
// by keyword/value associations.  Most use of \clsnm{KeyVal} need not be
// concerned with this.
class KeyVal: public VRefCount {
    // these classes need to directly access the key_value member
    friend class AggregateKeyVal;
    friend class PrefixKeyVal;
  public:
    enum {MaxKeywordLength = 256};
    enum KeyValError { OK, HasNoValue, WrongType,
                       UnknownKeyword, OperationFailed };
  private:
    KeyValError errcod;
    // do not allow a copy constructor or assignment
    KeyVal(const KeyVal&);
    void operator=(const KeyVal&);
  protected:
    int verbose_;

    KeyVal();

    //. Set the current error condition.
    void seterror(KeyValError err);
    void seterror(KeyValValue::KeyValValueError err);

    //. Ultimately called by \srccd{exists}.
    virtual int    key_exists(const char*) = 0;
    //. Ultimately called by \srccd{count}.
    virtual int    key_count(const char* =0);
    //. Ultimately called by \srccd{value}.
    virtual RefKeyValValue key_value(const char*,
                                     const KeyValValue& def) = 0;
    //. Ultimately called by \srccd{booleanvalue}.
    virtual int    key_booleanvalue(const char*,const KeyValValue& def);
    //. Ultimately called by \srccd{doublevalue}.
    virtual double key_doublevalue(const char* key,const KeyValValue& def);
    //. Ultimately called by \srccd{floatvalue}.
    virtual float  key_floatvalue(const char* key,const KeyValValue& def);
    //. Ultimately called by \srccd{charvalue}.
    virtual char   key_charvalue(const char* key,const KeyValValue& def);
    //. Ultimately called by \srccd{intvalue}.
    virtual int    key_intvalue(const char* key,const KeyValValue& def);
    //. Ultimately called by \srccd{pcharvalue}.
    virtual char*  key_pcharvalue(const char* key,const KeyValValue& def);
    //. Ultimately called by \srccd{describedclassvalue}.
    virtual RefDescribedClass key_describedclassvalue(const char* key,
                                                      const KeyValValue& def);

  public:
    virtual ~KeyVal();

    // For nonindexed things.   If a subclass defines one of these,
    // then the overloaded functions will be hidden.  The key_... functions
    // should be overridden instead.

    //. This takes as its only argument a keyword.
    // Returns 1 if the keyword has a value and 0 otherwise.
    int    exists(const char*);
    //. If the value of a keyword is an array, then return its length.
    // If no arguments are given then the top level will be checked to
    // see if it is an array and, if so, the number of elements will be
    // counted.
    int    count(const char* =0);
    //. Return the value associated with the keyword.
    RefKeyValValue value(const char* = 0,
                         const KeyValValue& def=KeyValValue());
    //. Returns the boolean value (0 = false, 1 = true) of \vrbl{key}.
    int    booleanvalue(const char* key = 0,
                        const KeyValValue& def=KeyValValueboolean());
    //. Returns the double value of \vrbl{key}.
    double doublevalue(const char* key = 0,
                       const KeyValValue& def=KeyValValuedouble());
    //. Returns the float value of \vrbl{key}.
    float  floatvalue(const char* key = 0,
                      const KeyValValue& def=KeyValValuefloat());
    //. Returns the char value of \vrbl{key}.
    char   charvalue(const char* key = 0,
                     const KeyValValue& def=KeyValValuechar());
    //. Returns the int value of \vrbl{key}.
    int    intvalue(const char* key = 0,
                    const KeyValValue& def=KeyValValueint());
    //. Returns a copy of the string representation of the \vrbl{key}'s
    // value. Storage for the copy is obtained with new.
    char*  pcharvalue(const char* key = 0,
                      const KeyValValue& def=KeyValValuepchar());
    //. Returns a reference to an object of type DescribedClass
    // (@pxref{The DescribedClass Class}).
    RefDescribedClass describedclassvalue(const char* key = 0,
                     const KeyValValue& def=KeyValValueRefDescribedClass());

    // For vectors:
    int    exists(const char*,int);
    int    count(const char*,int);
    int    booleanvalue(const char*,int,
                        const KeyValValue& def=KeyValValueboolean());
    double doublevalue(const char* key,int,
                       const KeyValValue& def=KeyValValuedouble());
    float  floatvalue(const char* key,int,
                      const KeyValValue& def=KeyValValuefloat());
    char   charvalue(const char* key,int,
                     const KeyValValue& def=KeyValValuechar());
    int    intvalue(const char* key,int,
                    const KeyValValue& def=KeyValValueint());
    char*  pcharvalue(const char* key,int,
                      const KeyValValue& def=KeyValValuepchar());
    RefDescribedClass describedclassvalue(const char* key,int,
                     const KeyValValue& def=KeyValValueRefDescribedClass());

    int    exists(int i);
    int    count(int i);
    int    booleanvalue(int i,
                        const KeyValValue& def=KeyValValueboolean());
    double doublevalue(int i,
                       const KeyValValue& def=KeyValValuedouble());
    float  floatvalue(int i,
                      const KeyValValue& def=KeyValValuefloat());
    char   charvalue(int i,
                     const KeyValValue& def=KeyValValuechar());
    int    intvalue(int i,
                    const KeyValValue& def=KeyValValueint());
    char*  pcharvalue(int i,
                      const KeyValValue& def=KeyValValuepchar());
    RefDescribedClass describedclassvalue(int i,
                     const KeyValValue& def=KeyValValueRefDescribedClass());

    // For arrays:
    int    exists(const char*,int,int);
    int    count(const char*,int,int);
    int    booleanvalue(const char*,int,int,
                        const KeyValValue& def=KeyValValueboolean());
    double doublevalue(const char* key,int,int,
                       const KeyValValue& def=KeyValValuedouble());
    float  floatvalue(const char* key,int,int,
                      const KeyValValue& def=KeyValValuefloat());
    char   charvalue(const char* key,int,int,
                     const KeyValValue& def=KeyValValuechar());
    int    intvalue(const char* key,int,int,
                    const KeyValValue& def=KeyValValueint());
    char*  pcharvalue(const char* key,int,int,
                      const KeyValValue& def=KeyValValuepchar());
    RefDescribedClass describedclassvalue(const char* key,int,int,
                     const KeyValValue& def=KeyValValueRefDescribedClass());

    int    exists(int i,int j);
    int    count(int i,int j);
    int    booleanvalue(int i,int j,
                        const KeyValValue& def=KeyValValueboolean());
    double doublevalue(int i,int j,
                       const KeyValValue& def=KeyValValuedouble());
    float  floatvalue(int i,int j,
                      const KeyValValue& def=KeyValValuefloat());
    char   charvalue(int i,int j,
                     const KeyValValue& def=KeyValValuechar());
    int    intvalue(int i,int j,
                    const KeyValValue& def=KeyValValueint());
    char*  pcharvalue(int i,int j,
                      const KeyValValue& def=KeyValValuepchar());
    RefDescribedClass describedclassvalue(int i,int j,
                     const KeyValValue& def=KeyValValueRefDescribedClass());

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

    //. Return the current error condition.
    KeyValError error();
    //. Return a textual representation of \vrbl{err}.  The current error
    // condition will be used if the argument is omitted.
    char*  errormsg(KeyValError err);
    //. Return a textual representation of the current error.
    char*  errormsg();

    virtual void errortrace(ostream&fp=cerr);
    virtual void dump(ostream&fp=cerr);

    //. Print keywords that were never looked at, if possible.
    virtual void print_unseen(ostream&fp=cout);
    //. Return 1 if there were unseen keywords, 0 if there are
    //none, or -1 this keyval doesn't keep track of unseen
    //keywords.
    virtual int have_unseen();

    //. Control printing of assignments.
    void verbose(int v) { verbose_ = v; }
    int verbose() const { return verbose_; }
};

REF_dec(KeyVal);

// this class allows keyval associations to be set up by the program,
// rather than determined by an external file
class AssignedKeyVal: public KeyVal {
  private:
    AVLMap<KeyValKeyword,RefKeyValValue> _map;
    // do not allow a copy constructor or assignment
    AssignedKeyVal(const AssignedKeyVal&);
    void operator=(const AssignedKeyVal&);
  protected:
    int    key_exists(const char*);
    RefKeyValValue key_value(const char*,
                             const KeyValValue& def);
  public:
    AssignedKeyVal();
    ~AssignedKeyVal();

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
    AVLMap<KeyValKeyword,RefKeyValValue> _map;
    // do not allow a copy constructor or assignment
    StringKeyVal(const StringKeyVal&);
    void operator=(const StringKeyVal&);
  protected:
    StringKeyVal();
    int    key_exists(const char*);
    RefKeyValValue key_value(const char*,
                             const KeyValValue& def);
  public:
    virtual ~StringKeyVal();
    virtual const char* stringvalue(const char *) = 0;
    // returns the name of the exact class the object at the keyword
    virtual const char* classname(const char*);
    // returns a string which is the actual keyword if some sort
    // of variable substitution takes place (needed to make multiple
    // references to the same object work in input files)
    virtual const char* truekeyword(const char*);

    virtual void errortrace(ostream&fp=cerr);
    virtual void dump(ostream&fp=cerr);
};

class AggregateKeyVal : public KeyVal {
  private:
    enum { MaxKeyVal = 4 };
    RefKeyVal kv[MaxKeyVal];
    RefKeyVal getkeyval(const char*key);
    // do not allow a copy constructor or assignment
    AggregateKeyVal(const AggregateKeyVal&);
    void operator=(const AggregateKeyVal&);
  protected:
    int    key_exists(const char*);
    RefKeyValValue key_value(const char*,
                             const KeyValValue& def);
  public:
    AggregateKeyVal(const RefKeyVal&);
    AggregateKeyVal(const RefKeyVal&,const RefKeyVal&);
    AggregateKeyVal(const RefKeyVal&,const RefKeyVal&,const RefKeyVal&);
    AggregateKeyVal(const RefKeyVal&,const RefKeyVal&,const RefKeyVal&,
                    const RefKeyVal&);
    ~AggregateKeyVal();
    void errortrace(ostream&fp=cerr);
    void dump(ostream&fp=cerr);
};

class PrefixKeyVal : public KeyVal {
  private:
    char* prefix;
    RefKeyVal keyval;
    void setup(const char*,int,int,int,int,int);
    int getnewprefixkey(const char*key,char*newkey);
    // do not allow a copy constructor or assignment
    PrefixKeyVal(const PrefixKeyVal&);
    void operator=(const PrefixKeyVal&);
    int    key_exists(const char*);
    RefKeyValValue key_value(const char*,
                             const KeyValValue& def);
  public:
    PrefixKeyVal(const RefKeyVal&,int);
    PrefixKeyVal(const RefKeyVal&,int,int);
    PrefixKeyVal(const RefKeyVal&,int,int,int);
    PrefixKeyVal(const RefKeyVal&,int,int,int,int);
    PrefixKeyVal(const RefKeyVal&,const char*);
    PrefixKeyVal(const RefKeyVal&,const char*,int);
    PrefixKeyVal(const RefKeyVal&,const char*,int,int);
    PrefixKeyVal(const RefKeyVal&,const char*,int,int,int);
    PrefixKeyVal(const RefKeyVal&,const char*,int,int,int,int);
    // old CTOR syntax (use the above instead)
    PrefixKeyVal(const char*,const RefKeyVal&);
    PrefixKeyVal(const char*,const RefKeyVal&,int);
    PrefixKeyVal(const char*,const RefKeyVal&,int,int);
    PrefixKeyVal(const char*,const RefKeyVal&,int,int,int);
    PrefixKeyVal(const char*,const RefKeyVal&,int,int,int,int);
    ~PrefixKeyVal();
    void errortrace(ostream&fp=cerr);
    void dump(ostream&fp=cerr);
};

class IPV2;
//. The \clsnm{ParsedKeyVal} class converts textual information into
//keyword/value assocations.  The parsing is done with an \clsnm{IPV2}
//object, which is the encapsulation of the ipv2 C library.
class ParsedKeyVal : public StringKeyVal {
  private:
    int nfile;
    char**file;
    int nfp;
    IPV2* ipv2;
    // do not allow a copy constructor or assignment
    ParsedKeyVal(const ParsedKeyVal&);
    void operator=(const ParsedKeyVal&);
  public:
    //. Create an empty \clsnm{ParsedKeyVal}.
    ParsedKeyVal();
    //. Parse the given input file.
    ParsedKeyVal(const char*);
    //. Read input from the given \clsnm{istream}.
    ParsedKeyVal(istream&);
    //. Use the given \clsnm{IPV2*} object.  The new \clsnm{ParsedKeyVal}
    //takes wnership of the passed \clsnm{IPV2} object.
    ParsedKeyVal(IPV2*);
    //. This ctor is given a string which is used to form keywords
    // that are sought in the keyval argument.  The associated values
    // are used to construct file names that are used to initialize
    // the parsedkeyval.  The keywords sought are string'dir' for the
    // directory prefix and string'files' for an array of file names.
    ParsedKeyVal(const char*,const RefKeyVal&);
    //. Cleanup, deleting the \clsnm{IPV2} object.
    ~ParsedKeyVal();

    //. This is like the \srccd{ParsedKeyVal(const char*,const RefKeyVal&)}
    // ctor, but writes the contents of the files to the given ostream.
    static void cat_files(const char*,const RefKeyVal&,ostream &o);

    //. Read input data from the given filename, stream, or string.
    void read(const char*);
    void read(istream&);
    void parse_string(const char *);

    //. Overrides of parent members.
    const char* stringvalue(const char*);
    const char* classname(const char*);
    const char* truekeyword(const char*);
    void errortrace(ostream&fp=cerr);
    void dump(ostream&fp=cerr);
    void print_unseen(ostream&fp=cout);
    int have_unseen();
};

REF_dec(ParsedKeyVal);

#endif /* _KeyVal_h */

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
