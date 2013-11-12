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

#include <iostream>
#include <string>
#include <map>
#include <vector>

#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include <util/class/class.h>
#include <util/misc/scexception.h>
#include <util/keyval/keyvalval.h>

namespace sc {

/**
 * @ingroup CoreKeyVal
 The KeyVal class is designed to simplify the process of allowing
 a user to specify keyword/value associations to a C++ program.  A
 flexible input style and ease of use for the programmer is achieved with
 this method.  Keywords are represented by null terminated character arrays.
 The keywords are organized hierarchially, in a manner similar to the way
 that many file systems are organized.  One character is special,
 ":", which is used to separate the various hierarchial labels,
 which are referred to as "segments", in the keyword.

 A convention for specifying arrays is provided by KeyVal.  Each
 index of the array is given by appending a segment containing the
 character representation of the index.  Thus, "array:3:4" would be
 a the keyword corresponding to fourth row and fifth column of
 "array", since indexing starts at zero.

 To allow the KeyVal class to have associations that can represent
 data for classes, the keyword can be associated with a class as well as
 a value.  This permits polymorphic data to be unambiguously represented
 by keyword/value associations.  Most use of KeyVal need not be
 concerned with this.
*/
class KeyVal: public RefCount {
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

    /// Set the current error condition.
    void seterror(KeyValError err);
    /// Set the current error condition.
    void seterror(KeyValValue::KeyValValueError err);

    /// Ultimately called by exists.
    virtual int    key_exists(const char*) = 0;
    /// Ultimately called by count.
    virtual int    key_count(const char* =0);
    /// Ultimately called by value.
    virtual Ref<KeyValValue> key_value(const char*,
                                     const KeyValValue& def) = 0;
    /// Ultimately called by booleanvalue.
    virtual int    key_booleanvalue(const char*,const KeyValValue& def);
    /// Ultimately called by doublevalue.
    virtual double key_doublevalue(const char* key,const KeyValValue& def);
    /// Ultimately called by floatvalue.
    virtual float  key_floatvalue(const char* key,const KeyValValue& def);
    /// Ultimately called by charvalue.
    virtual char   key_charvalue(const char* key,const KeyValValue& def);
    /// Ultimately called by intvalue.
    virtual int    key_intvalue(const char* key,const KeyValValue& def);
    /// Ultimately called by longvalue.
    virtual long   key_longvalue(const char* key,const KeyValValue& def);
    /// Ultimately called by sizevalue.
    virtual size_t key_sizevalue(const char* key,const KeyValValue& def);
    /// Ultimately called by pcharvalue.
    DEPRECATED virtual char*  key_pcharvalue(const char* key,const KeyValValue& def);
    /// Ultimately called by stringvalue.
    virtual std::string key_stringvalue(const char* key,
                                        const KeyValValue& def);
    /// Ultimately called by describedclassvalue.
    virtual Ref<DescribedClass> key_describedclassvalue(const char* key,
                                                      const KeyValValue& def);

  public:
    virtual ~KeyVal();

    // For nonindexed things.   If a subclass defines one of these,
    // then the overloaded functions will be hidden.  The key_... functions
    // should be overridden instead.

    /** This takes as its only argument a keyword.
        Returns 1 if the keyword has a value and 0 otherwise. */
    int    exists(const char*);
    /** If the value of a keyword is an array, then return its length.
        If no arguments are given then the top level will be checked to
        see if it is an array and, if so, the number of elements will be
        counted. */
    int    count(const char* =0);
    /// Return the value associated with the keyword.
    Ref<KeyValValue> value(const char* key = 0,
                         const KeyValValue& def=KeyValValue());
    /// Returns the boolean value (0 = false, 1 = true) of key.
    int    booleanvalue(const char* key = 0,
                        const KeyValValue& def=KeyValValueboolean());
    /// Returns the double value of key.
    double doublevalue(const char* key = 0,
                       const KeyValValue& def=KeyValValuedouble());
    /// Returns the float value of key.
    float  floatvalue(const char* key = 0,
                      const KeyValValue& def=KeyValValuefloat());
    /// Returns the char value of key.
    char   charvalue(const char* key = 0,
                     const KeyValValue& def=KeyValValuechar());
    /// Returns the int value of key.
    int    intvalue(const char* key = 0,
                    const KeyValValue& def=KeyValValueint());
    /// Returns the long value of key.
    long   longvalue(const char* key = 0,
                     const KeyValValue& def=KeyValValuelong());
    /// Returns the size_t value of key.
    size_t sizevalue(const char* key = 0,
                     const KeyValValue& def=KeyValValuesize());
    /** Returns a copy of the string representation of the key's
        value. Storage for the copy is obtained with new.
        This is deprecated--use stringvalue instead. */
    DEPRECATED char*  pcharvalue(const char* key = 0,
                      const KeyValValue& def=KeyValValuestring());
    /** Returns a string representation of the key's value. */
    std::string stringvalue(const char* key = 0,
                            const KeyValValue& def=KeyValValuestring());
    /// Returns a reference to an object of type DescribedClass.
    Ref<DescribedClass> describedclassvalue(const char* key = 0,
                     const KeyValValue& def=KeyValValueRefDescribedClass());
    /** Returns the name of the exact class of the object at the keyword.
        If no classname is assigned, or this KeyVal does not support classes,
        then 0 is returned. */
    virtual const char* classname(const char*);

    /// Returns a reference to an object of type DescribedClass using the top level keywords of this KeyVal
    /// @param[in] classname specifies the class name, must be derived from DescribedClass
    /// @return Ref<> to the object, null if constructor not successful
    virtual Ref<DescribedClass> describedclass(const char* classname);

    /** @name Reading Vectors.
        These members correspond to the above members, but take
        an additional integer argument, i, which is a vector index.
        This is equivalent to getting a value for a keyword named
        "<i>key</i>:<i>i</i>".  The routines that do not take
        key arguments get the value for the keyword named "<i>i</i>".
     */
    //@{
    int    exists(const char* key,int i);
    int    count(const char* key,int i);
    int    booleanvalue(const char* key,int i,
                        const KeyValValue& def=KeyValValueboolean());
    double doublevalue(const char* key,int i,
                       const KeyValValue& def=KeyValValuedouble());
    float  floatvalue(const char* key,int i,
                      const KeyValValue& def=KeyValValuefloat());
    char   charvalue(const char* key,int i,
                     const KeyValValue& def=KeyValValuechar());
    int    intvalue(const char* key,int i,
                    const KeyValValue& def=KeyValValueint());
    long   longvalue(const char* key,int i,
                     const KeyValValue& def=KeyValValuelong());
    size_t sizevalue(const char* key,int i,
                     const KeyValValue& def=KeyValValuesize());
    DEPRECATED char*  pcharvalue(const char* key,int i,
                      const KeyValValue& def=KeyValValuestring());
    std::string stringvalue(const char* key,int i,
                            const KeyValValue& def=KeyValValuestring());
    Ref<DescribedClass> describedclassvalue(const char* key,int,
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
    long   longvalue(int i,
                     const KeyValValue& def=KeyValValuelong());
    size_t sizevalue(int i,
                     const KeyValValue& def=KeyValValuesize());
    DEPRECATED char*  pcharvalue(int i,
                      const KeyValValue& def=KeyValValuestring());
    std::string stringvalue(int i,
                            const KeyValValue& def=KeyValValuestring());
    Ref<DescribedClass> describedclassvalue(int i,
                     const KeyValValue& def=KeyValValueRefDescribedClass());
    //@}

    /** @name Reading 2D Arrays.
        These members correspond to the above members, but take additional
        integer arguments, i and j, which is an array index.  This is
        equivalent to getting a value for a keyword named
        "<i>key</i>:<i>i</i>:<i>j</i>".  The routines that do not take key
        arguments get the value for the keyword named "<i>i</i>:<i>j</i>".  */
    //@{
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
    long   longvalue(const char* key,int,int,
                     const KeyValValue& def=KeyValValuelong());
    size_t sizevalue(const char* key,int,int,
                     const KeyValValue& def=KeyValValuesize());
    DEPRECATED char*  pcharvalue(const char* key,int,int,
                      const KeyValValue& def=KeyValValuestring());
    std::string stringvalue(const char* key,int,int,
                            const KeyValValue& def=KeyValValuestring());
    Ref<DescribedClass> describedclassvalue(const char* key,int,int,
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
    long   longvalue(int i,int j,
                     const KeyValValue& def=KeyValValuelong());
    size_t sizevalue(int i,int j,
                     const KeyValValue& def=KeyValValuesize());
    DEPRECATED char*  pcharvalue(int i,int j,
                      const KeyValValue& def=KeyValValuestring());
    std::string stringvalue(int i,int j,
                            const KeyValValue& def=KeyValValuestring());
    Ref<DescribedClass> describedclassvalue(int i,int j,
                     const KeyValValue& def=KeyValValueRefDescribedClass());
    //@}

    /** @name Reading 3D Arrays.
        These members correspond to the above members, but can be used
        to read in arrays with more than two dimensions.  The nindex
        argument is the number of indices in the array.  It is followed
        by an int giving the value of each index.  */
    //@{
    int    Va_exists(const char* key,int nindex,...);
    int    Va_count(const char* key,int nindex,...);
    int    Va_booleanvalue(const char* key,int nindex,...);
    double Va_doublevalue(const char* key,int nindex,...);
    float  Va_floatvalue(const char* key,int nindex,...);
    char   Va_charvalue(const char* key,int nindex,...);
    int    Va_intvalue(const char* key,int nindex,...);
    long   Va_longvalue(const char* key,int nindex,...);
    size_t Va_sizevalue(const char* key,int nindex,...);
    DEPRECATED char*  Va_pcharvalue(const char* key,int nindex,...);
    std::string Va_stringvalue(const char* key,int nindex,...);
    Ref<DescribedClass> Va_describedclassvalue(const char* key,int nindex,...);
    //@}

    /// Return the current error condition.
    KeyValError error();
    /// Return a textual representation of err.
    const char*  errormsg(KeyValError err);
    /// Return a textual representation of the current error.
    const char*  errormsg();
    /// Write a message to fp describing the error.
    virtual void errortrace(std::ostream&fp=ExEnv::err0());
    /// Write a message to fp describing the error.
    virtual void dump(std::ostream&fp=ExEnv::err0());

    /// Print keywords that were never looked at, if possible.
    virtual void print_unseen(std::ostream&fp=ExEnv::out0());
    /** Return 1 if there were unseen keywords, 0 if there are
        none, or -1 this keyval doesn't keep track of unseen
        keywords. */
    virtual int have_unseen();

    /// Control printing of assignments.
    void verbose(int v) { verbose_ = v; }
    /// Returns nonzero if assignments are printed.
    int verbose() const { return verbose_; }
};



/** @ingroup CoreKeyVal
 *  This class allows keyval associations to be set up by the program,
    rather than determined by an external file. */
class AssignedKeyVal: public KeyVal {
  private:
    typedef std::map<std::string,Ref<KeyValValue> > _map_t;
    _map_t _map;
    // do not allow a copy constructor or assignment
    AssignedKeyVal(const AssignedKeyVal&);
    void operator=(const AssignedKeyVal&);
  protected:
    int    key_exists(const char*);
    Ref<KeyValValue> key_value(const char*,
                             const KeyValValue& def);
  public:
    AssignedKeyVal();
    ~AssignedKeyVal();

    /** @name Assignments.
        Each of this routines assigns key to val.  */
    //@{
    void assign(const char* key, const Ref<KeyValValue>& val);
    void assign(const char* key, double val);
    void assignboolean(const char* key, int val);
    void assign(const char* key, float val);
    void assign(const char* key, char val);
    void assign(const char* key, int val);
    void assign(const char* key, long val);
    void assign(const char* key, const char* val);
    void assign(const char* key, const std::string& val);
    void assign(const char* key, const Ref<DescribedClass>& val);
    //@}

    const char* classname(const char*);

    /// Erase all of the stored assignments.
    void clear();

    /// @param os output stream object
    /// Prints contents to @a os
    void print(std::ostream& os = ExEnv::out0()) const;

    template <typename ValueType>
    static Ref<AssignedKeyVal> instance(const char* key, const ValueType& value) {
      Ref<AssignedKeyVal> result;
      result->assign(key,value);
      return result;
    }
};

/** @ingroup CoreKeyVal
 *  StringKeyVal is a base class for KeyVal implementations
    that store all values in a string format.  These are
    converted to other data types through KeyValValue.
*/
class StringKeyVal: public KeyVal {
  private:
    // once a described class is found it is kept here so
    // multiple references to it return the same instance
    std::map<std::string,Ref<KeyValValue> > _map;
    // do not allow a copy constructor or assignment
    StringKeyVal(const StringKeyVal&);
    void operator=(const StringKeyVal&);
  protected:
    StringKeyVal();
    int    key_exists(const char*);
    Ref<KeyValValue> key_value(const char*,
                             const KeyValValue& def);
  public:
    virtual ~StringKeyVal();
    /// Returns the string representation of the value assigned to key.
    virtual std::string stringrep(const char *key) = 0;
    /** Returns the name of the exact class of the object at the keyword.
        If no classname is assigned then 0 is returned. */
    virtual const char* classname(const char*);
    /** Returns a string which is the actual keyword if some sort
        of variable substitution takes place (needed to make multiple
        references to the same object work in input files). */
    virtual const char* truekeyword(const char*);

    /** @name Debugging.
        See the parent class documentation for descriptions of these functions.
    */
    //@{
    virtual void errortrace(std::ostream&fp=ExEnv::err0());
    virtual void dump(std::ostream&fp=ExEnv::err0());
    //@}
};

/** @ingroup CoreKeyVal
 *  This takes several KeyVal objects and makes them look like
    one KeyVal object.  When a key is sought first KeyVal, then
    the next, and so on is searched until the keyword is found.
*/
class AggregateKeyVal : public KeyVal {
  private:
    enum { MaxKeyVal = 4 };
    Ref<KeyVal> kv[MaxKeyVal];
    Ref<KeyVal> getkeyval(const char*key);
    // do not allow a copy constructor or assignment
    AggregateKeyVal(const AggregateKeyVal&);
    void operator=(const AggregateKeyVal&);
  protected:
    int    key_exists(const char*);
    Ref<KeyValValue> key_value(const char*,
                             const KeyValValue& def);
  public:
    /** @name Constructors.
        These contructors create an AggregateKeyVal that is formed from
        several other KeyVal objects.  The search order is keyval1,
        keyval2, and so on.  All KeyVal objects including and after the
        first null KeyVal will be ignored.
    */
    //@{
    AggregateKeyVal(const Ref<KeyVal>& keyval1);
    AggregateKeyVal(const Ref<KeyVal>& keyval1,const Ref<KeyVal>& keyval2);
    AggregateKeyVal(const Ref<KeyVal>& keyval1,const Ref<KeyVal>& keyval2,
                    const Ref<KeyVal>& keyval3);
    AggregateKeyVal(const Ref<KeyVal>& keyval1,const Ref<KeyVal>& keyval2,
                    const Ref<KeyVal>& keyval3, const Ref<KeyVal>& keyval4);
    //@}
    ~AggregateKeyVal();

    const char* classname(const char*);
    void errortrace(std::ostream&fp=ExEnv::err0());
    void dump(std::ostream&fp=ExEnv::err0());
};

/** @ingroup CoreKeyVal
 *  PrefixKeyVal is a KeyVal that searches a different KeyVal using
    modified keys.  This is convenient for reading keys grouped together
    with a common prefix.  Consider the following code:
    <pre>
    sc::Ref<sc::KeyVal> keyval = new sc::PrefixKeyVal(original_keyval,"A");
    int r = keyval->intvalue("x");
    </pre>
    This code will assign to r the value associated with "x" in keyval.
    keyval will search for "x" by searching for "A::x" in original_keyval.

    This class is important for implementing constructors that take
    KeyVal arguments.  When an object is being constructed from a KeyVal,
    it may contain another object that must be constructed from a KeyVal.
    In order to let the sub-object read the correct keywords from the
    KeyVal, without knowledge of the containing objects keyword prefix,
    a PrefixKeyVal can be constructed.  For example, the code
    \code
    class A: public DescribedClass {
       double f0_;
      public:
       A(const Ref<KeyVal> &keyval): f0_(keyval->doublevalue("f0")) {}
    }
    class B: public DescribedClass {
       double f1_;
       Ref<A> a_;
      public:
       B(const Ref<KeyVal> &keyval):
         f1_(keyval->doublevalue("f1")),
         a_(new PrefixKeyVal(keyval,"a"))
       {}
    };
    \endcode
    can be used to read ParsedKeyVal input that looks like
    <pre>
    b\<B>: (
      f1 = 1.0
      a\<A>: (
        f0 = 2.0
      )
    )
    </pre>
 */
class PrefixKeyVal : public KeyVal {
  private:
    char* prefix;
    Ref<KeyVal> keyval;
    void setup(const char*,int,int,int,int,int);
    int getnewprefixkey(const char*key,char*newkey);
    // do not allow a copy constructor or assignment
    PrefixKeyVal(const PrefixKeyVal&);
    void operator=(const PrefixKeyVal&);
    int    key_exists(const char*);
    Ref<KeyValValue> key_value(const char*,
                             const KeyValValue& def);
  public:
    /** @name Constructors.
        Construct a PrefixKeyVal, using the given prefix and indices. */
    //@{
    PrefixKeyVal(const Ref<KeyVal>&,int i);
    PrefixKeyVal(const Ref<KeyVal>&,int i,int j);
    PrefixKeyVal(const Ref<KeyVal>&,int i,int j,int k);
    PrefixKeyVal(const Ref<KeyVal>&,int i,int j,int k,int l);
    PrefixKeyVal(const Ref<KeyVal>&,const char*prefix);
    PrefixKeyVal(const Ref<KeyVal>&,const char*prefix,int i);
    PrefixKeyVal(const Ref<KeyVal>&,const char*prefix,int i,int j);
    PrefixKeyVal(const Ref<KeyVal>&,const char*prefix,int i,int j,int k);
    PrefixKeyVal(const Ref<KeyVal>&,const char*prefix,int i,int j,int k,int l);
    //@}
    ~PrefixKeyVal();

    const char* classname(const char*);
    void errortrace(std::ostream&fp=ExEnv::err0());
    void dump(std::ostream&fp=ExEnv::err0());
};

class IPV2;
/** @ingroup CoreKeyVal
 *  Converts textual information into keyword/value assocations.  The
    parsing is done with an IPV2 object.  The \ref keyval for more
    information on the input format.  */
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
    /// Create an empty ParsedKeyVal.
    ParsedKeyVal();
    /// Parse the given input file.
    ParsedKeyVal(const char*file);
    /// Read input from s.
    ParsedKeyVal(std::istream&s);
    /** Use the given IPV2* object.  The new ParsedKeyVal
        takes wnership of the passed IPV2 object. */
    ParsedKeyVal(IPV2*);
    /** This ctor is given a string which is used to form keywords
        that are sought in the keyval argument.  The associated values
        are used to construct file names that are used to initialize
        the ParsedKeyVal.  The keywords sought are string'dir' for the
        directory prefix and string'files' for an array of file names. */
    ParsedKeyVal(const char*,const Ref<KeyVal>&);
    /// Cleanup, deleting the IPV2 object.
    ~ParsedKeyVal();

    /** This is like the ParsedKeyVal(const char*,const Ref<KeyVal>&)
        ctor, but writes the contents of the files to the given ostream. */
    static void cat_files(const char*,const Ref<KeyVal>&,std::ostream &o);

    /// Read input data from the given filename
    void read(const char*);
    /// Read input data from the given filename
    void read(const std::string &);
    /// Read input data from the given stream.
    void read(std::istream&);
    /// Read input data from the given string.
    void parse_string(const std::string&);

    /** @name Overrides of parent members.
        See parent class documentation. */
    //@{
    std::string stringrep(const char*);
    const char* classname(const char*);
    const char* truekeyword(const char*);
    void errortrace(std::ostream&fp=ExEnv::err0());
    void dump(std::ostream&fp=ExEnv::err0());
    void print_unseen(std::ostream&fp=ExEnv::out0());
    int have_unseen();
    //@}
};

namespace detail {
  /// GetValue(keyval, key, i) grabs the value corresponding to key
  template <typename T> struct GetValue;
  template <> struct GetValue<bool> {
    static bool eval(const Ref<KeyVal>& kv, const char* key, int i) {
      return kv->booleanvalue(key, i);
    }
    static bool eval(const Ref<KeyVal>& kv, const char* key, int i, int j) {
      return kv->booleanvalue(key, i, j);
    }
  };
  template <> struct GetValue<double> {
    static double eval(const Ref<KeyVal>& kv, const char* key, int i) {
      return kv->doublevalue(key, i);
    }
    static double eval(const Ref<KeyVal>& kv, const char* key, int i, int j) {
      return kv->doublevalue(key, i, j);
    }
  };
  template <> struct GetValue<float> {
    static float eval(const Ref<KeyVal>& kv, const char* key, int i) {
      return kv->floatvalue(key, i);
    }
    static float eval(const Ref<KeyVal>& kv, const char* key, int i, int j) {
      return kv->floatvalue(key, i, j);
    }
  };
  template <> struct GetValue<int> {
    static int eval(const Ref<KeyVal>& kv, const char* key, int i) {
      return kv->intvalue(key, i);
    }
    static int eval(const Ref<KeyVal>& kv, const char* key, int i, int j) {
      return kv->intvalue(key, i, j);
    }
  };
  template <> struct GetValue<long> {
    static long eval(const Ref<KeyVal>& kv, const char* key, int i) {
      return kv->longvalue(key, i);
    }
    static long eval(const Ref<KeyVal>& kv, const char* key, int i, int j) {
      return kv->longvalue(key, i, j);
    }
  };
  template <> struct GetValue<std::size_t> {
    static std::size_t eval(const Ref<KeyVal>& kv, const char* key, int i) {
      return kv->sizevalue(key, i);
    }
    static std::size_t eval(const Ref<KeyVal>& kv, const char* key, int i, int j) {
      return kv->sizevalue(key, i, j);
    }
  };
  template <> struct GetValue<char> {
    static char eval(const Ref<KeyVal>& kv, const char* key, int i) {
      return kv->charvalue(key, i);
    }
    static char eval(const Ref<KeyVal>& kv, const char* key, int i, int j) {
      return kv->charvalue(key, i, j);
    }
  };
  template <> struct GetValue<std::string> {
    static std::string eval(const Ref<KeyVal>& kv, const char* key, int i) {
      return kv->stringvalue(key, i);
    }
    static std::string eval(const Ref<KeyVal>& kv, const char* key, int i, int j) {
      return kv->stringvalue(key, i, j);
    }
  };
}

/** @ingroup CoreKeyVal
 * Provides convenient way to fill standard containers from KeyVal
 */
class Keyword {
  public:
    Keyword(const Ref<KeyVal>& kv,
            const std::string& key) : kv_(kv), key_(key) {}

    const Ref<KeyVal>& keyval() const { return kv_; }
    const std::string& key() const { return key_; }

    /// fills up std::vector using KeyVal object of the following form
    /// key = [ value0 value1 value2 .... ]
    template <typename Value, typename Alloc>
    Keyword& operator>>(std::vector<Value,Alloc>& vec) {
      const std::size_t n = kv_->count(key_.c_str());
      vec.resize(n);
      for(std::size_t i=0; i<n; ++i) {
        vec[i] = detail::GetValue<Value>::eval(kv_, key_.c_str(), i);
      }
      return *this;
    }

    /// fills up std::set using KeyVal object of the following form
    /// key = [ key0 key1 key2 ..... ]
    template <typename Key, typename Compare, typename Alloc>
    Keyword& operator>>(std::set<Key, Compare, Alloc>& c) {
      const std::size_t n = kv_->count(key_.c_str());
      for(std::size_t i=0; i<n; ++i) {
        c.insert( detail::GetValue<Key>::eval(kv_, key_.c_str(), i) );
      }
      return *this;
    }

    /// fills up std::map using KeyVal object of the following form
    /// key = [ [key0 value0] [key1 value1] [key2 value2] ....]
    template <typename Key, typename Data, typename Compare, typename Alloc>
    Keyword& operator>>(std::map<Key, Data, Compare, Alloc>& c) {
      const std::size_t n = kv_->count(key_.c_str());
      for(std::size_t i=0; i<n; ++i) {
        if (kv_->count(key_.c_str(),i) != 2) {
          std::ostringstream oss;
          oss << key_ << ":" << i;
          throw sc::InputError("invalid std::map specification in KeyVal",
                               __FILE__, __LINE__,
                               oss.str().c_str());
        }
        Key k = detail::GetValue<Key>::eval(kv_, key_.c_str(), i, 0);
        Data d = detail::GetValue<Data>::eval(kv_, key_.c_str(), i, 1);
        c.insert(std::make_pair(k,d));
      }
      return *this;
    }

  private:
    Ref<KeyVal> kv_;
    std::string key_;
};

}

#endif /* _KeyVal_h */

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
