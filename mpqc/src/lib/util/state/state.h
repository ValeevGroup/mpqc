
#ifndef _libqc_state_h
#define _libqc_state_h

#include <stdio.h>
#include <stdlib.h>

#include <util/class/class.h>
#include <util/state/qc_xdr.h>
#include <util/container/array.h>

class StateIn;
class StateOut;

class SavableState: virtual public DescribedClass {
#   define CLASSNAME SavableState
#   include <util/class/classda.h>
  protected:
    SavableState();
  public:
    virtual ~SavableState();

    //// save functions

    // Save an object with unkown exact type.  Duplicated save_states
    // for an object result in a single object being restored.
    void save_state(StateOut&);

    // Save an object with known exact type.
    // This is non virtual and checks to make sure that
    // the exact type == object type.  Each derived class should provide
    // its own version.
    void save_object_state(StateOut&);

    // Save the virtual bases for the object.
    // This must be done in the same order that the ctor
    // initializes the virtual bases.  This does not include
    // the DescribedClass and SavableState virtual base classes.
    // This must be implemented by the user if the class has other virtual
    // bases.  (These virtual bases must come after SavableState.)
    virtual void save_vbase_state(StateOut&);

    // Save the base classes (with save_data_state) and the members
    // in the same order that the StateIn CTOR initializes them.
    // This must be implemented by the derived class if the class has
    // data.
    virtual void save_data_state(StateOut&);

    //// restore functions

    // This restores objects with unknown exact type.
    static SavableState* restore_state(StateIn&);

    // Each derived class version of this CTOR handles the restore
    // corresponding to save_object_state, save_vbase_state, and
    // save_data_state listed above.  All derived class StateIn&
    // constructors must invoke the SavableState(StateIn&) constructor.
  protected:
    SavableState(StateIn&,const ClassDesc&);
  };
DescribedClass_REF_dec(SavableState);

////////////////////////////////////////////////////////////////////

class StateDataPtrSet;
class StateDataNumSet;
class ClassDescPintMap;

class StateOut: virtual public DescribedClass {
#   define CLASSNAME StateOut
#   include <util/class/classda.h>
  protected:
    StateDataPtrSet* ps_;
    ClassDescPintMap* _classidmap;
    int _nextclassid;
    virtual int put_array_void(const void*,int);
    virtual void putparents(const ClassDesc*);
  public:
    StateOut();
    virtual ~StateOut();
    virtual int putpointer(void*);
    virtual int putstring(char*);
    virtual int put_version(const ClassDesc*);
    virtual int put(const ClassDesc*);
    virtual int put(char r);
    virtual int put(int r);
    virtual int put(float r);
    virtual int put(double r);
    virtual int put(char*,int);
    virtual int put(int*,int);
    virtual int put(float*,int);
    virtual int put(double*,int);
    virtual int put_array_char(const char*p,int size);
    virtual int put_array_int(const int*p,int size);
    virtual int put_array_float(const float*p,int size);
    virtual int put_array_double(const double*p,int size);
    //virtual int put(SavableState&ss);
    //virtual int put(SavableState*);
  };
DescribedClass_REF_dec(StateOut);

class StateIn: virtual public DescribedClass {
#   define CLASSNAME StateIn
#   include <util/class/classda.h>
  protected:
    StateDataNumSet* ps_;
    int _nextobject;
    ArraysetClassDescP _cd;
    Arraysetint _version;
    virtual int get_array_void(void*,int);
  public:
    StateIn();
    virtual ~StateIn();
    virtual int  getpointer(void**);
    virtual void havepointer(int,void*);

    // a call to nextobject followed by havepointer(int) is equiv
    // to havepointer(int,void**);
    virtual void nextobject(int);
    virtual void havepointer(void*);
    
    virtual int getstring(char*&);
    virtual int get_version(const ClassDesc*);
    virtual int get(const ClassDesc**);
    virtual int get(char&r);
    virtual int get(int&r);
    virtual int get(float&r);
    virtual int get(double&r);
    virtual int get(char*&);
    virtual int get(int*&);
    virtual int get(float*&);
    virtual int get(double*&);
    virtual int get_array_char(char*p,int size);
    virtual int get_array_int(int*p,int size);
    virtual int get_array_float(float*p,int size);
    virtual int get_array_double(double*p,int size);
    //virtual int get(SavableState&ss);
  };
DescribedClass_REF_dec(StateIn);

////////////////////////////////////////////////////////////////////

class StateOutFile: virtual public StateOut {
  protected:
    int opened_;
    FILE* fp_;
  public:
    StateOutFile();
    StateOutFile(FILE*);
    StateOutFile(const char *, const char * = "w");
    ~StateOutFile();

    virtual void flush();
    virtual void close();
    virtual void rewind();
  };

class StateInFile: virtual public StateIn {
  protected:
    int opened_;
    FILE* fp_;
  public:
    StateInFile();
    StateInFile(FILE*);
    StateInFile(const char *, const char * ="r");
    ~StateInFile();

    virtual void flush();
    virtual void close();
    virtual void rewind();
  };

////////////////////////////////////////////////////////////////////

class StateOutText: public StateOutFile {
  protected:
    void newline();
    void comment(const char*,...);
    void start_array();
    void end_array();
    int putpointer(void*);
    void putparents(const ClassDesc*);
  public:
    StateOutText();
    StateOutText(FILE*);
    StateOutText(const char *, const char * = "w");
    ~StateOutText();
    int putstring(char*);
    int put_array_char(const char*,int);
    int put_array_int(const int*,int);
    int put_array_float(const float*,int);
    int put_array_double(const double*,int);
    int put(const ClassDesc*);
    int put(char r);
    int put(int r);
    int put(float r);
    int put(double r);
    int put(char*,int);
    int put(int*,int);
    int put(float*,int);
    int put(double*,int);
  };

class StateInText: public StateInFile {
  protected:
    int _newlines;
    
    void read(char*);
    void read(int&);
    void read(float&);
    void read(double&);
    void newline();
    void comment();
    void start_array();
    void end_array();
    int  getpointer(void**);

    void abort();
  public:
    StateInText();
    StateInText(FILE*);
    StateInText(const char *, const char * = "r");
    ~StateInText();
    int getstring(char*&);
    int get_array_char(char*,int);
    int get_array_int(int*,int);
    int get_array_float(float*,int);
    int get_array_double(double*,int);
    int get(const ClassDesc**);
    int get(char&r);
    int get(int&r);
    int get(float&r);
    int get(double&r);
    int get(char*&);
    int get(int*&);
    int get(float*&);
    int get(double*&);
  };

////////////////////////////////////////////////////////////////////

class StateOutBin: public StateOutFile {
  protected:
    int put_array_void(const void*,int);
  public:
    StateOutBin();
    StateOutBin(FILE*);
    StateOutBin(const char *, const char * = "w");
    ~StateOutBin();
  };

class StateInBin: public StateInFile {
  protected:
    int get_array_void(void*,int);
  public:
    StateInBin();
    StateInBin(FILE*);
    StateInBin(const char *, const char * = "r");
    ~StateInBin();
  };

////////////////////////////////////////////////////////////////////

class StateOutXDR : virtual public StateOut, public QCXDR {
  protected:
  public:
    StateOutXDR();
    ~StateOutXDR();
    int put_array_char(const char*,int);
    int put_array_int(const int*,int);
    int put_array_float(const float*,int);
    int put_array_double(const double*,int);
  };

class StateOutBinXDR : public StateOutBin,
                       public StateOutXDR
{
  protected:
    //this is needed for a mips-sgi-irix4 gcc 2.5.2 bug
    inline int put_array_void(const void*v,int i)
    { StateOutBin::put_array_void(v,i); }
  public:
    StateOutBinXDR();
    StateOutBinXDR(FILE*);
    StateOutBinXDR(const char *, const char * = "w");
    ~StateOutBinXDR();
};

class StateInXDR : virtual public StateIn, public QCXDR {
  protected:
  public:
    StateInXDR();
    ~StateInXDR();
    int get_array_char(char*,int);
    int get_array_int(int*,int);
    int get_array_float(float*,int);
    int get_array_double(double*,int);
  };

class StateInBinXDR : public StateInBin,
                      public StateInXDR
{
  protected:
    //this is needed for a mips-sgi-irix4 gcc 2.5.2 bug
    inline int get_array_void(void*v,int i)
    { StateInBin::get_array_void(v,i); }
  public:
    StateInBinXDR();
    StateInBinXDR(FILE*);
    StateInBinXDR(const char *, const char * = "r");
    ~StateInBinXDR();
  };

#endif
