
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_misc_compute_h
#define _util_misc_compute_h

#include <stdio.h>
#include <util/container/set.h>

class Result;

typedef Result* ResultP;

ARRAY_dec(ResultP);
SET_dec(ResultP);

class Compute
{
   friend class Result;
  private:
    SetResultP _results;
    void add(Result*);

    // Prohibit copy
    Compute(const Compute&) {};

  protected:
    // Recompute at least the results that have compute true
    // and are not already computed.  This should only be called
    // by Result's members.
    virtual void compute() = 0;
  public:
    Compute();
    virtual ~Compute();

    // Marks all results as being out of date.  Any subsequent access
    // to results will cause Compute::compute() to be called.
    void obsolete();
};

// Usually Result_dec(Type) will be used to create a result that has
// a particular datum associated with it, however simple Result's can
// also be declared to keep track of datum's for which it is awkward
// to use Result_dec.
class Result
{
  private:
    int _compute;
    int _computed;
    Compute* _c;
  protected:
    // This make sure that the datum is up to date.  If it is not then
    // Compute::compute() will be called.
    void update();
  protected:
    Result(StateIn&,Compute*);
    virtual void save_data_state(StateOut&);
  public:
    Result(Compute*c);
    virtual ~Result();
    int& compute();
    int compute(int c);
    int& computed();
    int needed();
};

// This is like result but the accuracy with which a result was computed
// as well as the desired accuracy are stored.  A computed_ datum always
// has an actual accuracy greater than or equal to the computed accuracy.
class AccResult: public Result
{
  private:
    double _actual_accuracy;
    double _desired_accuracy;
  protected:
    AccResult(StateIn&,Compute*);
    virtual void save_data_state(StateOut&);
  public:
    AccResult(Compute*c);
    ~AccResult();
    double actual_accuracy();
    double desired_accuracy();
    void set_desired_accuracy(double);
    void set_actual_accuracy(double);
};

// This provide access to the result.  Before the result is
// returned Result::update() is called to make sure that the
// datum is up to date.
#define Result_inline_dec(T)						      \
class Result ## T: public Result {					      \
  private:								      \
    T _result;								      \
  public:								      \
    inline Result ## T (Compute*c):Result(c) {};			      \
    inline operator T&() { update(); return _result; };			      \
    inline T* operator ->() { update(); return &_result; };		      \
    inline T& result() { update(); return _result; };			      \
    inline T& result_noupdate() { return _result; };			      \
    inline void operator=(const T& a) { _result = a; };			      \
}

#define _Result_dec(PREFIX,BASE,T,PUBLIC)				      \
class PREFIX ## T: public BASE {					      \
  private:								      \
    T _result;								      \
  public:								      \
    PREFIX ## T (Compute*c);						      \
    operator T&();							      \
    T* operator ->();							      \
    T& result();							      \
    T& result_noupdate();						      \
    void operator=(const T& a);						      \
    PUBLIC								      \
}

#define Result_dec(T) _Result_dec(Result,Result,T,)
#define AccResult_dec(T) _Result_dec(AccResult,AccResult,T,)
#define SSAccResult_dec(T) _Result_dec(AccResult,AccResult,T,\
                                    void save_data_state(StateOut&);\
                                    AccResult ## T (StateIn&,Compute*);)

#define _Result_def(PREFIX,BASE,T,OTHER)				      \
PREFIX ## T:: PREFIX ## T (Compute*c): BASE (c) {};			      \
PREFIX ## T::operator T&() { update(); return _result; };		      \
T* PREFIX ## T::operator ->() { update(); return &_result; };		      \
T& PREFIX ## T::result() { update(); return _result; };			      \
T& PREFIX ## T::result_noupdate() { return _result; };			      \
void PREFIX ## T::operator=(const T& a) { _result = a; };			      \
OTHER

#define Result_def(T) _Result_def(Result,Result,T,)
#define AccResult_def(T) _Result_def(AccResult,AccResult,T,)
#define SSAccResult_def(T) _Result_def(AccResult,AccResult,T,		      \
void AccResult ## T::save_data_state(StateOut&s)			      \
{									      \
  AccResult::save_data_state(s);					      \
  _result.save_data_state(s);						      \
}									      \
AccResult ## T::AccResult ## T(StateIn&s,Compute*c):			      \
  AccResult(s,c) COMPUTE_H_COMMA  _result(s){})

#define COMPUTE_H_COMMA ,

// these are for builtin types like int, double, etc since the ARM says
// operator->() must return a class object
#define Result_inline_dec_nc(T)						      \
class Result ## T: public Result {					      \
  private:								      \
    T _result;								      \
  public:								      \
    inline Result ## T (Compute*c):Result(c) {};			      \
    inline operator T&() { update(); return _result; };			      \
    inline T* pointer() { update(); return &_result; };		              \
    inline T& result() { update(); return _result; };			      \
    inline T& result_noupdate() { return _result; };			      \
    inline void operator=(const T& a) { _result = a; };			      \
}

#define _Result_dec_nc(PREFIX,BASE,T,PUBLIC)				      \
class PREFIX ## T: public BASE {					      \
  private:								      \
    T _result;								      \
  public:								      \
    PREFIX ## T (Compute*c);						      \
    operator T&();							      \
    T* pointer();							      \
    T& result();							      \
    T& result_noupdate();						      \
    void operator=(const T& a);						      \
    PUBLIC								      \
}

#define Result_dec_nc(T) _Result_dec_nc(Result,Result,T,)
#define AccResult_dec_nc(T) _Result_dec_nc(AccResult,AccResult,T,)
#define SSAccResult_dec_nc(T) _Result_dec_nc(AccResult,AccResult,T,\
                                       void save_data_state(StateOut&);\
                                       AccResult ## T (StateIn&,Compute*);)

#define _Result_def_nc(PREFIX,BASE,T,OTHER)				      \
PREFIX ## T::PREFIX ## T (Compute*c): BASE (c) {};			      \
PREFIX ## T::operator T&() { update(); return _result; };		      \
T* PREFIX ## T::pointer() { update(); return &_result; };		      \
T& PREFIX ## T::result() { update(); return _result; };			      \
T& PREFIX ## T::result_noupdate() { return _result; };			      \
void PREFIX ## T::operator=(const T& a) { _result = a; };		      \
OTHER

#define Result_def_nc(T) _Result_def_nc(Result,Result,T,)
#define AccResult_def_nc(T) _Result_def_nc(AccResult,AccResult,T,)
#define SSAccResult_def_nc(T) _Result_def_nc(AccResult,AccResult,T,	      \
void AccResult ## T::save_data_state(StateOut&s)			      \
{									      \
  AccResult::save_data_state(s);					      \
  s.put(_result);							      \
}									      \
AccResult ## T::AccResult ## T(StateIn&s,Compute*c):			      \
  AccResult(s,c){s.get(_result);}					      \
)

#ifdef INLINE_FUNCTIONS
#include <util/misc/compute_i.h>
#undef Result_def
#undef Result_dec
#define Result_def(T)
#define Result_dec(T) Result_inline_dec(T)

#undef Result_def_nc
#undef Result_dec_nc
#define Result_def_nc(T)
#define Result_dec_nc(T) Result_inline_dec_nc(T)
#endif

// Results for some common types
Result_dec_nc(int);
Result_dec_nc(float);
Result_dec_nc(double);
SSAccResult_dec_nc(double);

#endif
