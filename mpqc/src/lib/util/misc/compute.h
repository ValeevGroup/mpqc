
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
  public:
    Result(Compute*c);
    int& compute();
    int compute(int c);
    int& computed();
    int needed();
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
    inline void operator=(T& a) { _result = a; };			      \
}

#define Result_dec(T)							      \
class Result ## T: public Result {					      \
  private:								      \
    T _result;								      \
  public:								      \
    Result ## T (Compute*c);						      \
    operator T&();							      \
    T* operator ->();							      \
    T& result();							      \
    T& result_noupdate();						      \
    void operator=(T& a);						      \
}

#define Result_def(T)							      \
Result ## T::Result ## T (Compute*c):Result(c) {};			      \
Result ## T::operator T&() { update(); return _result; };		      \
T* Result ## T::operator ->() { update(); return &_result; };		      \
T& Result ## T::result() { update(); return _result; };			      \
T& Result ## T::result_noupdate() { return _result; };			      \
void Result ## T::operator=(T& a) { _result = a; };

#ifdef INLINE_FUNCTIONS
#include <util/misc/compute_i.h>
#undef Result_def
#undef Result_dec
#define Result_def(T)
#define Result_dec(T) Result_inline_dec(T)
#endif

// Results for some common types
Result_dec(int);
Result_dec(float);
Result_dec(double);

#endif
