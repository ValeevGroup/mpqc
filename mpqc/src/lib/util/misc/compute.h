
#ifndef _util_container_result_h
#define _util_container_result_h

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
    inline Result(Compute*c):_c(c) { c->add(this); };
    inline int& compute() { return _compute; };
    inline int compute(int c) { int r = _compute; _compute = c; return r; };
    inline int& computed() { return _computed; };
};

// This provide access to the result.  Before the result is
// returned Result::update() is called to make sure that the
// datum is up to date.
#define Result_dec(T)							      \
class Result ## T: public Result {					      \
  private:								      \
    T _result;								      \
  public:								      \
    inline Result ## T (Compute*c):Result(c) {};			      \
    inline operator T&() { update(); return _result; };			      \
    inline T* operator ->() { update(); return &_result; };		      \
    inline T& result() { update(); return _result; };			      \
    inline T& result_noupdate() { return _result; };			      \
}

// Results for some common type
Result_dec(int);
Result_dec(float);
Result_dec(double);

#endif
