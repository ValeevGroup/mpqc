
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_misc_compute_h
#define _util_misc_compute_h

#include <stdio.h>
#include <util/container/set.h>
#include <util/state/state.h>

class ResultBase;
class StateIn;
class StateOut;

typedef ResultBase* ResultBaseP;

ARRAY_dec(ResultBaseP);
SET_dec(ResultBaseP);

class Compute
{
   friend class ResultBase;
  private:
    SetResultBaseP _results;
    void add(ResultBase*);

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
class ResultBase
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
    ResultBase(StateIn&,Compute*);
    ResultBase(const ResultBase&,Compute*);
    virtual void save_data_state(StateOut&);
    ResultBase& operator=(const ResultBase&);
  public:
    ResultBase(Compute*c);
    virtual ~ResultBase();
    int& compute() { return _compute; }
    int compute(int c) { int r = _compute; _compute = c; return r; }
    int& computed() { return _computed; }
    int needed() { return _compute && (!_computed); }
};

// This is like result but the accuracy with which a result was computed
// as well as the desired accuracy are stored.  A computed_ datum always
// has an actual accuracy greater than or equal to the computed accuracy.
class AccResultBase: public ResultBase
{
  private:
    double _actual_accuracy;
    double _desired_accuracy;
  protected:
    AccResultBase(StateIn&,Compute*);
    AccResultBase(const AccResultBase&,Compute*);
    virtual void save_data_state(StateOut&);
    AccResultBase& operator=(const AccResultBase&);
  public:
    AccResultBase(Compute*c);
    ~AccResultBase();
    double actual_accuracy() const;
    double desired_accuracy() const;
    void set_desired_accuracy(double);
    void set_actual_accuracy(double);
};

#include <util/misc/comptmpl.h>

typedef Result<int> Resultint;
typedef Result<double> Resultdouble;
typedef NCAccResult<double> AccResultdouble;

#endif
