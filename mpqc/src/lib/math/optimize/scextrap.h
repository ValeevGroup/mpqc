
#ifndef _math_optimize_scextrap_h
#define _math_optimize_scextrap_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/class/class.h>
#include <util/keyval/keyval.h>

DescribedClass_REF_fwddec(SCExtrapData);
class SCExtrapData: public DescribedClass {
#   define CLASSNAME SCExtrapData
#   include <util/class/classda.h>
  public:
    ~SCExtrapData();
    virtual SCExtrapData* copy() = 0;
    virtual void zero() = 0;
    virtual void accumulate_scaled(double scale, const RefSCExtrapData&) = 0;
};
DescribedClass_REF_dec(SCExtrapData);

DescribedClass_REF_fwddec(SCExtrapError);
class SCExtrapError: public DescribedClass {
#   define CLASSNAME SCExtrapError
#   include <util/class/classda.h>
  public:
    ~SCExtrapError();
    virtual double error() = 0;
    virtual double scalar_product(const RefSCExtrapError&) = 0;
};
DescribedClass_REF_dec(SCExtrapError);

class SelfConsistentExtrapolation: public DescribedClass {
#   define CLASSNAME SelfConsistentExtrapolation
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classda.h>
  private:
    double error_;
    int errorset_;
    double tolerance_;
  protected:
    void set_error(double e) { error_ = e; errorset_ = 1; }
  public:
    SelfConsistentExtrapolation();
    SelfConsistentExtrapolation(const RefKeyVal&);
    ~SelfConsistentExtrapolation();

    void set_tolerance(double t) { tolerance_ = t; }
    double tolerance() { return tolerance_; }
    double error() { return error_; }

    int converged() { return errorset_? error_ <= tolerance_ : 0; }

    // Makes a copy of data and returns the extrapolation in
    // data.  A reference to error is saved so a copy must
    // be given to extrapolate if error could be changed.
    virtual int extrapolate(const RefSCExtrapData& data,
                            const RefSCExtrapError& error) = 0;
};
DescribedClass_REF_dec(SelfConsistentExtrapolation);

#endif
