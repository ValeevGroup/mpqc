
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _math_scmat_dim_h
#define _math_scmat_dim_h

#include <util/keyval/keyval.h>
#include <util/state/state.h>

class RefSCDimension;
class SCBlockInfo: public SavableState {
#   define CLASSNAME SCBlockInfo
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int n_;
    int nblocks_;
    int *start_;
    int *size_;
    RefSCDimension *subdims_;
    void init_start();
  public:
    SCBlockInfo(int n_, int nblocks = 0, const int *blocksizes_ = 0);
    SCBlockInfo(StateIn&);
    SCBlockInfo(const RefKeyVal& keyval);
    ~SCBlockInfo();
    void save_data_state(StateOut&);
    int equiv(SCBlockInfo *);
    int nelem() const { return n_; }
    int nblock() const { return nblocks_; }
    int start(int i) const { return start_[i]; }
    int size(int i) const { return size_[i]; }
    int fence(int i) const { return start_[i] + size_[i]; }
    void elem_to_block(int i, int &block, int &offset);
    RefSCDimension subdim(int i);
    void set_subdim(int i, const RefSCDimension &dim);

    void print(ostream&o=cout);
};
SavableState_REF_dec(SCBlockInfo);

//texi The @code{SCDimension} class is used to determine the size and
// blocking of matrices.
class SCDimension: public SavableState {
#   define CLASSNAME SCDimension
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    char *name_;
    int n_;
    RefSCBlockInfo blocks_;
    SCDimension(const char* name = 0);
  public:
    //texi Create a dimension with an optional name.  The
    // name is a copy of the @code{'\0'} terminated string @var{name}.
    SCDimension(int n, const char* name = 0);
    SCDimension(const RefSCBlockInfo&, const char *name = 0);
    SCDimension(int n, int nblocks, const int *blocksizes = 0,
                const char* name = 0);
    SCDimension(const RefKeyVal&);
    SCDimension(StateIn&s);
    ~SCDimension();
    void save_data_state(StateOut&);

    //texi Test to see if two dimensions are equivalent.
    int equiv(const SCDimension*) const;
    
    //texi Return the dimension.
    int n() const { return n_; }
    //texi Return the name of the dimension.  If no name was given
    // to the constructor, then return @code{0}.
    const char* name() const { return name_; }

    RefSCBlockInfo blocks() { return blocks_; }

    void print(ostream&o=cout);
};

DCRef_declare(SCDimension);
SSRef_declare(SCDimension);

//texi
//  The @code{RefSCDimension} class is a smart pointer to an
//   @code{SCDimension} specialization.
class RefSCDimension: public SSRefSCDimension {
    // standard overrides
  public:
    //texi Initializes the dimension pointer to @code{0}.  The
    // reference must be initialized before it is used.
    RefSCDimension();
    //texi Make this and @var{d} refer to the same @code{SCDimension}.
    RefSCDimension(const RefSCDimension& d);
    //texi Make this refer to @var{d}.
    RefSCDimension(SCDimension *d);

    RefSCDimension(const DCRefBase&);
    ~RefSCDimension();
    //texi Make this refer to @var{d}.
    RefSCDimension& operator=(SCDimension* d);

    RefSCDimension& operator=(const DCRefBase & c);
    //texi Make this and @var{d} refer to the same @code{SCDimension}.
    RefSCDimension& operator=(const RefSCDimension & d);

    // dimension specific functions
  public:
    //texi Return the dimension.
    operator int() const;
    int n() const;

    void print(ostream&o=cout);
};

#endif
