
#ifndef _chemistry_qc_integral_integral_h
#define _chemistry_qc_integral_integral_h

#include <math/newmat7/newmat.h>
#include <chemistry/qc/basis/basis.h>

class OneBodyInt;

// OneBodyIntShellIter* obii = ...;
// OneBodyInt* obi = ...;
// for (obii=obi; obii; obii.next()) {
//    for (int i = 0; i<obii.i_len(); i++) {
//      for (int j = 0; i<obii.j_len(); j++) {
//        double integral = obii[i][j];
//        ...
//        }
//      }
//    }
// for a range of shells:
// (obii=obi).set_range(i_start,i_length,j_start,j_length)

class OneBodyIntShellIter
{
 private:
  void init();
 protected:
  int i_;
  int j_;
  int i_len_;
  int j_len_;
  int i_function_;
  int j_function_;
  
  int i_start_;
  int i_end_;
  int j_start_;
  int j_end_;
  
  double* buffer_;
  OneBodyInt* obi;
 public:
  OneBodyIntShellIter(OneBodyInt*);
  
  inline int i() { return i_; }
  inline int j() { return j_; }
  inline int i_len() { return i_len_; }
  inline int j_len() { return j_len_; }
  inline double operator()(int i,int j) { return val_by_shell_bf(i,j); }
  inline double val_by_shell_bf(int i,int j)
  {
    return buffer_[i * j_len_ + j];
  }
  inline double val_by_overall_bf(int i,int j)
  {
    return buffer_[(i-i_function_) * j_len_ + (j-j_function_)];
  }
  inline int i_function() { return i_function_; }
  inline int j_function() { return j_function_; }
  inline int i_function_fence() { return i_function_ + i_len_; }
  inline int j_function_fence() { return j_function_ + j_len_; }
  inline double* values() { return buffer_; }
  void set_range(int,int,int,int);
  
  virtual void start() = 0;
  virtual operator int() = 0;
  
  virtual void next() = 0;
  
  virtual ~OneBodyIntShellIter();
};

class OneBodyIntRedundantShellIter:
  public OneBodyIntShellIter
{
 private:
 public:
  OneBodyIntRedundantShellIter(OneBodyInt*obi);
  ~OneBodyIntRedundantShellIter();
  
  void next();
  void start();
  operator int();
};

class OneBodyIntNonredundantNonrepeatedShellIter:
  public OneBodyIntShellIter
{
 private:
  int diagonal;
 public:
  OneBodyIntNonredundantNonrepeatedShellIter(OneBodyInt*obi);
  ~OneBodyIntNonredundantNonrepeatedShellIter();
  
  void set_range(int,int,int,int);
  void next();
  void start();
  operator int();
};

class OneBodyIntNonredundantRestShellIter:
  public OneBodyIntShellIter
{
 private:
  int diagonal;
 public:
  OneBodyIntNonredundantRestShellIter(OneBodyInt*obi);
  ~OneBodyIntNonredundantRestShellIter();
  
  void set_range(int,int,int,int);
  void next();
  void start();
  operator int();
};

class OneBodyInt
{
 private:
  const GaussianBasisSet* bs;
 public:
  OneBodyInt(const GaussianBasisSet*b);
  
  virtual int nbasis();
  virtual int nshell();
  virtual GaussianBasisSet& basis();
  //virtual double integral(int,int);
  virtual void compute_shell(int,int,double*) = 0;
  virtual void compute(SymmetricMatrix&);
  virtual ~OneBodyInt();
  virtual OneBodyIntShellIter* get_nonredundant_nonrepeated_shell_iter();
  virtual OneBodyIntShellIter* get_nonredundant_rest_shell_iter();
  virtual OneBodyIntShellIter* get_redundant_shell_iter();
};

class TestOneBodyInt:
  public OneBodyInt
{
 private:
 public:
  TestOneBodyInt(const GaussianBasisSet*bs);
  ~TestOneBodyInt();
  void compute_shell(int,int,double*);
};

#ifdef ONE_BODY_SPECIALIZATIONS /* Mon Feb  1 11:07:42 1993 */
/**/class GaussianOverlapInt:
/**/  public OneBodyInt
/**/  {
/**/  private:
/**/  public:
/**/    GaussianOverlapInt(GaussianBasisSet*);
/**/    void compute_shell(int,int,double*);
/**/  };
/**/
/**/class GaussianKineticInt:
/**/  public OneBodyInt
/**/  {
/**/  private:
/**/  public:
/**/    GaussianKineticInt(GaussianBasisSet*);
/**/    void compute_shell(int,int,double*);
/**/  };
/**/
/**/class GaussianPointChargeInt:
/**/  public OneBodyInt
/**/  {
/**/  private:
/**/    PointSet<double>* ps;
/**/  public:
/**/    GaussianPointChargeInt(PointSet<double>*,GaussianBasisSet*);
/**/    GaussianPointChargeInt(Molecule*,GaussianBasisSet*);
/**/    void compute_shell(int,int,double*);
/**/  };
#endif /* ONE_BODY_SPECIALIZATIONS Mon Feb  1 11:07:42 1993 */

#ifdef TWO_BODY /* Mon Feb  1 11:02:30 1993 */
/**/class TwoBodyInt
/**/  {
/**/  private:
/**/    GaussianBasisSet* bs;
/**/  public:
/**/    TwoBodyInt(GaussianBasisSet*);
/**/    int nbasis();
/**/    virtual double integral(int,int,int,int) = 0;
/**/  };
/**/
/**/class GaussianElectronRepulsionInt:
/**/  public TwoBodyInt
/**/  {
/**/  private:
/**/  public:
/**/    GaussianElectronRepulsionInt(GaussianBasisSet*);
/**/  };
#endif /* TWO_BODY Mon Feb  1 11:02:30 1993 */

#endif
