
#ifndef _math_scmat_matrix_h
#define _math_scmat_matrix_h

// the mother of all matrices
class SCMatrix: virtual public SavableState {
  public:
    SCMatrix();
    virtual ~SCMatrix();

    virtual int nrow() = 0;
    virtual int ncol() = 0;

    /////////////////////////////////////////////////////////////////
    // access to elements
    // read only access to elements
    virtual const double& element(int i,int j) const = 0;
    // read/write access to elements
    virtual double& element(int i,int j) = 0;
    // other ways to access the elements
    inline const double& operator()(int i,int j) const {return element(i,j);};
    inline double& operator()(int i,int j) {return element(i,j);};

    /////////////////////////////////////////////////////////////////
    // access with generic block iterators
    virtual void first(SCMatrixBlock&);
    virtual void next(SCMatrixBlock&) const;
    // read only access
    virtual void first_ro(SCMatrixBlock&) const;

    /////////////////////////////////////////////////////////////////
    // access with rectangular and special block iterators
    virtual void first(SCMatrixRectBlock&);
    virtual void next(SCMatrixRectBlock&);
    virtual void first(SCMatrixSpecBlock&);
    virtual void next(SCMatrixSpecBlock&);
};
DescribedClass_REF_dec(SCMatrix);

class LocalSCMatrix: public SCMatrix {
  private:
    int _dim;
    double* _data;
  public:
    /////////////////////////////////////////////////////////////////
    // functions for only this class
    LocalSCMatrix(int dim=0);
    ~LocalSCMatrix();
    int dim();
    int resize(int dim);
    const double& canonical(int i,int j) const;
    double& canonical(int i,int j);

    /////////////////////////////////////////////////////////////////
    // overriding functions for SCMatrix
    int nrow();
    int ncol();

    const double& element(int i,int j) const;
    double& element(int i,int j);
};

#endif
