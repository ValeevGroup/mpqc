
/* symm.h -- definition of the symmetrized internal coordinate class
 *
 *      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
 *      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
 *      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
 *      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
 *      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
 *      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
 *
 *  Author:
 *      E. T. Seidl
 *      Bldg. 12A, Rm. 2033
 *      Computer Systems Laboratory
 *      Division of Computer Research and Technology
 *      National Institutes of Health
 *      Bethesda, Maryland 20892
 *      Internet: seidl@alw.nih.gov
 *      February, 1993
 */

#ifndef _libQC_symm_h
#define _libQC_symm_h

#include <stdio.h>
#include <iostream.h>

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <math/nihmatrix/nihmatrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/simple.h>
#include <chemistry/molecule/simpleQCList.h>

#define V_BASE virtual public SavableState

//////////////////////////////////////////////////////////////////////////

class SymmCoList;
class SymmCoPtr;
class SimpleCoList;
class RefSimpleCo;

// class SymmCCount {
//   friend class SymmCoPtr;
// 
//   private:
//     int count;
//   protected:
//     SymmCCount(): count(0) {}
//   };

//////////////////////////////////////////////////////////////////////////

class RefSimpleCoList;
class RefSymmCoList;

class SymmCo : V_BASE {
#   define CLASSNAME SymmCo
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    char *label_;
    int n_simple;
    RefSimpleCo *simple;
    double *coeff_;
    double val;

  public:
    SymmCo();
    SymmCo(const char*,int);
    SymmCo(const char*,RefSimpleCo);
    SymmCo(KeyVal&);
    SymmCo(RefSimpleCoList,RefKeyVal,const char*,int=0);
    SymmCo(const SymmCo&);
    ~SymmCo();

    void init();

    inline const char * label() const { return label_; }
    inline int nsimple() const { return n_simple; }
    //inline SimpleCoList * simples() { return simple; }
    inline RefSimpleCo& get_simple(int i) { return simple[i]; }
    inline double value() const { return val; }

    virtual SymmCo& operator=(const SymmCo&);
    int operator==(SymmCo&);
    inline int operator!=(SymmCo&u) { return !(*this==u); }

    void calc_intco(Molecule&,double* =0);
    void normalize();

    void print(ostream&, const char* =" ") const;
    void print(FILE* =stdout, const char* =" ") const;

    void save_data_state(StateOut&);
    SymmCo(StateIn&);

    friend RefSymmCoList Geom_form_symm(Molecule&,RefSimpleCoList =0,int =1,
                                       RefSymmCoList fixed=0);
    friend RefSymmCoList Geom_form_symm(Molecule&,RefSymmCoList,int =1,
                                       RefSymmCoList fixed=0);

    friend DMatrix Geom_form_hessian(Molecule&,RefSimpleCoList =0);
    friend DMatrix Geom_form_hessian(Molecule&,RefSimpleCoList,RefSymmCoList);

    friend DMatrix Geom_form_K(Molecule&,RefSimpleCoList =0,int =1);
    friend DMatrix Geom_form_K(Molecule&,RefSymmCoList,int =1,
                               RefSymmCoList fixed=0,DMatrix*Kfixed=0);

    friend ostream& operator<<(ostream&,SymmCo&);
  };
DescribedClass_REF_dec(SymmCo);

/////////////////////////////////////////////////////////////////////////

RefSymmCoList Geom_read_symm(RefKeyVal,const char*,RefSimpleCoList =0);

void Geom_add_symm(RefKeyVal,const char*,RefSymmCoList,RefSimpleCoList =0);
DMatrix Geom_make_bmat(RefSymmCoList,Molecule&);
RefSymmCoList Geom_symm_from_simple(RefSimpleCoList,int=0);


void Geom_calc_simples(RefSymmCoList,Molecule&);
void Geom_normalize(RefSymmCoList);

void Geom_print_pretty(RefSymmCoList);
void Geom_print_pretty(ostream&,RefSymmCoList);

#undef V_BASE

#endif /* _intco_symm_h */
