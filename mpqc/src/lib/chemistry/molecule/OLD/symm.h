
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

#define V_BASE virtual public DescribedClass, virtual public SavableState

//////////////////////////////////////////////////////////////////////////

class SymmCoList;
class SymmCoPtr;
class SimpleCoList;
class SimpleCoPtr;

class SymmCCount {
  friend class SymmCoPtr;

  private:
    int count;
  protected:
    SymmCCount(): count(0) {}
  };

//////////////////////////////////////////////////////////////////////////

class SymmCo : public SymmCCount, V_BASE {
DescribedClass_DECLARE(SymmCo)
SavableState_DECLARE(SymmCo)
  private:
    char *label_;
    int n_simple;
    SimpleCoPtr *simple;
    double *coeff_;
    double val;

  public:
    SymmCo();
    SymmCo(const char*,int);
    SymmCo(const char*,SimpleCoPtr&);
    SymmCo(SimpleCoList*,KeyVal*,const char*,int=0);
    SymmCo(const SymmCo&);
    ~SymmCo();

    void init();

    inline const char * label() const { return label_; }
    inline int nsimple() const { return n_simple; }
    //inline SimpleCoList * simples() { return simple; }
    inline SimpleCoPtr& get_simple(int i) { return simple[i]; }
    inline double value() const { return val; }

    virtual SymmCo& operator=(const SymmCo&);
    int operator==(SymmCo&);
    inline int operator!=(SymmCo&u) { return !(*this==u); }

    void calc_intco(Molecule&,double* =0);
    void normalize();

    void print(ostream&, const char* =" ") const;
    void print(FILE* =stdout, const char* =" ") const;

    void save_data_state(StateOut&);
    void restore_data_state(int,StateIn&);

    friend SymmCoList * Geom_form_symm(Molecule&,SimpleCoList* =0,int =1);
    friend SymmCoList * Geom_form_symm(Molecule&,SymmCoList*,int =1);

    friend DMatrix Geom_form_hessian(Molecule&,SimpleCoList* =0);
    friend DMatrix Geom_form_hessian(Molecule&,SimpleCoList*,SymmCoList*);

    friend DMatrix Geom_form_K(Molecule&,SimpleCoList* =0,int =1);
    friend DMatrix Geom_form_K(Molecule&,SymmCoList*,int =1);

    friend ostream& operator<<(ostream&,SymmCo&);
  };

/////////////////////////////////////////////////////////////////////////

SymmCoList * Geom_read_symm(KeyVal*,const char*,SimpleCoList* =0);

void Geom_add_symm(KeyVal*,const char*,SymmCoList*,SimpleCoList* =0);
DMatrix Geom_make_bmat(SymmCoList*,Molecule&);
SymmCoList * Geom_symm_from_simple(SimpleCoList*,int=0);


void Geom_calc_simples(SymmCoList*,Molecule&);
void Geom_normalize(SymmCoList*);

void Geom_print_pretty(SymmCoList*);
void Geom_print_pretty(ostream&,SymmCoList*);

#undef V_BASE

#endif /* _intco_symm_h */
