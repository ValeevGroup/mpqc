
/* pointgrp.h -- definition of the point group classes
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
 *      June, 1993
 */

#ifndef _math_symmetry_pointgrp_h
#define _math_symmetry_pointgrp_h

#include <stdio.h>
#include <iostream.h>

#include <math/nihmatrix/nihmatrix.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>


class CharacterTable;

class IrreducibleRepresentation {
  friend class CharacterTable;

  private:
    int g;        // the order of the group
    int degen;    // the degeneracy of the irrep
    int nrot_;    // the number of rotations in this irrep
    int ntrans_;  // the number of translations in this irrep
    char *symb;   // mulliken symbol for this irrep
    double *rep;  // the characters for this irrep

    void init();

  public:
    IrreducibleRepresentation();
    IrreducibleRepresentation(int,int,const char*);
    IrreducibleRepresentation(const IrreducibleRepresentation&);
    ~IrreducibleRepresentation();

    IrreducibleRepresentation& operator=(const IrreducibleRepresentation&);

    inline int order() const { return g; }
    inline int degeneracy() const { return degen; }
    inline int nrot() const { return nrot_; }
    inline int ntrans() const { return ntrans_; }
    inline const char * symbol() const { return symb; }
    inline double character(int i) const { return rep[i]; }

    inline double operator[](int i) const { return rep[i]; }

    void print(FILE* =stdout, const char * =" ");
};

/////////////////////////////////////////////////////////////

class CharacterTable {
  public:
    enum pgroups {C1, CS, CI, CN, CNV, CNH, DN, DND, DNH, SN, T, TH, TD, O,
                  OH, I, IH};

  private:
    int g;                               // the order of the point group
    int nt;                              // order of the princ rot axis
    pgroups pg;                          // the class of the point group
    int nirrep_;                         // the number of irreps in this pg
    IrreducibleRepresentation *gamma_;   // an array of irreps
    DMatrix *symop;                      // the matrices describing sym ops
    char *symb;                          // the Schoenflies symbol for the pg

    int parse_symbol();
    int make_table();

  public:
    CharacterTable();
    CharacterTable(const char*);
    CharacterTable(const CharacterTable&);
    ~CharacterTable();

    CharacterTable& operator=(const CharacterTable&);

    inline int nirrep() const { return nirrep_; }
    inline int order() const { return g; }
    inline const char * symbol() const { return symb; }
    inline IrreducibleRepresentation& gamma(int i) { return gamma_[i]; }
    inline DMatrix& symm_operation(int i) { return symop[i]; }

    inline IrreducibleRepresentation& operator[](int i) { return gamma_[i]; }

    void print(FILE* =stdout, const char * =" ");
};

/////////////////////////////////////////////////////////////

class PointGroup
  : virtual public DescribedClass, virtual public SavableState {
#   define CLASSNAME PointGroup
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    char *symb;

  public:
    PointGroup();
    PointGroup(const char*);
    PointGroup(KeyVal&);
    PointGroup(StateIn&);
    PointGroup(PointGroup&);
    ~PointGroup();

    PointGroup& operator=(PointGroup&);
    
    inline CharacterTable char_table() const { return CharacterTable(symb); }
    inline const char * symbol() const { return symb; }

    void save_data_state(StateOut& so);
    void restore_data_state(int, StateIn& si);
};

#endif
