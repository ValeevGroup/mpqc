
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _math_symmetry_pointgrp_h
#define _math_symmetry_pointgrp_h

#include <stdio.h>

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <math/topology/point.h>

////////////////////////////////////////////////////////////////////

//texi
// The @code{SymmetryOperation} class provides a 3 x 3 matrix representation
// of a symmetry operation, such as a rotation or reflection.
class SymmetryOperation {
  private:
    double d[3][3];

  public:
    SymmetryOperation();
    ~SymmetryOperation();

    //texi returns the trace of the transformation matrix
    double trace() const { return d[0][0]+d[1][1]+d[2][2]; }

    //texi returns the i'th row of the transformation matrix
    double* operator[](int i) { return d[i]; }

    //texi const version of the above
    const double* operator[](int i) const { return d[i]; }

    //texi returns a reference to the (i,j)th element of the transformation
    // matrix
    double& operator()(int i, int j) { return d[i][j]; }

    //texi const version of the above
    const double operator()(int i, int j) const { return d[i][j]; }

    //texi zero out the symop
    void zero() { memset(d,0,sizeof(double)*9); }

    //texi This operates on @code{this} with r (i.e. return r * @code{this}).
    SymmetryOperation operate(const SymmetryOperation& r) const;

    //texi This performs the transform r * @code{this} * r~
    SymmetryOperation sim_transform(const SymmetryOperation& r) const;
    
    //texi Set equal to a unit matrix
    void unit() { zero(); d[0][0] = d[1][1] = d[2][2] = 1.0; }

    //texi Set equal to E
    void E() { unit(); }
    
    //texi Set equal to an inversion
    void i() { zero(); d[0][0] = d[1][1] = d[2][2] = -1.0; }

    //texi Set equal to reflection in xy plane
    void sigma_h() { unit(); d[2][2] = -1.0; }

    //texi Set equal to reflection in xz plane
    void sigma_xz() { unit(); d[1][1] = -1.0; }

    //texi Set equal to reflection in yz plane
    void sigma_yz() { unit(); d[0][0] = -1.0; }

    //texi Set equal to a clockwise rotation by 2pi/n
    void rotation(int n);
    void rotation(double theta);
    
    //texi Set equal to C2 about the x axis
    void c2_x() { i(); d[0][0] = 1.0; }

    //texi Set equal to C2 about the x axis
    void c2_y() { i(); d[1][1] = 1.0; }

    //texi print the matrix 
    void print(FILE* =stdout) const;
};

////////////////////////////////////////////////////////////////////

//texi
// The @code{SymRep} class provides an n dimensional matrix representation
// of a symmetry operation, such as a rotation or reflection.  The trace
// of a @code{SymRep} can be used as the character for that symmetry operation.
// d is hardwired to 5x5 since the H irrep in Ih is 5 dimensional.
class SymRep {
  private:
    int n;
    double d[5][5];

  public:
    SymRep(int =0);
    SymRep(const SymmetryOperation&);
    ~SymRep();

    //texi Cast to a SymmetryOperation.
    operator SymmetryOperation() const;
    
    //texi returns the trace of the transformation matrix
    inline double trace() const;

    //texi set the dimension of d
    void set_dim(int i) { n=i; }
    
    //texi returns the i'th row of the transformation matrix
    double* operator[](int i) { return d[i]; }
    //texi const version of the above
    const double* operator[](int i) const { return d[i]; }

    //texi returns a reference to the (i,j)th element of the transformation
    // matrix
    double& operator()(int i, int j) { return d[i][j]; }
    //texi const version of the above
    const double operator()(int i, int j) const { return d[i][j]; }

    //texi zero out the symop
    void zero() { memset(d,0,sizeof(double)*25); }

    //texi This operates on @code{this} with r (i.e. return r * @code{this}).
    SymRep operate(const SymRep& r) const;

    //texi This performs the transform r * @code{this} * r~
    SymRep sim_transform(const SymRep& r) const;
    
    //texi Set equal to a unit matrix
    void unit() {
      zero(); d[0][0] = d[1][1] = d[2][2] = d[3][3] = d[4][4] = 1.0;
    }
    
    //texi Set equal to the identity
    void E() { unit(); }
    
    //texi Set equal to an inversion
    void i() { zero(); d[0][0] = d[1][1] = d[2][2] = d[3][3] = d[4][4] = -1.0;}

    //texi Set equal to reflection in xy plane
    void sigma_h();

    //texi Set equal to reflection in xz plane
    void sigma_xz();

    //texi Set equal to reflection in yz plane
    void sigma_yz();

    //texi Set equal to a clockwise rotation by 2pi/n
    void rotation(int n);
    void rotation(double theta);
    
    //texi Set equal to C2 about the x axis
    void c2_x();

    //texi Set equal to C2 about the x axis
    void c2_y();

    //texi print the matrix 
    void print(FILE* =stdout) const;
};

inline double
SymRep::trace() const
{
  double r=0;
  for (int i=0; i < n; i++)
    r += d[i][i];
  return r;
}

////////////////////////////////////////////////////////////////////


class CharacterTable;

//texi
// The @code{IrreducibleRepresentation} class provides information associated
// with a particular irreducible representation of a point group.  This
// includes the Mulliken symbol for the irrep, the degeneracy of the irrep,
// the characters which represent the irrep, and the number of translations
// and rotations in the irrep.  The order of the point group is also provided
// (this is equal to the number of characters in an irrep). 
class IrreducibleRepresentation {
  friend class CharacterTable;

  private:
    int g;         // the order of the group
    int degen;     // the degeneracy of the irrep
    int nrot_;     // the number of rotations in this irrep
    int ntrans_;   // the number of translations in this irrep
    int complex_;  // true if this irrep has a complex representation
    char *symb;    // mulliken symbol for this irrep

    SymRep *rep;   // representation matrices for the symops

  public:
    IrreducibleRepresentation();
    IrreducibleRepresentation(const IrreducibleRepresentation&);
    //texi This constructor takes as arguments the order of the point group,
    // the degeneracy of the irrep, and the Mulliken symbol of the irrep.
    // The Mulliken symbol is copied internally.
    IrreducibleRepresentation(int,int,const char*);

    ~IrreducibleRepresentation();

    IrreducibleRepresentation& operator=(const IrreducibleRepresentation&);

    //texi Initialize the order, degeneracy, and Mulliken symbol of the irrep.
    void init(int =0, int =0, const char* =0);
    
    //texi Returns the order of the group.
    int order() const { return g; }

    //texi Returns the degeneracy of the irrep.
    int degeneracy() const { return degen; }

    //texi Returns the number of projection operators for the irrep.
    int nproj() const { return degen*degen; }

    //texi Returns the number of rotations associated with the irrep.
    int nrot() const { return nrot_; }

    //texi Returns the number of translations associated with the irrep.
    int ntrans() const { return ntrans_; }

    //texi Returns the Mulliken symbol for the irrep.
    const char * symbol() const { return symb; }

    //texi
    // Returns the character for the i'th symmetry operation of the point
    // group.
    double character(int i) const {
      return complex_ ? 0.5*rep[i].trace() : rep[i].trace();
    }

    //texi This is equivalent to the @b{character} member.
    double operator[](int i) const { return character(i); }

    //texi
    // Returns the character for the d'th contribution to the i'th symmetry
    // operation of the point group.
    double p(int d, int i) const {
      int dc=d/degen; int dr=d%degen;
      return rep[i](dr,dc);
    }

    double operator()(int d, int i) const { return p(d,i); }

    //texi
    // This prints the irrep to the given file, or stdout if none is given.
    // The second argument is an optional string of spaces to offset by.
    void print(FILE* =stdout, const char * =" ");
};

/////////////////////////////////////////////////////////////
//texi
// The @code{CharacterTable} class provides a workable character table for
// all of the non-cubic point groups.  While I have tried to match the
// ordering in Cotton's book, I don't guarantee that it is always followed.
// It shouldn't matter anyway.  Also note that I don't lump symmetry operations
// of the same class together.  For example, in C3v there are two distinct
// C3 rotations and 3 distinct reflections, each with a separate character.
// Thus symop has 6 elements rather than the 3 you'll find in most published
// character tables.
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
    SymmetryOperation *symop;            // the matrices describing sym ops
    char *symb;                          // the Schoenflies symbol for the pg

    //texi this determines what type of point group we're dealing with
    int parse_symbol();
    //texi this fills in the irrep and symop arrays.
    int make_table();

    // these create the character tables for the cubic groups
    void t();
    void th();
    void td();
    void o();
    void oh();
    void i();
    void ih();

  public:
    CharacterTable();
    //texi This constructor takes the Schoenflies symbol of a point group as
    // input.
    CharacterTable(const char*);
    //texi This is like the above, but it also takes a reference to a
    // @code{SymmetryOperation} which is the frame of reference.  All symmetry
    // operations are transformed to this frame of reference.
    CharacterTable(const char*,const SymmetryOperation&);
    CharacterTable(const CharacterTable&);
    ~CharacterTable();

    CharacterTable& operator=(const CharacterTable&);

    //texi Returns the number of irreps.
    int nirrep() const { return nirrep_; }
    //texi Returns the order of the point group
    int order() const { return g; }
    //texi Returns the Schoenflies symbol for the point group
    const char * symbol() const { return symb; }
    //texi Returns the i'th irrep.
    IrreducibleRepresentation& gamma(int i) { return gamma_[i]; }
    //texi Returns the i'th symmetry operation.
    SymmetryOperation& symm_operation(int i) { return symop[i]; }

    //texi Shorthand for @code{gamma}.
    IrreducibleRepresentation& operator[](int i) { return gamma_[i]; }

    //texi Cn, Cnh, and Sn point groups have complex representations.
    // This function returns 1 if the point group has a complex representation,
    // 0 otherwise.
    int complex() const { if(pg==CN || pg==SN || pg==CNH) return 1; return 0; }

    //texi 
    // This prints the irrep to the given file, or stdout if none is given.
    // The second argument is an optional string of spaces to offset by.
    void print(FILE* =stdout, const char * =" ");
};

/////////////////////////////////////////////////////////////

//texi
// The @code{PointGroup} class is really a place holder for a
// @code{CharacterTable}.  It contains a string representation of the
// Schoenflies symbol of a point group, a frame of reference for the
// symmetry operation transformation matrices, and a point of origin.
// The origin may one day disappear...I once thought it would be useful.
// Maybe I'll remember why some day.
//
// @code{PointGroup} is the only class in libsymmetry which is a
// @code{SavableState}.  I did this to save space...it takes less than a
// second to generate a character table from the Schoenflies symbol, but
// the character table for Ih will occupy a great deal of memory.  Best to
// free up that memory when it's not needed.
class PointGroup: public SavableState {
#   define CLASSNAME PointGroup
#   define HAVE_CTOR
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    char *symb;
    SymmetryOperation frame;
    Point origin_;

  public:
    PointGroup();
    //texi This constructor takes a string containing the Schoenflies symbol
    // of the point group as its only argument.
    PointGroup(const char*);
    //texi Like the above, but this constructor also takes a frame of reference
    // as an argument.
    PointGroup(const char*,SymmetryOperation&);
    //texi Like the above, but this constructor also takes a point of origin
    // as an argument.
    PointGroup(const char*,SymmetryOperation&,Point&);
    //texi The KeyVal constructor (@ref{The PointGroup KeyVal Constructor}).
    PointGroup(const RefKeyVal&);

    PointGroup(StateIn&);
    PointGroup(PointGroup&);
    ~PointGroup();

    PointGroup& operator=(PointGroup&);

    //texi Returns the @code{CharacterTable} for this point group.
    CharacterTable char_table() const;
    //texi Returns the Schoenflies symbol for this point group.
    const char * symbol() const { return symb; }
    //texi Returns the frame of reference for this point group.
    SymmetryOperation& symm_frame() { return frame; }
    //texi A const version of the above
    const SymmetryOperation& symm_frame() const { return frame; }
    //texi Returns the origin of the symmetry frame.
    Point& origin() { return origin_; }
    //texi A const version of the above.
    const Point& origin() const { return origin_; }

    //texi Sets (or resets) the Schoenflies symbol.
    void set_symbol(const char*);

    void save_data_state(StateOut& so);
};

#endif
