
#ifndef _util_misc_identity_h
#define _util_misc_identity_h

#ifdef __GNUG__
#pragma interface
#endif

// Identifier's are used to distinguish and order objects.  On
// many architectures a pointer to the object will suffice, but
// the C++ standard only guarantees that this works for two pointers
// pointing within the same structure or array.  Classes need to
// inherit from Identity to use this mechanism.  Identity,
// Identifier, and the shorthand boolean operations may have to
// be modified for some different architectures.

class Identity;
class Identifier {
  private:
    const void* id;
  public:
    Identifier(): id(0) {}
    Identifier(const Identity* i): id((void*)i) {}
    Identifier(const Identifier& i): id(i.id) {}
    void operator = (const Identifier& i) { id = i.id; }
    ~Identifier() {}
    
    operator < (const Identifier&i) const { return id < i.id; }
    operator > (const Identifier&i) const { return id > i.id; }
    operator == (const Identifier&i) const { return id == i.id; }
    operator <= (const Identifier&i) const { return id <= i.id; }
    operator >= (const Identifier&i) const { return id >= i.id; }
    operator != (const Identifier&i) const { return id != i.id; }
};

class Identity {
  public:
    virtual ~Identity();
    Identifier identifier() { return this; }
};

// shorthand boolean operation for pointer arguments
inline int lt(const Identity*i, const Identity*j) { return i < j; }
inline int gt(const Identity*i, const Identity*j) { return i > j; }
inline int le(const Identity*i, const Identity*j) { return i <= j; }
inline int ge(const Identity*i, const Identity*j) { return i >= j; }
inline int eq(const Identity*i, const Identity*j) { return i == j; }
inline int ne(const Identity*i, const Identity*j) { return i != j; }
inline int cmp(const Identity*i, const Identity*j)
{
  return (i==j)?0:((i<j)?-1:1);
}

#endif
