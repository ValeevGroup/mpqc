
#ifndef _util_container_set_h
#define _util_container_set_h

#include <Pix.h>
#include <stdlib.h>
#include <util/container/array.h>

#define SET_dec(Type)							      \
class Set ## Type							      \
{									      \
 protected:								      \
  int nelement;								      \
  int array_length;							      \
  enum { array_incr=100 };						      \
  Array ## Type element;						      \
  int range_check(int i) const;						      \
  Pix index_to_pix(int i) const;					      \
  int pix_to_index_nocheck(Pix i) const;				      \
  int pix_to_index(Pix i) const;					      \
 public:								      \
  Set ## Type ();							      \
  Set ## Type (const Set ## Type & s);					      \
  Set ## Type& operator = (const Set ## Type & s);			      \
  ~Set ## Type ();							      \
  int length() const;							      \
  void clear();								      \
  Pix add( Type & e);							      \
  Type & operator()(Pix i) const;					      \
  Set ## Type & operator += (const Set ## Type&s);			      \
  Pix seek( Type &item);						      \
  int contains( Type &item);						      \
  void del( Type &item);						      \
  int owns(Pix i);							      \
  Pix first();								      \
  void next(Pix&i);							      \
}

#define SET_def(Type)							      \
Pix Set ## Type :: index_to_pix(int i) const				      \
{ return (Pix)(range_check(i)+1); }					      \
int Set ## Type :: pix_to_index_nocheck(Pix i) const			      \
{ return ((int)i)-1; }							      \
int Set ## Type :: pix_to_index(Pix i) const				      \
{ return range_check(((int)i)-1); }					      \
int Set ## Type :: length() const { return nelement; };			      \
int Set ## Type :: range_check(int i) const				      \
{									      \
  if ((i<0) || (i >= nelement)) {					      \
      fprintf(stderr,"Set::range_check(%d): nelement=%d\n",i,nelement);	      \
      abort();								      \
    }									      \
  return i;								      \
}									      \
Set ## Type :: Set ## Type ():nelement(0),element(array_incr) {}	      \
Set ## Type :: Set ## Type (const Set ## Type & s)			      \
{									      \
  this->operator = (s);							      \
}									      \
Set ## Type& Set ## Type :: operator = (const Set ## Type & s)		      \
{									      \
  nelement = s.nelement;						      \
  element.set_length(s.element.length());				      \
  for (int i=0; i<nelement; i++) element[i] = s.element[i];		      \
  return *this;								      \
}									      \
Set ## Type :: ~ Set ## Type () { clear(); }				      \
void Set ## Type :: clear()						      \
{									      \
  element.clear();							      \
  nelement = 0;								      \
}									      \
Pix Set ## Type :: add( Type & e)					      \
{									      \
  int i;								      \
  for (i=0; i<nelement; i++) {						      \
      if (e == element[i]) {						      \
          return index_to_pix(i);					      \
        }								      \
    }									      \
  if (nelement == element.length()) {					      \
      Array ## Type tmpelement(element.length() + array_incr);		      \
      for (i=0; i<nelement; i++) {					      \
          tmpelement[i] = element[i];					      \
        }								      \
      element = tmpelement;						      \
    }									      \
  element[nelement] = e;						      \
  nelement++;								      \
  return index_to_pix(nelement-1);					      \
}									      \
void Set ## Type :: del( Type & e)					      \
{									      \
  int i;								      \
  for (i=0; i<nelement; i++) {						      \
      if (e == element[i]) {						      \
          nelement--;							      \
          for (int j=i; j<nelement; j++) {				      \
              element[j] = element[j+1];				      \
            }								      \
          element[j] = 0;						      \
          break;							      \
        }								      \
    }									      \
}									      \
Type & Set ## Type :: operator()(Pix i) const				      \
{									      \
  if (pix_to_index(i)<0 || pix_to_index(i)>=nelement) {			      \
      fprintf(stderr,"Set::operator() out of range: %d (nelement = %d)\n",    \
              pix_to_index(i),nelement);				      \
      abort();								      \
    };									      \
  return element[pix_to_index(i)];					      \
}									      \
Set ## Type & Set ## Type :: operator += (const Set ## Type&s)		      \
{									      \
  for (int i=0; i<s.nelement; i++) {					      \
      add(s.element[i]);						      \
    }									      \
  return *this;								      \
}									      \
Pix Set ## Type :: seek( Type &item)					      \
{									      \
  for (int i=0; i<nelement; i++) {					      \
      if (item == element[i]) return index_to_pix(i);			      \
    }									      \
  return 0;								      \
}									      \
int Set ## Type :: contains( Type &item) { return (int) seek(item); };	      \
int Set ## Type :: owns(Pix i)						      \
{									      \
  if (pix_to_index(i) >=0 && pix_to_index(i) < nelement) return 1;	      \
  else return 0;							      \
}									      \
Pix Set ## Type :: first()						      \
{									      \
  if (nelement) return index_to_pix(0); else return 0;			      \
}									      \
void Set ## Type :: next(Pix&i)						      \
{									      \
  if (pix_to_index_nocheck(i) < nelement - 1)				      \
    i = index_to_pix(pix_to_index(i)+1);				      \
  else i = 0;								      \
}

// This class combines array capabilities with set capabilities.
// It is basically a set with the iseek and operator[](int i) members
// added.  At the moment, it requires that Set use an array internally.
// When Set is improved, Arrayset must be updated.
#define ARRAYSET_dec(Type)						      \
class Arrayset ## Type : public Set ## Type				      \
{									      \
 public:								      \
  Arrayset ## Type ();							      \
  Arrayset ## Type (const Arrayset ## Type&);				      \
  Type & operator[](int i);						      \
  const Type & operator[](int i) const;					      \
  int iseek( Type &item);						      \
  ~Arrayset ## Type();							      \
}

#define ARRAYSET_def(Type)						      \
Arrayset ## Type :: Arrayset ## Type ()					      \
{									      \
}									      \
Arrayset ## Type :: Arrayset ## Type (const Arrayset ## Type&a):	      \
Set ## Type(a)								      \
{									      \
}									      \
Type & Arrayset ## Type :: operator[](int i)				      \
{									      \
  if (i<0 || i>=nelement) {						      \
      fprintf(stderr,"Arrayset::operator[] out of range: %d (nelement = %d)\n", \
              i,nelement);						      \
      abort();								      \
    };									      \
  return element[i];							      \
}									      \
const Type & Arrayset ## Type :: operator[](int i) const		      \
{									      \
  if (i<0 || i>=nelement) {						      \
      fprintf(stderr,"Arrayset::operator[] out of range: %d (nelement = %d)\n", \
              i,nelement);						      \
      abort();								      \
    };									      \
  return element[i];							      \
}									      \
int Arrayset ## Type :: iseek( Type &item)				      \
{									      \
  for (int i=0; i<nelement; i++) {					      \
      if (item == element[i]) return i;					      \
    }									      \
  return -1;								      \
}									      \
Arrayset ## Type ::~Arrayset ## Type ()					      \
{									      \
}

// declare sets for the basic types
SET_dec(int);
ARRAYSET_dec(int);
ARRAY_dec(Arraysetint);
SET_dec(double);
ARRAYSET_dec(double);

#endif
