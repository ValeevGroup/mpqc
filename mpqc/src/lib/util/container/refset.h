
#ifndef _refset_h
#define _refset_h

#include <Pix.h>
#include <util/ref/ref.h>
#include <stdlib.h>

#define REFARRAY(Type)							      \
class RefArray ## Type							      \
{									      \
 private:								      \
  int _length;								      \
  Ref ## Type * _array;							      \
 public:								      \
  RefArray ## Type():_array(0),_length(0) {}				      \
  RefArray ## Type(int size):_array(0),_length(0) { set_length(size); }	      \
  ~RefArray ## Type() { clear(); }					      \
  void set_length(int size)						      \
  {									      \
    if (_array) delete[] _array;					      \
    if (size) _array = new Ref ## Type [ size ];			      \
    else _array = 0;							      \
    _length = size;							      \
  }									      \
  inline int length() const { return _length; };			      \
  RefArray ## Type& operator = (RefArray ## Type & s)			      \
  {									      \
    if (_array) delete[] _array;					      \
    _length = s._length;						      \
    if (_length) _array = new Ref ## Type [ _length ];			      \
    else _array = 0; \
    for (int i=0; i<_length; i++) { \
        _array[i] = s._array[i]; \
      } \
    return (*this); \
  }									      \
  inline void clear() { set_length(0); }				      \
  Ref ## Type& operator[] (int i) const \
  { \
    if (i<0 || i>=_length) { \
        fprintf(stderr,"RefArray::operator[] out of range: %d (nelement = %d)\n",\
                i,_length); \
        abort(); \
      }; \
    return _array[i]; \
  }	      \
}

// the old RefArray operator []
//inline Ref ## Type& operator[] (int i) const { return _array[i]; }

// This class combines array capabilities with set capabilities.
// It is basically a set with the iseek and operator[](int i) members
// added.  At the moment, it requires that RefSet use an array internally.
// When RefSet is improved, RefArraySet must be updated.
#define REFARRAYSET(Type) \
class RefArraySet ## Type : public RefSet ## Type \
{ \
 public: \
   Ref ## Type & operator[](int i) const \
     { \
       if (i<0 || i>=nelement) { \
          fprintf(stderr,"RefArraySet::operator[] out of range: %d (nelement = %d)\n",\
                  i,nelement); \
          abort(); \
         }; \
       return element[i]; \
     };   \
    int iseek ( Ref ## Type &item) const				      \
    {									      \
      for (int i=0; i<nelement; i++) {					      \
          if (item == element[i]) return i;		      \
        }								      \
      return 0;								      \
    };									      \
}

#define REFSET(Type)							      \
class  RefSet ## Type 							      \
{									      \
 protected:								      \
  int nelement;								      \
  enum { array_incr=100 }; \
  RefArray ## Type element;						      \
  int range_check(int i) const \
    { \
      if ((i<0) || (i >= nelement)) { \
          fprintf(stderr,"RefSet::range_check(%d): nelement=%d\n",i,nelement);\
          abort(); \
        } \
      return i; \
    }; \
  inline Pix index_to_pix(int i) const {  return (Pix)(range_check(i)+1); }   \
  inline int pix_to_index_nocheck(Pix i) const { return ((int)i)-1; } \
  inline int pix_to_index(Pix i) const { return range_check(((int)i)-1); } \
 public:								    \
   RefSet ## Type ():nelement(0),element(array_incr) {} \
   RefSet ## Type (RefSet ## Type & s) \
     { \
       this->operator = (s); \
     }; \
   RefSet ## Type& operator = (RefSet ## Type & s) \
     { \
       nelement = s.nelement; \
       element.set_length(s.element.length()); \
       for (int i=0; i<nelement; i++) element[i] = s.element[i]; \
       return *this; \
     }; \
  ~ RefSet ## Type () { clear(); }					      \
  void print(FILE*fp=stdout) \
    { \
      int i; \
      fprintf(fp,"RefSet: printing:\n"); \
      for (i=0; i<nelement; i++) { \
          element[i]->print(fp); \
	} \
      fprintf(fp,"RefSet: done printing\n"); \
    } \
  int length() const { return nelement; };				      \
  void clear()								      \
    {									      \
      element.clear();							      \
      nelement = 0;							      \
    };									      \
  Pix add( Ref ## Type & e)						      \
    {									      \
      int i;								      \
      for (i=0; i<nelement; i++) {					      \
         if (e == element[i]) { \
             return index_to_pix(i);			      \
           } \
        }								      \
      if (nelement == element.length()) { \
          RefArray ## Type tmpelement(element.length() + array_incr); \
          for (i=0; i<nelement; i++) { \
              tmpelement[i] = element[i];		      \
             } \
          element = tmpelement; \
	} \
      element[nelement] = e; \
      nelement++;							      \
      return index_to_pix(nelement-1);					      \
    };									      \
   Ref ## Type & operator()(Pix i) const \
     { \
       if (pix_to_index(i)<0 || pix_to_index(i)>=nelement) { \
          fprintf(stderr,"Ref::operator() out of range: %d (nelement = %d)\n",\
                  pix_to_index(i),nelement); \
          abort(); \
         }; \
       return element[pix_to_index(i)]; \
     };   \
  RefSet ## Type & operator += (RefSet ## Type&s) \
    { \
      for (int i=0; i<s.nelement; i++) { \
          add(s.element[i]); \
        } \
    return *this; \
    } \
  RefSet ## Type & operator -= (RefSet ## Type&s) \
    { \
      for (int i=0; i<s.nelement; i++) { \
          del(s.element[i]); \
        } \
    return *this; \
    } \
  void del( Ref ## Type &item) \
    { \
      Pix pix = seek(item); \
      if (!pix) return; \
      for (int i=pix_to_index(pix); i<nelement-1; i++) { \
          element[i] = element[i+1]; \
        } \
      nelement--; \
     } \
  Pix seek( Ref ## Type &item)						      \
    {									      \
      for (int i=0; i<nelement; i++) {					      \
          if (item == element[i]) return index_to_pix(i);		      \
        }								      \
      return 0;								      \
    };									      \
  int contains( Ref ## Type &item) { return (int) seek(item); };	      \
  int owns(Pix i)							      \
    {									      \
      if (pix_to_index(i) >=0 && pix_to_index(i) < nelement) return 1;	      \
      else return 0;							      \
    };									      \
  Pix first()								      \
    {									      \
      if (nelement) return index_to_pix(0); else return 0;		      \
    };									      \
  void next(Pix&i)							      \
    {									      \
      if (pix_to_index_nocheck(i) < nelement - 1) \
        i = index_to_pix(pix_to_index(i)+1);    \
      else i = 0;							      \
    };									      \
}

// WARNING: The following code is wrong, see above for corrections
// template <class Type>
// class RefSet
// {
//  private:
//   int nelement;
//   Ref<Type>** element;
//   inline Pix index_to_pix(int i) { return (Pix)(i+1); }
//   inline int pix_to_index(Pix i) { return ((int)i)-1; }
//  public:
//   RefSet():nelement(0),element(0) {};
//   ~RefSet() { clear(); }
//   int length() { return nelement; };
//   void clear()
//     {
//       for (int i=0; i<nelement; i++) delete element[i];
//       delete[] element;
//       nelement = 0;
//       element = 0;
//     };
//   Pix add(Ref<Type>& e)
//     {
//       int i;
//       for (i=0; i<nelement; i++) {
//          if (e == *element[i]) return index_to_pix(i);
//         }
//       Ref<Type>** tmpelement = new Ref<Type>*[nelement+1];
//       for (i=0; i<nelement; i++) tmpelement[i] = element[i];
//       if (element) delete[] element;
//       element = tmpelement;
//       element[nelement] = new Ref<Type>(e);
//       nelement++;
//       return index_to_pix(nelement-1);
//     };
//   Ref<Type>& operator()(Pix i) { return *element[pix_to_index(i)]; };
//   Pix seek(Ref<Type>&item)
//     {
//       for (int i=0; i<nelement; i++) {
//           if (item == *element[i]) return index_to_pix(i);
//         }
//       return 0;
//     };
//   int contains(Ref<Type>&item) { return (int) seek(item); };
//   int owns(Pix i)
//     {
//       if (pix_to_index(i) >=0 && pix_to_index(i) < nelement) return 1;
//       else return 0;
//     };
//   Pix first()
//     {
//       if (nelement) return index_to_pix(1); else return 0;
//     };
//   void next(Pix&i)
//     {
//       if (pix_to_index(i) < nelement) i = index_to_pix(pix_to_index(i)+1);
//       else i = 0;
//     };
// };

#endif
