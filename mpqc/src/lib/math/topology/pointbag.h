
#ifndef _libQC_pointbag_h
#define _libQC_pointbag_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/topology/point.h>

//This is not working gcc 2.3.3 for obscure reasons
// #include "ptrset.h"
// template <class Type>
// class PointBag
// {
//  private:
//   PtrSet<PointBagElem<Type>> impl;
//  public:
//   PointBag() {};
//   ~PointBag() {
//       for (Pix i=first(); i!=0; next(i)) {
// 	  delete impl(i);
// 	}
//     };
//   inline length() { return impl.length(); };
//   add(Point p, Type x) {
//       impl.add(new PointBagElem<Type>(p,x));
//     };
//   inline int owns(Pix i) { return impl.owns(i); }
//   inline Type& operator()(Pix i) { return impl(i).object(); }
//   inline Point& point(Pix i) { return impl(i).point(); }
//   inline int first() { return impl.first(); }
//   inline void next(Pix&i) { impl.next(i); }
// };

#include <util/container/voidptrSet.h>
#include <util/container/voidptrAVLSet.h>

// This is not working for even more obscure reasons:
// template <class Type>
// class PointBagElem
// {
//  private:
//   Type obj;
//   Point pnt;
//  public:
//   PointBagElem(Point&p, Type&t):pnt(p),obj(t) {};
//   inline Point& point() { return pnt; };
//   inline Type& object() { return obj; };
// };
// 
// template <class Type>
// class PointBag
// {
//  private:
//   VoidPtrAVLSet impl;
//  public:
//   PointBag() {};
//   ~PointBag() {
//       for (Pix i=first(); i!=0; next(i)) {
// 	  delete (PointBagElem<Type>*) impl(i);
// 	}
//     };
//   inline length() { return impl.length(); };
//   add(Point p, Type x) {
//       impl.add((void*)new PointBagElem<Type>(p,x));
//     };
//   inline int owns(Pix i) { return impl.owns(i); };
//   inline Type& operator()(Pix i) {
//       return ((PointBagElem<Type>*)impl(i))->object();
//     };
//   inline Type& get(Pix i) { return operator()(i); }
//   inline Point& point(Pix i) { return ((PointBagElem<Type>*)impl(i))->point(); };
//   inline Pix first() { return impl.first(); };
//   inline void next(Pix&i) { impl.next(i); };
//};

// Overcome bugs in gcc 2.3.3 by getting rid of templates for now:

class PointBagElem_double
{
  private:
    double obj;
    Point pnt;
  public:
    PointBagElem_double(const Point&p, const double&t);
    PointBagElem_double(const PointBagElem_double&);
    ~PointBagElem_double();
    Point& point();
    double& object();
};

class PointBag_double
{
  private:
    VoidPtrAVLSet impl;
  public:
    PointBag_double();
    PointBag_double(PointBag_double&);
    ~PointBag_double();
    int length();
    void add(const Point& p, double x);
    void add(const PointBagElem_double&);
    int owns(Pix i);
    double& operator()(Pix i);
    double& get(Pix i);
    Point& point(Pix i);
    Pix first();
    void next(Pix&i);
};

#ifdef INLINE_FUNCTIONS
#include <math/topology/pntbag_i.h>
#endif

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
