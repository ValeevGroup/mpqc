
#ifndef _libQC_pointset_h
#define _libQC_pointset_h

#include <math/topology/point.h>
#include <util/container/ptrset.h>

template <class Type>
class PointSetElem<Type>
{
 private:
  Type obj;
  Point pnt;
 public:
  PointSetElem(Point p, Type t):pnt(p),obj(t) {};
  Point& point() { return pnt; };
  Type& object() { return obj; };
};

template <class Type>
class PointSet:
{
 private:
  PtrSet<PointSetElem<Type>> impl;
 public:
  PointSet() {};
  ~PointSet() {};
  length() { return impl.length(); }
  empty() { return 
};

#endif
