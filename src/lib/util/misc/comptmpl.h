//
// comptmpl.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <util/misc/compute.h>

namespace sc {

/** Result are members of Compute specializations that keep track of
  whether or not a particular result should be computed or if it has
  already been computed.  For non-class template parameters, use
  NCResult. */
template <class T>
class Result: public ResultInfo {
  private:
    T _result;
  public:
    Result(Compute*c):ResultInfo(c) {};
    Result(const Result<T> &r, Compute*c):ResultInfo(c)
    { _result=r._result; }
    operator T&() { update(); return _result; };
    T* operator ->() { update(); return &_result; };
    T& result() { update(); return _result; };
    T& result_noupdate() { return _result; };
    const T& result_noupdate() const { return _result; };
    void operator=(const T& a) { _result = a; }
    void operator=(const Result<T> &r)
       { ResultInfo::operator=(r); _result = r._result; };
};

/** This is similar to Result, but can be used with
    non-class types. */
template <class T>
class NCResult: public ResultInfo {
  private:
    T _result;
  public:
    NCResult(Compute*c):ResultInfo(c) {};
    NCResult(const NCResult<T> &r, Compute*c):ResultInfo(c)
    { _result=r._result; }
    operator T&() { update(); return _result; };
    T& result() { update(); return _result; };
    T& result_noupdate() { return _result; };
    const T& result_noupdate() const { return _result; };
    void operator=(const T& a) { _result = a; }
    void operator=(const NCResult<T> &r)
       { ResultInfo::operator=(r); _result = r._result; };
};

/** This associates a result datum with an accuracy.  If the result datum
    is to be saved or restored use SSAccResult. */
template <class T>
class AccResult: public AccResultInfo {
  private:
    T _result;
  public:
    AccResult(Compute*c):AccResultInfo(c) {};
    AccResult(const AccResult<T> &r, Compute*c):AccResultInfo(c)
    { _result=r._result; }
    operator T&() { update(); return _result; };
    T* operator ->() { update(); return &_result; };
    T& result() { update(); return _result; };
    T& result_noupdate() { return _result; };
    const T& result_noupdate() const { return _result; };
    void operator=(const T& a) { _result = a; }
    void operator=(const AccResult<T> &r)
       { AccResultInfo::operator=(r); _result = r._result; };
    void restore_state(StateIn&s) {
      AccResultInfo::restore_state(s);
    }
    void save_data_state(StateOut&s)
    {
      AccResultInfo::save_data_state(s);
    }
    AccResult(StateIn&s,Compute*c): AccResultInfo(s,c) {}
};

/** This associates a result datum with an accuracy.  The
    datum must be a SavableState and will be saved and restored. */
template <class T>
class SSAccResult: public AccResultInfo {
  private:
    T _result;
  public:
    SSAccResult(Compute*c):AccResultInfo(c) {};
    SSAccResult(const SSAccResult<T> &r, Compute*c):AccResultInfo(c)
    { _result=r._result; }
    operator T&() { update(); return _result; };
    T* operator ->() { update(); return &_result; };
    T& result() { update(); return _result; };
    T& result_noupdate() { return _result; };
    const T& result_noupdate() const { return _result; };
    void operator=(const T& a) { _result = a; }
    void operator=(const SSAccResult<T> &r)
       { AccResultInfo::operator=(r); _result = r._result; };
    void restore_state(StateIn&s) {
      AccResultInfo::restore_state(s);
      _result.restore_state(s);
    }
    void save_data_state(StateOut&s)
    {
      AccResultInfo::save_data_state(s);
      _result.save_data_state(s);
    }
    SSAccResult(StateIn&s,Compute*c): AccResultInfo(s,c), _result(s) {}
};

/** This associates a result non-class datum with an accuracy. */
template <class T>
class NCAccResult: public AccResultInfo {
  private:
    T _result;
  public:
    NCAccResult(Compute*c):AccResultInfo(c) {};
    NCAccResult(const NCAccResult<T> &r, Compute*c):AccResultInfo(c)
    { _result=r._result; }
    operator T&() { update(); return _result; };
    T& result() { update(); return _result; };
    T& result_noupdate() { return _result; };
    const T& result_noupdate() const { return _result; };
    void operator=(const T& a) { _result = a; }
    void operator=(const NCAccResult<T> &r)
       { AccResultInfo::operator=(r); _result = r._result; };
    void restore_state(StateIn&s) {
      AccResultInfo::restore_state(s);
      s.get(_result);
    }
    void save_data_state(StateOut&s)
    {
      AccResultInfo::save_data_state(s);
      s.put(_result);
    }
    NCAccResult(StateIn&s,Compute*c): AccResultInfo(s,c) {s.get(_result);}
};

}

// ///////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
