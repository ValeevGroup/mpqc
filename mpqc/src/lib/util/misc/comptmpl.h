
template <class T>
class Result: public ResultBase {
  private:
    T _result;
  public:
    Result(Compute*c):ResultBase(c) {};
    Result(const Result<T> &r, Compute*c):ResultBase(c)
    { _result=r._result; }
    operator T&() { update(); return _result; };
    operator T() { update(); return _result; };
    T* operator ->() { update(); return &_result; };
    T& result() { update(); return _result; };
    T& result_noupdate() { return _result; };
    const T& result_noupdate() const { return _result; };
    void operator=(const T& a) { _result = a; }
    void operator=(const Result<T> &r)
       { Result::operator=(r); _result = r._result; };
};

// Accuracy result with any type.  The result datum is not saved
// or restored.
template <class T>
class AccResult: public AccResultBase {
  private:
    T _result;
  public:
    AccResult(Compute*c):AccResultBase(c) {};
    AccResult(const AccResult<T> &r, Compute*c):AccResultBase(c)
    { _result=r._result; }
    operator T&() { update(); return _result; };
    operator T() { update(); return _result; };
    T* operator ->() { update(); return &_result; };
    T& result() { update(); return _result; };
    T& result_noupdate() { return _result; };
    const T& result_noupdate() const { return _result; };
    void operator=(const T& a) { _result = a; }
    void operator=(const AccResult<T> &r)
       { AccResultBase::operator=(r); _result = r._result; };
    void save_data_state(StateOut&s)
    {
      AccResultBase::save_data_state(s);
    }
    AccResult(StateIn&s,Compute*c): AccResultBase(s,c) {}
};

// Accuracy Result with SavableState type
template <class T>
class SSAccResult: public AccResultBase {
  private:
    T _result;
  public:
    SSAccResult(Compute*c):AccResultBase(c) {};
    SSAccResult(const SSAccResult<T> &r, Compute*c):AccResultBase(c)
    { _result=r._result; }
    operator T&() { update(); return _result; };
    operator T() { update(); return _result; };
    T* operator ->() { update(); return &_result; };
    T& result() { update(); return _result; };
    T& result_noupdate() { return _result; };
    const T& result_noupdate() const { return _result; };
    void operator=(const T& a) { _result = a; }
    void operator=(const SSAccResult<T> &r)
       { AccResultBase::operator=(r); _result = r._result; };
    void save_data_state(StateOut&s)
    {
      AccResultBase::save_data_state(s);
      _result.save_data_state(s);
    }
    SSAccResult(StateIn&s,Compute*c): AccResultBase(s,c), _result(s) {}
};

// Accuracy Result with non-class type
template <class T>
class NCAccResult: public AccResultBase {
  private:
    T _result;
  public:
    NCAccResult(Compute*c):AccResultBase(c) {};
    NCAccResult(const NCAccResult<T> &r, Compute*c):AccResultBase(c)
    { _result=r._result; }
    operator T&() { update(); return _result; };
    operator T() { update(); return _result; };
    T* operator ->() { update(); return &_result; };
    T& result() { update(); return _result; };
    T& result_noupdate() { return _result; };
    const T& result_noupdate() const { return _result; };
    void operator=(const T& a) { _result = a; }
    void operator=(const NCAccResult<T> &r)
       { AccResultBase::operator=(r); _result = r._result; };
    void save_data_state(StateOut&s)
    {
      AccResultBase::save_data_state(s);
      s.put(_result);
    }
    NCAccResult(StateIn&s,Compute*c): AccResultBase(s,c) {s.get(_result);}
};
