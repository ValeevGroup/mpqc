
#ifdef __GNUG__
#pragma implementation
#endif

template <class Type>
class Array {
  protected:
    int _length;
    int _managed;
    Type * _array;
  public:
    Array():_length(0),_array(0) {}
    Array(const Array<Type>&a):_length(0),_array(0) { operator=(a); }
    Array(Type* data,int size):_length(size),_managed(0),_array(data){}
    Array(int size):_length(0),_array(0) { set_length(size); }
    ~Array() { clear(); }
    int length() const { return _length; };
    void clear() { set_length(0); }
    void set_length(int size) {
        if (_array && _managed) delete[] _array;
        _managed = 1;
        if (size) _array = new Type [ size ];
        else _array = 0;
        _length = size;
      }
    void reset_length(int size) {
        Type*tmp=_array;
        if (size) _array = new Type [ size ];
        else _array = 0;
        int maxi;
        if (size < _length) maxi = size;
        else maxi = _length;
        for (int i=0; i<maxi; i++) _array[i] = tmp[i];
        if (_managed && tmp) delete[] tmp;
        _managed = 1;
        _length = size;
      }
    Array<Type>& operator = (const Array<Type> & s) {
        if (_managed && _array) delete[] _array;
        _managed = 1;
        _length = s._length;
        if (_length) _array = new Type [ _length ];
        else _array = 0;
        for (int i=0; i<_length; i++) {
            _array[i] = s._array[i];
          }
        return (*this);
      }
    Type& operator[] (int i) const {
        if (i<0 || i>=_length) {
            fprintf(stderr,"Array::operator[](%d) out of range (%d)\n",
                    i,_length);
            abort();
          };
        return _array[i];
      }
    Type& operator() (int i) const {
        if (i<0 || i>=_length) {
            fprintf(stderr,"Array::operator()(%d) out of range (%d)\n",
                    i,_length);
            abort();
          };
        return _array[i];
      }
};
