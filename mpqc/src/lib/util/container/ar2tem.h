
template <class Type>
class Array2 {
  protected:
    int _length0;
    int _length1;
    int _managed;
    Type * _array;
  public:
    Array2():_length0(0),_length1(0),_array(0) {}
    Array2(const Array2<Type> &a): _length0(0),_length1(0),_array(0) {
        operator=(a);
      }
    Array2(Type* data,int size0,int size1):
      _length0(size0),_length1(size1),_managed(0),_array(data) {}
    Array2(int size0,int size1): _length0(0),_length1(0),_array(0) {
        set_lengths(size0,size1);
      }
    ~Array2() { clear(); }
    int length0() const { return _length0; };
    int length1() const { return _length1; };
    void clear() { set_lengths(0,0); }
    void set_lengths(int size0,int size1) {
        if (_managed && _array) delete[] _array;
        _managed = 1;
        if (size0*size1) _array = new Type [ size0*size1 ];
        else _array = 0;
        _length0 = size0;
        _length1 = size1;
      }
    Array2<Type>& operator = (const Array2<Type> & s) {
        if (_managed && _array) delete[] _array;
        _managed = 1;
        _length0 = s._length0;
        _length1 = s._length1;
        if (_length0*_length1) _array = new Type [ _length0*_length1 ];
        else _array = 0;
        for (int i=0; i<_length0*_length1; i++) {
            _array[i] = s._array[i];
          }
        return (*this);
      }
    Type& operator() (int i,int j) {
        if (i<0 || i>=_length0 || j<0 || j>=_length1) {
            fprintf(stderr,"Array2::operator()(%d,%d): out of range (%d,%d)\n",
                    i,j,_length0,_length1);
            abort();
          };
        return _array[i*_length1+j];
      }
    const Type& operator() (int i,int j) const {
        if (i<0 || i>=_length0 || j<0 || j>=_length1) {
            fprintf(stderr,"Array2::operator()(%d,%d): out of range (%d,%d)\n",
                    i,j,_length0,_length1);
            abort();
          };
        return _array[i*_length1+j];
      }
};
