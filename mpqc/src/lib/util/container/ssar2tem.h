
#ifdef __GNUC__
#pragma interface
#endif

template <class Type>
class SSBArray2: public Array2<Type> {
  public:
    SSBArray2() {}
    SSBArray2(const Array2<Type> &a): Array2<Type>(a) {}
    SSBArray2(Type* data,int size0,int size1): Array2<Type>(data,size0,size1){}
    SSBArray2(int size0,int size1): Array2<Type>(size0,size1) {}
    SSBArray2(StateIn&s) {
        s.get(_length0);
        s.get(_length1);
        if (_length0&&_length1) s.get(_array);
      }
    void save_object_state(StateOut&s) {
        s.put(_length0);
        s.put(_length1);
        if (_length0&&_length1) s.put(_array,_length0*_length1);
      }
};
