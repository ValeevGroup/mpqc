
#ifdef __GNUC__
#pragma interface
#endif

template <class Type>
class SSBArray: public Array<Type> {
  public:
    SSBArray() {}
    SSBArray(const Array<Type>&a): Array<Type>(a) {}
    SSBArray(Type* data,int size): Array<Type>(data,size) {}
    SSBArray(int size): Array<Type>(size) {}
    SSBArray(StateIn&s) {
      s.get(_length);
      if (_length) s.get(_array);
      _managed=1;
    }
    void save_object_state(StateOut&s) {
        s.put(_length);
        if (_length) s.put(_array,_length);
      }
};
