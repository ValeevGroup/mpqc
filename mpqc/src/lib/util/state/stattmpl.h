
#ifdef __GNUG__
#pragma interface
#endif

template <class T>
class SSRef: public DCRef<T>, public SSRefBase {
  public:
    SSRef() {}
    SSRef (const SSRef<T> & o): DCRef<T> (o) {}
    SSRef (T * o): DCRef<T> (o) {}
    SSRef (const DCRefBase&o): DCRef<T> (o) {}
    ~SSRef () {}
    SSRef<T>& operator=(T* cr) {
        DCRef<T>::operator=(cr); return *this; }
    SSRef<T>& operator=(const DCRefBase & c) {
        DCRef<T>::operator=(c); return *this; }
    SSRef<T>& operator=(const SSRef<T>&c) {
        DCRef<T>::operator=(c); return *this; }
    SSRef (StateIn&s) { restore_state(s); }
    SavableState *sspointer() { return p; }
    void restore_state(StateIn&si) {
        SavableState* ss = restore_ss(si);
        T* t = T::castdown(ss);
        check_castdown_result((void*)t,ss);
        assign_pointer(t);
      };
};
