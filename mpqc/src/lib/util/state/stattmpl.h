
#ifdef __GNUG__
#pragma interface
#endif

template <class T>
class SSRef: public DCRef<T> {
  public:
    SSRef() {}
    SSRef (const SSRef<T> & o): DCRef<T> (o) {}
    SSRef (T * o): DCRef<T> (o) {}
    SSRef (const RefDescribedClassBase&o): DCRef<T> (o) {}
    ~SSRef () {}
    SSRef<T>& operator=(T* cr) {
        DCRef<T>::operator=(cr); return *this; }
    SSRef<T>& operator=(const RefDescribedClassBase & c) {
        DCRef<T>::operator=(c); return *this; }
    SSRef<T>& operator=(const SSRef<T>&c) {
        DCRef<T>::operator=(c); return *this; }
    SSRef (StateIn&s) { restore_state(s); }
    void save_data_state(StateOut&s) { save_state(s); }
    void save_state(StateOut&so) {
        if (so.putpointer(pointer())) {
            so.put(pointer()->class_desc());
            so.have_classdesc();
            pointer()->save_vbase_state(so);
            pointer()->save_data_state(so);
          }
      }
    void restore_state(StateIn&si) {
        SavableState* ss;
        int objnum = si.getpointer((void**)&ss);
        if (objnum) {
            const ClassDesc* cd;
            si.get(&cd);
            si.nextobject(objnum);
            si.have_classdesc();
            DescribedClass* dc = cd->create(si);
            ss = SavableState::castdown(dc);
          }
        T* t = T::castdown(ss);
        if (!t && ss) {
            fprintf(stderr,
                    "SSRef::restore_state() got type \"%s\"\n",
                    ss->class_name());
            abort();
          }
        assign_pointer(t);
      };
};
