
#ifndef _libqc_state_debug_h
#define _libqc_state_debug_h

#define SavableState_DECLARE_debug(classname) \
  public: \
    static classname* restore_state(StateIn&); \
    void save_parent_state(StateOut&); \
    void restore_parent_state(StateIn&); \
  private:

#define SavableState_IMPL_debug(classname) \
    classname* classname::restore_state(StateIn&si) \
      { \
      ClassDesc* cd; \
      int version; \
      DescribedClass* dc; \
      classname* ss = 0; \
      int objnum = si.getpointer((void**)&dc); \
      printf("objnum = %d, pointer = 0x%x\n",objnum,(void*)dc); \
      if (objnum) { \
        si.get(&cd,version); \
        if (!cd) { \
          fprintf(stderr,#classname "::restore_state: cannot restore\n"); \
          exit(1); \
          } \
        dc = cd->create_described_class(); \
        if (!dc) { \
          fprintf(stderr,#classname "::restore_state: create failed\n"); \
          exit(1); \
          } \
        ss = classname::castdown(dc); \
        if (!ss) { \
          fprintf(stderr,#classname "::restore_state: castdown failed(1)\n"); \
          exit(1); \
          } \
        si.havepointer(objnum,(void*)dc); \
        ss->restore_vbase_state(version,si); \
        ss->restore_data_state(version,si); \
        } \
      else if (dc) { \
        ss = classname::castdown(dc); \
        if (!ss) { \
          fprintf(stderr,#classname "::restore_state: castdown failed(2)\n"); \
          exit(1); \
          } \
        } \
      return ss; \
      } \
    void classname::save_parent_state(StateOut&so) \
      { \
      so.put(classname::class_desc()->version()); \
      classname::save_data_state(so); \
      } \
    void classname::restore_parent_state(StateIn&si) \
      { \
      int version; \
      si.get(version); \
      classname::restore_data_state(version,si); \
      }

#endif
