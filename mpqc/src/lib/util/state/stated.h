
#ifndef CLASSNAME
#define CLASSNAME you_forgot_to_define_CLASSNAME
#endif

// this should be included in the middle of a class declaration

public:
void save_object_state(StateOut&);
static CLASSNAME* restore_state(StateIn&);
private:
