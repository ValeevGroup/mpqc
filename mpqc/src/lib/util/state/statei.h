
#ifndef CLASSNAME
#define CLASSNAME you_forgot_to_define_CLASSNAME
#endif

#define stringize_(arg) # arg
#define stringize(arg) stringize_(arg)

void
CLASSNAME::save_object_state(StateOut&so)
{
  if (class_desc() != static_class_desc()) {
      fprintf(stderr, "Warning:" stringize(CLASSNAME) "::save_object_state:");
      fprintf(stderr," exact type not known -- object not saved\n");
      return;
    }
  // save the version info
  //so.put_version(static_class_desc());
  save_vbase_state(so);
  save_data_state(so);
}

CLASSNAME*
CLASSNAME::restore_state(StateIn&si)
{
#ifdef __GNUC__
  return CLASSNAME::castdown(SavableState::restore_state(si));
#else
  return CLASSNAME::castdown(att_hack_job(si));
#endif
}

#undef stringize
#undef stringize_
