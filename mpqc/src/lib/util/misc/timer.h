
#ifndef _util_misc_timer_h
#define _util_misc_timer_h

#ifdef __cplusplus
extern "C" {
#endif

  void tim_enter(const char *);
  void tim_exit(const char *);
  void tim_change(const char *);
  void tim_print(int);

#ifdef __cplusplus
}
#endif

#endif /* _util_misc_timer_h */
