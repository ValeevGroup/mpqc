
#ifndef _util_keyval_ipv2_scan_h
#define _util_keyval_ipv2_scan_h

struct ip_string_list_struct {
  char *string;
  struct ip_string_list_struct *p;
  };
typedef struct ip_string_list_struct ip_string_list_t;

#endif
