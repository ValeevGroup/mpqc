
#ifndef _chemistry_qc_force_read_h
#define _chemistry_qc_force_read_h

#ifdef __cplusplus
void
dmt_force_osscf_keyval_init(KeyVal*,FILE*);
void
dmt_force_csscf_keyval_init(KeyVal*,FILE*);

extern "C" {
int
dmt_force_read_and_bcast_boolean(KeyVal*keyval,FILE*,char*name,int*boolval);
int
dmt_force_read_and_bcast_int(KeyVal*keyval,FILE*,char*name,int*intval);
int
dmt_force_read_string(KeyVal*keyval,FILE*,char*name,char**val);
}
#else
int
dmt_force_read_and_bcast_boolean(void*keyval,FILE*,char*name,int*boolval);
int
dmt_force_read_and_bcast_int(void*keyval,FILE*,char*name,int*intval);
int
dmt_force_read_string(void*keyval,FILE*,char*name,char**val);
#endif

#endif
