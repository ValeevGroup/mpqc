
#ifndef _chemistry_qc_intv2_storage_h
#define _chemistry_qc_intv2_storage_h

#ifdef __cplusplus
extern "C" {
#endif

    void int_storage(int);
    void int_reduce_storage_threshold();
    void int_done_storage();
    int int_have_stored_integral(int,int,int,int,int,int,int);
    void int_store_integral(int,int,int,int,int,int,int,int);

#ifdef __cplusplus
}
#include <chemistry/qc/intv3/storage.h>
#endif

#endif
