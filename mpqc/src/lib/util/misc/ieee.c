
#ifdef SGI
#include <sys/fpu.h>
ieee_trap_errors()
{
  union fpc_csr fc;

  fc.fc_word = get_fpc_csr();
  fc.fc_struct.en_divide0 = 1;
  fc.fc_struct.en_invalid = 1;
  fc.fc_struct.en_overflow = 1;
  set_fpc_csr(fc.fc_word);
  }

#else
ieee_trap_errors()
{
}
#endif
