
#ifndef _math_scmat_disthql_h
#define _math_scmat_disthql_h

#include <util/group/message.h>

namespace sc {

void dist_diagonalize(int n, int m, double *a, double *d, double *v,
                      const Ref<MessageGrp> &);

}

#endif
