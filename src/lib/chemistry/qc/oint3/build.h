#ifndef _chemistry_qc_oint3_build_h
#define _chemistry_qc_oint3_build_h

#include <mpqc_config.h>
#include <chemistry/qc/intv3/array.h>

namespace sc {

#define MG 3

#define DECLARE_BUILD(ii,j,k,l) \
  int i ## ii ## j ## k ## l ();\
  int i ## ii ## j ## k ## l ## eAB ()

class BuildIntV3 {
  public:
    double int_v_ooze;
    double int_v_zeta12;
    double int_v_zeta34;
    double int_v_oo2zeta12;
    double int_v_oo2zeta34;
    double int_v_W0;
    double int_v_W1;
    double int_v_W2;
    double int_v_p120;
    double int_v_p121;
    double int_v_p122;
    double int_v_p340;
    double int_v_p341;
    double int_v_p342;
    double int_v_r10;
    double int_v_r11;
    double int_v_r12;
    double int_v_r20;
    double int_v_r21;
    double int_v_r22;
    double int_v_r30;
    double int_v_r31;
    double int_v_r32;
    double int_v_r40;
    double int_v_r41;
    double int_v_r42;
    double int_v_k12;
    double int_v_k34;
    IntV3Arraydoublep3 int_v_list;
  public:
    BuildIntV3();
    ~BuildIntV3();

    int impossible_integral();

#if (MG == 1) || (MG == 2) || (MG == 3) || (MG == 4)
    DECLARE_BUILD(0,1,0,0);
    DECLARE_BUILD(0,1,0,1);
    DECLARE_BUILD(0,1,1,1);
    DECLARE_BUILD(1,1,0,0);
    DECLARE_BUILD(1,1,1,1);
#endif

#if (MG == 2) || (MG == 3) || (MG == 4)
    DECLARE_BUILD(0,2,0,0);
    DECLARE_BUILD(0,2,0,1);
    DECLARE_BUILD(0,2,0,2);
    DECLARE_BUILD(0,2,1,1);
    DECLARE_BUILD(0,2,1,2);
    DECLARE_BUILD(0,2,2,2);
    DECLARE_BUILD(1,2,0,0);
    DECLARE_BUILD(1,2,0,1);
    DECLARE_BUILD(1,2,1,1);
    DECLARE_BUILD(1,2,1,2);
    DECLARE_BUILD(1,2,2,2);
    DECLARE_BUILD(2,2,0,0);
    DECLARE_BUILD(2,2,0,1);
    DECLARE_BUILD(2,2,1,1);
    DECLARE_BUILD(2,2,2,2);
#endif

#if (MG == 3) || (MG == 4)
    DECLARE_BUILD(0,3,0,0);
    DECLARE_BUILD(0,3,0,1);
    DECLARE_BUILD(0,3,0,2);
    DECLARE_BUILD(0,3,0,3);
    DECLARE_BUILD(0,3,1,1);
    DECLARE_BUILD(0,3,1,2);
    DECLARE_BUILD(0,3,1,3);
    DECLARE_BUILD(0,3,2,2);
    DECLARE_BUILD(0,3,2,3);
    DECLARE_BUILD(0,3,3,3);
    DECLARE_BUILD(1,3,0,0);
    DECLARE_BUILD(1,3,0,1);
    DECLARE_BUILD(1,3,0,2);
    DECLARE_BUILD(1,3,1,1);
    DECLARE_BUILD(1,3,1,2);
    DECLARE_BUILD(1,3,1,3);
    DECLARE_BUILD(1,3,2,2);
    DECLARE_BUILD(1,3,2,3);
    DECLARE_BUILD(1,3,3,3);
    DECLARE_BUILD(2,3,0,0);
    DECLARE_BUILD(2,3,0,1);
    DECLARE_BUILD(2,3,0,2);
    DECLARE_BUILD(2,3,1,1);
    DECLARE_BUILD(2,3,1,2);
    DECLARE_BUILD(2,3,2,2);
    DECLARE_BUILD(2,3,2,3);
    DECLARE_BUILD(2,3,3,3);
    DECLARE_BUILD(3,3,0,0);
    DECLARE_BUILD(3,3,0,1);
    DECLARE_BUILD(3,3,0,2);
    DECLARE_BUILD(3,3,1,1);
    DECLARE_BUILD(3,3,1,2);
    DECLARE_BUILD(3,3,2,2);
    DECLARE_BUILD(3,3,3,3);
#endif
};

}

#endif
