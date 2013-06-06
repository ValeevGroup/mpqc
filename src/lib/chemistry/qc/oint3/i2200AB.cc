#include <chemistry/qc/oint3/build.h>
int sc::BuildIntV3::i2200eAB(){
/* the cost is 21 */
double t1;
double t2;
double t3;
double t4;
double t5;
double t6;
double t7;
double t8;
double t9;
t1=int_v_zeta34*int_v_ooze;
t2=int_v_oo2zeta12*t1;
t1=(-1)*t2;
double***RESTRICT int_v_list0=int_v_list(0);
double**RESTRICT int_v_list00=int_v_list0[0];
double*RESTRICT int_v_list001=int_v_list00[1];
t2=int_v_list001[0]*t1;
double*RESTRICT int_v_list000=int_v_list00[0];
t1=int_v_list000[0]*int_v_oo2zeta12;
t3=t1+t2;
t1=int_v_W0-int_v_p120;
double*RESTRICT int_v_list002=int_v_list00[2];
t2=int_v_list002[0]*t1;
t4=t1*t2;
t1=t4+t3;
double***RESTRICT int_v_list2=int_v_list(2);
double**RESTRICT int_v_list20=int_v_list2[0];
double*RESTRICT int_v_list200=int_v_list20[0];
int_v_list200[5]=t1;
t4=int_v_W2-int_v_p122;
t5=t4*t2;
int_v_list200[4]=t5;
t6=int_v_W1-int_v_p121;
t7=t2*t6;
int_v_list200[3]=t7;
t2=int_v_list002[0]*t4;
t8=t4*t2;
t2=t3+t8;
int_v_list200[2]=t2;
t8=int_v_list002[0]*t6;
t9=t4*t8;
int_v_list200[1]=t9;
t4=t6*t8;
t6=t3+t4;
int_v_list200[0]=t6;
return 1;}
