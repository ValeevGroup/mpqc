#include <chemistry/qc/oint3/build.h>
int sc::BuildIntV3::i1111eAB(){
/* the cost is 32 */
double t1;
double t2;
double t3;
double t4;
double t5;
double t6;
double t7;
double t8;
double t9;
double t10;
double t11;
double t12;
t1=0.5*int_v_ooze;
double***int_v_list0=int_v_list(0);
double**int_v_list00=int_v_list0[0];
double*int_v_list001=int_v_list00[1];
t2=t1*int_v_list001[0];
t1=int_v_W0-int_v_p340;
double*int_v_list002=int_v_list00[2];
t3=t1*int_v_list002[0];
t1=int_v_p340-int_v_r30;
t4=t1*int_v_list001[0];
t1=t4+t3;
t3=int_v_W0-int_v_p120;
t4=t3*t1;
t5=t4+t2;
double***int_v_list1=int_v_list(1);
double**int_v_list11=int_v_list1[1];
double*int_v_list110=int_v_list11[0];
int_v_list110[8]=t5;
t4=int_v_W2-int_v_p342;
t6=t4*int_v_list002[0];
t4=int_v_p342-int_v_r32;
t7=t4*int_v_list001[0];
t4=t7+t6;
t6=t3*t4;
int_v_list110[7]=t6;
t7=int_v_W1-int_v_p341;
t8=t7*int_v_list002[0];
t7=int_v_p341-int_v_r31;
t9=t7*int_v_list001[0];
t7=t9+t8;
t8=t3*t7;
int_v_list110[6]=t8;
t3=int_v_W2-int_v_p122;
t9=t3*t1;
int_v_list110[5]=t9;
t10=t3*t4;
t11=t2+t10;
int_v_list110[4]=t11;
t10=t3*t7;
int_v_list110[3]=t10;
t3=int_v_W1-int_v_p121;
t12=t1*t3;
int_v_list110[2]=t12;
t1=t3*t4;
int_v_list110[1]=t1;
t4=t7*t3;
t3=t2+t4;
int_v_list110[0]=t3;
return 1;}
