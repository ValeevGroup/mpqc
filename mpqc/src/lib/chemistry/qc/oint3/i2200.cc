#include <chemistry/qc/oint3/build.h>
int sc::BuildIntV3::i2200(){
/* the cost is 51 */
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
double t13;
t1=int_v_zeta34*int_v_ooze;
t2=int_v_oo2zeta12*t1;
t1=(-1)*t2;
double***int_v_list0=int_v_list(0);
double**int_v_list00=int_v_list0[0];
double*int_v_list001=int_v_list00[1];
t2=t1*int_v_list001[0];
double*int_v_list000=int_v_list00[0];
t1=int_v_oo2zeta12*int_v_list000[0];
t3=t1+t2;
t1=int_v_W0-int_v_p120;
double*int_v_list002=int_v_list00[2];
t2=t1*int_v_list002[0];
t4=int_v_p120-int_v_r10;
t5=t4*int_v_list001[0];
t6=t5+t2;
t2=t1*t6;
t5=t2+t3;
t2=t1*int_v_list001[0];
t1=t4*int_v_list000[0];
t7=t1+t2;
t1=t4*t7;
t2=t1+t5;
double***int_v_list2=int_v_list(2);
double**int_v_list20=int_v_list2[0];
double*int_v_list200=int_v_list20[0];
int_v_list200[5]=t2;
t1=int_v_W2-int_v_p122;
t4=t1*t6;
t5=int_v_p122-int_v_r12;
t8=t5*t7;
t9=t8+t4;
int_v_list200[4]=t9;
t4=int_v_W1-int_v_p121;
t8=t6*t4;
t6=int_v_p121-int_v_r11;
t10=t7*t6;
t7=t10+t8;
int_v_list200[3]=t7;
t8=t1*int_v_list002[0];
t10=t5*int_v_list001[0];
t11=t10+t8;
t8=t1*t11;
t10=t3+t8;
t8=t1*int_v_list001[0];
t11=t5*int_v_list000[0];
t12=t11+t8;
t8=t5*t12;
t11=t8+t10;
int_v_list200[2]=t11;
t8=int_v_list002[0]*t4;
t10=t6*int_v_list001[0];
t12=t10+t8;
t8=t1*t12;
t1=int_v_list001[0]*t4;
t10=t6*int_v_list000[0];
t13=t10+t1;
t1=t5*t13;
t5=t1+t8;
int_v_list200[1]=t5;
t1=t4*t12;
t4=t3+t1;
t1=t6*t13;
t3=t1+t4;
int_v_list200[0]=t3;
return 1;}
