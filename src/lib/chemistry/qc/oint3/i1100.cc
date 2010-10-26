#include <chemistry/qc/oint3/build.h>
int sc::BuildIntV3::i1100(){
/* the cost is 15 */
double t1;
double t2;
double t3;
double t4;
double t5;
t1=int_v_W0-int_v_p120;
double***restrictxx int_v_list0=int_v_list(0);
double**restrictxx int_v_list00=int_v_list0[0];
double*restrictxx int_v_list001=int_v_list00[1];
t2=int_v_list001[0]*t1;
t1=int_v_p120-int_v_r10;
double*restrictxx int_v_list000=int_v_list00[0];
t3=int_v_list000[0]*t1;
t1=t3+t2;
double***restrictxx int_v_list1=int_v_list(1);
double**restrictxx int_v_list10=int_v_list1[0];
double*restrictxx int_v_list100=int_v_list10[0];
int_v_list100[2]=t1;
t2=int_v_W2-int_v_p122;
t3=int_v_list001[0]*t2;
t2=int_v_p122-int_v_r12;
t4=int_v_list000[0]*t2;
t2=t4+t3;
int_v_list100[1]=t2;
t3=int_v_W1-int_v_p121;
t4=t3*int_v_list001[0];
t3=int_v_p121-int_v_r11;
t5=t3*int_v_list000[0];
t3=t5+t4;
int_v_list100[0]=t3;
return 1;}
