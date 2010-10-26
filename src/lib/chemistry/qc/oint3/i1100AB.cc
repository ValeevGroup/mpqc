#include <chemistry/qc/oint3/build.h>
int sc::BuildIntV3::i1100eAB(){
/* the cost is 6 */
double t1;
double t2;
double t3;
double t4;
t1=int_v_W0-int_v_p120;
double***restrictxx int_v_list0=int_v_list(0);
double**restrictxx int_v_list00=int_v_list0[0];
double*restrictxx int_v_list001=int_v_list00[1];
t2=int_v_list001[0]*t1;
double***restrictxx int_v_list1=int_v_list(1);
double**restrictxx int_v_list10=int_v_list1[0];
double*restrictxx int_v_list100=int_v_list10[0];
int_v_list100[2]=t2;
t1=int_v_W2-int_v_p122;
t3=t1*int_v_list001[0];
int_v_list100[1]=t3;
t1=int_v_W1-int_v_p121;
t4=t1*int_v_list001[0];
int_v_list100[0]=t4;
return 1;}
