#include <chemistry/qc/oint3/build.h>
int BuildIntV3::i1100eAB(){
/* the cost is 6 */
double t1;
double t2;
double t3;
double t4;
t1=int_v_W0-int_v_p120;
t2=int_v_list(0,0,1)[0]*t1;
int_v_list(1,0,0)[2]=t2;
t1=int_v_W2-int_v_p122;
t3=t1*int_v_list(0,0,1)[0];
int_v_list(1,0,0)[1]=t3;
t1=int_v_W1-int_v_p121;
t4=t1*int_v_list(0,0,1)[0];
int_v_list(1,0,0)[0]=t4;
return 1;}
