#include <chemistry/qc/oint3/build.h>
int BuildIntV3::i0100(){
/* the cost is 15 */
double t1;
double t2;
double t3;
double t4;
double t5;
t1=int_v_W0-int_v_p120;
t2=int_v_list.dp[0][0][1][0]*t1;
t1=int_v_p120-int_v_r10;
t3=int_v_list.dp[0][0][0][0]*t1;
t1=t3+t2;
int_v_list.dp[1][0][0][2]=t1;
t2=int_v_W2-int_v_p122;
t3=int_v_list.dp[0][0][1][0]*t2;
t2=int_v_p122-int_v_r12;
t4=int_v_list.dp[0][0][0][0]*t2;
t2=t4+t3;
int_v_list.dp[1][0][0][1]=t2;
t3=int_v_W1-int_v_p121;
t4=t3*int_v_list.dp[0][0][1][0];
t3=int_v_p121-int_v_r11;
t5=t3*int_v_list.dp[0][0][0][0];
t3=t5+t4;
int_v_list.dp[1][0][0][0]=t3;
return 1;}
