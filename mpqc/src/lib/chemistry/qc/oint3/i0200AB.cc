#include <chemistry/qc/oint3/build.h>
int BuildIntV3::i0200eAB(){
/* the cost is 24 */
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
t1=int_v_zeta34*int_v_ooze;
t2=int_v_oo2zeta12*t1;
t1=(-1)*t2;
t2=int_v_list.dp[0][0][1][0]*t1;
t1=int_v_list.dp[0][0][0][0]*int_v_oo2zeta12;
t3=t1+t2;
t1=int_v_W0-int_v_p120;
t2=int_v_list.dp[0][0][2][0]*t1;
t4=t1*t2;
t5=t4+t3;
int_v_list.dp[2][0][0][5]=t5;
t4=int_v_W2-int_v_p122;
t6=t4*t2;
int_v_list.dp[2][0][0][4]=t6;
t7=int_v_W1-int_v_p121;
t8=t2*t7;
int_v_list.dp[2][0][0][3]=t8;
t2=int_v_list.dp[0][0][2][0]*t4;
t9=t4*t2;
t2=t3+t9;
int_v_list.dp[2][0][0][2]=t2;
t9=int_v_list.dp[0][0][2][0]*t7;
t10=t4*t9;
int_v_list.dp[2][0][0][1]=t10;
t11=t7*t9;
t9=t3+t11;
int_v_list.dp[2][0][0][0]=t9;
t3=int_v_list.dp[0][0][1][0]*t1;
int_v_list.dp[1][0][0][2]=t3;
t1=int_v_list.dp[0][0][1][0]*t4;
int_v_list.dp[1][0][0][1]=t1;
t4=int_v_list.dp[0][0][1][0]*t7;
int_v_list.dp[1][0][0][0]=t4;
return 1;}
