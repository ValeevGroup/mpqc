#include <chemistry/qc/oint3/build.h>
int BuildIntV3::i1300eAB(){
/* the cost is 75 */
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
double t14;
double t15;
double t16;
double t17;
double t18;
double t19;
double t20;
double t21;
double t22;
t1=int_v_W0-int_v_p120;
t2=t1*int_v_list(0,0,2)[0];
t3=int_v_ooze*2;
t4=int_v_zeta34*t3;
t3=int_v_oo2zeta12*t4;
t4=(-1)*t3;
t3=t4*t2;
t5=t1*int_v_list(0,0,1)[0];
int_v_list(1,0,0)[2]=t5;
t6=int_v_oo2zeta12*2;
t7=t6*t5;
t8=t7+t3;
t3=int_v_zeta34*int_v_ooze;
t7=int_v_oo2zeta12*t3;
t3=(-1)*t7;
t7=t3*int_v_list(0,0,2)[0];
t9=int_v_oo2zeta12*int_v_list(0,0,1)[0];
t10=t9+t7;
t7=t1*int_v_list(0,0,3)[0];
t9=t1*t7;
t11=t9+t10;
t9=t1*t11;
t12=t9+t8;
int_v_list(3,0,0)[9]=t12;
t8=int_v_W2-int_v_p122;
t9=t8*t11;
int_v_list(3,0,0)[8]=t9;
t13=int_v_W1-int_v_p121;
t14=t11*t13;
int_v_list(3,0,0)[7]=t14;
t11=t3*t2;
t15=int_v_oo2zeta12*t5;
t5=t15+t11;
t11=t8*t7;
t15=t8*t11;
t11=t15+t5;
int_v_list(3,0,0)[6]=t11;
t15=t13*t7;
t7=t8*t15;
int_v_list(3,0,0)[5]=t7;
t16=t13*t15;
t15=t5+t16;
int_v_list(3,0,0)[4]=t15;
t5=t8*int_v_list(0,0,2)[0];
t16=t4*t5;
t17=t8*int_v_list(0,0,1)[0];
int_v_list(1,0,0)[1]=t17;
t18=t6*t17;
t17=t18+t16;
t16=t8*int_v_list(0,0,3)[0];
t18=t8*t16;
t16=t10+t18;
t18=t8*t16;
t16=t18+t17;
int_v_list(3,0,0)[3]=t16;
t17=t13*int_v_list(0,0,2)[0];
t18=t3*t17;
t19=t13*int_v_list(0,0,1)[0];
int_v_list(1,0,0)[0]=t19;
t20=int_v_oo2zeta12*t19;
t21=t20+t18;
t18=t13*int_v_list(0,0,3)[0];
t20=t8*t18;
t22=t8*t20;
t20=t22+t21;
int_v_list(3,0,0)[2]=t20;
t21=t13*t18;
t18=t10+t21;
t10=t8*t18;
int_v_list(3,0,0)[1]=t10;
t21=t4*t17;
t4=t6*t19;
t6=t4+t21;
t4=t13*t18;
t18=t4+t6;
int_v_list(3,0,0)[0]=t18;
t4=t3*int_v_list(0,0,1)[0];
t3=int_v_oo2zeta12*int_v_list(0,0,0)[0];
t6=t3+t4;
t3=t1*t2;
t1=t3+t6;
int_v_list(2,0,0)[5]=t1;
t3=t8*t2;
int_v_list(2,0,0)[4]=t3;
t4=t13*t2;
int_v_list(2,0,0)[3]=t4;
t2=t8*t5;
t5=t6+t2;
int_v_list(2,0,0)[2]=t5;
t2=t8*t17;
int_v_list(2,0,0)[1]=t2;
t8=t13*t17;
t13=t6+t8;
int_v_list(2,0,0)[0]=t13;
return 1;}
