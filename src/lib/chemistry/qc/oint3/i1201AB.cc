#include <chemistry/qc/oint3/build.h>
int sc::BuildIntV3::i1201eAB(){
/* the cost is 132 */
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
double t23;
double t24;
double t25;
double t26;
double t27;
double t28;
double t29;
double t30;
double t31;
double t32;
double t33;
double t34;
double t35;
double t36;
double t37;
double t38;
t1=int_v_W0-int_v_p120;
double***RESTRICT int_v_list0=int_v_list(0);
double**RESTRICT int_v_list00=int_v_list0[0];
double*RESTRICT int_v_list002=int_v_list00[2];
t2=t1*int_v_list002[0];
t3=0.5*int_v_ooze;
t4=t3*t2;
t5=int_v_W0-int_v_p340;
t6=t5*int_v_list002[0];
t7=int_v_p340-int_v_r30;
double*RESTRICT int_v_list001=int_v_list00[1];
t8=t7*int_v_list001[0];
t9=t8+t6;
t6=int_v_zeta34*int_v_ooze;
t8=int_v_oo2zeta12*t6;
t6=(-1)*t8;
t8=t6*t9;
t10=t8+t4;
t11=t5*int_v_list001[0];
double*RESTRICT int_v_list000=int_v_list00[0];
t12=t7*int_v_list000[0];
t13=t12+t11;
t11=int_v_oo2zeta12*t13;
t12=t11+t10;
t10=t3*int_v_list002[0];
double*RESTRICT int_v_list003=int_v_list00[3];
t13=t5*int_v_list003[0];
t5=t7*int_v_list002[0];
t7=t5+t13;
t5=t1*t7;
t13=t5+t10;
t5=t1*t13;
t14=t5+t12;
double***RESTRICT int_v_list2=int_v_list(2);
double**RESTRICT int_v_list21=int_v_list2[1];
double*RESTRICT int_v_list210=int_v_list21[0];
int_v_list210[17]=t14;
t5=int_v_W2-int_v_p342;
t12=t5*int_v_list002[0];
t15=int_v_p342-int_v_r32;
t16=t15*int_v_list001[0];
t17=t16+t12;
t12=t6*t17;
t16=t5*int_v_list001[0];
t18=t15*int_v_list000[0];
t19=t18+t16;
t16=int_v_oo2zeta12*t19;
t18=t16+t12;
t19=t5*int_v_list003[0];
t5=t15*int_v_list002[0];
t15=t5+t19;
t5=t1*t15;
t19=t1*t5;
t20=t19+t18;
int_v_list210[16]=t20;
t19=int_v_W1-int_v_p341;
t21=t19*int_v_list002[0];
t22=int_v_p341-int_v_r31;
t23=t22*int_v_list001[0];
t24=t23+t21;
t21=t6*t24;
t23=t19*int_v_list001[0];
t25=t22*int_v_list000[0];
t26=t25+t23;
t23=int_v_oo2zeta12*t26;
t25=t23+t21;
t26=t19*int_v_list003[0];
t19=t22*int_v_list002[0];
t22=t19+t26;
t19=t1*t22;
t26=t1*t19;
t27=t26+t25;
int_v_list210[15]=t27;
t26=int_v_W2-int_v_p122;
t28=t26*t13;
int_v_list210[14]=t28;
t29=t26*t5;
t30=t4+t29;
int_v_list210[13]=t30;
t29=t26*t19;
int_v_list210[12]=t29;
t31=int_v_W1-int_v_p121;
t32=t13*t31;
int_v_list210[11]=t32;
t13=t31*t5;
int_v_list210[10]=t13;
t5=t31*t19;
t19=t4+t5;
int_v_list210[9]=t19;
t4=t11+t8;
t5=t26*t7;
t8=t26*t5;
t5=t8+t4;
int_v_list210[8]=t5;
t8=t26*int_v_list002[0];
t11=t3*t8;
t33=t12+t11;
t11=t16+t33;
t12=t26*t15;
t16=t10+t12;
t12=t26*t16;
t16=t12+t11;
int_v_list210[7]=t16;
t11=t26*t22;
t12=t26*t11;
t11=t25+t12;
int_v_list210[6]=t11;
t12=t31*t7;
t7=t26*t12;
int_v_list210[5]=t7;
t25=t31*t15;
t15=t26*t25;
t33=t31*int_v_list002[0];
t34=t3*t33;
t35=t34+t15;
int_v_list210[4]=t35;
t15=t31*t22;
t22=t10+t15;
t10=t26*t22;
int_v_list210[3]=t10;
t15=t31*t12;
t12=t4+t15;
int_v_list210[2]=t12;
t4=t31*t25;
t15=t18+t4;
int_v_list210[1]=t15;
t4=t21+t34;
t18=t23+t4;
t4=t31*t22;
t21=t4+t18;
int_v_list210[0]=t21;
t4=t6*int_v_list001[0];
t6=int_v_oo2zeta12*int_v_list000[0];
t18=t6+t4;
t4=t1*t2;
t6=t4+t18;
double**RESTRICT int_v_list20=int_v_list2[0];
double*RESTRICT int_v_list200=int_v_list20[0];
int_v_list200[5]=t6;
t4=t26*t2;
int_v_list200[4]=t4;
t22=t31*t2;
int_v_list200[3]=t22;
t2=t26*t8;
t8=t18+t2;
int_v_list200[2]=t8;
t2=t26*t33;
int_v_list200[1]=t2;
t23=t31*t33;
t25=t18+t23;
int_v_list200[0]=t25;
t18=t3*int_v_list001[0];
t3=t1*t9;
t23=t3+t18;
double***RESTRICT int_v_list1=int_v_list(1);
double**RESTRICT int_v_list11=int_v_list1[1];
double*RESTRICT int_v_list110=int_v_list11[0];
int_v_list110[8]=t23;
t3=t1*t17;
int_v_list110[7]=t3;
t33=t1*t24;
int_v_list110[6]=t33;
t34=t26*t9;
int_v_list110[5]=t34;
t36=t26*t17;
t37=t18+t36;
int_v_list110[4]=t37;
t36=t26*t24;
int_v_list110[3]=t36;
t38=t31*t9;
int_v_list110[2]=t38;
t9=t31*t17;
int_v_list110[1]=t9;
t17=t31*t24;
t24=t18+t17;
int_v_list110[0]=t24;
t17=t1*int_v_list001[0];
double**RESTRICT int_v_list10=int_v_list1[0];
double*RESTRICT int_v_list100=int_v_list10[0];
int_v_list100[2]=t17;
t1=t26*int_v_list001[0];
int_v_list100[1]=t1;
t18=t31*int_v_list001[0];
int_v_list100[0]=t18;
return 1;}