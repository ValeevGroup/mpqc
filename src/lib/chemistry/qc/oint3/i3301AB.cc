#include <chemistry/qc/oint3/build.h>
int sc::BuildIntV3::i3301eAB(){
/* the cost is 295 */
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
double t39;
double t40;
double t41;
double t42;
double t43;
double t44;
double t45;
double t46;
double t47;
double t48;
double t49;
double t50;
double t51;
double t52;
double t53;
double t54;
double t55;
double t56;
t1=int_v_zeta34*int_v_ooze;
t2=int_v_oo2zeta12*t1;
t1=(-1)*t2;
double***RESTRICT int_v_list0=int_v_list(0);
double**RESTRICT int_v_list00=int_v_list0[0];
double*RESTRICT int_v_list002=int_v_list00[2];
t2=t1*int_v_list002[0];
double*RESTRICT int_v_list001=int_v_list00[1];
t3=int_v_oo2zeta12*int_v_list001[0];
t4=t3+t2;
t2=int_v_W0-int_v_p120;
double*RESTRICT int_v_list003=int_v_list00[3];
t3=t2*int_v_list003[0];
t5=t2*t3;
t6=t5+t4;
t5=0.5*int_v_ooze;
t7=t5*t6;
t8=t5*int_v_list002[0];
t9=int_v_W0-int_v_p340;
t10=t9*int_v_list003[0];
t11=int_v_p340-int_v_r30;
t12=t11*int_v_list002[0];
t13=t12+t10;
t10=t2*t13;
t12=t10+t8;
t10=2*int_v_ooze;
t14=int_v_zeta34*t10;
t10=int_v_oo2zeta12*t14;
t14=(-1)*t10;
t10=t14*t12;
t15=t10+t7;
t10=t5*int_v_list001[0];
t16=t9*int_v_list002[0];
t17=t11*int_v_list001[0];
t18=t17+t16;
t16=t2*t18;
t17=t16+t10;
t16=int_v_oo2zeta12*2;
t19=t16*t17;
t20=t19+t15;
t15=t5*t3;
t19=t1*t13;
t21=t19+t15;
t22=int_v_oo2zeta12*t18;
t23=t22+t21;
t21=t5*int_v_list003[0];
double*RESTRICT int_v_list004=int_v_list00[4];
t24=t9*int_v_list004[0];
t9=t11*int_v_list003[0];
t11=t9+t24;
t9=t2*t11;
t24=t9+t21;
t9=t2*t24;
t25=t9+t23;
t9=t2*t25;
t23=t9+t20;
double***RESTRICT int_v_list3=int_v_list(3);
double**RESTRICT int_v_list31=int_v_list3[1];
double*RESTRICT int_v_list310=int_v_list31[0];
int_v_list310[29]=t23;
t9=int_v_W2-int_v_p342;
t20=t9*int_v_list003[0];
t26=int_v_p342-int_v_r32;
t27=t26*int_v_list002[0];
t28=t27+t20;
t20=t2*t28;
t27=t14*t20;
t29=t9*int_v_list002[0];
t30=t26*int_v_list001[0];
t31=t30+t29;
t29=t2*t31;
t30=t16*t29;
t32=t30+t27;
t27=t1*t28;
t30=int_v_oo2zeta12*t31;
t33=t30+t27;
t34=t9*int_v_list004[0];
t9=t26*int_v_list003[0];
t26=t9+t34;
t9=t2*t26;
t34=t2*t9;
t35=t34+t33;
t34=t2*t35;
t36=t34+t32;
int_v_list310[28]=t36;
t32=int_v_W1-int_v_p341;
t34=t32*int_v_list003[0];
t37=int_v_p341-int_v_r31;
t38=t37*int_v_list002[0];
t39=t38+t34;
t34=t2*t39;
t38=t14*t34;
t40=t32*int_v_list002[0];
t41=t37*int_v_list001[0];
t42=t41+t40;
t40=t2*t42;
t41=t16*t40;
t43=t41+t38;
t38=t1*t39;
t41=int_v_oo2zeta12*t42;
t44=t41+t38;
t45=t32*int_v_list004[0];
t32=t37*int_v_list003[0];
t37=t32+t45;
t32=t2*t37;
t45=t2*t32;
t46=t45+t44;
t45=t2*t46;
t47=t45+t43;
int_v_list310[27]=t47;
t43=int_v_W2-int_v_p122;
t45=t43*t25;
int_v_list310[26]=t45;
t48=t43*t35;
t49=t7+t48;
int_v_list310[25]=t49;
t48=t43*t46;
int_v_list310[24]=t48;
t50=int_v_W1-int_v_p121;
t51=t25*t50;
int_v_list310[23]=t51;
t25=t50*t35;
int_v_list310[22]=t25;
t35=t50*t46;
t46=t7+t35;
int_v_list310[21]=t46;
t7=t1*t12;
t12=int_v_oo2zeta12*t17;
t17=t12+t7;
t7=t43*t24;
t12=t43*t7;
t7=t12+t17;
int_v_list310[20]=t7;
t12=t43*t3;
t35=t5*t12;
t52=t1*t20;
t20=t52+t35;
t35=int_v_oo2zeta12*t29;
t29=t35+t20;
t20=t43*t9;
t53=t15+t20;
t20=t43*t53;
t53=t20+t29;
int_v_list310[19]=t53;
t20=t1*t34;
t29=int_v_oo2zeta12*t40;
t34=t29+t20;
t40=t43*t32;
t54=t43*t40;
t40=t54+t34;
int_v_list310[18]=t40;
t34=t50*t24;
t24=t43*t34;
int_v_list310[17]=t24;
t54=t50*t3;
t3=t5*t54;
t55=t50*t9;
t9=t43*t55;
t56=t9+t3;
int_v_list310[16]=t56;
t9=t50*t32;
t32=t15+t9;
t9=t43*t32;
int_v_list310[15]=t9;
t15=t50*t34;
t34=t17+t15;
int_v_list310[14]=t34;
t15=t35+t52;
t17=t50*t55;
t35=t17+t15;
int_v_list310[13]=t35;
t15=t20+t3;
t3=t29+t15;
t15=t50*t32;
t17=t15+t3;
int_v_list310[12]=t17;
t3=t43*t13;
t15=t14*t3;
t3=t43*t18;
t20=t16*t3;
t3=t20+t15;
t15=t22+t19;
t19=t43*t11;
t20=t43*t19;
t19=t20+t15;
t20=t43*t19;
t19=t20+t3;
int_v_list310[11]=t19;
t3=t43*t28;
t20=t8+t3;
t3=t14*t20;
t20=t43*int_v_list003[0];
t22=t43*t20;
t29=t4+t22;
t22=t5*t29;
t32=t22+t3;
t3=t43*t31;
t22=t10+t3;
t3=t16*t22;
t22=t3+t32;
t3=t5*t20;
t20=t27+t3;
t3=t30+t20;
t20=t43*t26;
t27=t21+t20;
t20=t43*t27;
t27=t20+t3;
t3=t43*t27;
t20=t3+t22;
int_v_list310[10]=t20;
t3=t43*t39;
t22=t14*t3;
t3=t43*t42;
t27=t16*t3;
t3=t27+t22;
t22=t43*t37;
t27=t43*t22;
t22=t44+t27;
t27=t43*t22;
t22=t27+t3;
int_v_list310[9]=t22;
t3=t50*t13;
t13=t1*t3;
t27=t50*t18;
t18=int_v_oo2zeta12*t27;
t30=t18+t13;
t13=t50*t11;
t11=t43*t13;
t18=t43*t11;
t11=t18+t30;
int_v_list310[8]=t11;
t18=t50*t28;
t28=t1*t18;
t30=t50*int_v_list003[0];
t32=t43*t30;
t44=t5*t32;
t52=t44+t28;
t28=t50*t31;
t31=int_v_oo2zeta12*t28;
t44=t31+t52;
t31=t50*t26;
t26=t43*t31;
t52=t5*t30;
t55=t52+t26;
t26=t43*t55;
t55=t26+t44;
int_v_list310[7]=t55;
t26=t50*t39;
t39=t8+t26;
t8=t1*t39;
t26=t50*t42;
t42=t10+t26;
t10=int_v_oo2zeta12*t42;
t26=t10+t8;
t8=t50*t37;
t10=t21+t8;
t8=t43*t10;
t21=t43*t8;
t8=t21+t26;
int_v_list310[6]=t8;
t21=t50*t13;
t13=t15+t21;
t15=t43*t13;
int_v_list310[5]=t15;
t21=t50*t31;
t26=t33+t21;
t21=t43*t26;
t31=t50*t30;
t30=t4+t31;
t4=t5*t30;
t5=t4+t21;
int_v_list310[4]=t5;
t21=t38+t52;
t31=t41+t21;
t21=t50*t10;
t10=t21+t31;
t21=t43*t10;
int_v_list310[3]=t21;
t31=t14*t3;
t3=t16*t27;
t27=t3+t31;
t3=t50*t13;
t13=t3+t27;
int_v_list310[2]=t13;
t3=t14*t18;
t18=t16*t28;
t27=t18+t3;
t3=t50*t26;
t18=t3+t27;
int_v_list310[1]=t18;
t3=t14*t39;
t26=t4+t3;
t3=t16*t42;
t4=t3+t26;
t3=t50*t10;
t10=t3+t4;
int_v_list310[0]=t10;
t3=t2*int_v_list002[0];
t4=t14*t3;
t26=t2*int_v_list001[0];
t27=t16*t26;
t28=t27+t4;
t4=t2*t6;
t2=t4+t28;
double**RESTRICT int_v_list30=int_v_list3[0];
double*RESTRICT int_v_list300=int_v_list30[0];
int_v_list300[9]=t2;
t4=t43*t6;
int_v_list300[8]=t4;
t27=t50*t6;
int_v_list300[7]=t27;
t6=t1*t3;
t3=int_v_oo2zeta12*t26;
t26=t3+t6;
t3=t43*t12;
t6=t3+t26;
int_v_list300[6]=t6;
t3=t43*t54;
int_v_list300[5]=t3;
t12=t50*t54;
t28=t26+t12;
int_v_list300[4]=t28;
t12=t43*int_v_list002[0];
t26=t14*t12;
t12=t43*int_v_list001[0];
t31=t16*t12;
t12=t31+t26;
t26=t43*t29;
t29=t26+t12;
int_v_list300[3]=t29;
t12=t50*int_v_list002[0];
t26=t1*t12;
t1=t50*int_v_list001[0];
t31=int_v_oo2zeta12*t1;
t33=t31+t26;
t26=t43*t32;
t31=t26+t33;
int_v_list300[2]=t31;
t26=t43*t30;
int_v_list300[1]=t26;
t32=t14*t12;
t12=t16*t1;
t1=t12+t32;
t12=t50*t30;
t14=t12+t1;
int_v_list300[0]=t14;
return 1;}