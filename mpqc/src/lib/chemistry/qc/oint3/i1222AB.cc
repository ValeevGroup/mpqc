#include <chemistry/qc/oint3/build.h>
int sc::BuildIntV3::i1222eAB(){
/* the cost is 318 */
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
double t57;
double t58;
double t59;
double t60;
double t61;
double t62;
double t63;
double t64;
double t65;
t1=0.5*int_v_ooze;
double***int_v_list0=int_v_list(0);
double**int_v_list00=int_v_list0[0];
double*int_v_list002=int_v_list00[2];
t2=t1*int_v_list002[0];
t3=int_v_W0-int_v_p340;
double*int_v_list003=int_v_list00[3];
t4=t3*int_v_list003[0];
t5=int_v_p340-int_v_r30;
t6=t5*int_v_list002[0];
t7=t6+t4;
t4=int_v_W0-int_v_p120;
t6=t4*t7;
t8=t6+t2;
t6=int_v_ooze*2;
t9=0.5*t6;
t6=t9*t8;
t10=int_v_zeta12*int_v_ooze;
t11=int_v_oo2zeta34*t10;
t10=t11*(-1);
t11=t10*int_v_list002[0];
double*int_v_list001=int_v_list00[1];
t12=int_v_oo2zeta34*int_v_list001[0];
t13=t12+t11;
t11=t3*t7;
t12=t11+t13;
t11=t3*int_v_list002[0];
t14=t5*int_v_list001[0];
t15=t14+t11;
t11=t5*t15;
t14=t11+t12;
t11=int_v_zeta34*int_v_ooze;
t12=int_v_oo2zeta12*t11;
t11=(-1)*t12;
t12=t11*t14;
t16=t12+t6;
t6=t10*int_v_list001[0];
double*int_v_list000=int_v_list00[0];
t17=int_v_oo2zeta34*int_v_list000[0];
t18=t17+t6;
t6=t3*t15;
t17=t6+t18;
t6=t3*int_v_list001[0];
t19=t5*int_v_list000[0];
t20=t19+t6;
t6=t5*t20;
t19=t6+t17;
t6=int_v_oo2zeta12*t19;
t17=t6+t16;
t16=t9*t7;
t19=t10*int_v_list003[0];
t10=int_v_oo2zeta34*int_v_list002[0];
t21=t10+t19;
double*int_v_list004=int_v_list00[4];
t10=t3*int_v_list004[0];
t19=t5*int_v_list003[0];
t22=t19+t10;
t10=t3*t22;
t3=t10+t21;
t10=t5*t7;
t5=t10+t3;
t3=t4*t5;
t10=t3+t16;
t3=t4*t10;
t16=t3+t17;
double***int_v_list2=int_v_list(2);
double**int_v_list22=int_v_list2[2];
double*int_v_list220=int_v_list22[0];
int_v_list220[35]=t16;
t3=int_v_W2-int_v_p342;
t17=t3*int_v_list003[0];
t19=int_v_p342-int_v_r32;
t23=t19*int_v_list002[0];
t24=t23+t17;
t17=t4*t24;
t23=t1*t17;
t25=t3*t7;
t26=t19*t15;
t27=t26+t25;
t25=t11*t27;
t26=t25+t23;
t28=t3*t15;
t29=t19*t20;
t30=t29+t28;
t28=int_v_oo2zeta12*t30;
t29=t28+t26;
t26=t1*t24;
t30=t3*t22;
t31=t19*t7;
t32=t31+t30;
t30=t4*t32;
t31=t30+t26;
t30=t4*t31;
t33=t30+t29;
int_v_list220[34]=t33;
t29=int_v_W1-int_v_p341;
t30=t29*int_v_list003[0];
t34=int_v_p341-int_v_r31;
t35=t34*int_v_list002[0];
t36=t35+t30;
t30=t4*t36;
t35=t1*t30;
t37=t29*t7;
t38=t34*t15;
t39=t38+t37;
t37=t11*t39;
t38=t37+t35;
t40=t29*t15;
t41=t34*t20;
t20=t41+t40;
t40=int_v_oo2zeta12*t20;
t20=t40+t38;
t38=t1*t36;
t41=t29*t22;
t22=t34*t7;
t42=t22+t41;
t22=t4*t42;
t41=t22+t38;
t22=t4*t41;
t43=t22+t20;
int_v_list220[33]=t43;
t20=t3*t24;
t22=t13+t20;
t20=t3*int_v_list002[0];
t44=t19*int_v_list001[0];
t45=t44+t20;
t20=t19*t45;
t44=t20+t22;
t20=t11*t44;
t22=t3*t45;
t46=t18+t22;
t22=t3*int_v_list001[0];
t47=t19*int_v_list000[0];
t48=t47+t22;
t22=t19*t48;
t47=t22+t46;
t22=int_v_oo2zeta12*t47;
t46=t22+t20;
t47=t3*int_v_list004[0];
t48=t19*int_v_list003[0];
t49=t48+t47;
t47=t3*t49;
t48=t21+t47;
t47=t19*t24;
t49=t47+t48;
t47=t4*t49;
t48=t4*t47;
t50=t48+t46;
int_v_list220[32]=t50;
t48=t3*t36;
t51=t34*int_v_list001[0];
t52=t29*int_v_list002[0];
t53=t52+t51;
t51=t19*t53;
t52=t51+t48;
t48=t11*t52;
t51=t3*t53;
t54=t29*int_v_list001[0];
t55=t34*int_v_list000[0];
t56=t55+t54;
t54=t19*t56;
t55=t54+t51;
t51=int_v_oo2zeta12*t55;
t54=t51+t48;
t55=t29*int_v_list004[0];
t57=t34*int_v_list003[0];
t58=t57+t55;
t55=t3*t58;
t3=t19*t36;
t19=t3+t55;
t3=t4*t19;
t55=t4*t3;
t57=t55+t54;
int_v_list220[31]=t57;
t54=t29*t36;
t55=t13+t54;
t13=t34*t53;
t54=t13+t55;
t13=t11*t54;
t11=t29*t53;
t55=t18+t11;
t11=t34*t56;
t18=t11+t55;
t11=int_v_oo2zeta12*t18;
t18=t11+t13;
t55=t29*t58;
t29=t21+t55;
t21=t34*t36;
t34=t21+t29;
t21=t4*t34;
t29=t4*t21;
t55=t29+t18;
int_v_list220[30]=t55;
t29=int_v_W2-int_v_p122;
t56=t29*t10;
int_v_list220[29]=t56;
t58=t1*t8;
t8=t29*t31;
t59=t8+t58;
int_v_list220[28]=t59;
t8=t29*t41;
int_v_list220[27]=t8;
t60=t9*t17;
t17=t29*t47;
t61=t17+t60;
int_v_list220[26]=t61;
t17=t29*t3;
t60=t35+t17;
int_v_list220[25]=t60;
t17=t29*t21;
int_v_list220[24]=t17;
t35=int_v_W1-int_v_p121;
t62=t10*t35;
int_v_list220[23]=t62;
t10=t35*t31;
int_v_list220[22]=t10;
t31=t35*t41;
t41=t58+t31;
int_v_list220[21]=t41;
t31=t35*t47;
int_v_list220[20]=t31;
t47=t35*t3;
t3=t23+t47;
int_v_list220[19]=t3;
t23=t9*t30;
t30=t35*t21;
t21=t30+t23;
int_v_list220[18]=t21;
t23=t6+t12;
t6=t29*t5;
t12=t29*t6;
t6=t12+t23;
int_v_list220[17]=t6;
t12=t29*t7;
t30=t1*t12;
t12=t25+t30;
t30=t28+t12;
t12=t29*t32;
t47=t1*t7;
t58=t47+t12;
t12=t29*t58;
t58=t12+t30;
int_v_list220[16]=t58;
t12=t40+t37;
t30=t29*t42;
t63=t29*t30;
t30=t63+t12;
int_v_list220[15]=t30;
t12=t29*t24;
t63=t2+t12;
t12=t9*t63;
t63=t20+t12;
t12=t22+t63;
t20=t9*t24;
t22=t29*t49;
t63=t22+t20;
t20=t29*t63;
t22=t20+t12;
int_v_list220[14]=t22;
t12=t29*t36;
t20=t1*t12;
t12=t48+t20;
t20=t51+t12;
t12=t29*t19;
t63=t38+t12;
t12=t29*t63;
t38=t12+t20;
int_v_list220[13]=t38;
t12=t29*t34;
t20=t29*t12;
t12=t18+t20;
int_v_list220[12]=t12;
t18=t35*t5;
t5=t29*t18;
int_v_list220[11]=t5;
t20=t35*t32;
t32=t29*t20;
t63=t35*t7;
t7=t1*t63;
t63=t7+t32;
int_v_list220[10]=t63;
t32=t35*t42;
t42=t47+t32;
t32=t29*t42;
int_v_list220[9]=t32;
t47=t35*t24;
t24=t9*t47;
t64=t35*t49;
t49=t29*t64;
t65=t49+t24;
int_v_list220[8]=t65;
t24=t35*t36;
t49=t2+t24;
t2=t1*t49;
t24=t35*t19;
t19=t26+t24;
t24=t29*t19;
t26=t24+t2;
int_v_list220[7]=t26;
t2=t9*t36;
t24=t35*t34;
t34=t24+t2;
t2=t29*t34;
int_v_list220[6]=t2;
t24=t35*t18;
t18=t23+t24;
int_v_list220[5]=t18;
t23=t28+t25;
t24=t35*t20;
t20=t24+t23;
int_v_list220[4]=t20;
t23=t37+t7;
t7=t40+t23;
t23=t35*t42;
t24=t23+t7;
int_v_list220[3]=t24;
t7=t35*t64;
t23=t46+t7;
int_v_list220[2]=t23;
t7=t1*t47;
t25=t48+t7;
t7=t51+t25;
t25=t35*t19;
t19=t25+t7;
int_v_list220[1]=t19;
t7=t9*t49;
t25=t13+t7;
t7=t11+t25;
t11=t35*t34;
t13=t11+t7;
int_v_list220[0]=t13;
t7=t9*t15;
t11=t4*t14;
t25=t11+t7;
double***int_v_list1=int_v_list(1);
double**int_v_list12=int_v_list1[2];
double*int_v_list120=int_v_list12[0];
int_v_list120[17]=t25;
t7=t1*t45;
t11=t4*t27;
t28=t11+t7;
int_v_list120[16]=t28;
t11=t1*t53;
t34=t4*t39;
t36=t34+t11;
int_v_list120[15]=t36;
t34=t4*t44;
int_v_list120[14]=t34;
t37=t4*t52;
int_v_list120[13]=t37;
t40=t4*t54;
int_v_list120[12]=t40;
t4=t29*t14;
int_v_list120[11]=t4;
t42=t29*t27;
t46=t1*t15;
t1=t46+t42;
int_v_list120[10]=t1;
t15=t29*t39;
int_v_list120[9]=t15;
t42=t9*t45;
t45=t29*t44;
t47=t45+t42;
int_v_list120[8]=t47;
t42=t29*t52;
t45=t11+t42;
int_v_list120[7]=t45;
t11=t29*t54;
int_v_list120[6]=t11;
t29=t35*t14;
int_v_list120[5]=t29;
t14=t35*t27;
int_v_list120[4]=t14;
t27=t35*t39;
t39=t46+t27;
int_v_list120[3]=t39;
t27=t35*t44;
int_v_list120[2]=t27;
t42=t35*t52;
t44=t7+t42;
int_v_list120[1]=t44;
t7=t9*t53;
t9=t35*t54;
t35=t9+t7;
int_v_list120[0]=t35;
return 1;}
