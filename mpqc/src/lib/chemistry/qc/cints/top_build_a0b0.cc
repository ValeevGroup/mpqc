
#include <chemistry/qc/cints/int2jf.h>

void
TwoBodyIntJF::top_build_p000(iclass *H_Cl, iclass *V_Cl) /* type = 2 */
{
  int m;
  double u00 = U[0][0];
  double u01 = U[0][1];
  double u02 = U[0][2];
  double u40 = U[4][0];
  double u41 = U[4][1];
  double u42 = U[4][2];
  double fm, fmp1;
  double *vp = H_Cl->Val;

  m = H_Cl->m;
  fm = F[m];
  fmp1 = F[m+1];

  *(vp) += u00*fm + u40*fmp1;
  *(vp+1) += u01*fm + u41*fmp1;
  *(vp+2) += u02*fm + u42*fmp1;
}

void
TwoBodyIntJF::top_build_00p0(iclass *H_Cl, iclass *V_Cl) /* type = 3 */
{
  double *u2 = &(U[2][0]);
  double *u5 = &(U[5][0]);
  double *vp = H_Cl->Val;
  int m, mp1;
  int i;

  m = H_Cl->m;
  mp1 = m+1;

  for(i=0;i<3;i++)
    *(vp+i) += u2[i]*F[m] + u5[i]*F[mp1];
   
}

void
TwoBodyIntJF::top_build_d000(iclass *H_Cl, iclass *V_Cl) /* type = 4 */
{
  double *I0, *I1;
  double K1;

  I0 = build_p000(V_Cl, H_Cl->operands[0]);
  I1 = build_p000(V_Cl, H_Cl->operands[1]);

  K1 = oo2z*(F[H_Cl->m]-poz*F[H_Cl->m+1]);

  H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0] + K1;
  H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1];
  H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2];
  H_Cl->Val[3] += U[0][1]*I0[1] + U[4][1]*I1[1] + K1;
  H_Cl->Val[4] += U[0][1]*I0[2] + U[4][1]*I1[2];
  H_Cl->Val[5] += U[0][2]*I0[2] + U[4][2]*I1[2] + K1;

}

void
TwoBodyIntJF::top_build_00d0(iclass *H_Cl, iclass *V_Cl) /* type = 5 */
{
  double *I0, *I1;
  double K1;

  I0 = build_00p0(V_Cl, H_Cl->operands[0]);
  I1 = build_00p0(V_Cl, H_Cl->operands[1]);

  K1 = oo2n*(F[H_Cl->m]-pon*F[H_Cl->m+1]);

  H_Cl->Val[0] += U[2][0]*I0[0] + U[5][0]*I1[0] + K1;
  H_Cl->Val[1] += U[2][0]*I0[1] + U[5][0]*I1[1];
  H_Cl->Val[2] += U[2][0]*I0[2] + U[5][0]*I1[2];
  H_Cl->Val[3] += U[2][1]*I0[1] + U[5][1]*I1[1] + K1;
  H_Cl->Val[4] += U[2][1]*I0[2] + U[5][1]*I1[2];
  H_Cl->Val[5] += U[2][2]*I0[2] + U[5][2]*I1[2] + K1;

}


void
TwoBodyIntJF::top_build_p0p0(iclass *H_Cl, iclass *V_Cl) /* type = 6 */
{
  double *I0, *I1;
  double K1;

  I0 = build_00p0(V_Cl, H_Cl->operands[0]);
  I1 = build_00p0(V_Cl, H_Cl->operands[1]);
  K1 = oo2zn*F[H_Cl->m+1];

  H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0] + K1;
  H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1];
  H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2];
  H_Cl->Val[3] += U[0][1]*I0[0] + U[4][1]*I1[0];
  H_Cl->Val[4] += U[0][1]*I0[1] + U[4][1]*I1[1] + K1;
  H_Cl->Val[5] += U[0][1]*I0[2] + U[4][1]*I1[2];
  H_Cl->Val[6] += U[0][2]*I0[0] + U[4][2]*I1[0];
  H_Cl->Val[7] += U[0][2]*I0[1] + U[4][2]*I1[1];
  H_Cl->Val[8] += U[0][2]*I0[2] + U[4][2]*I1[2] + K1;

}

void
TwoBodyIntJF::top_build_f000(iclass *H_Cl, iclass *V_Cl) /* type = 7 */
{
  double *I0, *I1, *I2, *I3;
  double K1 = 2*oo2z;

  I0 = build_d000(V_Cl, H_Cl->operands[0]);
  I1 = build_d000(V_Cl, H_Cl->operands[1]);
  I2 = build_p000(V_Cl, H_Cl->operands[2]);
  I3 = build_p000(V_Cl, H_Cl->operands[3]);

  H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0] + K1*(I2[0]-poz*I3[0]);
  H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1] + oo2z*(I2[1]-poz*I3[1]);
  H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2] + oo2z*(I2[2]-poz*I3[2]);
  H_Cl->Val[3] += U[0][0]*I0[3] + U[4][0]*I1[3];
  H_Cl->Val[4] += U[0][0]*I0[4] + U[4][0]*I1[4];
  H_Cl->Val[5] += U[0][0]*I0[5] + U[4][0]*I1[5];
  H_Cl->Val[6] += U[0][1]*I0[3] + U[4][1]*I1[3] + K1*(I2[1]-poz*I3[1]);
  H_Cl->Val[7] += U[0][1]*I0[4] + U[4][1]*I1[4] + oo2z*(I2[2]-poz*I3[2]);
  H_Cl->Val[8] += U[0][1]*I0[5] + U[4][1]*I1[5];
  H_Cl->Val[9] += U[0][2]*I0[5] + U[4][2]*I1[5] + K1*(I2[2]-poz*I3[2]);

}

void
TwoBodyIntJF::top_build_00f0(iclass *H_Cl, iclass *V_Cl) /* type = 8 */
{
  double *I0, *I1, *I2, *I3;
  double K1 = 2*oo2n;

  I0 = build_00d0(V_Cl, H_Cl->operands[0]);
  I1 = build_00d0(V_Cl, H_Cl->operands[1]);
  I2 = build_00p0(V_Cl, H_Cl->operands[2]);
  I3 = build_00p0(V_Cl, H_Cl->operands[3]);

  H_Cl->Val[0] += U[2][0]*I0[0] + U[5][0]*I1[0] + K1*(I2[0]-pon*I3[0]);
  H_Cl->Val[1] += U[2][0]*I0[1] + U[5][0]*I1[1] + oo2n*(I2[1]-pon*I3[1]);
  H_Cl->Val[2] += U[2][0]*I0[2] + U[5][0]*I1[2] + oo2n*(I2[2]-pon*I3[2]);
  H_Cl->Val[3] += U[2][0]*I0[3] + U[5][0]*I1[3];
  H_Cl->Val[4] += U[2][0]*I0[4] + U[5][0]*I1[4];
  H_Cl->Val[5] += U[2][0]*I0[5] + U[5][0]*I1[5];
  H_Cl->Val[6] += U[2][1]*I0[3] + U[5][1]*I1[3] + K1*(I2[1]-pon*I3[1]);
  H_Cl->Val[7] += U[2][1]*I0[4] + U[5][1]*I1[4] + oo2n*(I2[2]-pon*I3[2]);
  H_Cl->Val[8] += U[2][1]*I0[5] + U[5][1]*I1[5];
  H_Cl->Val[9] += U[2][2]*I0[5] + U[5][2]*I1[5] + K1*(I2[2]-pon*I3[2]);

}

void
TwoBodyIntJF::top_build_d0p0(iclass *H_Cl, iclass *V_Cl) /* type = 9 */
{
  double *I0, *I1, *I2, *I3, *I4;

  I0 = build_p0p0(V_Cl, H_Cl->operands[0]);
  I1 = build_p0p0(V_Cl, H_Cl->operands[1]);
  I2 = build_00p0(V_Cl, H_Cl->operands[2]);
  I3 = build_00p0(V_Cl, H_Cl->operands[3]);
  I4 = build_p000(V_Cl, H_Cl->operands[4]);

H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0] + oo2z*(I2[0] - poz*I3[0]) + oo2zn*I4[0];
H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1] + oo2z*(I2[1] - poz*I3[1]);
H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2] + oo2z*(I2[2] - poz*I3[2]);
H_Cl->Val[3] += U[0][0]*I0[3] + U[4][0]*I1[3] + oo2zn*I4[1];
H_Cl->Val[4] += U[0][0]*I0[4] + U[4][0]*I1[4];
H_Cl->Val[5] += U[0][0]*I0[5] + U[4][0]*I1[5];
H_Cl->Val[6] += U[0][0]*I0[6] + U[4][0]*I1[6] + oo2zn*I4[2];
H_Cl->Val[7] += U[0][0]*I0[7] + U[4][0]*I1[7];
H_Cl->Val[8] += U[0][0]*I0[8] + U[4][0]*I1[8];
H_Cl->Val[9] += U[0][1]*I0[3] + U[4][1]*I1[3] + oo2z*(I2[0] - poz*I3[0]);
H_Cl->Val[10] += U[0][1]*I0[4] + U[4][1]*I1[4] + oo2z*(I2[1] - poz*I3[1]) + oo2zn*I4[1];
H_Cl->Val[11] += U[0][1]*I0[5] + U[4][1]*I1[5] + oo2z*(I2[2] - poz*I3[2]);
H_Cl->Val[12] += U[0][1]*I0[6] + U[4][1]*I1[6];
H_Cl->Val[13] += U[0][1]*I0[7] + U[4][1]*I1[7] + oo2zn*I4[2];
H_Cl->Val[14] += U[0][1]*I0[8] + U[4][1]*I1[8];
H_Cl->Val[15] += U[0][2]*I0[6] + U[4][2]*I1[6] + oo2z*(I2[0] - poz*I3[0]);
H_Cl->Val[16] += U[0][2]*I0[7] + U[4][2]*I1[7] + oo2z*(I2[1] - poz*I3[1]);
H_Cl->Val[17] += U[0][2]*I0[8] + U[4][2]*I1[8] + oo2z*(I2[2] - poz*I3[2]) + oo2zn*I4[2];
}

void
TwoBodyIntJF::top_build_p0d0(iclass *H_Cl, iclass *V_Cl) /* type = 10 */
{
  double *I0, *I1, *I4;
  double K1 = 2.0*oo2zn;

  I0 = build_00d0(V_Cl, H_Cl->operands[0]);
  I1 = build_00d0(V_Cl, H_Cl->operands[1]);
  I4 = build_p000(V_Cl, H_Cl->operands[4]);

H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0] + K1*I4[0];
H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1] + oo2zn*I4[1];
H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2] + oo2zn*I4[2];
H_Cl->Val[3] += U[0][0]*I0[3] + U[4][0]*I1[3];
H_Cl->Val[4] += U[0][0]*I0[4] + U[4][0]*I1[4];
H_Cl->Val[5] += U[0][0]*I0[5] + U[4][0]*I1[5];
H_Cl->Val[6] += U[0][1]*I0[0] + U[4][1]*I1[0];
H_Cl->Val[7] += U[0][1]*I0[1] + U[4][1]*I1[1] + oo2zn*I4[0];
H_Cl->Val[8] += U[0][1]*I0[2] + U[4][1]*I1[2];
H_Cl->Val[9] += U[0][1]*I0[3] + U[4][1]*I1[3] + K1*I4[1];
H_Cl->Val[10] += U[0][1]*I0[4] + U[4][1]*I1[4] + oo2zn*I4[2];
H_Cl->Val[11] += U[0][1]*I0[5] + U[4][1]*I1[5];
H_Cl->Val[12] += U[0][2]*I0[0] + U[4][2]*I1[0];
H_Cl->Val[13] += U[0][2]*I0[1] + U[4][2]*I1[1];
H_Cl->Val[14] += U[0][2]*I0[2] + U[4][2]*I1[2] + oo2zn*I4[0];
H_Cl->Val[15] += U[0][2]*I0[3] + U[4][2]*I1[3];
H_Cl->Val[16] += U[0][2]*I0[4] + U[4][2]*I1[4] + oo2zn*I4[1];
H_Cl->Val[17] += U[0][2]*I0[5] + U[4][2]*I1[5] + K1*I4[2];

}

void
TwoBodyIntJF::top_build_f0p0(iclass *H_Cl, iclass *V_Cl) /* type = 11 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double K1 = 2.0*oo2z;

  I0 = build_d0p0(V_Cl, H_Cl->operands[0]);
  I1 = build_d0p0(V_Cl, H_Cl->operands[1]);
  I2 = build_p0p0(V_Cl, H_Cl->operands[2]);
  I3 = build_p0p0(V_Cl, H_Cl->operands[3]);
  I4 = build_d000(V_Cl, H_Cl->operands[4]);

H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0] + K1*(I2[0] - poz*I3[0]) +
         oo2zn*I4[0];
H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1] + K1*(I2[1] - poz*I3[1]);
H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2] + K1*(I2[2] - poz*I3[2]);
H_Cl->Val[3] += U[0][0]*I0[3] + U[4][0]*I1[3] + oo2z*(I2[3] - poz*I3[3]) +
         oo2zn*I4[1];
H_Cl->Val[4] += U[0][0]*I0[4] + U[4][0]*I1[4] + oo2z*(I2[4] - poz*I3[4]);
H_Cl->Val[5] += U[0][0]*I0[5] + U[4][0]*I1[5] + oo2z*(I2[5] - poz*I3[5]);
H_Cl->Val[6] += U[0][0]*I0[6] + U[4][0]*I1[6] + oo2z*(I2[6] - poz*I3[6]) +
         oo2zn*I4[2];
H_Cl->Val[7] += U[0][0]*I0[7] + U[4][0]*I1[7] + oo2z*(I2[7] - poz*I3[7]);
H_Cl->Val[8] += U[0][0]*I0[8] + U[4][0]*I1[8] + oo2z*(I2[8] - poz*I3[8]);
H_Cl->Val[9] += U[0][0]*I0[9] + U[4][0]*I1[9] + oo2zn*I4[3];
H_Cl->Val[10] += U[0][0]*I0[10] + U[4][0]*I1[10];
H_Cl->Val[11] += U[0][0]*I0[11] + U[4][0]*I1[11];
H_Cl->Val[12] += U[0][0]*I0[12] + U[4][0]*I1[12] + oo2zn*I4[4];
H_Cl->Val[13] += U[0][0]*I0[13] + U[4][0]*I1[13];
H_Cl->Val[14] += U[0][0]*I0[14] + U[4][0]*I1[14];
H_Cl->Val[15] += U[0][0]*I0[15] + U[4][0]*I1[15] + oo2zn*I4[5];
H_Cl->Val[16] += U[0][0]*I0[16] + U[4][0]*I1[16];
H_Cl->Val[17] += U[0][0]*I0[17] + U[4][0]*I1[17];
H_Cl->Val[18] += U[0][1]*I0[9] + U[4][1]*I1[9] + K1*(I2[3] - poz*I3[3]);
H_Cl->Val[19] += U[0][1]*I0[10] + U[4][1]*I1[10] + K1*(I2[4] - poz*I3[4]) +
         oo2zn*I4[3];
H_Cl->Val[20] += U[0][1]*I0[11] + U[4][1]*I1[11] + K1*(I2[5] - poz*I3[5]);
H_Cl->Val[21] += U[0][1]*I0[12] + U[4][1]*I1[12] + oo2z*(I2[6] - poz*I3[6]);
H_Cl->Val[22] += U[0][1]*I0[13] + U[4][1]*I1[13] + oo2z*(I2[7] - poz*I3[7]) +
         oo2zn*I4[4];
H_Cl->Val[23] += U[0][1]*I0[14] + U[4][1]*I1[14] + oo2z*(I2[8] - poz*I3[8]);
H_Cl->Val[24] += U[0][1]*I0[15] + U[4][1]*I1[15];
H_Cl->Val[25] += U[0][1]*I0[16] + U[4][1]*I1[16] + oo2zn*I4[5];
H_Cl->Val[26] += U[0][1]*I0[17] + U[4][1]*I1[17];
H_Cl->Val[27] += U[0][2]*I0[15] + U[4][2]*I1[15] + K1*(I2[6] - poz*I3[6]);
H_Cl->Val[28] += U[0][2]*I0[16] + U[4][2]*I1[16] + K1*(I2[7] - poz*I3[7]);
H_Cl->Val[29] += U[0][2]*I0[17] + U[4][2]*I1[17] + K1*(I2[8] - poz*I3[8]) +
         oo2zn*I4[5];
}

void
TwoBodyIntJF::top_build_p0f0(iclass *H_Cl, iclass *V_Cl) /* type = 12 */
{
  double *I0, *I1, *I4;
  double K1 = 2.0*oo2zn;
  double K2 = 3.0*oo2zn;

  I0 = build_00f0(V_Cl, H_Cl->operands[0]);
  I1 = build_00f0(V_Cl, H_Cl->operands[1]);
  I4 = build_00d0(V_Cl, H_Cl->operands[4]);

H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0] + K2*I4[0];
H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1] + K1*I4[1];
H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2] + K1*I4[2];
H_Cl->Val[3] += U[0][0]*I0[3] + U[4][0]*I1[3] + oo2zn*I4[3];
H_Cl->Val[4] += U[0][0]*I0[4] + U[4][0]*I1[4] + oo2zn*I4[4];
H_Cl->Val[5] += U[0][0]*I0[5] + U[4][0]*I1[5] + oo2zn*I4[5];
H_Cl->Val[6] += U[0][0]*I0[6] + U[4][0]*I1[6];
H_Cl->Val[7] += U[0][0]*I0[7] + U[4][0]*I1[7];
H_Cl->Val[8] += U[0][0]*I0[8] + U[4][0]*I1[8];
H_Cl->Val[9] += U[0][0]*I0[9] + U[4][0]*I1[9];
H_Cl->Val[10] += U[0][1]*I0[0] + U[4][1]*I1[0];
H_Cl->Val[11] += U[0][1]*I0[1] + U[4][1]*I1[1] + oo2zn*I4[0];
H_Cl->Val[12] += U[0][1]*I0[2] + U[4][1]*I1[2];
H_Cl->Val[13] += U[0][1]*I0[3] + U[4][1]*I1[3] + K1*I4[1];
H_Cl->Val[14] += U[0][1]*I0[4] + U[4][1]*I1[4] + oo2zn*I4[2];
H_Cl->Val[15] += U[0][1]*I0[5] + U[4][1]*I1[5];
H_Cl->Val[16] += U[0][1]*I0[6] + U[4][1]*I1[6] + K2*I4[3];
H_Cl->Val[17] += U[0][1]*I0[7] + U[4][1]*I1[7] + K1*I4[4];
H_Cl->Val[18] += U[0][1]*I0[8] + U[4][1]*I1[8] + oo2zn*I4[5];
H_Cl->Val[19] += U[0][1]*I0[9] + U[4][1]*I1[9];
H_Cl->Val[20] += U[0][2]*I0[0] + U[4][2]*I1[0];
H_Cl->Val[21] += U[0][2]*I0[1] + U[4][2]*I1[1];
H_Cl->Val[22] += U[0][2]*I0[2] + U[4][2]*I1[2] + oo2zn*I4[0];
H_Cl->Val[23] += U[0][2]*I0[3] + U[4][2]*I1[3];
H_Cl->Val[24] += U[0][2]*I0[4] + U[4][2]*I1[4] + oo2zn*I4[1];
H_Cl->Val[25] += U[0][2]*I0[5] + U[4][2]*I1[5] + K1*I4[2];
H_Cl->Val[26] += U[0][2]*I0[6] + U[4][2]*I1[6];
H_Cl->Val[27] += U[0][2]*I0[7] + U[4][2]*I1[7] + oo2zn*I4[3];
H_Cl->Val[28] += U[0][2]*I0[8] + U[4][2]*I1[8] + K1*I4[4];
H_Cl->Val[29] += U[0][2]*I0[9] + U[4][2]*I1[9] + K2*I4[5];
}

void
TwoBodyIntJF::top_build_d0d0(iclass *H_Cl, iclass *V_Cl) /* type = 13 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double K1 = 2.0*oo2zn;

  I0 = build_p0d0(V_Cl, H_Cl->operands[0]);
  I1 = build_p0d0(V_Cl, H_Cl->operands[1]);
  I2 = build_00d0(V_Cl, H_Cl->operands[2]);
  I3 = build_00d0(V_Cl, H_Cl->operands[3]);
  I4 = build_p0p0(V_Cl, H_Cl->operands[4]);

H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0] + oo2z*(I2[0] - poz*I3[0]) +
         K1*I4[0];
H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1] + oo2z*(I2[1] - poz*I3[1]) +
         oo2zn*I4[1];
H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2] + oo2z*(I2[2] - poz*I3[2]) +
         oo2zn*I4[2];
H_Cl->Val[3] += U[0][0]*I0[3] + U[4][0]*I1[3] + oo2z*(I2[3] - poz*I3[3]);
H_Cl->Val[4] += U[0][0]*I0[4] + U[4][0]*I1[4] + oo2z*(I2[4] - poz*I3[4]);
H_Cl->Val[5] += U[0][0]*I0[5] + U[4][0]*I1[5] + oo2z*(I2[5] - poz*I3[5]);
H_Cl->Val[6] += U[0][0]*I0[6] + U[4][0]*I1[6] + K1*I4[3];
H_Cl->Val[7] += U[0][0]*I0[7] + U[4][0]*I1[7] + oo2zn*I4[4];
H_Cl->Val[8] += U[0][0]*I0[8] + U[4][0]*I1[8] + oo2zn*I4[5];
H_Cl->Val[9] += U[0][0]*I0[9] + U[4][0]*I1[9];
H_Cl->Val[10] += U[0][0]*I0[10] + U[4][0]*I1[10];
H_Cl->Val[11] += U[0][0]*I0[11] + U[4][0]*I1[11];
H_Cl->Val[12] += U[0][0]*I0[12] + U[4][0]*I1[12] + K1*I4[6];
H_Cl->Val[13] += U[0][0]*I0[13] + U[4][0]*I1[13] + oo2zn*I4[7];
H_Cl->Val[14] += U[0][0]*I0[14] + U[4][0]*I1[14] + oo2zn*I4[8];
H_Cl->Val[15] += U[0][0]*I0[15] + U[4][0]*I1[15];
H_Cl->Val[16] += U[0][0]*I0[16] + U[4][0]*I1[16];
H_Cl->Val[17] += U[0][0]*I0[17] + U[4][0]*I1[17];
H_Cl->Val[18] += U[0][1]*I0[6] + U[4][1]*I1[6] + oo2z*(I2[0] - poz*I3[0]);
H_Cl->Val[19] += U[0][1]*I0[7] + U[4][1]*I1[7] + oo2z*(I2[1] - poz*I3[1]) +
         oo2zn*I4[3];
H_Cl->Val[20] += U[0][1]*I0[8] + U[4][1]*I1[8] + oo2z*(I2[2] - poz*I3[2]);
H_Cl->Val[21] += U[0][1]*I0[9] + U[4][1]*I1[9] + oo2z*(I2[3] - poz*I3[3]) +
         K1*I4[4];
H_Cl->Val[22] += U[0][1]*I0[10] + U[4][1]*I1[10] + oo2z*(I2[4] - poz*I3[4]) +
         oo2zn*I4[5];
H_Cl->Val[23] += U[0][1]*I0[11] + U[4][1]*I1[11] + oo2z*(I2[5] - poz*I3[5]);
H_Cl->Val[24] += U[0][1]*I0[12] + U[4][1]*I1[12];
H_Cl->Val[25] += U[0][1]*I0[13] + U[4][1]*I1[13] + oo2zn*I4[6];
H_Cl->Val[26] += U[0][1]*I0[14] + U[4][1]*I1[14];
H_Cl->Val[27] += U[0][1]*I0[15] + U[4][1]*I1[15] + K1*I4[7];
H_Cl->Val[28] += U[0][1]*I0[16] + U[4][1]*I1[16] + oo2zn*I4[8];
H_Cl->Val[29] += U[0][1]*I0[17] + U[4][1]*I1[17];
H_Cl->Val[30] += U[0][2]*I0[12] + U[4][2]*I1[12] + oo2z*(I2[0] - poz*I3[0]);
H_Cl->Val[31] += U[0][2]*I0[13] + U[4][2]*I1[13] + oo2z*(I2[1] - poz*I3[1]);
H_Cl->Val[32] += U[0][2]*I0[14] + U[4][2]*I1[14] + oo2z*(I2[2] - poz*I3[2]) +
         oo2zn*I4[6];
H_Cl->Val[33] += U[0][2]*I0[15] + U[4][2]*I1[15] + oo2z*(I2[3] - poz*I3[3]);
H_Cl->Val[34] += U[0][2]*I0[16] + U[4][2]*I1[16] + oo2z*(I2[4] - poz*I3[4]) +
         oo2zn*I4[7];
H_Cl->Val[35] += U[0][2]*I0[17] + U[4][2]*I1[17] + oo2z*(I2[5] - poz*I3[5]) +
         K1*I4[8];
}


void
TwoBodyIntJF::top_build_f0d0(iclass *H_Cl, iclass *V_Cl) /* type = 14 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2zn = 2.0*oo2zn;
  double twooo2z = 2.0*oo2z;

  I0 = build_d0d0(V_Cl, H_Cl->operands[0]);
  I1 = build_d0d0(V_Cl, H_Cl->operands[1]);
  I2 = build_p0d0(V_Cl, H_Cl->operands[2]);
  I3 = build_p0d0(V_Cl, H_Cl->operands[3]);
  I4 = build_d0p0(V_Cl, H_Cl->operands[4]);


H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0]
           + (twooo2z)*(I2[0] - (poz)*I3[0])
           + (twooo2zn)*I4[0];
H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1]
           + (twooo2z)*(I2[1] - (poz)*I3[1])
           + (oo2zn)*I4[1];
H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2]
           + (twooo2z)*(I2[2] - (poz)*I3[2])
           + (oo2zn)*I4[2];
H_Cl->Val[3] += U[0][0]*I0[3] + U[4][0]*I1[3]
           + (twooo2z)*(I2[3] - (poz)*I3[3]);
H_Cl->Val[4] += U[0][0]*I0[4] + U[4][0]*I1[4]
           + (twooo2z)*(I2[4] - (poz)*I3[4]);
H_Cl->Val[5] += U[0][0]*I0[5] + U[4][0]*I1[5]
           + (twooo2z)*(I2[5] - (poz)*I3[5]);
H_Cl->Val[6] += U[0][0]*I0[6] + U[4][0]*I1[6]
           + (oo2z)*(I2[6] - (poz)*I3[6])
           + (twooo2zn)*I4[3];
H_Cl->Val[7] += U[0][0]*I0[7] + U[4][0]*I1[7]
           + (oo2z)*(I2[7] - (poz)*I3[7])
           + (oo2zn)*I4[4];
H_Cl->Val[8] += U[0][0]*I0[8] + U[4][0]*I1[8]
           + (oo2z)*(I2[8] - (poz)*I3[8])
           + (oo2zn)*I4[5];
H_Cl->Val[9] += U[0][0]*I0[9] + U[4][0]*I1[9]
           + (oo2z)*(I2[9] - (poz)*I3[9]);
H_Cl->Val[10] += U[0][0]*I0[10] + U[4][0]*I1[10]
           + (oo2z)*(I2[10] - (poz)*I3[10]);
H_Cl->Val[11] += U[0][0]*I0[11] + U[4][0]*I1[11]
           + (oo2z)*(I2[11] - (poz)*I3[11]);
H_Cl->Val[12] += U[0][0]*I0[12] + U[4][0]*I1[12]
           + (oo2z)*(I2[12] - (poz)*I3[12])
           + (twooo2zn)*I4[6];
H_Cl->Val[13] += U[0][0]*I0[13] + U[4][0]*I1[13]
           + (oo2z)*(I2[13] - (poz)*I3[13])
           + (oo2zn)*I4[7];
H_Cl->Val[14] += U[0][0]*I0[14] + U[4][0]*I1[14]
           + (oo2z)*(I2[14] - (poz)*I3[14])
           + (oo2zn)*I4[8];
H_Cl->Val[15] += U[0][0]*I0[15] + U[4][0]*I1[15]
           + (oo2z)*(I2[15] - (poz)*I3[15]);
H_Cl->Val[16] += U[0][0]*I0[16] + U[4][0]*I1[16]
           + (oo2z)*(I2[16] - (poz)*I3[16]);
H_Cl->Val[17] += U[0][0]*I0[17] + U[4][0]*I1[17]
           + (oo2z)*(I2[17] - (poz)*I3[17]);
H_Cl->Val[18] += U[0][0]*I0[18] + U[4][0]*I1[18]
           + (twooo2zn)*I4[9];
H_Cl->Val[19] += U[0][0]*I0[19] + U[4][0]*I1[19]
           + (oo2zn)*I4[10];
H_Cl->Val[20] += U[0][0]*I0[20] + U[4][0]*I1[20]
           + (oo2zn)*I4[11];
H_Cl->Val[21] += U[0][0]*I0[21] + U[4][0]*I1[21];
H_Cl->Val[22] += U[0][0]*I0[22] + U[4][0]*I1[22];
H_Cl->Val[23] += U[0][0]*I0[23] + U[4][0]*I1[23];
H_Cl->Val[24] += U[0][0]*I0[24] + U[4][0]*I1[24]
           + (twooo2zn)*I4[12];
H_Cl->Val[25] += U[0][0]*I0[25] + U[4][0]*I1[25]
           + (oo2zn)*I4[13];
H_Cl->Val[26] += U[0][0]*I0[26] + U[4][0]*I1[26]
           + (oo2zn)*I4[14];
H_Cl->Val[27] += U[0][0]*I0[27] + U[4][0]*I1[27];
H_Cl->Val[28] += U[0][0]*I0[28] + U[4][0]*I1[28];
H_Cl->Val[29] += U[0][0]*I0[29] + U[4][0]*I1[29];
H_Cl->Val[30] += U[0][0]*I0[30] + U[4][0]*I1[30]
           + (twooo2zn)*I4[15];
H_Cl->Val[31] += U[0][0]*I0[31] + U[4][0]*I1[31]
           + (oo2zn)*I4[16];
H_Cl->Val[32] += U[0][0]*I0[32] + U[4][0]*I1[32]
           + (oo2zn)*I4[17];
H_Cl->Val[33] += U[0][0]*I0[33] + U[4][0]*I1[33];
H_Cl->Val[34] += U[0][0]*I0[34] + U[4][0]*I1[34];
H_Cl->Val[35] += U[0][0]*I0[35] + U[4][0]*I1[35];
H_Cl->Val[36] += U[0][1]*I0[18] + U[4][1]*I1[18]
           + (twooo2z)*(I2[6] - (poz)*I3[6]);
H_Cl->Val[37] += U[0][1]*I0[19] + U[4][1]*I1[19]
           + (twooo2z)*(I2[7] - (poz)*I3[7])
           + (oo2zn)*I4[9];
H_Cl->Val[38] += U[0][1]*I0[20] + U[4][1]*I1[20]
           + (twooo2z)*(I2[8] - (poz)*I3[8]);
H_Cl->Val[39] += U[0][1]*I0[21] + U[4][1]*I1[21]
           + (twooo2z)*(I2[9] - (poz)*I3[9])
           + (twooo2zn)*I4[10];
H_Cl->Val[40] += U[0][1]*I0[22] + U[4][1]*I1[22]
           + (twooo2z)*(I2[10] - (poz)*I3[10])
           + (oo2zn)*I4[11];
H_Cl->Val[41] += U[0][1]*I0[23] + U[4][1]*I1[23]
           + (twooo2z)*(I2[11] - (poz)*I3[11]);
H_Cl->Val[42] += U[0][1]*I0[24] + U[4][1]*I1[24]
           + (oo2z)*(I2[12] - (poz)*I3[12]);
H_Cl->Val[43] += U[0][1]*I0[25] + U[4][1]*I1[25]
           + (oo2z)*(I2[13] - (poz)*I3[13])
           + (oo2zn)*I4[12];
H_Cl->Val[44] += U[0][1]*I0[26] + U[4][1]*I1[26]
           + (oo2z)*(I2[14] - (poz)*I3[14]);
H_Cl->Val[45] += U[0][1]*I0[27] + U[4][1]*I1[27]
           + (oo2z)*(I2[15] - (poz)*I3[15])
           + (twooo2zn)*I4[13];
H_Cl->Val[46] += U[0][1]*I0[28] + U[4][1]*I1[28]
           + (oo2z)*(I2[16] - (poz)*I3[16])
           + (oo2zn)*I4[14];
H_Cl->Val[47] += U[0][1]*I0[29] + U[4][1]*I1[29]
           + (oo2z)*(I2[17] - (poz)*I3[17]);
H_Cl->Val[48] += U[0][1]*I0[30] + U[4][1]*I1[30];
H_Cl->Val[49] += U[0][1]*I0[31] + U[4][1]*I1[31]
           + (oo2zn)*I4[15];
H_Cl->Val[50] += U[0][1]*I0[32] + U[4][1]*I1[32];
H_Cl->Val[51] += U[0][1]*I0[33] + U[4][1]*I1[33]
           + (twooo2zn)*I4[16];
H_Cl->Val[52] += U[0][1]*I0[34] + U[4][1]*I1[34]
           + (oo2zn)*I4[17];
H_Cl->Val[53] += U[0][1]*I0[35] + U[4][1]*I1[35];
H_Cl->Val[54] += U[0][2]*I0[30] + U[4][2]*I1[30]
           + (twooo2z)*(I2[12] - (poz)*I3[12]);
H_Cl->Val[55] += U[0][2]*I0[31] + U[4][2]*I1[31]
           + (twooo2z)*(I2[13] - (poz)*I3[13]);
H_Cl->Val[56] += U[0][2]*I0[32] + U[4][2]*I1[32]
           + (twooo2z)*(I2[14] - (poz)*I3[14])
           + (oo2zn)*I4[15];
H_Cl->Val[57] += U[0][2]*I0[33] + U[4][2]*I1[33]
           + (twooo2z)*(I2[15] - (poz)*I3[15]);
H_Cl->Val[58] += U[0][2]*I0[34] + U[4][2]*I1[34]
           + (twooo2z)*(I2[16] - (poz)*I3[16])
           + (oo2zn)*I4[16];
H_Cl->Val[59] += U[0][2]*I0[35] + U[4][2]*I1[35]
           + (twooo2z)*(I2[17] - (poz)*I3[17])
           + (twooo2zn)*I4[17];
}

void
TwoBodyIntJF::top_build_d0f0(iclass *H_Cl, iclass *V_Cl) /* type = 15 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2zn = 2.0*oo2zn;
  double threeoo2zn = 3.0*oo2zn;
  double twooo2z = 2.0*oo2z;

  I0 = build_p0f0(V_Cl, H_Cl->operands[0]);
  I1 = build_p0f0(V_Cl, H_Cl->operands[1]);
  I2 = build_00f0(V_Cl, H_Cl->operands[2]);
  I3 = build_00f0(V_Cl, H_Cl->operands[3]);
  I4 = build_p0d0(V_Cl, H_Cl->operands[4]);
 

H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0]
           + (oo2z)*(I2[0] - (poz)*I3[0])
           + (threeoo2zn)*I4[0];
H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1]
           + (oo2z)*(I2[1] - (poz)*I3[1])
           + (twooo2zn)*I4[1];
H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2]
           + (oo2z)*(I2[2] - (poz)*I3[2])
           + (twooo2zn)*I4[2];
H_Cl->Val[3] += U[0][0]*I0[3] + U[4][0]*I1[3]
           + (oo2z)*(I2[3] - (poz)*I3[3])
           + (oo2zn)*I4[3];
H_Cl->Val[4] += U[0][0]*I0[4] + U[4][0]*I1[4]
           + (oo2z)*(I2[4] - (poz)*I3[4])
           + (oo2zn)*I4[4];
H_Cl->Val[5] += U[0][0]*I0[5] + U[4][0]*I1[5]
           + (oo2z)*(I2[5] - (poz)*I3[5])
           + (oo2zn)*I4[5];
H_Cl->Val[6] += U[0][0]*I0[6] + U[4][0]*I1[6]
           + (oo2z)*(I2[6] - (poz)*I3[6]);
H_Cl->Val[7] += U[0][0]*I0[7] + U[4][0]*I1[7]
           + (oo2z)*(I2[7] - (poz)*I3[7]);
H_Cl->Val[8] += U[0][0]*I0[8] + U[4][0]*I1[8]
           + (oo2z)*(I2[8] - (poz)*I3[8]);
H_Cl->Val[9] += U[0][0]*I0[9] + U[4][0]*I1[9]
           + (oo2z)*(I2[9] - (poz)*I3[9]);
H_Cl->Val[10] += U[0][0]*I0[10] + U[4][0]*I1[10]
           + (threeoo2zn)*I4[6];
H_Cl->Val[11] += U[0][0]*I0[11] + U[4][0]*I1[11]
           + (twooo2zn)*I4[7];
H_Cl->Val[12] += U[0][0]*I0[12] + U[4][0]*I1[12]
           + (twooo2zn)*I4[8];
H_Cl->Val[13] += U[0][0]*I0[13] + U[4][0]*I1[13]
           + (oo2zn)*I4[9];
H_Cl->Val[14] += U[0][0]*I0[14] + U[4][0]*I1[14]
           + (oo2zn)*I4[10];
H_Cl->Val[15] += U[0][0]*I0[15] + U[4][0]*I1[15]
           + (oo2zn)*I4[11];
H_Cl->Val[16] += U[0][0]*I0[16] + U[4][0]*I1[16];
H_Cl->Val[17] += U[0][0]*I0[17] + U[4][0]*I1[17];
H_Cl->Val[18] += U[0][0]*I0[18] + U[4][0]*I1[18];
H_Cl->Val[19] += U[0][0]*I0[19] + U[4][0]*I1[19];
H_Cl->Val[20] += U[0][0]*I0[20] + U[4][0]*I1[20]
           + (threeoo2zn)*I4[12];
H_Cl->Val[21] += U[0][0]*I0[21] + U[4][0]*I1[21]
           + (twooo2zn)*I4[13];
H_Cl->Val[22] += U[0][0]*I0[22] + U[4][0]*I1[22]
           + (twooo2zn)*I4[14];
H_Cl->Val[23] += U[0][0]*I0[23] + U[4][0]*I1[23]
           + (oo2zn)*I4[15];
H_Cl->Val[24] += U[0][0]*I0[24] + U[4][0]*I1[24]
           + (oo2zn)*I4[16];
H_Cl->Val[25] += U[0][0]*I0[25] + U[4][0]*I1[25]
           + (oo2zn)*I4[17];
H_Cl->Val[26] += U[0][0]*I0[26] + U[4][0]*I1[26];
H_Cl->Val[27] += U[0][0]*I0[27] + U[4][0]*I1[27];
H_Cl->Val[28] += U[0][0]*I0[28] + U[4][0]*I1[28];
H_Cl->Val[29] += U[0][0]*I0[29] + U[4][0]*I1[29];
H_Cl->Val[30] += U[0][1]*I0[10] + U[4][1]*I1[10]
           + (oo2z)*(I2[0] - (poz)*I3[0]);
H_Cl->Val[31] += U[0][1]*I0[11] + U[4][1]*I1[11]
           + (oo2z)*(I2[1] - (poz)*I3[1])
           + (oo2zn)*I4[6];
H_Cl->Val[32] += U[0][1]*I0[12] + U[4][1]*I1[12]
           + (oo2z)*(I2[2] - (poz)*I3[2]);
H_Cl->Val[33] += U[0][1]*I0[13] + U[4][1]*I1[13]
           + (oo2z)*(I2[3] - (poz)*I3[3])
           + (twooo2zn)*I4[7];
H_Cl->Val[34] += U[0][1]*I0[14] + U[4][1]*I1[14]
           + (oo2z)*(I2[4] - (poz)*I3[4])
           + (oo2zn)*I4[8];
H_Cl->Val[35] += U[0][1]*I0[15] + U[4][1]*I1[15]
           + (oo2z)*(I2[5] - (poz)*I3[5]);
H_Cl->Val[36] += U[0][1]*I0[16] + U[4][1]*I1[16]
           + (oo2z)*(I2[6] - (poz)*I3[6])
           + (threeoo2zn)*I4[9];
H_Cl->Val[37] += U[0][1]*I0[17] + U[4][1]*I1[17]
           + (oo2z)*(I2[7] - (poz)*I3[7])
           + (twooo2zn)*I4[10];
H_Cl->Val[38] += U[0][1]*I0[18] + U[4][1]*I1[18]
           + (oo2z)*(I2[8] - (poz)*I3[8])
           + (oo2zn)*I4[11];
H_Cl->Val[39] += U[0][1]*I0[19] + U[4][1]*I1[19]
           + (oo2z)*(I2[9] - (poz)*I3[9]);
H_Cl->Val[40] += U[0][1]*I0[20] + U[4][1]*I1[20];
H_Cl->Val[41] += U[0][1]*I0[21] + U[4][1]*I1[21]
           + (oo2zn)*I4[12];
H_Cl->Val[42] += U[0][1]*I0[22] + U[4][1]*I1[22];
H_Cl->Val[43] += U[0][1]*I0[23] + U[4][1]*I1[23]
           + (twooo2zn)*I4[13];
H_Cl->Val[44] += U[0][1]*I0[24] + U[4][1]*I1[24]
           + (oo2zn)*I4[14];
H_Cl->Val[45] += U[0][1]*I0[25] + U[4][1]*I1[25];
H_Cl->Val[46] += U[0][1]*I0[26] + U[4][1]*I1[26]
           + (threeoo2zn)*I4[15];
H_Cl->Val[47] += U[0][1]*I0[27] + U[4][1]*I1[27]
           + (twooo2zn)*I4[16];
H_Cl->Val[48] += U[0][1]*I0[28] + U[4][1]*I1[28]
           + (oo2zn)*I4[17];
H_Cl->Val[49] += U[0][1]*I0[29] + U[4][1]*I1[29];
H_Cl->Val[50] += U[0][2]*I0[20] + U[4][2]*I1[20]
           + (oo2z)*(I2[0] - (poz)*I3[0]);
H_Cl->Val[51] += U[0][2]*I0[21] + U[4][2]*I1[21]
           + (oo2z)*(I2[1] - (poz)*I3[1]);
H_Cl->Val[52] += U[0][2]*I0[22] + U[4][2]*I1[22]
           + (oo2z)*(I2[2] - (poz)*I3[2])
           + (oo2zn)*I4[12];
H_Cl->Val[53] += U[0][2]*I0[23] + U[4][2]*I1[23]
           + (oo2z)*(I2[3] - (poz)*I3[3]);
H_Cl->Val[54] += U[0][2]*I0[24] + U[4][2]*I1[24]
           + (oo2z)*(I2[4] - (poz)*I3[4])
           + (oo2zn)*I4[13];
H_Cl->Val[55] += U[0][2]*I0[25] + U[4][2]*I1[25]
           + (oo2z)*(I2[5] - (poz)*I3[5])
           + (twooo2zn)*I4[14];
H_Cl->Val[56] += U[0][2]*I0[26] + U[4][2]*I1[26]
           + (oo2z)*(I2[6] - (poz)*I3[6]);
H_Cl->Val[57] += U[0][2]*I0[27] + U[4][2]*I1[27]
           + (oo2z)*(I2[7] - (poz)*I3[7])
           + (oo2zn)*I4[15];
H_Cl->Val[58] += U[0][2]*I0[28] + U[4][2]*I1[28]
           + (oo2z)*(I2[8] - (poz)*I3[8])
           + (twooo2zn)*I4[16];
H_Cl->Val[59] += U[0][2]*I0[29] + U[4][2]*I1[29]
           + (oo2z)*(I2[9] - (poz)*I3[9])
           + (threeoo2zn)*I4[17];
}

void
TwoBodyIntJF::top_build_f0f0(iclass *H_Cl, iclass *V_Cl) /* type = 16 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2zn = 2.0*oo2zn;
  double threeoo2zn = 3.0*oo2zn;
  double twooo2z = 2.0*oo2z;

  I0 = build_d0f0(V_Cl, H_Cl->operands[0]);
  I1 = build_d0f0(V_Cl, H_Cl->operands[1]);
  I2 = build_p0f0(V_Cl, H_Cl->operands[2]);
  I3 = build_p0f0(V_Cl, H_Cl->operands[3]);
  I4 = build_d0d0(V_Cl, H_Cl->operands[4]);


H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0]
           + (twooo2z)*(I2[0] - (poz)*I3[0])
           + (threeoo2zn)*I4[0];
H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1]
           + (twooo2z)*(I2[1] - (poz)*I3[1])
           + (twooo2zn)*I4[1];
H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2]
           + (twooo2z)*(I2[2] - (poz)*I3[2])
           + (twooo2zn)*I4[2];
H_Cl->Val[3] += U[0][0]*I0[3] + U[4][0]*I1[3]
           + (twooo2z)*(I2[3] - (poz)*I3[3])
           + (oo2zn)*I4[3];
H_Cl->Val[4] += U[0][0]*I0[4] + U[4][0]*I1[4]
           + (twooo2z)*(I2[4] - (poz)*I3[4])
           + (oo2zn)*I4[4];
H_Cl->Val[5] += U[0][0]*I0[5] + U[4][0]*I1[5]
           + (twooo2z)*(I2[5] - (poz)*I3[5])
           + (oo2zn)*I4[5];
H_Cl->Val[6] += U[0][0]*I0[6] + U[4][0]*I1[6]
           + (twooo2z)*(I2[6] - (poz)*I3[6]);
H_Cl->Val[7] += U[0][0]*I0[7] + U[4][0]*I1[7]
           + (twooo2z)*(I2[7] - (poz)*I3[7]);
H_Cl->Val[8] += U[0][0]*I0[8] + U[4][0]*I1[8]
           + (twooo2z)*(I2[8] - (poz)*I3[8]);
H_Cl->Val[9] += U[0][0]*I0[9] + U[4][0]*I1[9]
           + (twooo2z)*(I2[9] - (poz)*I3[9]);
H_Cl->Val[10] += U[0][0]*I0[10] + U[4][0]*I1[10]
           + (oo2z)*(I2[10] - (poz)*I3[10])
           + (threeoo2zn)*I4[6];
H_Cl->Val[11] += U[0][0]*I0[11] + U[4][0]*I1[11]
           + (oo2z)*(I2[11] - (poz)*I3[11])
           + (twooo2zn)*I4[7];
H_Cl->Val[12] += U[0][0]*I0[12] + U[4][0]*I1[12]
           + (oo2z)*(I2[12] - (poz)*I3[12])
           + (twooo2zn)*I4[8];
H_Cl->Val[13] += U[0][0]*I0[13] + U[4][0]*I1[13]
           + (oo2z)*(I2[13] - (poz)*I3[13])
           + (oo2zn)*I4[9];
H_Cl->Val[14] += U[0][0]*I0[14] + U[4][0]*I1[14]
           + (oo2z)*(I2[14] - (poz)*I3[14])
           + (oo2zn)*I4[10];
H_Cl->Val[15] += U[0][0]*I0[15] + U[4][0]*I1[15]
           + (oo2z)*(I2[15] - (poz)*I3[15])
           + (oo2zn)*I4[11];
H_Cl->Val[16] += U[0][0]*I0[16] + U[4][0]*I1[16]
           + (oo2z)*(I2[16] - (poz)*I3[16]);
H_Cl->Val[17] += U[0][0]*I0[17] + U[4][0]*I1[17]
           + (oo2z)*(I2[17] - (poz)*I3[17]);
H_Cl->Val[18] += U[0][0]*I0[18] + U[4][0]*I1[18]
           + (oo2z)*(I2[18] - (poz)*I3[18]);
H_Cl->Val[19] += U[0][0]*I0[19] + U[4][0]*I1[19]
           + (oo2z)*(I2[19] - (poz)*I3[19]);
H_Cl->Val[20] += U[0][0]*I0[20] + U[4][0]*I1[20]
           + (oo2z)*(I2[20] - (poz)*I3[20])
           + (threeoo2zn)*I4[12];
H_Cl->Val[21] += U[0][0]*I0[21] + U[4][0]*I1[21]
           + (oo2z)*(I2[21] - (poz)*I3[21])
           + (twooo2zn)*I4[13];
H_Cl->Val[22] += U[0][0]*I0[22] + U[4][0]*I1[22]
           + (oo2z)*(I2[22] - (poz)*I3[22])
           + (twooo2zn)*I4[14];
H_Cl->Val[23] += U[0][0]*I0[23] + U[4][0]*I1[23]
           + (oo2z)*(I2[23] - (poz)*I3[23])
           + (oo2zn)*I4[15];
H_Cl->Val[24] += U[0][0]*I0[24] + U[4][0]*I1[24]
           + (oo2z)*(I2[24] - (poz)*I3[24])
           + (oo2zn)*I4[16];
H_Cl->Val[25] += U[0][0]*I0[25] + U[4][0]*I1[25]
           + (oo2z)*(I2[25] - (poz)*I3[25])
           + (oo2zn)*I4[17];
H_Cl->Val[26] += U[0][0]*I0[26] + U[4][0]*I1[26]
           + (oo2z)*(I2[26] - (poz)*I3[26]);
H_Cl->Val[27] += U[0][0]*I0[27] + U[4][0]*I1[27]
           + (oo2z)*(I2[27] - (poz)*I3[27]);
H_Cl->Val[28] += U[0][0]*I0[28] + U[4][0]*I1[28]
           + (oo2z)*(I2[28] - (poz)*I3[28]);
H_Cl->Val[29] += U[0][0]*I0[29] + U[4][0]*I1[29]
           + (oo2z)*(I2[29] - (poz)*I3[29]);
H_Cl->Val[30] += U[0][0]*I0[30] + U[4][0]*I1[30]
           + (threeoo2zn)*I4[18];
H_Cl->Val[31] += U[0][0]*I0[31] + U[4][0]*I1[31]
           + (twooo2zn)*I4[19];
H_Cl->Val[32] += U[0][0]*I0[32] + U[4][0]*I1[32]
           + (twooo2zn)*I4[20];
H_Cl->Val[33] += U[0][0]*I0[33] + U[4][0]*I1[33]
           + (oo2zn)*I4[21];
H_Cl->Val[34] += U[0][0]*I0[34] + U[4][0]*I1[34]
           + (oo2zn)*I4[22];
H_Cl->Val[35] += U[0][0]*I0[35] + U[4][0]*I1[35]
           + (oo2zn)*I4[23];
H_Cl->Val[36] += U[0][0]*I0[36] + U[4][0]*I1[36];
H_Cl->Val[37] += U[0][0]*I0[37] + U[4][0]*I1[37];
H_Cl->Val[38] += U[0][0]*I0[38] + U[4][0]*I1[38];
H_Cl->Val[39] += U[0][0]*I0[39] + U[4][0]*I1[39];
H_Cl->Val[40] += U[0][0]*I0[40] + U[4][0]*I1[40]
           + (threeoo2zn)*I4[24];
H_Cl->Val[41] += U[0][0]*I0[41] + U[4][0]*I1[41]
           + (twooo2zn)*I4[25];
H_Cl->Val[42] += U[0][0]*I0[42] + U[4][0]*I1[42]
           + (twooo2zn)*I4[26];
H_Cl->Val[43] += U[0][0]*I0[43] + U[4][0]*I1[43]
           + (oo2zn)*I4[27];
H_Cl->Val[44] += U[0][0]*I0[44] + U[4][0]*I1[44]
           + (oo2zn)*I4[28];
H_Cl->Val[45] += U[0][0]*I0[45] + U[4][0]*I1[45]
           + (oo2zn)*I4[29];
H_Cl->Val[46] += U[0][0]*I0[46] + U[4][0]*I1[46];
H_Cl->Val[47] += U[0][0]*I0[47] + U[4][0]*I1[47];
H_Cl->Val[48] += U[0][0]*I0[48] + U[4][0]*I1[48];
H_Cl->Val[49] += U[0][0]*I0[49] + U[4][0]*I1[49];
H_Cl->Val[50] += U[0][0]*I0[50] + U[4][0]*I1[50]
           + (threeoo2zn)*I4[30];
H_Cl->Val[51] += U[0][0]*I0[51] + U[4][0]*I1[51]
           + (twooo2zn)*I4[31];
H_Cl->Val[52] += U[0][0]*I0[52] + U[4][0]*I1[52]
           + (twooo2zn)*I4[32];
H_Cl->Val[53] += U[0][0]*I0[53] + U[4][0]*I1[53]
           + (oo2zn)*I4[33];
H_Cl->Val[54] += U[0][0]*I0[54] + U[4][0]*I1[54]
           + (oo2zn)*I4[34];
H_Cl->Val[55] += U[0][0]*I0[55] + U[4][0]*I1[55]
           + (oo2zn)*I4[35];
H_Cl->Val[56] += U[0][0]*I0[56] + U[4][0]*I1[56];
H_Cl->Val[57] += U[0][0]*I0[57] + U[4][0]*I1[57];
H_Cl->Val[58] += U[0][0]*I0[58] + U[4][0]*I1[58];
H_Cl->Val[59] += U[0][0]*I0[59] + U[4][0]*I1[59];
H_Cl->Val[60] += U[0][1]*I0[30] + U[4][1]*I1[30]
           + (twooo2z)*(I2[10] - (poz)*I3[10]);
H_Cl->Val[61] += U[0][1]*I0[31] + U[4][1]*I1[31]
           + (twooo2z)*(I2[11] - (poz)*I3[11])
           + (oo2zn)*I4[18];
H_Cl->Val[62] += U[0][1]*I0[32] + U[4][1]*I1[32]
           + (twooo2z)*(I2[12] - (poz)*I3[12]);
H_Cl->Val[63] += U[0][1]*I0[33] + U[4][1]*I1[33]
           + (twooo2z)*(I2[13] - (poz)*I3[13])
           + (twooo2zn)*I4[19];
H_Cl->Val[64] += U[0][1]*I0[34] + U[4][1]*I1[34]
           + (twooo2z)*(I2[14] - (poz)*I3[14])
           + (oo2zn)*I4[20];
H_Cl->Val[65] += U[0][1]*I0[35] + U[4][1]*I1[35]
           + (twooo2z)*(I2[15] - (poz)*I3[15]);
H_Cl->Val[66] += U[0][1]*I0[36] + U[4][1]*I1[36]
           + (twooo2z)*(I2[16] - (poz)*I3[16])
           + (threeoo2zn)*I4[21];
H_Cl->Val[67] += U[0][1]*I0[37] + U[4][1]*I1[37]
           + (twooo2z)*(I2[17] - (poz)*I3[17])
           + (twooo2zn)*I4[22];
H_Cl->Val[68] += U[0][1]*I0[38] + U[4][1]*I1[38]
           + (twooo2z)*(I2[18] - (poz)*I3[18])
           + (oo2zn)*I4[23];
H_Cl->Val[69] += U[0][1]*I0[39] + U[4][1]*I1[39]
           + (twooo2z)*(I2[19] - (poz)*I3[19]);
H_Cl->Val[70] += U[0][1]*I0[40] + U[4][1]*I1[40]
           + (oo2z)*(I2[20] - (poz)*I3[20]);
H_Cl->Val[71] += U[0][1]*I0[41] + U[4][1]*I1[41]
           + (oo2z)*(I2[21] - (poz)*I3[21])
           + (oo2zn)*I4[24];
H_Cl->Val[72] += U[0][1]*I0[42] + U[4][1]*I1[42]
           + (oo2z)*(I2[22] - (poz)*I3[22]);
H_Cl->Val[73] += U[0][1]*I0[43] + U[4][1]*I1[43]
           + (oo2z)*(I2[23] - (poz)*I3[23])
           + (twooo2zn)*I4[25];
H_Cl->Val[74] += U[0][1]*I0[44] + U[4][1]*I1[44]
           + (oo2z)*(I2[24] - (poz)*I3[24])
           + (oo2zn)*I4[26];
H_Cl->Val[75] += U[0][1]*I0[45] + U[4][1]*I1[45]
           + (oo2z)*(I2[25] - (poz)*I3[25]);
H_Cl->Val[76] += U[0][1]*I0[46] + U[4][1]*I1[46]
           + (oo2z)*(I2[26] - (poz)*I3[26])
           + (threeoo2zn)*I4[27];
H_Cl->Val[77] += U[0][1]*I0[47] + U[4][1]*I1[47]
           + (oo2z)*(I2[27] - (poz)*I3[27])
           + (twooo2zn)*I4[28];
H_Cl->Val[78] += U[0][1]*I0[48] + U[4][1]*I1[48]
           + (oo2z)*(I2[28] - (poz)*I3[28])
           + (oo2zn)*I4[29];
H_Cl->Val[79] += U[0][1]*I0[49] + U[4][1]*I1[49]
           + (oo2z)*(I2[29] - (poz)*I3[29]);
H_Cl->Val[80] += U[0][1]*I0[50] + U[4][1]*I1[50];
H_Cl->Val[81] += U[0][1]*I0[51] + U[4][1]*I1[51]
           + (oo2zn)*I4[30];
H_Cl->Val[82] += U[0][1]*I0[52] + U[4][1]*I1[52];
H_Cl->Val[83] += U[0][1]*I0[53] + U[4][1]*I1[53]
           + (twooo2zn)*I4[31];
H_Cl->Val[84] += U[0][1]*I0[54] + U[4][1]*I1[54]
           + (oo2zn)*I4[32];
H_Cl->Val[85] += U[0][1]*I0[55] + U[4][1]*I1[55];
H_Cl->Val[86] += U[0][1]*I0[56] + U[4][1]*I1[56]
           + (threeoo2zn)*I4[33];
H_Cl->Val[87] += U[0][1]*I0[57] + U[4][1]*I1[57]
           + (twooo2zn)*I4[34];
H_Cl->Val[88] += U[0][1]*I0[58] + U[4][1]*I1[58]
           + (oo2zn)*I4[35];
H_Cl->Val[89] += U[0][1]*I0[59] + U[4][1]*I1[59];
H_Cl->Val[90] += U[0][2]*I0[50] + U[4][2]*I1[50]
           + (twooo2z)*(I2[20] - (poz)*I3[20]);
H_Cl->Val[91] += U[0][2]*I0[51] + U[4][2]*I1[51]
           + (twooo2z)*(I2[21] - (poz)*I3[21]);
H_Cl->Val[92] += U[0][2]*I0[52] + U[4][2]*I1[52]
           + (twooo2z)*(I2[22] - (poz)*I3[22])
           + (oo2zn)*I4[30];
H_Cl->Val[93] += U[0][2]*I0[53] + U[4][2]*I1[53]
           + (twooo2z)*(I2[23] - (poz)*I3[23]);
H_Cl->Val[94] += U[0][2]*I0[54] + U[4][2]*I1[54]
           + (twooo2z)*(I2[24] - (poz)*I3[24])
           + (oo2zn)*I4[31];
H_Cl->Val[95] += U[0][2]*I0[55] + U[4][2]*I1[55]
           + (twooo2z)*(I2[25] - (poz)*I3[25])
           + (twooo2zn)*I4[32];
H_Cl->Val[96] += U[0][2]*I0[56] + U[4][2]*I1[56]
           + (twooo2z)*(I2[26] - (poz)*I3[26]);
H_Cl->Val[97] += U[0][2]*I0[57] + U[4][2]*I1[57]
           + (twooo2z)*(I2[27] - (poz)*I3[27])
           + (oo2zn)*I4[33];
H_Cl->Val[98] += U[0][2]*I0[58] + U[4][2]*I1[58]
           + (twooo2z)*(I2[28] - (poz)*I3[28])
           + (twooo2zn)*I4[34];
H_Cl->Val[99] += U[0][2]*I0[59] + U[4][2]*I1[59]
           + (twooo2z)*(I2[29] - (poz)*I3[29])
           + (threeoo2zn)*I4[35];
}

void
TwoBodyIntJF::top_build_g000(iclass *H_Cl, iclass *V_Cl) /* type = 17 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2z = 2.0*oo2z;
  double threeoo2z = 3.0*oo2z;

  I0 = build_f000(V_Cl, H_Cl->operands[0]);
  I1 = build_f000(V_Cl, H_Cl->operands[1]);
  I2 = build_d000(V_Cl, H_Cl->operands[2]);
  I3 = build_d000(V_Cl, H_Cl->operands[3]);


H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0]
           + (threeoo2z)*(I2[0] - (poz)*I3[0]);
H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1]
           + (twooo2z)*(I2[1] - (poz)*I3[1]);
H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2]
           + (twooo2z)*(I2[2] - (poz)*I3[2]);
H_Cl->Val[3] += U[0][0]*I0[3] + U[4][0]*I1[3]
           + (oo2z)*(I2[3] - (poz)*I3[3]);
H_Cl->Val[4] += U[0][0]*I0[4] + U[4][0]*I1[4]
           + (oo2z)*(I2[4] - (poz)*I3[4]);
H_Cl->Val[5] += U[0][0]*I0[5] + U[4][0]*I1[5]
           + (oo2z)*(I2[5] - (poz)*I3[5]);
H_Cl->Val[6] += U[0][0]*I0[6] + U[4][0]*I1[6];
H_Cl->Val[7] += U[0][0]*I0[7] + U[4][0]*I1[7];
H_Cl->Val[8] += U[0][0]*I0[8] + U[4][0]*I1[8];
H_Cl->Val[9] += U[0][0]*I0[9] + U[4][0]*I1[9];
H_Cl->Val[10] += U[0][1]*I0[6] + U[4][1]*I1[6]
           + (threeoo2z)*(I2[3] - (poz)*I3[3]);
H_Cl->Val[11] += U[0][1]*I0[7] + U[4][1]*I1[7]
           + (twooo2z)*(I2[4] - (poz)*I3[4]);
H_Cl->Val[12] += U[0][1]*I0[8] + U[4][1]*I1[8]
           + (oo2z)*(I2[5] - (poz)*I3[5]);
H_Cl->Val[13] += U[0][1]*I0[9] + U[4][1]*I1[9];
H_Cl->Val[14] += U[0][2]*I0[9] + U[4][2]*I1[9]
           + (threeoo2z)*(I2[5] - (poz)*I3[5]);
}

void
TwoBodyIntJF::top_build_00g0(iclass *H_Cl, iclass *V_Cl) /* type = 18 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2n = 2.0*oo2n;
  double threeoo2n = 3.0*oo2n;

  I0 = build_00f0(V_Cl, H_Cl->operands[0]);
  I1 = build_00f0(V_Cl, H_Cl->operands[1]);
  I2 = build_00d0(V_Cl, H_Cl->operands[2]);
  I3 = build_00d0(V_Cl, H_Cl->operands[3]);


H_Cl->Val[0] += U[2][0]*I0[0] + U[5][0]*I1[0]
           + (threeoo2n)*(I2[0] - (pon)*I3[0]);
H_Cl->Val[1] += U[2][0]*I0[1] + U[5][0]*I1[1]
           + (twooo2n)*(I2[1] - (pon)*I3[1]);
H_Cl->Val[2] += U[2][0]*I0[2] + U[5][0]*I1[2]
           + (twooo2n)*(I2[2] - (pon)*I3[2]);
H_Cl->Val[3] += U[2][0]*I0[3] + U[5][0]*I1[3]
           + (oo2n)*(I2[3] - (pon)*I3[3]);
H_Cl->Val[4] += U[2][0]*I0[4] + U[5][0]*I1[4]
           + (oo2n)*(I2[4] - (pon)*I3[4]);
H_Cl->Val[5] += U[2][0]*I0[5] + U[5][0]*I1[5]
           + (oo2n)*(I2[5] - (pon)*I3[5]);
H_Cl->Val[6] += U[2][0]*I0[6] + U[5][0]*I1[6];
H_Cl->Val[7] += U[2][0]*I0[7] + U[5][0]*I1[7];
H_Cl->Val[8] += U[2][0]*I0[8] + U[5][0]*I1[8];
H_Cl->Val[9] += U[2][0]*I0[9] + U[5][0]*I1[9];
H_Cl->Val[10] += U[2][1]*I0[6] + U[5][1]*I1[6]
           + (threeoo2n)*(I2[3] - (pon)*I3[3]);
H_Cl->Val[11] += U[2][1]*I0[7] + U[5][1]*I1[7]
           + (twooo2n)*(I2[4] - (pon)*I3[4]);
H_Cl->Val[12] += U[2][1]*I0[8] + U[5][1]*I1[8]
           + (oo2n)*(I2[5] - (pon)*I3[5]);
H_Cl->Val[13] += U[2][1]*I0[9] + U[5][1]*I1[9];
H_Cl->Val[14] += U[2][2]*I0[9] + U[5][2]*I1[9]
           + (threeoo2n)*(I2[5] - (pon)*I3[5]);
}

void
TwoBodyIntJF::top_build_g0p0(iclass *H_Cl, iclass *V_Cl) /* type = 19 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2z = 2.0*oo2z;
  double threeoo2z = 3.0*oo2z;

  I0 = build_f0p0(V_Cl, H_Cl->operands[0]);
  I1 = build_f0p0(V_Cl, H_Cl->operands[1]);
  I2 = build_d0p0(V_Cl, H_Cl->operands[2]);
  I3 = build_d0p0(V_Cl, H_Cl->operands[3]);
  I4 = build_f000(V_Cl, H_Cl->operands[4]);


H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0]
           + (threeoo2z)*(I2[0] - (poz)*I3[0])
           + (oo2zn)*I4[0];
H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1]
           + (threeoo2z)*(I2[1] - (poz)*I3[1]);
H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2]
           + (threeoo2z)*(I2[2] - (poz)*I3[2]);
H_Cl->Val[3] += U[0][0]*I0[3] + U[4][0]*I1[3]
           + (twooo2z)*(I2[3] - (poz)*I3[3])
           + (oo2zn)*I4[1];
H_Cl->Val[4] += U[0][0]*I0[4] + U[4][0]*I1[4]
           + (twooo2z)*(I2[4] - (poz)*I3[4]);
H_Cl->Val[5] += U[0][0]*I0[5] + U[4][0]*I1[5]
           + (twooo2z)*(I2[5] - (poz)*I3[5]);
H_Cl->Val[6] += U[0][0]*I0[6] + U[4][0]*I1[6]
           + (twooo2z)*(I2[6] - (poz)*I3[6])
           + (oo2zn)*I4[2];
H_Cl->Val[7] += U[0][0]*I0[7] + U[4][0]*I1[7]
           + (twooo2z)*(I2[7] - (poz)*I3[7]);
H_Cl->Val[8] += U[0][0]*I0[8] + U[4][0]*I1[8]
           + (twooo2z)*(I2[8] - (poz)*I3[8]);
H_Cl->Val[9] += U[0][0]*I0[9] + U[4][0]*I1[9]
           + (oo2z)*(I2[9] - (poz)*I3[9])
           + (oo2zn)*I4[3];
H_Cl->Val[10] += U[0][0]*I0[10] + U[4][0]*I1[10]
           + (oo2z)*(I2[10] - (poz)*I3[10]);
H_Cl->Val[11] += U[0][0]*I0[11] + U[4][0]*I1[11]
           + (oo2z)*(I2[11] - (poz)*I3[11]);
H_Cl->Val[12] += U[0][0]*I0[12] + U[4][0]*I1[12]
           + (oo2z)*(I2[12] - (poz)*I3[12])
           + (oo2zn)*I4[4];
H_Cl->Val[13] += U[0][0]*I0[13] + U[4][0]*I1[13]
           + (oo2z)*(I2[13] - (poz)*I3[13]);
H_Cl->Val[14] += U[0][0]*I0[14] + U[4][0]*I1[14]
           + (oo2z)*(I2[14] - (poz)*I3[14]);
H_Cl->Val[15] += U[0][0]*I0[15] + U[4][0]*I1[15]
           + (oo2z)*(I2[15] - (poz)*I3[15])
           + (oo2zn)*I4[5];
H_Cl->Val[16] += U[0][0]*I0[16] + U[4][0]*I1[16]
           + (oo2z)*(I2[16] - (poz)*I3[16]);
H_Cl->Val[17] += U[0][0]*I0[17] + U[4][0]*I1[17]
           + (oo2z)*(I2[17] - (poz)*I3[17]);
H_Cl->Val[18] += U[0][0]*I0[18] + U[4][0]*I1[18]
           + (oo2zn)*I4[6];
H_Cl->Val[19] += U[0][0]*I0[19] + U[4][0]*I1[19];
H_Cl->Val[20] += U[0][0]*I0[20] + U[4][0]*I1[20];
H_Cl->Val[21] += U[0][0]*I0[21] + U[4][0]*I1[21]
           + (oo2zn)*I4[7];
H_Cl->Val[22] += U[0][0]*I0[22] + U[4][0]*I1[22];
H_Cl->Val[23] += U[0][0]*I0[23] + U[4][0]*I1[23];
H_Cl->Val[24] += U[0][0]*I0[24] + U[4][0]*I1[24]
           + (oo2zn)*I4[8];
H_Cl->Val[25] += U[0][0]*I0[25] + U[4][0]*I1[25];
H_Cl->Val[26] += U[0][0]*I0[26] + U[4][0]*I1[26];
H_Cl->Val[27] += U[0][0]*I0[27] + U[4][0]*I1[27]
           + (oo2zn)*I4[9];
H_Cl->Val[28] += U[0][0]*I0[28] + U[4][0]*I1[28];
H_Cl->Val[29] += U[0][0]*I0[29] + U[4][0]*I1[29];
H_Cl->Val[30] += U[0][1]*I0[18] + U[4][1]*I1[18]
           + (threeoo2z)*(I2[9] - (poz)*I3[9]);
H_Cl->Val[31] += U[0][1]*I0[19] + U[4][1]*I1[19]
           + (threeoo2z)*(I2[10] - (poz)*I3[10])
           + (oo2zn)*I4[6];
H_Cl->Val[32] += U[0][1]*I0[20] + U[4][1]*I1[20]
           + (threeoo2z)*(I2[11] - (poz)*I3[11]);
H_Cl->Val[33] += U[0][1]*I0[21] + U[4][1]*I1[21]
           + (twooo2z)*(I2[12] - (poz)*I3[12]);
H_Cl->Val[34] += U[0][1]*I0[22] + U[4][1]*I1[22]
           + (twooo2z)*(I2[13] - (poz)*I3[13])
           + (oo2zn)*I4[7];
H_Cl->Val[35] += U[0][1]*I0[23] + U[4][1]*I1[23]
           + (twooo2z)*(I2[14] - (poz)*I3[14]);
H_Cl->Val[36] += U[0][1]*I0[24] + U[4][1]*I1[24]
           + (oo2z)*(I2[15] - (poz)*I3[15]);
H_Cl->Val[37] += U[0][1]*I0[25] + U[4][1]*I1[25]
           + (oo2z)*(I2[16] - (poz)*I3[16])
           + (oo2zn)*I4[8];
H_Cl->Val[38] += U[0][1]*I0[26] + U[4][1]*I1[26]
           + (oo2z)*(I2[17] - (poz)*I3[17]);
H_Cl->Val[39] += U[0][1]*I0[27] + U[4][1]*I1[27];
H_Cl->Val[40] += U[0][1]*I0[28] + U[4][1]*I1[28]
           + (oo2zn)*I4[9];
H_Cl->Val[41] += U[0][1]*I0[29] + U[4][1]*I1[29];
H_Cl->Val[42] += U[0][2]*I0[27] + U[4][2]*I1[27]
           + (threeoo2z)*(I2[15] - (poz)*I3[15]);
H_Cl->Val[43] += U[0][2]*I0[28] + U[4][2]*I1[28]
           + (threeoo2z)*(I2[16] - (poz)*I3[16]);
H_Cl->Val[44] += U[0][2]*I0[29] + U[4][2]*I1[29]
           + (threeoo2z)*(I2[17] - (poz)*I3[17])
           + (oo2zn)*I4[9];
}

void
TwoBodyIntJF::top_build_p0g0(iclass *H_Cl, iclass *V_Cl) /* type = 20 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2zn = 2.0*oo2zn;
  double threeoo2zn = 3.0*oo2zn;
  double fouroo2zn = 4.0*oo2zn;

  I0 = build_00g0(V_Cl, H_Cl->operands[0]);
  I1 = build_00g0(V_Cl, H_Cl->operands[1]);
  I4 = build_00f0(V_Cl, H_Cl->operands[4]);


H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0]
           + (fouroo2zn)*I4[0];
H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1]
           + (threeoo2zn)*I4[1];
H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2]
           + (threeoo2zn)*I4[2];
H_Cl->Val[3] += U[0][0]*I0[3] + U[4][0]*I1[3]
           + (twooo2zn)*I4[3];
H_Cl->Val[4] += U[0][0]*I0[4] + U[4][0]*I1[4]
           + (twooo2zn)*I4[4];
H_Cl->Val[5] += U[0][0]*I0[5] + U[4][0]*I1[5]
           + (twooo2zn)*I4[5];
H_Cl->Val[6] += U[0][0]*I0[6] + U[4][0]*I1[6]
           + (oo2zn)*I4[6];
H_Cl->Val[7] += U[0][0]*I0[7] + U[4][0]*I1[7]
           + (oo2zn)*I4[7];
H_Cl->Val[8] += U[0][0]*I0[8] + U[4][0]*I1[8]
           + (oo2zn)*I4[8];
H_Cl->Val[9] += U[0][0]*I0[9] + U[4][0]*I1[9]
           + (oo2zn)*I4[9];
H_Cl->Val[10] += U[0][0]*I0[10] + U[4][0]*I1[10];
H_Cl->Val[11] += U[0][0]*I0[11] + U[4][0]*I1[11];
H_Cl->Val[12] += U[0][0]*I0[12] + U[4][0]*I1[12];
H_Cl->Val[13] += U[0][0]*I0[13] + U[4][0]*I1[13];
H_Cl->Val[14] += U[0][0]*I0[14] + U[4][0]*I1[14];
H_Cl->Val[15] += U[0][1]*I0[0] + U[4][1]*I1[0];
H_Cl->Val[16] += U[0][1]*I0[1] + U[4][1]*I1[1]
           + (oo2zn)*I4[0];
H_Cl->Val[17] += U[0][1]*I0[2] + U[4][1]*I1[2];
H_Cl->Val[18] += U[0][1]*I0[3] + U[4][1]*I1[3]
           + (twooo2zn)*I4[1];
H_Cl->Val[19] += U[0][1]*I0[4] + U[4][1]*I1[4]
           + (oo2zn)*I4[2];
H_Cl->Val[20] += U[0][1]*I0[5] + U[4][1]*I1[5];
H_Cl->Val[21] += U[0][1]*I0[6] + U[4][1]*I1[6]
           + (threeoo2zn)*I4[3];
H_Cl->Val[22] += U[0][1]*I0[7] + U[4][1]*I1[7]
           + (twooo2zn)*I4[4];
H_Cl->Val[23] += U[0][1]*I0[8] + U[4][1]*I1[8]
           + (oo2zn)*I4[5];
H_Cl->Val[24] += U[0][1]*I0[9] + U[4][1]*I1[9];
H_Cl->Val[25] += U[0][1]*I0[10] + U[4][1]*I1[10]
           + (fouroo2zn)*I4[6];
H_Cl->Val[26] += U[0][1]*I0[11] + U[4][1]*I1[11]
           + (threeoo2zn)*I4[7];
H_Cl->Val[27] += U[0][1]*I0[12] + U[4][1]*I1[12]
           + (twooo2zn)*I4[8];
H_Cl->Val[28] += U[0][1]*I0[13] + U[4][1]*I1[13]
           + (oo2zn)*I4[9];
H_Cl->Val[29] += U[0][1]*I0[14] + U[4][1]*I1[14];
H_Cl->Val[30] += U[0][2]*I0[0] + U[4][2]*I1[0];
H_Cl->Val[31] += U[0][2]*I0[1] + U[4][2]*I1[1];
H_Cl->Val[32] += U[0][2]*I0[2] + U[4][2]*I1[2]
           + (oo2zn)*I4[0];
H_Cl->Val[33] += U[0][2]*I0[3] + U[4][2]*I1[3];
H_Cl->Val[34] += U[0][2]*I0[4] + U[4][2]*I1[4]
           + (oo2zn)*I4[1];
H_Cl->Val[35] += U[0][2]*I0[5] + U[4][2]*I1[5]
           + (twooo2zn)*I4[2];
H_Cl->Val[36] += U[0][2]*I0[6] + U[4][2]*I1[6];
H_Cl->Val[37] += U[0][2]*I0[7] + U[4][2]*I1[7]
           + (oo2zn)*I4[3];
H_Cl->Val[38] += U[0][2]*I0[8] + U[4][2]*I1[8]
           + (twooo2zn)*I4[4];
H_Cl->Val[39] += U[0][2]*I0[9] + U[4][2]*I1[9]
           + (threeoo2zn)*I4[5];
H_Cl->Val[40] += U[0][2]*I0[10] + U[4][2]*I1[10];
H_Cl->Val[41] += U[0][2]*I0[11] + U[4][2]*I1[11]
           + (oo2zn)*I4[6];
H_Cl->Val[42] += U[0][2]*I0[12] + U[4][2]*I1[12]
           + (twooo2zn)*I4[7];
H_Cl->Val[43] += U[0][2]*I0[13] + U[4][2]*I1[13]
           + (threeoo2zn)*I4[8];
H_Cl->Val[44] += U[0][2]*I0[14] + U[4][2]*I1[14]
           + (fouroo2zn)*I4[9];
}

void
TwoBodyIntJF::top_build_d0g0(iclass *H_Cl, iclass *V_Cl) /* type = 21 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2zn = 2.0*oo2zn;
  double threeoo2zn = 3.0*oo2zn;
  double fouroo2zn = 4.0*oo2zn;

  I0 = build_p0g0(V_Cl, H_Cl->operands[0]);
  I1 = build_p0g0(V_Cl, H_Cl->operands[1]);
  I2 = build_00g0(V_Cl, H_Cl->operands[2]);
  I3 = build_00g0(V_Cl, H_Cl->operands[3]);
  I4 = build_p0f0(V_Cl, H_Cl->operands[4]);

H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0]
           + (oo2z)*(I2[0] - (poz)*I3[0])
           + (fouroo2zn)*I4[0];
H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1]
           + (oo2z)*(I2[1] - (poz)*I3[1])
           + (threeoo2zn)*I4[1];
H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2]
           + (oo2z)*(I2[2] - (poz)*I3[2])
           + (threeoo2zn)*I4[2];
H_Cl->Val[3] += U[0][0]*I0[3] + U[4][0]*I1[3]
           + (oo2z)*(I2[3] - (poz)*I3[3])
           + (twooo2zn)*I4[3];
H_Cl->Val[4] += U[0][0]*I0[4] + U[4][0]*I1[4]
           + (oo2z)*(I2[4] - (poz)*I3[4])
           + (twooo2zn)*I4[4];
H_Cl->Val[5] += U[0][0]*I0[5] + U[4][0]*I1[5]
           + (oo2z)*(I2[5] - (poz)*I3[5])
           + (twooo2zn)*I4[5];
H_Cl->Val[6] += U[0][0]*I0[6] + U[4][0]*I1[6]
           + (oo2z)*(I2[6] - (poz)*I3[6])
           + (oo2zn)*I4[6];
H_Cl->Val[7] += U[0][0]*I0[7] + U[4][0]*I1[7]
           + (oo2z)*(I2[7] - (poz)*I3[7])
           + (oo2zn)*I4[7];
H_Cl->Val[8] += U[0][0]*I0[8] + U[4][0]*I1[8]
           + (oo2z)*(I2[8] - (poz)*I3[8])
           + (oo2zn)*I4[8];
H_Cl->Val[9] += U[0][0]*I0[9] + U[4][0]*I1[9]
           + (oo2z)*(I2[9] - (poz)*I3[9])
           + (oo2zn)*I4[9];
H_Cl->Val[10] += U[0][0]*I0[10] + U[4][0]*I1[10]
           + (oo2z)*(I2[10] - (poz)*I3[10]);
H_Cl->Val[11] += U[0][0]*I0[11] + U[4][0]*I1[11]
           + (oo2z)*(I2[11] - (poz)*I3[11]);
H_Cl->Val[12] += U[0][0]*I0[12] + U[4][0]*I1[12]
           + (oo2z)*(I2[12] - (poz)*I3[12]);
H_Cl->Val[13] += U[0][0]*I0[13] + U[4][0]*I1[13]
           + (oo2z)*(I2[13] - (poz)*I3[13]);
H_Cl->Val[14] += U[0][0]*I0[14] + U[4][0]*I1[14]
           + (oo2z)*(I2[14] - (poz)*I3[14]);
H_Cl->Val[15] += U[0][0]*I0[15] + U[4][0]*I1[15]
           + (fouroo2zn)*I4[10];
H_Cl->Val[16] += U[0][0]*I0[16] + U[4][0]*I1[16]
           + (threeoo2zn)*I4[11];
H_Cl->Val[17] += U[0][0]*I0[17] + U[4][0]*I1[17]
           + (threeoo2zn)*I4[12];
H_Cl->Val[18] += U[0][0]*I0[18] + U[4][0]*I1[18]
           + (twooo2zn)*I4[13];
H_Cl->Val[19] += U[0][0]*I0[19] + U[4][0]*I1[19]
           + (twooo2zn)*I4[14];
H_Cl->Val[20] += U[0][0]*I0[20] + U[4][0]*I1[20]
           + (twooo2zn)*I4[15];
H_Cl->Val[21] += U[0][0]*I0[21] + U[4][0]*I1[21]
           + (oo2zn)*I4[16];
H_Cl->Val[22] += U[0][0]*I0[22] + U[4][0]*I1[22]
           + (oo2zn)*I4[17];
H_Cl->Val[23] += U[0][0]*I0[23] + U[4][0]*I1[23]
           + (oo2zn)*I4[18];
H_Cl->Val[24] += U[0][0]*I0[24] + U[4][0]*I1[24]
           + (oo2zn)*I4[19];
H_Cl->Val[25] += U[0][0]*I0[25] + U[4][0]*I1[25];
H_Cl->Val[26] += U[0][0]*I0[26] + U[4][0]*I1[26];
H_Cl->Val[27] += U[0][0]*I0[27] + U[4][0]*I1[27];
H_Cl->Val[28] += U[0][0]*I0[28] + U[4][0]*I1[28];
H_Cl->Val[29] += U[0][0]*I0[29] + U[4][0]*I1[29];
H_Cl->Val[30] += U[0][0]*I0[30] + U[4][0]*I1[30]
           + (fouroo2zn)*I4[20];
H_Cl->Val[31] += U[0][0]*I0[31] + U[4][0]*I1[31]
           + (threeoo2zn)*I4[21];
H_Cl->Val[32] += U[0][0]*I0[32] + U[4][0]*I1[32]
           + (threeoo2zn)*I4[22];
H_Cl->Val[33] += U[0][0]*I0[33] + U[4][0]*I1[33]
           + (twooo2zn)*I4[23];
H_Cl->Val[34] += U[0][0]*I0[34] + U[4][0]*I1[34]
           + (twooo2zn)*I4[24];
H_Cl->Val[35] += U[0][0]*I0[35] + U[4][0]*I1[35]
           + (twooo2zn)*I4[25];
H_Cl->Val[36] += U[0][0]*I0[36] + U[4][0]*I1[36]
           + (oo2zn)*I4[26];
H_Cl->Val[37] += U[0][0]*I0[37] + U[4][0]*I1[37]
           + (oo2zn)*I4[27];
H_Cl->Val[38] += U[0][0]*I0[38] + U[4][0]*I1[38]
           + (oo2zn)*I4[28];
H_Cl->Val[39] += U[0][0]*I0[39] + U[4][0]*I1[39]
           + (oo2zn)*I4[29];
H_Cl->Val[40] += U[0][0]*I0[40] + U[4][0]*I1[40];
H_Cl->Val[41] += U[0][0]*I0[41] + U[4][0]*I1[41];
H_Cl->Val[42] += U[0][0]*I0[42] + U[4][0]*I1[42];
H_Cl->Val[43] += U[0][0]*I0[43] + U[4][0]*I1[43];
H_Cl->Val[44] += U[0][0]*I0[44] + U[4][0]*I1[44];
H_Cl->Val[45] += U[0][1]*I0[15] + U[4][1]*I1[15]
           + (oo2z)*(I2[0] - (poz)*I3[0]);
H_Cl->Val[46] += U[0][1]*I0[16] + U[4][1]*I1[16]
           + (oo2z)*(I2[1] - (poz)*I3[1])
           + (oo2zn)*I4[10];
H_Cl->Val[47] += U[0][1]*I0[17] + U[4][1]*I1[17]
           + (oo2z)*(I2[2] - (poz)*I3[2]);
H_Cl->Val[48] += U[0][1]*I0[18] + U[4][1]*I1[18]
           + (oo2z)*(I2[3] - (poz)*I3[3])
           + (twooo2zn)*I4[11];
H_Cl->Val[49] += U[0][1]*I0[19] + U[4][1]*I1[19]
           + (oo2z)*(I2[4] - (poz)*I3[4])
           + (oo2zn)*I4[12];
H_Cl->Val[50] += U[0][1]*I0[20] + U[4][1]*I1[20]
           + (oo2z)*(I2[5] - (poz)*I3[5]);
H_Cl->Val[51] += U[0][1]*I0[21] + U[4][1]*I1[21]
           + (oo2z)*(I2[6] - (poz)*I3[6])
           + (threeoo2zn)*I4[13];
H_Cl->Val[52] += U[0][1]*I0[22] + U[4][1]*I1[22]
           + (oo2z)*(I2[7] - (poz)*I3[7])
           + (twooo2zn)*I4[14];
H_Cl->Val[53] += U[0][1]*I0[23] + U[4][1]*I1[23]
           + (oo2z)*(I2[8] - (poz)*I3[8])
           + (oo2zn)*I4[15];
H_Cl->Val[54] += U[0][1]*I0[24] + U[4][1]*I1[24]
           + (oo2z)*(I2[9] - (poz)*I3[9]);
H_Cl->Val[55] += U[0][1]*I0[25] + U[4][1]*I1[25]
           + (oo2z)*(I2[10] - (poz)*I3[10])
           + (fouroo2zn)*I4[16];
H_Cl->Val[56] += U[0][1]*I0[26] + U[4][1]*I1[26]
           + (oo2z)*(I2[11] - (poz)*I3[11])
           + (threeoo2zn)*I4[17];
H_Cl->Val[57] += U[0][1]*I0[27] + U[4][1]*I1[27]
           + (oo2z)*(I2[12] - (poz)*I3[12])
           + (twooo2zn)*I4[18];
H_Cl->Val[58] += U[0][1]*I0[28] + U[4][1]*I1[28]
           + (oo2z)*(I2[13] - (poz)*I3[13])
           + (oo2zn)*I4[19];
H_Cl->Val[59] += U[0][1]*I0[29] + U[4][1]*I1[29]
           + (oo2z)*(I2[14] - (poz)*I3[14]);
H_Cl->Val[60] += U[0][1]*I0[30] + U[4][1]*I1[30];
H_Cl->Val[61] += U[0][1]*I0[31] + U[4][1]*I1[31]
           + (oo2zn)*I4[20];
H_Cl->Val[62] += U[0][1]*I0[32] + U[4][1]*I1[32];
H_Cl->Val[63] += U[0][1]*I0[33] + U[4][1]*I1[33]
           + (twooo2zn)*I4[21];
H_Cl->Val[64] += U[0][1]*I0[34] + U[4][1]*I1[34]
           + (oo2zn)*I4[22];
H_Cl->Val[65] += U[0][1]*I0[35] + U[4][1]*I1[35];
H_Cl->Val[66] += U[0][1]*I0[36] + U[4][1]*I1[36]
           + (threeoo2zn)*I4[23];
H_Cl->Val[67] += U[0][1]*I0[37] + U[4][1]*I1[37]
           + (twooo2zn)*I4[24];
H_Cl->Val[68] += U[0][1]*I0[38] + U[4][1]*I1[38]
           + (oo2zn)*I4[25];
H_Cl->Val[69] += U[0][1]*I0[39] + U[4][1]*I1[39];
H_Cl->Val[70] += U[0][1]*I0[40] + U[4][1]*I1[40]
           + (fouroo2zn)*I4[26];
H_Cl->Val[71] += U[0][1]*I0[41] + U[4][1]*I1[41]
           + (threeoo2zn)*I4[27];
H_Cl->Val[72] += U[0][1]*I0[42] + U[4][1]*I1[42]
           + (twooo2zn)*I4[28];
H_Cl->Val[73] += U[0][1]*I0[43] + U[4][1]*I1[43]
           + (oo2zn)*I4[29];
H_Cl->Val[74] += U[0][1]*I0[44] + U[4][1]*I1[44];
H_Cl->Val[75] += U[0][2]*I0[30] + U[4][2]*I1[30]
           + (oo2z)*(I2[0] - (poz)*I3[0]);
H_Cl->Val[76] += U[0][2]*I0[31] + U[4][2]*I1[31]
           + (oo2z)*(I2[1] - (poz)*I3[1]);
H_Cl->Val[77] += U[0][2]*I0[32] + U[4][2]*I1[32]
           + (oo2z)*(I2[2] - (poz)*I3[2])
           + (oo2zn)*I4[20];
H_Cl->Val[78] += U[0][2]*I0[33] + U[4][2]*I1[33]
           + (oo2z)*(I2[3] - (poz)*I3[3]);
H_Cl->Val[79] += U[0][2]*I0[34] + U[4][2]*I1[34]
           + (oo2z)*(I2[4] - (poz)*I3[4])
           + (oo2zn)*I4[21];
H_Cl->Val[80] += U[0][2]*I0[35] + U[4][2]*I1[35]
           + (oo2z)*(I2[5] - (poz)*I3[5])
           + (twooo2zn)*I4[22];
H_Cl->Val[81] += U[0][2]*I0[36] + U[4][2]*I1[36]
           + (oo2z)*(I2[6] - (poz)*I3[6]);
H_Cl->Val[82] += U[0][2]*I0[37] + U[4][2]*I1[37]
           + (oo2z)*(I2[7] - (poz)*I3[7])
           + (oo2zn)*I4[23];
H_Cl->Val[83] += U[0][2]*I0[38] + U[4][2]*I1[38]
           + (oo2z)*(I2[8] - (poz)*I3[8])
           + (twooo2zn)*I4[24];
H_Cl->Val[84] += U[0][2]*I0[39] + U[4][2]*I1[39]
           + (oo2z)*(I2[9] - (poz)*I3[9])
           + (threeoo2zn)*I4[25];
H_Cl->Val[85] += U[0][2]*I0[40] + U[4][2]*I1[40]
           + (oo2z)*(I2[10] - (poz)*I3[10]);
H_Cl->Val[86] += U[0][2]*I0[41] + U[4][2]*I1[41]
           + (oo2z)*(I2[11] - (poz)*I3[11])
           + (oo2zn)*I4[26];
H_Cl->Val[87] += U[0][2]*I0[42] + U[4][2]*I1[42]
           + (oo2z)*(I2[12] - (poz)*I3[12])
           + (twooo2zn)*I4[27];
H_Cl->Val[88] += U[0][2]*I0[43] + U[4][2]*I1[43]
           + (oo2z)*(I2[13] - (poz)*I3[13])
           + (threeoo2zn)*I4[28];
H_Cl->Val[89] += U[0][2]*I0[44] + U[4][2]*I1[44]
           + (oo2z)*(I2[14] - (poz)*I3[14])
           + (fouroo2zn)*I4[29];
}

void
TwoBodyIntJF::top_build_g0d0(iclass *H_Cl, iclass *V_Cl) /* type = 22 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2zn = 2.0*oo2zn;
  double twooo2z = 2.0*oo2z;
  double threeoo2z = 3.0*oo2z;

  I0 = build_f0d0(V_Cl, H_Cl->operands[0]);
  I1 = build_f0d0(V_Cl, H_Cl->operands[1]);
  I2 = build_d0d0(V_Cl, H_Cl->operands[2]);
  I3 = build_d0d0(V_Cl, H_Cl->operands[3]);
  I4 = build_f0p0(V_Cl, H_Cl->operands[4]);


H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0]
           + (threeoo2z)*(I2[0] - (poz)*I3[0])
           + (twooo2zn)*I4[0];
H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1]
           + (threeoo2z)*(I2[1] - (poz)*I3[1])
           + (oo2zn)*I4[1];
H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2]
           + (threeoo2z)*(I2[2] - (poz)*I3[2])
           + (oo2zn)*I4[2];
H_Cl->Val[3] += U[0][0]*I0[3] + U[4][0]*I1[3]
           + (threeoo2z)*(I2[3] - (poz)*I3[3]);
H_Cl->Val[4] += U[0][0]*I0[4] + U[4][0]*I1[4]
           + (threeoo2z)*(I2[4] - (poz)*I3[4]);
H_Cl->Val[5] += U[0][0]*I0[5] + U[4][0]*I1[5]
           + (threeoo2z)*(I2[5] - (poz)*I3[5]);
H_Cl->Val[6] += U[0][0]*I0[6] + U[4][0]*I1[6]
           + (twooo2z)*(I2[6] - (poz)*I3[6])
           + (twooo2zn)*I4[3];
H_Cl->Val[7] += U[0][0]*I0[7] + U[4][0]*I1[7]
           + (twooo2z)*(I2[7] - (poz)*I3[7])
           + (oo2zn)*I4[4];
H_Cl->Val[8] += U[0][0]*I0[8] + U[4][0]*I1[8]
           + (twooo2z)*(I2[8] - (poz)*I3[8])
           + (oo2zn)*I4[5];
H_Cl->Val[9] += U[0][0]*I0[9] + U[4][0]*I1[9]
           + (twooo2z)*(I2[9] - (poz)*I3[9]);
H_Cl->Val[10] += U[0][0]*I0[10] + U[4][0]*I1[10]
           + (twooo2z)*(I2[10] - (poz)*I3[10]);
H_Cl->Val[11] += U[0][0]*I0[11] + U[4][0]*I1[11]
           + (twooo2z)*(I2[11] - (poz)*I3[11]);
H_Cl->Val[12] += U[0][0]*I0[12] + U[4][0]*I1[12]
           + (twooo2z)*(I2[12] - (poz)*I3[12])
           + (twooo2zn)*I4[6];
H_Cl->Val[13] += U[0][0]*I0[13] + U[4][0]*I1[13]
           + (twooo2z)*(I2[13] - (poz)*I3[13])
           + (oo2zn)*I4[7];
H_Cl->Val[14] += U[0][0]*I0[14] + U[4][0]*I1[14]
           + (twooo2z)*(I2[14] - (poz)*I3[14])
           + (oo2zn)*I4[8];
H_Cl->Val[15] += U[0][0]*I0[15] + U[4][0]*I1[15]
           + (twooo2z)*(I2[15] - (poz)*I3[15]);
H_Cl->Val[16] += U[0][0]*I0[16] + U[4][0]*I1[16]
           + (twooo2z)*(I2[16] - (poz)*I3[16]);
H_Cl->Val[17] += U[0][0]*I0[17] + U[4][0]*I1[17]
           + (twooo2z)*(I2[17] - (poz)*I3[17]);
H_Cl->Val[18] += U[0][0]*I0[18] + U[4][0]*I1[18]
           + (oo2z)*(I2[18] - (poz)*I3[18])
           + (twooo2zn)*I4[9];
H_Cl->Val[19] += U[0][0]*I0[19] + U[4][0]*I1[19]
           + (oo2z)*(I2[19] - (poz)*I3[19])
           + (oo2zn)*I4[10];
H_Cl->Val[20] += U[0][0]*I0[20] + U[4][0]*I1[20]
           + (oo2z)*(I2[20] - (poz)*I3[20])
           + (oo2zn)*I4[11];
H_Cl->Val[21] += U[0][0]*I0[21] + U[4][0]*I1[21]
           + (oo2z)*(I2[21] - (poz)*I3[21]);
H_Cl->Val[22] += U[0][0]*I0[22] + U[4][0]*I1[22]
           + (oo2z)*(I2[22] - (poz)*I3[22]);
H_Cl->Val[23] += U[0][0]*I0[23] + U[4][0]*I1[23]
           + (oo2z)*(I2[23] - (poz)*I3[23]);
H_Cl->Val[24] += U[0][0]*I0[24] + U[4][0]*I1[24]
           + (oo2z)*(I2[24] - (poz)*I3[24])
           + (twooo2zn)*I4[12];
H_Cl->Val[25] += U[0][0]*I0[25] + U[4][0]*I1[25]
           + (oo2z)*(I2[25] - (poz)*I3[25])
           + (oo2zn)*I4[13];
H_Cl->Val[26] += U[0][0]*I0[26] + U[4][0]*I1[26]
           + (oo2z)*(I2[26] - (poz)*I3[26])
           + (oo2zn)*I4[14];
H_Cl->Val[27] += U[0][0]*I0[27] + U[4][0]*I1[27]
           + (oo2z)*(I2[27] - (poz)*I3[27]);
H_Cl->Val[28] += U[0][0]*I0[28] + U[4][0]*I1[28]
           + (oo2z)*(I2[28] - (poz)*I3[28]);
H_Cl->Val[29] += U[0][0]*I0[29] + U[4][0]*I1[29]
           + (oo2z)*(I2[29] - (poz)*I3[29]);
H_Cl->Val[30] += U[0][0]*I0[30] + U[4][0]*I1[30]
           + (oo2z)*(I2[30] - (poz)*I3[30])
           + (twooo2zn)*I4[15];
H_Cl->Val[31] += U[0][0]*I0[31] + U[4][0]*I1[31]
           + (oo2z)*(I2[31] - (poz)*I3[31])
           + (oo2zn)*I4[16];
H_Cl->Val[32] += U[0][0]*I0[32] + U[4][0]*I1[32]
           + (oo2z)*(I2[32] - (poz)*I3[32])
           + (oo2zn)*I4[17];
H_Cl->Val[33] += U[0][0]*I0[33] + U[4][0]*I1[33]
           + (oo2z)*(I2[33] - (poz)*I3[33]);
H_Cl->Val[34] += U[0][0]*I0[34] + U[4][0]*I1[34]
           + (oo2z)*(I2[34] - (poz)*I3[34]);
H_Cl->Val[35] += U[0][0]*I0[35] + U[4][0]*I1[35]
           + (oo2z)*(I2[35] - (poz)*I3[35]);
H_Cl->Val[36] += U[0][0]*I0[36] + U[4][0]*I1[36]
           + (twooo2zn)*I4[18];
H_Cl->Val[37] += U[0][0]*I0[37] + U[4][0]*I1[37]
           + (oo2zn)*I4[19];
H_Cl->Val[38] += U[0][0]*I0[38] + U[4][0]*I1[38]
           + (oo2zn)*I4[20];
H_Cl->Val[39] += U[0][0]*I0[39] + U[4][0]*I1[39];
H_Cl->Val[40] += U[0][0]*I0[40] + U[4][0]*I1[40];
H_Cl->Val[41] += U[0][0]*I0[41] + U[4][0]*I1[41];
H_Cl->Val[42] += U[0][0]*I0[42] + U[4][0]*I1[42]
           + (twooo2zn)*I4[21];
H_Cl->Val[43] += U[0][0]*I0[43] + U[4][0]*I1[43]
           + (oo2zn)*I4[22];
H_Cl->Val[44] += U[0][0]*I0[44] + U[4][0]*I1[44]
           + (oo2zn)*I4[23];
H_Cl->Val[45] += U[0][0]*I0[45] + U[4][0]*I1[45];
H_Cl->Val[46] += U[0][0]*I0[46] + U[4][0]*I1[46];
H_Cl->Val[47] += U[0][0]*I0[47] + U[4][0]*I1[47];
H_Cl->Val[48] += U[0][0]*I0[48] + U[4][0]*I1[48]
           + (twooo2zn)*I4[24];
H_Cl->Val[49] += U[0][0]*I0[49] + U[4][0]*I1[49]
           + (oo2zn)*I4[25];
H_Cl->Val[50] += U[0][0]*I0[50] + U[4][0]*I1[50]
           + (oo2zn)*I4[26];
H_Cl->Val[51] += U[0][0]*I0[51] + U[4][0]*I1[51];
H_Cl->Val[52] += U[0][0]*I0[52] + U[4][0]*I1[52];
H_Cl->Val[53] += U[0][0]*I0[53] + U[4][0]*I1[53];
H_Cl->Val[54] += U[0][0]*I0[54] + U[4][0]*I1[54]
           + (twooo2zn)*I4[27];
H_Cl->Val[55] += U[0][0]*I0[55] + U[4][0]*I1[55]
           + (oo2zn)*I4[28];
H_Cl->Val[56] += U[0][0]*I0[56] + U[4][0]*I1[56]
           + (oo2zn)*I4[29];
H_Cl->Val[57] += U[0][0]*I0[57] + U[4][0]*I1[57];
H_Cl->Val[58] += U[0][0]*I0[58] + U[4][0]*I1[58];
H_Cl->Val[59] += U[0][0]*I0[59] + U[4][0]*I1[59];
H_Cl->Val[60] += U[0][1]*I0[36] + U[4][1]*I1[36]
           + (threeoo2z)*(I2[18] - (poz)*I3[18]);
H_Cl->Val[61] += U[0][1]*I0[37] + U[4][1]*I1[37]
           + (threeoo2z)*(I2[19] - (poz)*I3[19])
           + (oo2zn)*I4[18];
H_Cl->Val[62] += U[0][1]*I0[38] + U[4][1]*I1[38]
           + (threeoo2z)*(I2[20] - (poz)*I3[20]);
H_Cl->Val[63] += U[0][1]*I0[39] + U[4][1]*I1[39]
           + (threeoo2z)*(I2[21] - (poz)*I3[21])
           + (twooo2zn)*I4[19];
H_Cl->Val[64] += U[0][1]*I0[40] + U[4][1]*I1[40]
           + (threeoo2z)*(I2[22] - (poz)*I3[22])
           + (oo2zn)*I4[20];
H_Cl->Val[65] += U[0][1]*I0[41] + U[4][1]*I1[41]
           + (threeoo2z)*(I2[23] - (poz)*I3[23]);
H_Cl->Val[66] += U[0][1]*I0[42] + U[4][1]*I1[42]
           + (twooo2z)*(I2[24] - (poz)*I3[24]);
H_Cl->Val[67] += U[0][1]*I0[43] + U[4][1]*I1[43]
           + (twooo2z)*(I2[25] - (poz)*I3[25])
           + (oo2zn)*I4[21];
H_Cl->Val[68] += U[0][1]*I0[44] + U[4][1]*I1[44]
           + (twooo2z)*(I2[26] - (poz)*I3[26]);
H_Cl->Val[69] += U[0][1]*I0[45] + U[4][1]*I1[45]
           + (twooo2z)*(I2[27] - (poz)*I3[27])
           + (twooo2zn)*I4[22];
H_Cl->Val[70] += U[0][1]*I0[46] + U[4][1]*I1[46]
           + (twooo2z)*(I2[28] - (poz)*I3[28])
           + (oo2zn)*I4[23];
H_Cl->Val[71] += U[0][1]*I0[47] + U[4][1]*I1[47]
           + (twooo2z)*(I2[29] - (poz)*I3[29]);
H_Cl->Val[72] += U[0][1]*I0[48] + U[4][1]*I1[48]
           + (oo2z)*(I2[30] - (poz)*I3[30]);
H_Cl->Val[73] += U[0][1]*I0[49] + U[4][1]*I1[49]
           + (oo2z)*(I2[31] - (poz)*I3[31])
           + (oo2zn)*I4[24];
H_Cl->Val[74] += U[0][1]*I0[50] + U[4][1]*I1[50]
           + (oo2z)*(I2[32] - (poz)*I3[32]);
H_Cl->Val[75] += U[0][1]*I0[51] + U[4][1]*I1[51]
           + (oo2z)*(I2[33] - (poz)*I3[33])
           + (twooo2zn)*I4[25];
H_Cl->Val[76] += U[0][1]*I0[52] + U[4][1]*I1[52]
           + (oo2z)*(I2[34] - (poz)*I3[34])
           + (oo2zn)*I4[26];
H_Cl->Val[77] += U[0][1]*I0[53] + U[4][1]*I1[53]
           + (oo2z)*(I2[35] - (poz)*I3[35]);
H_Cl->Val[78] += U[0][1]*I0[54] + U[4][1]*I1[54];
H_Cl->Val[79] += U[0][1]*I0[55] + U[4][1]*I1[55]
           + (oo2zn)*I4[27];
H_Cl->Val[80] += U[0][1]*I0[56] + U[4][1]*I1[56];
H_Cl->Val[81] += U[0][1]*I0[57] + U[4][1]*I1[57]
           + (twooo2zn)*I4[28];
H_Cl->Val[82] += U[0][1]*I0[58] + U[4][1]*I1[58]
           + (oo2zn)*I4[29];
H_Cl->Val[83] += U[0][1]*I0[59] + U[4][1]*I1[59];
H_Cl->Val[84] += U[0][2]*I0[54] + U[4][2]*I1[54]
           + (threeoo2z)*(I2[30] - (poz)*I3[30]);
H_Cl->Val[85] += U[0][2]*I0[55] + U[4][2]*I1[55]
           + (threeoo2z)*(I2[31] - (poz)*I3[31]);
H_Cl->Val[86] += U[0][2]*I0[56] + U[4][2]*I1[56]
           + (threeoo2z)*(I2[32] - (poz)*I3[32])
           + (oo2zn)*I4[27];
H_Cl->Val[87] += U[0][2]*I0[57] + U[4][2]*I1[57]
           + (threeoo2z)*(I2[33] - (poz)*I3[33]);
H_Cl->Val[88] += U[0][2]*I0[58] + U[4][2]*I1[58]
           + (threeoo2z)*(I2[34] - (poz)*I3[34])
           + (oo2zn)*I4[28];
H_Cl->Val[89] += U[0][2]*I0[59] + U[4][2]*I1[59]
           + (threeoo2z)*(I2[35] - (poz)*I3[35])
           + (twooo2zn)*I4[29];
}

void
TwoBodyIntJF::top_build_g0f0(iclass *H_Cl, iclass *V_Cl) /* type = 23 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2zn = 2.0*oo2zn;
  double threeoo2zn = 3.0*oo2zn;
  double twooo2z = 2.0*oo2z;
  double threeoo2z = 3.0*oo2z;

  I0 = build_f0f0(V_Cl, H_Cl->operands[0]);
  I1 = build_f0f0(V_Cl, H_Cl->operands[1]);
  I2 = build_d0f0(V_Cl, H_Cl->operands[2]);
  I3 = build_d0f0(V_Cl, H_Cl->operands[3]);
  I4 = build_f0d0(V_Cl, H_Cl->operands[4]);


H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0]
           + (threeoo2z)*(I2[0] - (poz)*I3[0])
           + (threeoo2zn)*I4[0];
H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1]
           + (threeoo2z)*(I2[1] - (poz)*I3[1])
           + (twooo2zn)*I4[1];
H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2]
           + (threeoo2z)*(I2[2] - (poz)*I3[2])
           + (twooo2zn)*I4[2];
H_Cl->Val[3] += U[0][0]*I0[3] + U[4][0]*I1[3]
           + (threeoo2z)*(I2[3] - (poz)*I3[3])
           + (oo2zn)*I4[3];
H_Cl->Val[4] += U[0][0]*I0[4] + U[4][0]*I1[4]
           + (threeoo2z)*(I2[4] - (poz)*I3[4])
           + (oo2zn)*I4[4];
H_Cl->Val[5] += U[0][0]*I0[5] + U[4][0]*I1[5]
           + (threeoo2z)*(I2[5] - (poz)*I3[5])
           + (oo2zn)*I4[5];
H_Cl->Val[6] += U[0][0]*I0[6] + U[4][0]*I1[6]
           + (threeoo2z)*(I2[6] - (poz)*I3[6]);
H_Cl->Val[7] += U[0][0]*I0[7] + U[4][0]*I1[7]
           + (threeoo2z)*(I2[7] - (poz)*I3[7]);
H_Cl->Val[8] += U[0][0]*I0[8] + U[4][0]*I1[8]
           + (threeoo2z)*(I2[8] - (poz)*I3[8]);
H_Cl->Val[9] += U[0][0]*I0[9] + U[4][0]*I1[9]
           + (threeoo2z)*(I2[9] - (poz)*I3[9]);
H_Cl->Val[10] += U[0][0]*I0[10] + U[4][0]*I1[10]
           + (twooo2z)*(I2[10] - (poz)*I3[10])
           + (threeoo2zn)*I4[6];
H_Cl->Val[11] += U[0][0]*I0[11] + U[4][0]*I1[11]
           + (twooo2z)*(I2[11] - (poz)*I3[11])
           + (twooo2zn)*I4[7];
H_Cl->Val[12] += U[0][0]*I0[12] + U[4][0]*I1[12]
           + (twooo2z)*(I2[12] - (poz)*I3[12])
           + (twooo2zn)*I4[8];
H_Cl->Val[13] += U[0][0]*I0[13] + U[4][0]*I1[13]
           + (twooo2z)*(I2[13] - (poz)*I3[13])
           + (oo2zn)*I4[9];
H_Cl->Val[14] += U[0][0]*I0[14] + U[4][0]*I1[14]
           + (twooo2z)*(I2[14] - (poz)*I3[14])
           + (oo2zn)*I4[10];
H_Cl->Val[15] += U[0][0]*I0[15] + U[4][0]*I1[15]
           + (twooo2z)*(I2[15] - (poz)*I3[15])
           + (oo2zn)*I4[11];
H_Cl->Val[16] += U[0][0]*I0[16] + U[4][0]*I1[16]
           + (twooo2z)*(I2[16] - (poz)*I3[16]);
H_Cl->Val[17] += U[0][0]*I0[17] + U[4][0]*I1[17]
           + (twooo2z)*(I2[17] - (poz)*I3[17]);
H_Cl->Val[18] += U[0][0]*I0[18] + U[4][0]*I1[18]
           + (twooo2z)*(I2[18] - (poz)*I3[18]);
H_Cl->Val[19] += U[0][0]*I0[19] + U[4][0]*I1[19]
           + (twooo2z)*(I2[19] - (poz)*I3[19]);
H_Cl->Val[20] += U[0][0]*I0[20] + U[4][0]*I1[20]
           + (twooo2z)*(I2[20] - (poz)*I3[20])
           + (threeoo2zn)*I4[12];
H_Cl->Val[21] += U[0][0]*I0[21] + U[4][0]*I1[21]
           + (twooo2z)*(I2[21] - (poz)*I3[21])
           + (twooo2zn)*I4[13];
H_Cl->Val[22] += U[0][0]*I0[22] + U[4][0]*I1[22]
           + (twooo2z)*(I2[22] - (poz)*I3[22])
           + (twooo2zn)*I4[14];
H_Cl->Val[23] += U[0][0]*I0[23] + U[4][0]*I1[23]
           + (twooo2z)*(I2[23] - (poz)*I3[23])
           + (oo2zn)*I4[15];
H_Cl->Val[24] += U[0][0]*I0[24] + U[4][0]*I1[24]
           + (twooo2z)*(I2[24] - (poz)*I3[24])
           + (oo2zn)*I4[16];
H_Cl->Val[25] += U[0][0]*I0[25] + U[4][0]*I1[25]
           + (twooo2z)*(I2[25] - (poz)*I3[25])
           + (oo2zn)*I4[17];
H_Cl->Val[26] += U[0][0]*I0[26] + U[4][0]*I1[26]
           + (twooo2z)*(I2[26] - (poz)*I3[26]);
H_Cl->Val[27] += U[0][0]*I0[27] + U[4][0]*I1[27]
           + (twooo2z)*(I2[27] - (poz)*I3[27]);
H_Cl->Val[28] += U[0][0]*I0[28] + U[4][0]*I1[28]
           + (twooo2z)*(I2[28] - (poz)*I3[28]);
H_Cl->Val[29] += U[0][0]*I0[29] + U[4][0]*I1[29]
           + (twooo2z)*(I2[29] - (poz)*I3[29]);
H_Cl->Val[30] += U[0][0]*I0[30] + U[4][0]*I1[30]
           + (oo2z)*(I2[30] - (poz)*I3[30])
           + (threeoo2zn)*I4[18];
H_Cl->Val[31] += U[0][0]*I0[31] + U[4][0]*I1[31]
           + (oo2z)*(I2[31] - (poz)*I3[31])
           + (twooo2zn)*I4[19];
H_Cl->Val[32] += U[0][0]*I0[32] + U[4][0]*I1[32]
           + (oo2z)*(I2[32] - (poz)*I3[32])
           + (twooo2zn)*I4[20];
H_Cl->Val[33] += U[0][0]*I0[33] + U[4][0]*I1[33]
           + (oo2z)*(I2[33] - (poz)*I3[33])
           + (oo2zn)*I4[21];
H_Cl->Val[34] += U[0][0]*I0[34] + U[4][0]*I1[34]
           + (oo2z)*(I2[34] - (poz)*I3[34])
           + (oo2zn)*I4[22];
H_Cl->Val[35] += U[0][0]*I0[35] + U[4][0]*I1[35]
           + (oo2z)*(I2[35] - (poz)*I3[35])
           + (oo2zn)*I4[23];
H_Cl->Val[36] += U[0][0]*I0[36] + U[4][0]*I1[36]
           + (oo2z)*(I2[36] - (poz)*I3[36]);
H_Cl->Val[37] += U[0][0]*I0[37] + U[4][0]*I1[37]
           + (oo2z)*(I2[37] - (poz)*I3[37]);
H_Cl->Val[38] += U[0][0]*I0[38] + U[4][0]*I1[38]
           + (oo2z)*(I2[38] - (poz)*I3[38]);
H_Cl->Val[39] += U[0][0]*I0[39] + U[4][0]*I1[39]
           + (oo2z)*(I2[39] - (poz)*I3[39]);
H_Cl->Val[40] += U[0][0]*I0[40] + U[4][0]*I1[40]
           + (oo2z)*(I2[40] - (poz)*I3[40])
           + (threeoo2zn)*I4[24];
H_Cl->Val[41] += U[0][0]*I0[41] + U[4][0]*I1[41]
           + (oo2z)*(I2[41] - (poz)*I3[41])
           + (twooo2zn)*I4[25];
H_Cl->Val[42] += U[0][0]*I0[42] + U[4][0]*I1[42]
           + (oo2z)*(I2[42] - (poz)*I3[42])
           + (twooo2zn)*I4[26];
H_Cl->Val[43] += U[0][0]*I0[43] + U[4][0]*I1[43]
           + (oo2z)*(I2[43] - (poz)*I3[43])
           + (oo2zn)*I4[27];
H_Cl->Val[44] += U[0][0]*I0[44] + U[4][0]*I1[44]
           + (oo2z)*(I2[44] - (poz)*I3[44])
           + (oo2zn)*I4[28];
H_Cl->Val[45] += U[0][0]*I0[45] + U[4][0]*I1[45]
           + (oo2z)*(I2[45] - (poz)*I3[45])
           + (oo2zn)*I4[29];
H_Cl->Val[46] += U[0][0]*I0[46] + U[4][0]*I1[46]
           + (oo2z)*(I2[46] - (poz)*I3[46]);
H_Cl->Val[47] += U[0][0]*I0[47] + U[4][0]*I1[47]
           + (oo2z)*(I2[47] - (poz)*I3[47]);
H_Cl->Val[48] += U[0][0]*I0[48] + U[4][0]*I1[48]
           + (oo2z)*(I2[48] - (poz)*I3[48]);
H_Cl->Val[49] += U[0][0]*I0[49] + U[4][0]*I1[49]
           + (oo2z)*(I2[49] - (poz)*I3[49]);
H_Cl->Val[50] += U[0][0]*I0[50] + U[4][0]*I1[50]
           + (oo2z)*(I2[50] - (poz)*I3[50])
           + (threeoo2zn)*I4[30];
H_Cl->Val[51] += U[0][0]*I0[51] + U[4][0]*I1[51]
           + (oo2z)*(I2[51] - (poz)*I3[51])
           + (twooo2zn)*I4[31];
H_Cl->Val[52] += U[0][0]*I0[52] + U[4][0]*I1[52]
           + (oo2z)*(I2[52] - (poz)*I3[52])
           + (twooo2zn)*I4[32];
H_Cl->Val[53] += U[0][0]*I0[53] + U[4][0]*I1[53]
           + (oo2z)*(I2[53] - (poz)*I3[53])
           + (oo2zn)*I4[33];
H_Cl->Val[54] += U[0][0]*I0[54] + U[4][0]*I1[54]
           + (oo2z)*(I2[54] - (poz)*I3[54])
           + (oo2zn)*I4[34];
H_Cl->Val[55] += U[0][0]*I0[55] + U[4][0]*I1[55]
           + (oo2z)*(I2[55] - (poz)*I3[55])
           + (oo2zn)*I4[35];
H_Cl->Val[56] += U[0][0]*I0[56] + U[4][0]*I1[56]
           + (oo2z)*(I2[56] - (poz)*I3[56]);
H_Cl->Val[57] += U[0][0]*I0[57] + U[4][0]*I1[57]
           + (oo2z)*(I2[57] - (poz)*I3[57]);
H_Cl->Val[58] += U[0][0]*I0[58] + U[4][0]*I1[58]
           + (oo2z)*(I2[58] - (poz)*I3[58]);
H_Cl->Val[59] += U[0][0]*I0[59] + U[4][0]*I1[59]
           + (oo2z)*(I2[59] - (poz)*I3[59]);
H_Cl->Val[60] += U[0][0]*I0[60] + U[4][0]*I1[60]
           + (threeoo2zn)*I4[36];
H_Cl->Val[61] += U[0][0]*I0[61] + U[4][0]*I1[61]
           + (twooo2zn)*I4[37];
H_Cl->Val[62] += U[0][0]*I0[62] + U[4][0]*I1[62]
           + (twooo2zn)*I4[38];
H_Cl->Val[63] += U[0][0]*I0[63] + U[4][0]*I1[63]
           + (oo2zn)*I4[39];
H_Cl->Val[64] += U[0][0]*I0[64] + U[4][0]*I1[64]
           + (oo2zn)*I4[40];
H_Cl->Val[65] += U[0][0]*I0[65] + U[4][0]*I1[65]
           + (oo2zn)*I4[41];
H_Cl->Val[66] += U[0][0]*I0[66] + U[4][0]*I1[66];
H_Cl->Val[67] += U[0][0]*I0[67] + U[4][0]*I1[67];
H_Cl->Val[68] += U[0][0]*I0[68] + U[4][0]*I1[68];
H_Cl->Val[69] += U[0][0]*I0[69] + U[4][0]*I1[69];
H_Cl->Val[70] += U[0][0]*I0[70] + U[4][0]*I1[70]
           + (threeoo2zn)*I4[42];
H_Cl->Val[71] += U[0][0]*I0[71] + U[4][0]*I1[71]
           + (twooo2zn)*I4[43];
H_Cl->Val[72] += U[0][0]*I0[72] + U[4][0]*I1[72]
           + (twooo2zn)*I4[44];
H_Cl->Val[73] += U[0][0]*I0[73] + U[4][0]*I1[73]
           + (oo2zn)*I4[45];
H_Cl->Val[74] += U[0][0]*I0[74] + U[4][0]*I1[74]
           + (oo2zn)*I4[46];
H_Cl->Val[75] += U[0][0]*I0[75] + U[4][0]*I1[75]
           + (oo2zn)*I4[47];
H_Cl->Val[76] += U[0][0]*I0[76] + U[4][0]*I1[76];
H_Cl->Val[77] += U[0][0]*I0[77] + U[4][0]*I1[77];
H_Cl->Val[78] += U[0][0]*I0[78] + U[4][0]*I1[78];
H_Cl->Val[79] += U[0][0]*I0[79] + U[4][0]*I1[79];
H_Cl->Val[80] += U[0][0]*I0[80] + U[4][0]*I1[80]
           + (threeoo2zn)*I4[48];
H_Cl->Val[81] += U[0][0]*I0[81] + U[4][0]*I1[81]
           + (twooo2zn)*I4[49];
H_Cl->Val[82] += U[0][0]*I0[82] + U[4][0]*I1[82]
           + (twooo2zn)*I4[50];
H_Cl->Val[83] += U[0][0]*I0[83] + U[4][0]*I1[83]
           + (oo2zn)*I4[51];
H_Cl->Val[84] += U[0][0]*I0[84] + U[4][0]*I1[84]
           + (oo2zn)*I4[52];
H_Cl->Val[85] += U[0][0]*I0[85] + U[4][0]*I1[85]
           + (oo2zn)*I4[53];
H_Cl->Val[86] += U[0][0]*I0[86] + U[4][0]*I1[86];
H_Cl->Val[87] += U[0][0]*I0[87] + U[4][0]*I1[87];
H_Cl->Val[88] += U[0][0]*I0[88] + U[4][0]*I1[88];
H_Cl->Val[89] += U[0][0]*I0[89] + U[4][0]*I1[89];
H_Cl->Val[90] += U[0][0]*I0[90] + U[4][0]*I1[90]
           + (threeoo2zn)*I4[54];
H_Cl->Val[91] += U[0][0]*I0[91] + U[4][0]*I1[91]
           + (twooo2zn)*I4[55];
H_Cl->Val[92] += U[0][0]*I0[92] + U[4][0]*I1[92]
           + (twooo2zn)*I4[56];
H_Cl->Val[93] += U[0][0]*I0[93] + U[4][0]*I1[93]
           + (oo2zn)*I4[57];
H_Cl->Val[94] += U[0][0]*I0[94] + U[4][0]*I1[94]
           + (oo2zn)*I4[58];
H_Cl->Val[95] += U[0][0]*I0[95] + U[4][0]*I1[95]
           + (oo2zn)*I4[59];
H_Cl->Val[96] += U[0][0]*I0[96] + U[4][0]*I1[96];
H_Cl->Val[97] += U[0][0]*I0[97] + U[4][0]*I1[97];
H_Cl->Val[98] += U[0][0]*I0[98] + U[4][0]*I1[98];
H_Cl->Val[99] += U[0][0]*I0[99] + U[4][0]*I1[99];
H_Cl->Val[100] += U[0][1]*I0[60] + U[4][1]*I1[60]
           + (threeoo2z)*(I2[30] - (poz)*I3[30]);
H_Cl->Val[101] += U[0][1]*I0[61] + U[4][1]*I1[61]
           + (threeoo2z)*(I2[31] - (poz)*I3[31])
           + (oo2zn)*I4[36];
H_Cl->Val[102] += U[0][1]*I0[62] + U[4][1]*I1[62]
           + (threeoo2z)*(I2[32] - (poz)*I3[32]);
H_Cl->Val[103] += U[0][1]*I0[63] + U[4][1]*I1[63]
           + (threeoo2z)*(I2[33] - (poz)*I3[33])
           + (twooo2zn)*I4[37];
H_Cl->Val[104] += U[0][1]*I0[64] + U[4][1]*I1[64]
           + (threeoo2z)*(I2[34] - (poz)*I3[34])
           + (oo2zn)*I4[38];
H_Cl->Val[105] += U[0][1]*I0[65] + U[4][1]*I1[65]
           + (threeoo2z)*(I2[35] - (poz)*I3[35]);
H_Cl->Val[106] += U[0][1]*I0[66] + U[4][1]*I1[66]
           + (threeoo2z)*(I2[36] - (poz)*I3[36])
           + (threeoo2zn)*I4[39];
H_Cl->Val[107] += U[0][1]*I0[67] + U[4][1]*I1[67]
           + (threeoo2z)*(I2[37] - (poz)*I3[37])
           + (twooo2zn)*I4[40];
H_Cl->Val[108] += U[0][1]*I0[68] + U[4][1]*I1[68]
           + (threeoo2z)*(I2[38] - (poz)*I3[38])
           + (oo2zn)*I4[41];
H_Cl->Val[109] += U[0][1]*I0[69] + U[4][1]*I1[69]
           + (threeoo2z)*(I2[39] - (poz)*I3[39]);
H_Cl->Val[110] += U[0][1]*I0[70] + U[4][1]*I1[70]
           + (twooo2z)*(I2[40] - (poz)*I3[40]);
H_Cl->Val[111] += U[0][1]*I0[71] + U[4][1]*I1[71]
           + (twooo2z)*(I2[41] - (poz)*I3[41])
           + (oo2zn)*I4[42];
H_Cl->Val[112] += U[0][1]*I0[72] + U[4][1]*I1[72]
           + (twooo2z)*(I2[42] - (poz)*I3[42]);
H_Cl->Val[113] += U[0][1]*I0[73] + U[4][1]*I1[73]
           + (twooo2z)*(I2[43] - (poz)*I3[43])
           + (twooo2zn)*I4[43];
H_Cl->Val[114] += U[0][1]*I0[74] + U[4][1]*I1[74]
           + (twooo2z)*(I2[44] - (poz)*I3[44])
           + (oo2zn)*I4[44];
H_Cl->Val[115] += U[0][1]*I0[75] + U[4][1]*I1[75]
           + (twooo2z)*(I2[45] - (poz)*I3[45]);
H_Cl->Val[116] += U[0][1]*I0[76] + U[4][1]*I1[76]
           + (twooo2z)*(I2[46] - (poz)*I3[46])
           + (threeoo2zn)*I4[45];
H_Cl->Val[117] += U[0][1]*I0[77] + U[4][1]*I1[77]
           + (twooo2z)*(I2[47] - (poz)*I3[47])
           + (twooo2zn)*I4[46];
H_Cl->Val[118] += U[0][1]*I0[78] + U[4][1]*I1[78]
           + (twooo2z)*(I2[48] - (poz)*I3[48])
           + (oo2zn)*I4[47];
H_Cl->Val[119] += U[0][1]*I0[79] + U[4][1]*I1[79]
           + (twooo2z)*(I2[49] - (poz)*I3[49]);
H_Cl->Val[120] += U[0][1]*I0[80] + U[4][1]*I1[80]
           + (oo2z)*(I2[50] - (poz)*I3[50]);
H_Cl->Val[121] += U[0][1]*I0[81] + U[4][1]*I1[81]
           + (oo2z)*(I2[51] - (poz)*I3[51])
           + (oo2zn)*I4[48];
H_Cl->Val[122] += U[0][1]*I0[82] + U[4][1]*I1[82]
           + (oo2z)*(I2[52] - (poz)*I3[52]);
H_Cl->Val[123] += U[0][1]*I0[83] + U[4][1]*I1[83]
           + (oo2z)*(I2[53] - (poz)*I3[53])
           + (twooo2zn)*I4[49];
H_Cl->Val[124] += U[0][1]*I0[84] + U[4][1]*I1[84]
           + (oo2z)*(I2[54] - (poz)*I3[54])
           + (oo2zn)*I4[50];
H_Cl->Val[125] += U[0][1]*I0[85] + U[4][1]*I1[85]
           + (oo2z)*(I2[55] - (poz)*I3[55]);
H_Cl->Val[126] += U[0][1]*I0[86] + U[4][1]*I1[86]
           + (oo2z)*(I2[56] - (poz)*I3[56])
           + (threeoo2zn)*I4[51];
H_Cl->Val[127] += U[0][1]*I0[87] + U[4][1]*I1[87]
           + (oo2z)*(I2[57] - (poz)*I3[57])
           + (twooo2zn)*I4[52];
H_Cl->Val[128] += U[0][1]*I0[88] + U[4][1]*I1[88]
           + (oo2z)*(I2[58] - (poz)*I3[58])
           + (oo2zn)*I4[53];
H_Cl->Val[129] += U[0][1]*I0[89] + U[4][1]*I1[89]
           + (oo2z)*(I2[59] - (poz)*I3[59]);
H_Cl->Val[130] += U[0][1]*I0[90] + U[4][1]*I1[90];
H_Cl->Val[131] += U[0][1]*I0[91] + U[4][1]*I1[91]
           + (oo2zn)*I4[54];
H_Cl->Val[132] += U[0][1]*I0[92] + U[4][1]*I1[92];
H_Cl->Val[133] += U[0][1]*I0[93] + U[4][1]*I1[93]
           + (twooo2zn)*I4[55];
H_Cl->Val[134] += U[0][1]*I0[94] + U[4][1]*I1[94]
           + (oo2zn)*I4[56];
H_Cl->Val[135] += U[0][1]*I0[95] + U[4][1]*I1[95];
H_Cl->Val[136] += U[0][1]*I0[96] + U[4][1]*I1[96]
           + (threeoo2zn)*I4[57];
H_Cl->Val[137] += U[0][1]*I0[97] + U[4][1]*I1[97]
           + (twooo2zn)*I4[58];
H_Cl->Val[138] += U[0][1]*I0[98] + U[4][1]*I1[98]
           + (oo2zn)*I4[59];
H_Cl->Val[139] += U[0][1]*I0[99] + U[4][1]*I1[99];
H_Cl->Val[140] += U[0][2]*I0[90] + U[4][2]*I1[90]
           + (threeoo2z)*(I2[50] - (poz)*I3[50]);
H_Cl->Val[141] += U[0][2]*I0[91] + U[4][2]*I1[91]
           + (threeoo2z)*(I2[51] - (poz)*I3[51]);
H_Cl->Val[142] += U[0][2]*I0[92] + U[4][2]*I1[92]
           + (threeoo2z)*(I2[52] - (poz)*I3[52])
           + (oo2zn)*I4[54];
H_Cl->Val[143] += U[0][2]*I0[93] + U[4][2]*I1[93]
           + (threeoo2z)*(I2[53] - (poz)*I3[53]);
H_Cl->Val[144] += U[0][2]*I0[94] + U[4][2]*I1[94]
           + (threeoo2z)*(I2[54] - (poz)*I3[54])
           + (oo2zn)*I4[55];
H_Cl->Val[145] += U[0][2]*I0[95] + U[4][2]*I1[95]
           + (threeoo2z)*(I2[55] - (poz)*I3[55])
           + (twooo2zn)*I4[56];
H_Cl->Val[146] += U[0][2]*I0[96] + U[4][2]*I1[96]
           + (threeoo2z)*(I2[56] - (poz)*I3[56]);
H_Cl->Val[147] += U[0][2]*I0[97] + U[4][2]*I1[97]
           + (threeoo2z)*(I2[57] - (poz)*I3[57])
           + (oo2zn)*I4[57];
H_Cl->Val[148] += U[0][2]*I0[98] + U[4][2]*I1[98]
           + (threeoo2z)*(I2[58] - (poz)*I3[58])
           + (twooo2zn)*I4[58];
H_Cl->Val[149] += U[0][2]*I0[99] + U[4][2]*I1[99]
           + (threeoo2z)*(I2[59] - (poz)*I3[59])
           + (threeoo2zn)*I4[59];
}

void
TwoBodyIntJF::top_build_f0g0(iclass *H_Cl, iclass *V_Cl) /* type = 24 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2zn = 2.0*oo2zn;
  double threeoo2zn = 3.0*oo2zn;
  double fouroo2zn = 4.0*oo2zn;
  double twooo2z = 2.0*oo2z;

  I0 = build_d0g0(V_Cl, H_Cl->operands[0]);
  I1 = build_d0g0(V_Cl, H_Cl->operands[1]);
  I2 = build_p0g0(V_Cl, H_Cl->operands[2]);
  I3 = build_p0g0(V_Cl, H_Cl->operands[3]);
  I4 = build_d0f0(V_Cl, H_Cl->operands[4]);

H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0]
           + (twooo2z)*(I2[0] - (poz)*I3[0])
           + (fouroo2zn)*I4[0];
H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1]
           + (twooo2z)*(I2[1] - (poz)*I3[1])
           + (threeoo2zn)*I4[1];
H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2]
           + (twooo2z)*(I2[2] - (poz)*I3[2])
           + (threeoo2zn)*I4[2];
H_Cl->Val[3] += U[0][0]*I0[3] + U[4][0]*I1[3]
           + (twooo2z)*(I2[3] - (poz)*I3[3])
           + (twooo2zn)*I4[3];
H_Cl->Val[4] += U[0][0]*I0[4] + U[4][0]*I1[4]
           + (twooo2z)*(I2[4] - (poz)*I3[4])
           + (twooo2zn)*I4[4];
H_Cl->Val[5] += U[0][0]*I0[5] + U[4][0]*I1[5]
           + (twooo2z)*(I2[5] - (poz)*I3[5])
           + (twooo2zn)*I4[5];
H_Cl->Val[6] += U[0][0]*I0[6] + U[4][0]*I1[6]
           + (twooo2z)*(I2[6] - (poz)*I3[6])
           + (oo2zn)*I4[6];
H_Cl->Val[7] += U[0][0]*I0[7] + U[4][0]*I1[7]
           + (twooo2z)*(I2[7] - (poz)*I3[7])
           + (oo2zn)*I4[7];
H_Cl->Val[8] += U[0][0]*I0[8] + U[4][0]*I1[8]
           + (twooo2z)*(I2[8] - (poz)*I3[8])
           + (oo2zn)*I4[8];
H_Cl->Val[9] += U[0][0]*I0[9] + U[4][0]*I1[9]
           + (twooo2z)*(I2[9] - (poz)*I3[9])
           + (oo2zn)*I4[9];
H_Cl->Val[10] += U[0][0]*I0[10] + U[4][0]*I1[10]
           + (twooo2z)*(I2[10] - (poz)*I3[10]);
H_Cl->Val[11] += U[0][0]*I0[11] + U[4][0]*I1[11]
           + (twooo2z)*(I2[11] - (poz)*I3[11]);
H_Cl->Val[12] += U[0][0]*I0[12] + U[4][0]*I1[12]
           + (twooo2z)*(I2[12] - (poz)*I3[12]);
H_Cl->Val[13] += U[0][0]*I0[13] + U[4][0]*I1[13]
           + (twooo2z)*(I2[13] - (poz)*I3[13]);
H_Cl->Val[14] += U[0][0]*I0[14] + U[4][0]*I1[14]
           + (twooo2z)*(I2[14] - (poz)*I3[14]);
H_Cl->Val[15] += U[0][0]*I0[15] + U[4][0]*I1[15]
           + (oo2z)*(I2[15] - (poz)*I3[15])
           + (fouroo2zn)*I4[10];
H_Cl->Val[16] += U[0][0]*I0[16] + U[4][0]*I1[16]
           + (oo2z)*(I2[16] - (poz)*I3[16])
           + (threeoo2zn)*I4[11];
H_Cl->Val[17] += U[0][0]*I0[17] + U[4][0]*I1[17]
           + (oo2z)*(I2[17] - (poz)*I3[17])
           + (threeoo2zn)*I4[12];
H_Cl->Val[18] += U[0][0]*I0[18] + U[4][0]*I1[18]
           + (oo2z)*(I2[18] - (poz)*I3[18])
           + (twooo2zn)*I4[13];
H_Cl->Val[19] += U[0][0]*I0[19] + U[4][0]*I1[19]
           + (oo2z)*(I2[19] - (poz)*I3[19])
           + (twooo2zn)*I4[14];
H_Cl->Val[20] += U[0][0]*I0[20] + U[4][0]*I1[20]
           + (oo2z)*(I2[20] - (poz)*I3[20])
           + (twooo2zn)*I4[15];
H_Cl->Val[21] += U[0][0]*I0[21] + U[4][0]*I1[21]
           + (oo2z)*(I2[21] - (poz)*I3[21])
           + (oo2zn)*I4[16];
H_Cl->Val[22] += U[0][0]*I0[22] + U[4][0]*I1[22]
           + (oo2z)*(I2[22] - (poz)*I3[22])
           + (oo2zn)*I4[17];
H_Cl->Val[23] += U[0][0]*I0[23] + U[4][0]*I1[23]
           + (oo2z)*(I2[23] - (poz)*I3[23])
           + (oo2zn)*I4[18];
H_Cl->Val[24] += U[0][0]*I0[24] + U[4][0]*I1[24]
           + (oo2z)*(I2[24] - (poz)*I3[24])
           + (oo2zn)*I4[19];
H_Cl->Val[25] += U[0][0]*I0[25] + U[4][0]*I1[25]
           + (oo2z)*(I2[25] - (poz)*I3[25]);
H_Cl->Val[26] += U[0][0]*I0[26] + U[4][0]*I1[26]
           + (oo2z)*(I2[26] - (poz)*I3[26]);
H_Cl->Val[27] += U[0][0]*I0[27] + U[4][0]*I1[27]
           + (oo2z)*(I2[27] - (poz)*I3[27]);
H_Cl->Val[28] += U[0][0]*I0[28] + U[4][0]*I1[28]
           + (oo2z)*(I2[28] - (poz)*I3[28]);
H_Cl->Val[29] += U[0][0]*I0[29] + U[4][0]*I1[29]
           + (oo2z)*(I2[29] - (poz)*I3[29]);
H_Cl->Val[30] += U[0][0]*I0[30] + U[4][0]*I1[30]
           + (oo2z)*(I2[30] - (poz)*I3[30])
           + (fouroo2zn)*I4[20];
H_Cl->Val[31] += U[0][0]*I0[31] + U[4][0]*I1[31]
           + (oo2z)*(I2[31] - (poz)*I3[31])
           + (threeoo2zn)*I4[21];
H_Cl->Val[32] += U[0][0]*I0[32] + U[4][0]*I1[32]
           + (oo2z)*(I2[32] - (poz)*I3[32])
           + (threeoo2zn)*I4[22];
H_Cl->Val[33] += U[0][0]*I0[33] + U[4][0]*I1[33]
           + (oo2z)*(I2[33] - (poz)*I3[33])
           + (twooo2zn)*I4[23];
H_Cl->Val[34] += U[0][0]*I0[34] + U[4][0]*I1[34]
           + (oo2z)*(I2[34] - (poz)*I3[34])
           + (twooo2zn)*I4[24];
H_Cl->Val[35] += U[0][0]*I0[35] + U[4][0]*I1[35]
           + (oo2z)*(I2[35] - (poz)*I3[35])
           + (twooo2zn)*I4[25];
H_Cl->Val[36] += U[0][0]*I0[36] + U[4][0]*I1[36]
           + (oo2z)*(I2[36] - (poz)*I3[36])
           + (oo2zn)*I4[26];
H_Cl->Val[37] += U[0][0]*I0[37] + U[4][0]*I1[37]
           + (oo2z)*(I2[37] - (poz)*I3[37])
           + (oo2zn)*I4[27];
H_Cl->Val[38] += U[0][0]*I0[38] + U[4][0]*I1[38]
           + (oo2z)*(I2[38] - (poz)*I3[38])
           + (oo2zn)*I4[28];
H_Cl->Val[39] += U[0][0]*I0[39] + U[4][0]*I1[39]
           + (oo2z)*(I2[39] - (poz)*I3[39])
           + (oo2zn)*I4[29];
H_Cl->Val[40] += U[0][0]*I0[40] + U[4][0]*I1[40]
           + (oo2z)*(I2[40] - (poz)*I3[40]);
H_Cl->Val[41] += U[0][0]*I0[41] + U[4][0]*I1[41]
           + (oo2z)*(I2[41] - (poz)*I3[41]);
H_Cl->Val[42] += U[0][0]*I0[42] + U[4][0]*I1[42]
           + (oo2z)*(I2[42] - (poz)*I3[42]);
H_Cl->Val[43] += U[0][0]*I0[43] + U[4][0]*I1[43]
           + (oo2z)*(I2[43] - (poz)*I3[43]);
H_Cl->Val[44] += U[0][0]*I0[44] + U[4][0]*I1[44]
           + (oo2z)*(I2[44] - (poz)*I3[44]);
H_Cl->Val[45] += U[0][0]*I0[45] + U[4][0]*I1[45]
           + (fouroo2zn)*I4[30];
H_Cl->Val[46] += U[0][0]*I0[46] + U[4][0]*I1[46]
           + (threeoo2zn)*I4[31];
H_Cl->Val[47] += U[0][0]*I0[47] + U[4][0]*I1[47]
           + (threeoo2zn)*I4[32];
H_Cl->Val[48] += U[0][0]*I0[48] + U[4][0]*I1[48]
           + (twooo2zn)*I4[33];
H_Cl->Val[49] += U[0][0]*I0[49] + U[4][0]*I1[49]
           + (twooo2zn)*I4[34];
H_Cl->Val[50] += U[0][0]*I0[50] + U[4][0]*I1[50]
           + (twooo2zn)*I4[35];
H_Cl->Val[51] += U[0][0]*I0[51] + U[4][0]*I1[51]
           + (oo2zn)*I4[36];
H_Cl->Val[52] += U[0][0]*I0[52] + U[4][0]*I1[52]
           + (oo2zn)*I4[37];
H_Cl->Val[53] += U[0][0]*I0[53] + U[4][0]*I1[53]
           + (oo2zn)*I4[38];
H_Cl->Val[54] += U[0][0]*I0[54] + U[4][0]*I1[54]
           + (oo2zn)*I4[39];
H_Cl->Val[55] += U[0][0]*I0[55] + U[4][0]*I1[55];
H_Cl->Val[56] += U[0][0]*I0[56] + U[4][0]*I1[56];
H_Cl->Val[57] += U[0][0]*I0[57] + U[4][0]*I1[57];
H_Cl->Val[58] += U[0][0]*I0[58] + U[4][0]*I1[58];
H_Cl->Val[59] += U[0][0]*I0[59] + U[4][0]*I1[59];
H_Cl->Val[60] += U[0][0]*I0[60] + U[4][0]*I1[60]
           + (fouroo2zn)*I4[40];
H_Cl->Val[61] += U[0][0]*I0[61] + U[4][0]*I1[61]
           + (threeoo2zn)*I4[41];
H_Cl->Val[62] += U[0][0]*I0[62] + U[4][0]*I1[62]
           + (threeoo2zn)*I4[42];
H_Cl->Val[63] += U[0][0]*I0[63] + U[4][0]*I1[63]
           + (twooo2zn)*I4[43];
H_Cl->Val[64] += U[0][0]*I0[64] + U[4][0]*I1[64]
           + (twooo2zn)*I4[44];
H_Cl->Val[65] += U[0][0]*I0[65] + U[4][0]*I1[65]
           + (twooo2zn)*I4[45];
H_Cl->Val[66] += U[0][0]*I0[66] + U[4][0]*I1[66]
           + (oo2zn)*I4[46];
H_Cl->Val[67] += U[0][0]*I0[67] + U[4][0]*I1[67]
           + (oo2zn)*I4[47];
H_Cl->Val[68] += U[0][0]*I0[68] + U[4][0]*I1[68]
           + (oo2zn)*I4[48];
H_Cl->Val[69] += U[0][0]*I0[69] + U[4][0]*I1[69]
           + (oo2zn)*I4[49];
H_Cl->Val[70] += U[0][0]*I0[70] + U[4][0]*I1[70];
H_Cl->Val[71] += U[0][0]*I0[71] + U[4][0]*I1[71];
H_Cl->Val[72] += U[0][0]*I0[72] + U[4][0]*I1[72];
H_Cl->Val[73] += U[0][0]*I0[73] + U[4][0]*I1[73];
H_Cl->Val[74] += U[0][0]*I0[74] + U[4][0]*I1[74];
H_Cl->Val[75] += U[0][0]*I0[75] + U[4][0]*I1[75]
           + (fouroo2zn)*I4[50];
H_Cl->Val[76] += U[0][0]*I0[76] + U[4][0]*I1[76]
           + (threeoo2zn)*I4[51];
H_Cl->Val[77] += U[0][0]*I0[77] + U[4][0]*I1[77]
           + (threeoo2zn)*I4[52];
H_Cl->Val[78] += U[0][0]*I0[78] + U[4][0]*I1[78]
           + (twooo2zn)*I4[53];
H_Cl->Val[79] += U[0][0]*I0[79] + U[4][0]*I1[79]
           + (twooo2zn)*I4[54];
H_Cl->Val[80] += U[0][0]*I0[80] + U[4][0]*I1[80]
           + (twooo2zn)*I4[55];
H_Cl->Val[81] += U[0][0]*I0[81] + U[4][0]*I1[81]
           + (oo2zn)*I4[56];
H_Cl->Val[82] += U[0][0]*I0[82] + U[4][0]*I1[82]
           + (oo2zn)*I4[57];
H_Cl->Val[83] += U[0][0]*I0[83] + U[4][0]*I1[83]
           + (oo2zn)*I4[58];
H_Cl->Val[84] += U[0][0]*I0[84] + U[4][0]*I1[84]
           + (oo2zn)*I4[59];
H_Cl->Val[85] += U[0][0]*I0[85] + U[4][0]*I1[85];
H_Cl->Val[86] += U[0][0]*I0[86] + U[4][0]*I1[86];
H_Cl->Val[87] += U[0][0]*I0[87] + U[4][0]*I1[87];
H_Cl->Val[88] += U[0][0]*I0[88] + U[4][0]*I1[88];
H_Cl->Val[89] += U[0][0]*I0[89] + U[4][0]*I1[89];
H_Cl->Val[90] += U[0][1]*I0[45] + U[4][1]*I1[45]
           + (twooo2z)*(I2[15] - (poz)*I3[15]);
H_Cl->Val[91] += U[0][1]*I0[46] + U[4][1]*I1[46]
           + (twooo2z)*(I2[16] - (poz)*I3[16])
           + (oo2zn)*I4[30];
H_Cl->Val[92] += U[0][1]*I0[47] + U[4][1]*I1[47]
           + (twooo2z)*(I2[17] - (poz)*I3[17]);
H_Cl->Val[93] += U[0][1]*I0[48] + U[4][1]*I1[48]
           + (twooo2z)*(I2[18] - (poz)*I3[18])
           + (twooo2zn)*I4[31];
H_Cl->Val[94] += U[0][1]*I0[49] + U[4][1]*I1[49]
           + (twooo2z)*(I2[19] - (poz)*I3[19])
           + (oo2zn)*I4[32];
H_Cl->Val[95] += U[0][1]*I0[50] + U[4][1]*I1[50]
           + (twooo2z)*(I2[20] - (poz)*I3[20]);
H_Cl->Val[96] += U[0][1]*I0[51] + U[4][1]*I1[51]
           + (twooo2z)*(I2[21] - (poz)*I3[21])
           + (threeoo2zn)*I4[33];
H_Cl->Val[97] += U[0][1]*I0[52] + U[4][1]*I1[52]
           + (twooo2z)*(I2[22] - (poz)*I3[22])
           + (twooo2zn)*I4[34];
H_Cl->Val[98] += U[0][1]*I0[53] + U[4][1]*I1[53]
           + (twooo2z)*(I2[23] - (poz)*I3[23])
           + (oo2zn)*I4[35];
H_Cl->Val[99] += U[0][1]*I0[54] + U[4][1]*I1[54]
           + (twooo2z)*(I2[24] - (poz)*I3[24]);
H_Cl->Val[100] += U[0][1]*I0[55] + U[4][1]*I1[55]
           + (twooo2z)*(I2[25] - (poz)*I3[25])
           + (fouroo2zn)*I4[36];
H_Cl->Val[101] += U[0][1]*I0[56] + U[4][1]*I1[56]
           + (twooo2z)*(I2[26] - (poz)*I3[26])
           + (threeoo2zn)*I4[37];
H_Cl->Val[102] += U[0][1]*I0[57] + U[4][1]*I1[57]
           + (twooo2z)*(I2[27] - (poz)*I3[27])
           + (twooo2zn)*I4[38];
H_Cl->Val[103] += U[0][1]*I0[58] + U[4][1]*I1[58]
           + (twooo2z)*(I2[28] - (poz)*I3[28])
           + (oo2zn)*I4[39];
H_Cl->Val[104] += U[0][1]*I0[59] + U[4][1]*I1[59]
           + (twooo2z)*(I2[29] - (poz)*I3[29]);
H_Cl->Val[105] += U[0][1]*I0[60] + U[4][1]*I1[60]
           + (oo2z)*(I2[30] - (poz)*I3[30]);
H_Cl->Val[106] += U[0][1]*I0[61] + U[4][1]*I1[61]
           + (oo2z)*(I2[31] - (poz)*I3[31])
           + (oo2zn)*I4[40];
H_Cl->Val[107] += U[0][1]*I0[62] + U[4][1]*I1[62]
           + (oo2z)*(I2[32] - (poz)*I3[32]);
H_Cl->Val[108] += U[0][1]*I0[63] + U[4][1]*I1[63]
           + (oo2z)*(I2[33] - (poz)*I3[33])
           + (twooo2zn)*I4[41];
H_Cl->Val[109] += U[0][1]*I0[64] + U[4][1]*I1[64]
           + (oo2z)*(I2[34] - (poz)*I3[34])
           + (oo2zn)*I4[42];
H_Cl->Val[110] += U[0][1]*I0[65] + U[4][1]*I1[65]
           + (oo2z)*(I2[35] - (poz)*I3[35]);
H_Cl->Val[111] += U[0][1]*I0[66] + U[4][1]*I1[66]
           + (oo2z)*(I2[36] - (poz)*I3[36])
           + (threeoo2zn)*I4[43];
H_Cl->Val[112] += U[0][1]*I0[67] + U[4][1]*I1[67]
           + (oo2z)*(I2[37] - (poz)*I3[37])
           + (twooo2zn)*I4[44];
H_Cl->Val[113] += U[0][1]*I0[68] + U[4][1]*I1[68]
           + (oo2z)*(I2[38] - (poz)*I3[38])
           + (oo2zn)*I4[45];
H_Cl->Val[114] += U[0][1]*I0[69] + U[4][1]*I1[69]
           + (oo2z)*(I2[39] - (poz)*I3[39]);
H_Cl->Val[115] += U[0][1]*I0[70] + U[4][1]*I1[70]
           + (oo2z)*(I2[40] - (poz)*I3[40])
           + (fouroo2zn)*I4[46];
H_Cl->Val[116] += U[0][1]*I0[71] + U[4][1]*I1[71]
           + (oo2z)*(I2[41] - (poz)*I3[41])
           + (threeoo2zn)*I4[47];
H_Cl->Val[117] += U[0][1]*I0[72] + U[4][1]*I1[72]
           + (oo2z)*(I2[42] - (poz)*I3[42])
           + (twooo2zn)*I4[48];
H_Cl->Val[118] += U[0][1]*I0[73] + U[4][1]*I1[73]
           + (oo2z)*(I2[43] - (poz)*I3[43])
           + (oo2zn)*I4[49];
H_Cl->Val[119] += U[0][1]*I0[74] + U[4][1]*I1[74]
           + (oo2z)*(I2[44] - (poz)*I3[44]);
H_Cl->Val[120] += U[0][1]*I0[75] + U[4][1]*I1[75];
H_Cl->Val[121] += U[0][1]*I0[76] + U[4][1]*I1[76]
           + (oo2zn)*I4[50];
H_Cl->Val[122] += U[0][1]*I0[77] + U[4][1]*I1[77];
H_Cl->Val[123] += U[0][1]*I0[78] + U[4][1]*I1[78]
           + (twooo2zn)*I4[51];
H_Cl->Val[124] += U[0][1]*I0[79] + U[4][1]*I1[79]
           + (oo2zn)*I4[52];
H_Cl->Val[125] += U[0][1]*I0[80] + U[4][1]*I1[80];
H_Cl->Val[126] += U[0][1]*I0[81] + U[4][1]*I1[81]
           + (threeoo2zn)*I4[53];
H_Cl->Val[127] += U[0][1]*I0[82] + U[4][1]*I1[82]
           + (twooo2zn)*I4[54];
H_Cl->Val[128] += U[0][1]*I0[83] + U[4][1]*I1[83]
           + (oo2zn)*I4[55];
H_Cl->Val[129] += U[0][1]*I0[84] + U[4][1]*I1[84];
H_Cl->Val[130] += U[0][1]*I0[85] + U[4][1]*I1[85]
           + (fouroo2zn)*I4[56];
H_Cl->Val[131] += U[0][1]*I0[86] + U[4][1]*I1[86]
           + (threeoo2zn)*I4[57];
H_Cl->Val[132] += U[0][1]*I0[87] + U[4][1]*I1[87]
           + (twooo2zn)*I4[58];
H_Cl->Val[133] += U[0][1]*I0[88] + U[4][1]*I1[88]
           + (oo2zn)*I4[59];
H_Cl->Val[134] += U[0][1]*I0[89] + U[4][1]*I1[89];
H_Cl->Val[135] += U[0][2]*I0[75] + U[4][2]*I1[75]
           + (twooo2z)*(I2[30] - (poz)*I3[30]);
H_Cl->Val[136] += U[0][2]*I0[76] + U[4][2]*I1[76]
           + (twooo2z)*(I2[31] - (poz)*I3[31]);
H_Cl->Val[137] += U[0][2]*I0[77] + U[4][2]*I1[77]
           + (twooo2z)*(I2[32] - (poz)*I3[32])
           + (oo2zn)*I4[50];
H_Cl->Val[138] += U[0][2]*I0[78] + U[4][2]*I1[78]
           + (twooo2z)*(I2[33] - (poz)*I3[33]);
H_Cl->Val[139] += U[0][2]*I0[79] + U[4][2]*I1[79]
           + (twooo2z)*(I2[34] - (poz)*I3[34])
           + (oo2zn)*I4[51];
H_Cl->Val[140] += U[0][2]*I0[80] + U[4][2]*I1[80]
           + (twooo2z)*(I2[35] - (poz)*I3[35])
           + (twooo2zn)*I4[52];
H_Cl->Val[141] += U[0][2]*I0[81] + U[4][2]*I1[81]
           + (twooo2z)*(I2[36] - (poz)*I3[36]);
H_Cl->Val[142] += U[0][2]*I0[82] + U[4][2]*I1[82]
           + (twooo2z)*(I2[37] - (poz)*I3[37])
           + (oo2zn)*I4[53];
H_Cl->Val[143] += U[0][2]*I0[83] + U[4][2]*I1[83]
           + (twooo2z)*(I2[38] - (poz)*I3[38])
           + (twooo2zn)*I4[54];
H_Cl->Val[144] += U[0][2]*I0[84] + U[4][2]*I1[84]
           + (twooo2z)*(I2[39] - (poz)*I3[39])
           + (threeoo2zn)*I4[55];
H_Cl->Val[145] += U[0][2]*I0[85] + U[4][2]*I1[85]
           + (twooo2z)*(I2[40] - (poz)*I3[40]);
H_Cl->Val[146] += U[0][2]*I0[86] + U[4][2]*I1[86]
           + (twooo2z)*(I2[41] - (poz)*I3[41])
           + (oo2zn)*I4[56];
H_Cl->Val[147] += U[0][2]*I0[87] + U[4][2]*I1[87]
           + (twooo2z)*(I2[42] - (poz)*I3[42])
           + (twooo2zn)*I4[57];
H_Cl->Val[148] += U[0][2]*I0[88] + U[4][2]*I1[88]
           + (twooo2z)*(I2[43] - (poz)*I3[43])
           + (threeoo2zn)*I4[58];
H_Cl->Val[149] += U[0][2]*I0[89] + U[4][2]*I1[89]
           + (twooo2z)*(I2[44] - (poz)*I3[44])
           + (fouroo2zn)*I4[59];
}

void
TwoBodyIntJF::top_build_g0g0(iclass *H_Cl, iclass *V_Cl) /* type = 25 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2zn = 2.0*oo2zn;
  double threeoo2zn = 3.0*oo2zn;
  double fouroo2zn = 4.0*oo2zn;
  double twooo2z = 2.0*oo2z;
  double threeoo2z = 3.0*oo2z;

  I0 = build_f0g0(V_Cl, H_Cl->operands[0]);
  I1 = build_f0g0(V_Cl, H_Cl->operands[1]);
  I2 = build_d0g0(V_Cl, H_Cl->operands[2]);
  I3 = build_d0g0(V_Cl, H_Cl->operands[3]);
  I4 = build_f0f0(V_Cl, H_Cl->operands[4]);


H_Cl->Val[0] += U[0][0]*I0[0] + U[4][0]*I1[0]
           + (threeoo2z)*(I2[0] - (poz)*I3[0])
           + (fouroo2zn)*I4[0];
H_Cl->Val[1] += U[0][0]*I0[1] + U[4][0]*I1[1]
           + (threeoo2z)*(I2[1] - (poz)*I3[1])
           + (threeoo2zn)*I4[1];
H_Cl->Val[2] += U[0][0]*I0[2] + U[4][0]*I1[2]
           + (threeoo2z)*(I2[2] - (poz)*I3[2])
           + (threeoo2zn)*I4[2];
H_Cl->Val[3] += U[0][0]*I0[3] + U[4][0]*I1[3]
           + (threeoo2z)*(I2[3] - (poz)*I3[3])
           + (twooo2zn)*I4[3];
H_Cl->Val[4] += U[0][0]*I0[4] + U[4][0]*I1[4]
           + (threeoo2z)*(I2[4] - (poz)*I3[4])
           + (twooo2zn)*I4[4];
H_Cl->Val[5] += U[0][0]*I0[5] + U[4][0]*I1[5]
           + (threeoo2z)*(I2[5] - (poz)*I3[5])
           + (twooo2zn)*I4[5];
H_Cl->Val[6] += U[0][0]*I0[6] + U[4][0]*I1[6]
           + (threeoo2z)*(I2[6] - (poz)*I3[6])
           + (oo2zn)*I4[6];
H_Cl->Val[7] += U[0][0]*I0[7] + U[4][0]*I1[7]
           + (threeoo2z)*(I2[7] - (poz)*I3[7])
           + (oo2zn)*I4[7];
H_Cl->Val[8] += U[0][0]*I0[8] + U[4][0]*I1[8]
           + (threeoo2z)*(I2[8] - (poz)*I3[8])
           + (oo2zn)*I4[8];
H_Cl->Val[9] += U[0][0]*I0[9] + U[4][0]*I1[9]
           + (threeoo2z)*(I2[9] - (poz)*I3[9])
           + (oo2zn)*I4[9];
H_Cl->Val[10] += U[0][0]*I0[10] + U[4][0]*I1[10]
           + (threeoo2z)*(I2[10] - (poz)*I3[10]);
H_Cl->Val[11] += U[0][0]*I0[11] + U[4][0]*I1[11]
           + (threeoo2z)*(I2[11] - (poz)*I3[11]);
H_Cl->Val[12] += U[0][0]*I0[12] + U[4][0]*I1[12]
           + (threeoo2z)*(I2[12] - (poz)*I3[12]);
H_Cl->Val[13] += U[0][0]*I0[13] + U[4][0]*I1[13]
           + (threeoo2z)*(I2[13] - (poz)*I3[13]);
H_Cl->Val[14] += U[0][0]*I0[14] + U[4][0]*I1[14]
           + (threeoo2z)*(I2[14] - (poz)*I3[14]);
H_Cl->Val[15] += U[0][0]*I0[15] + U[4][0]*I1[15]
           + (twooo2z)*(I2[15] - (poz)*I3[15])
           + (fouroo2zn)*I4[10];
H_Cl->Val[16] += U[0][0]*I0[16] + U[4][0]*I1[16]
           + (twooo2z)*(I2[16] - (poz)*I3[16])
           + (threeoo2zn)*I4[11];
H_Cl->Val[17] += U[0][0]*I0[17] + U[4][0]*I1[17]
           + (twooo2z)*(I2[17] - (poz)*I3[17])
           + (threeoo2zn)*I4[12];
H_Cl->Val[18] += U[0][0]*I0[18] + U[4][0]*I1[18]
           + (twooo2z)*(I2[18] - (poz)*I3[18])
           + (twooo2zn)*I4[13];
H_Cl->Val[19] += U[0][0]*I0[19] + U[4][0]*I1[19]
           + (twooo2z)*(I2[19] - (poz)*I3[19])
           + (twooo2zn)*I4[14];
H_Cl->Val[20] += U[0][0]*I0[20] + U[4][0]*I1[20]
           + (twooo2z)*(I2[20] - (poz)*I3[20])
           + (twooo2zn)*I4[15];
H_Cl->Val[21] += U[0][0]*I0[21] + U[4][0]*I1[21]
           + (twooo2z)*(I2[21] - (poz)*I3[21])
           + (oo2zn)*I4[16];
H_Cl->Val[22] += U[0][0]*I0[22] + U[4][0]*I1[22]
           + (twooo2z)*(I2[22] - (poz)*I3[22])
           + (oo2zn)*I4[17];
H_Cl->Val[23] += U[0][0]*I0[23] + U[4][0]*I1[23]
           + (twooo2z)*(I2[23] - (poz)*I3[23])
           + (oo2zn)*I4[18];
H_Cl->Val[24] += U[0][0]*I0[24] + U[4][0]*I1[24]
           + (twooo2z)*(I2[24] - (poz)*I3[24])
           + (oo2zn)*I4[19];
H_Cl->Val[25] += U[0][0]*I0[25] + U[4][0]*I1[25]
           + (twooo2z)*(I2[25] - (poz)*I3[25]);
H_Cl->Val[26] += U[0][0]*I0[26] + U[4][0]*I1[26]
           + (twooo2z)*(I2[26] - (poz)*I3[26]);
H_Cl->Val[27] += U[0][0]*I0[27] + U[4][0]*I1[27]
           + (twooo2z)*(I2[27] - (poz)*I3[27]);
H_Cl->Val[28] += U[0][0]*I0[28] + U[4][0]*I1[28]
           + (twooo2z)*(I2[28] - (poz)*I3[28]);
H_Cl->Val[29] += U[0][0]*I0[29] + U[4][0]*I1[29]
           + (twooo2z)*(I2[29] - (poz)*I3[29]);
H_Cl->Val[30] += U[0][0]*I0[30] + U[4][0]*I1[30]
           + (twooo2z)*(I2[30] - (poz)*I3[30])
           + (fouroo2zn)*I4[20];
H_Cl->Val[31] += U[0][0]*I0[31] + U[4][0]*I1[31]
           + (twooo2z)*(I2[31] - (poz)*I3[31])
           + (threeoo2zn)*I4[21];
H_Cl->Val[32] += U[0][0]*I0[32] + U[4][0]*I1[32]
           + (twooo2z)*(I2[32] - (poz)*I3[32])
           + (threeoo2zn)*I4[22];
H_Cl->Val[33] += U[0][0]*I0[33] + U[4][0]*I1[33]
           + (twooo2z)*(I2[33] - (poz)*I3[33])
           + (twooo2zn)*I4[23];
H_Cl->Val[34] += U[0][0]*I0[34] + U[4][0]*I1[34]
           + (twooo2z)*(I2[34] - (poz)*I3[34])
           + (twooo2zn)*I4[24];
H_Cl->Val[35] += U[0][0]*I0[35] + U[4][0]*I1[35]
           + (twooo2z)*(I2[35] - (poz)*I3[35])
           + (twooo2zn)*I4[25];
H_Cl->Val[36] += U[0][0]*I0[36] + U[4][0]*I1[36]
           + (twooo2z)*(I2[36] - (poz)*I3[36])
           + (oo2zn)*I4[26];
H_Cl->Val[37] += U[0][0]*I0[37] + U[4][0]*I1[37]
           + (twooo2z)*(I2[37] - (poz)*I3[37])
           + (oo2zn)*I4[27];
H_Cl->Val[38] += U[0][0]*I0[38] + U[4][0]*I1[38]
           + (twooo2z)*(I2[38] - (poz)*I3[38])
           + (oo2zn)*I4[28];
H_Cl->Val[39] += U[0][0]*I0[39] + U[4][0]*I1[39]
           + (twooo2z)*(I2[39] - (poz)*I3[39])
           + (oo2zn)*I4[29];
H_Cl->Val[40] += U[0][0]*I0[40] + U[4][0]*I1[40]
           + (twooo2z)*(I2[40] - (poz)*I3[40]);
H_Cl->Val[41] += U[0][0]*I0[41] + U[4][0]*I1[41]
           + (twooo2z)*(I2[41] - (poz)*I3[41]);
H_Cl->Val[42] += U[0][0]*I0[42] + U[4][0]*I1[42]
           + (twooo2z)*(I2[42] - (poz)*I3[42]);
H_Cl->Val[43] += U[0][0]*I0[43] + U[4][0]*I1[43]
           + (twooo2z)*(I2[43] - (poz)*I3[43]);
H_Cl->Val[44] += U[0][0]*I0[44] + U[4][0]*I1[44]
           + (twooo2z)*(I2[44] - (poz)*I3[44]);
H_Cl->Val[45] += U[0][0]*I0[45] + U[4][0]*I1[45]
           + (oo2z)*(I2[45] - (poz)*I3[45])
           + (fouroo2zn)*I4[30];
H_Cl->Val[46] += U[0][0]*I0[46] + U[4][0]*I1[46]
           + (oo2z)*(I2[46] - (poz)*I3[46])
           + (threeoo2zn)*I4[31];
H_Cl->Val[47] += U[0][0]*I0[47] + U[4][0]*I1[47]
           + (oo2z)*(I2[47] - (poz)*I3[47])
           + (threeoo2zn)*I4[32];
H_Cl->Val[48] += U[0][0]*I0[48] + U[4][0]*I1[48]
           + (oo2z)*(I2[48] - (poz)*I3[48])
           + (twooo2zn)*I4[33];
H_Cl->Val[49] += U[0][0]*I0[49] + U[4][0]*I1[49]
           + (oo2z)*(I2[49] - (poz)*I3[49])
           + (twooo2zn)*I4[34];
H_Cl->Val[50] += U[0][0]*I0[50] + U[4][0]*I1[50]
           + (oo2z)*(I2[50] - (poz)*I3[50])
           + (twooo2zn)*I4[35];
H_Cl->Val[51] += U[0][0]*I0[51] + U[4][0]*I1[51]
           + (oo2z)*(I2[51] - (poz)*I3[51])
           + (oo2zn)*I4[36];
H_Cl->Val[52] += U[0][0]*I0[52] + U[4][0]*I1[52]
           + (oo2z)*(I2[52] - (poz)*I3[52])
           + (oo2zn)*I4[37];
H_Cl->Val[53] += U[0][0]*I0[53] + U[4][0]*I1[53]
           + (oo2z)*(I2[53] - (poz)*I3[53])
           + (oo2zn)*I4[38];
H_Cl->Val[54] += U[0][0]*I0[54] + U[4][0]*I1[54]
           + (oo2z)*(I2[54] - (poz)*I3[54])
           + (oo2zn)*I4[39];
H_Cl->Val[55] += U[0][0]*I0[55] + U[4][0]*I1[55]
           + (oo2z)*(I2[55] - (poz)*I3[55]);
H_Cl->Val[56] += U[0][0]*I0[56] + U[4][0]*I1[56]
           + (oo2z)*(I2[56] - (poz)*I3[56]);
H_Cl->Val[57] += U[0][0]*I0[57] + U[4][0]*I1[57]
           + (oo2z)*(I2[57] - (poz)*I3[57]);
H_Cl->Val[58] += U[0][0]*I0[58] + U[4][0]*I1[58]
           + (oo2z)*(I2[58] - (poz)*I3[58]);
H_Cl->Val[59] += U[0][0]*I0[59] + U[4][0]*I1[59]
           + (oo2z)*(I2[59] - (poz)*I3[59]);
H_Cl->Val[60] += U[0][0]*I0[60] + U[4][0]*I1[60]
           + (oo2z)*(I2[60] - (poz)*I3[60])
           + (fouroo2zn)*I4[40];
H_Cl->Val[61] += U[0][0]*I0[61] + U[4][0]*I1[61]
           + (oo2z)*(I2[61] - (poz)*I3[61])
           + (threeoo2zn)*I4[41];
H_Cl->Val[62] += U[0][0]*I0[62] + U[4][0]*I1[62]
           + (oo2z)*(I2[62] - (poz)*I3[62])
           + (threeoo2zn)*I4[42];
H_Cl->Val[63] += U[0][0]*I0[63] + U[4][0]*I1[63]
           + (oo2z)*(I2[63] - (poz)*I3[63])
           + (twooo2zn)*I4[43];
H_Cl->Val[64] += U[0][0]*I0[64] + U[4][0]*I1[64]
           + (oo2z)*(I2[64] - (poz)*I3[64])
           + (twooo2zn)*I4[44];
H_Cl->Val[65] += U[0][0]*I0[65] + U[4][0]*I1[65]
           + (oo2z)*(I2[65] - (poz)*I3[65])
           + (twooo2zn)*I4[45];
H_Cl->Val[66] += U[0][0]*I0[66] + U[4][0]*I1[66]
           + (oo2z)*(I2[66] - (poz)*I3[66])
           + (oo2zn)*I4[46];
H_Cl->Val[67] += U[0][0]*I0[67] + U[4][0]*I1[67]
           + (oo2z)*(I2[67] - (poz)*I3[67])
           + (oo2zn)*I4[47];
H_Cl->Val[68] += U[0][0]*I0[68] + U[4][0]*I1[68]
           + (oo2z)*(I2[68] - (poz)*I3[68])
           + (oo2zn)*I4[48];
H_Cl->Val[69] += U[0][0]*I0[69] + U[4][0]*I1[69]
           + (oo2z)*(I2[69] - (poz)*I3[69])
           + (oo2zn)*I4[49];
H_Cl->Val[70] += U[0][0]*I0[70] + U[4][0]*I1[70]
           + (oo2z)*(I2[70] - (poz)*I3[70]);
H_Cl->Val[71] += U[0][0]*I0[71] + U[4][0]*I1[71]
           + (oo2z)*(I2[71] - (poz)*I3[71]);
H_Cl->Val[72] += U[0][0]*I0[72] + U[4][0]*I1[72]
           + (oo2z)*(I2[72] - (poz)*I3[72]);
H_Cl->Val[73] += U[0][0]*I0[73] + U[4][0]*I1[73]
           + (oo2z)*(I2[73] - (poz)*I3[73]);
H_Cl->Val[74] += U[0][0]*I0[74] + U[4][0]*I1[74]
           + (oo2z)*(I2[74] - (poz)*I3[74]);
H_Cl->Val[75] += U[0][0]*I0[75] + U[4][0]*I1[75]
           + (oo2z)*(I2[75] - (poz)*I3[75])
           + (fouroo2zn)*I4[50];
H_Cl->Val[76] += U[0][0]*I0[76] + U[4][0]*I1[76]
           + (oo2z)*(I2[76] - (poz)*I3[76])
           + (threeoo2zn)*I4[51];
H_Cl->Val[77] += U[0][0]*I0[77] + U[4][0]*I1[77]
           + (oo2z)*(I2[77] - (poz)*I3[77])
           + (threeoo2zn)*I4[52];
H_Cl->Val[78] += U[0][0]*I0[78] + U[4][0]*I1[78]
           + (oo2z)*(I2[78] - (poz)*I3[78])
           + (twooo2zn)*I4[53];
H_Cl->Val[79] += U[0][0]*I0[79] + U[4][0]*I1[79]
           + (oo2z)*(I2[79] - (poz)*I3[79])
           + (twooo2zn)*I4[54];
H_Cl->Val[80] += U[0][0]*I0[80] + U[4][0]*I1[80]
           + (oo2z)*(I2[80] - (poz)*I3[80])
           + (twooo2zn)*I4[55];
H_Cl->Val[81] += U[0][0]*I0[81] + U[4][0]*I1[81]
           + (oo2z)*(I2[81] - (poz)*I3[81])
           + (oo2zn)*I4[56];
H_Cl->Val[82] += U[0][0]*I0[82] + U[4][0]*I1[82]
           + (oo2z)*(I2[82] - (poz)*I3[82])
           + (oo2zn)*I4[57];
H_Cl->Val[83] += U[0][0]*I0[83] + U[4][0]*I1[83]
           + (oo2z)*(I2[83] - (poz)*I3[83])
           + (oo2zn)*I4[58];
H_Cl->Val[84] += U[0][0]*I0[84] + U[4][0]*I1[84]
           + (oo2z)*(I2[84] - (poz)*I3[84])
           + (oo2zn)*I4[59];
H_Cl->Val[85] += U[0][0]*I0[85] + U[4][0]*I1[85]
           + (oo2z)*(I2[85] - (poz)*I3[85]);
H_Cl->Val[86] += U[0][0]*I0[86] + U[4][0]*I1[86]
           + (oo2z)*(I2[86] - (poz)*I3[86]);
H_Cl->Val[87] += U[0][0]*I0[87] + U[4][0]*I1[87]
           + (oo2z)*(I2[87] - (poz)*I3[87]);
H_Cl->Val[88] += U[0][0]*I0[88] + U[4][0]*I1[88]
           + (oo2z)*(I2[88] - (poz)*I3[88]);
H_Cl->Val[89] += U[0][0]*I0[89] + U[4][0]*I1[89]
           + (oo2z)*(I2[89] - (poz)*I3[89]);
H_Cl->Val[90] += U[0][0]*I0[90] + U[4][0]*I1[90]
           + (fouroo2zn)*I4[60];
H_Cl->Val[91] += U[0][0]*I0[91] + U[4][0]*I1[91]
           + (threeoo2zn)*I4[61];
H_Cl->Val[92] += U[0][0]*I0[92] + U[4][0]*I1[92]
           + (threeoo2zn)*I4[62];
H_Cl->Val[93] += U[0][0]*I0[93] + U[4][0]*I1[93]
           + (twooo2zn)*I4[63];
H_Cl->Val[94] += U[0][0]*I0[94] + U[4][0]*I1[94]
           + (twooo2zn)*I4[64];
H_Cl->Val[95] += U[0][0]*I0[95] + U[4][0]*I1[95]
           + (twooo2zn)*I4[65];
H_Cl->Val[96] += U[0][0]*I0[96] + U[4][0]*I1[96]
           + (oo2zn)*I4[66];
H_Cl->Val[97] += U[0][0]*I0[97] + U[4][0]*I1[97]
           + (oo2zn)*I4[67];
H_Cl->Val[98] += U[0][0]*I0[98] + U[4][0]*I1[98]
           + (oo2zn)*I4[68];
H_Cl->Val[99] += U[0][0]*I0[99] + U[4][0]*I1[99]
           + (oo2zn)*I4[69];
H_Cl->Val[100] += U[0][0]*I0[100] + U[4][0]*I1[100];
H_Cl->Val[101] += U[0][0]*I0[101] + U[4][0]*I1[101];
H_Cl->Val[102] += U[0][0]*I0[102] + U[4][0]*I1[102];
H_Cl->Val[103] += U[0][0]*I0[103] + U[4][0]*I1[103];
H_Cl->Val[104] += U[0][0]*I0[104] + U[4][0]*I1[104];
H_Cl->Val[105] += U[0][0]*I0[105] + U[4][0]*I1[105]
           + (fouroo2zn)*I4[70];
H_Cl->Val[106] += U[0][0]*I0[106] + U[4][0]*I1[106]
           + (threeoo2zn)*I4[71];
H_Cl->Val[107] += U[0][0]*I0[107] + U[4][0]*I1[107]
           + (threeoo2zn)*I4[72];
H_Cl->Val[108] += U[0][0]*I0[108] + U[4][0]*I1[108]
           + (twooo2zn)*I4[73];
H_Cl->Val[109] += U[0][0]*I0[109] + U[4][0]*I1[109]
           + (twooo2zn)*I4[74];
H_Cl->Val[110] += U[0][0]*I0[110] + U[4][0]*I1[110]
           + (twooo2zn)*I4[75];
H_Cl->Val[111] += U[0][0]*I0[111] + U[4][0]*I1[111]
           + (oo2zn)*I4[76];
H_Cl->Val[112] += U[0][0]*I0[112] + U[4][0]*I1[112]
           + (oo2zn)*I4[77];
H_Cl->Val[113] += U[0][0]*I0[113] + U[4][0]*I1[113]
           + (oo2zn)*I4[78];
H_Cl->Val[114] += U[0][0]*I0[114] + U[4][0]*I1[114]
           + (oo2zn)*I4[79];
H_Cl->Val[115] += U[0][0]*I0[115] + U[4][0]*I1[115];
H_Cl->Val[116] += U[0][0]*I0[116] + U[4][0]*I1[116];
H_Cl->Val[117] += U[0][0]*I0[117] + U[4][0]*I1[117];
H_Cl->Val[118] += U[0][0]*I0[118] + U[4][0]*I1[118];
H_Cl->Val[119] += U[0][0]*I0[119] + U[4][0]*I1[119];
H_Cl->Val[120] += U[0][0]*I0[120] + U[4][0]*I1[120]
           + (fouroo2zn)*I4[80];
H_Cl->Val[121] += U[0][0]*I0[121] + U[4][0]*I1[121]
           + (threeoo2zn)*I4[81];
H_Cl->Val[122] += U[0][0]*I0[122] + U[4][0]*I1[122]
           + (threeoo2zn)*I4[82];
H_Cl->Val[123] += U[0][0]*I0[123] + U[4][0]*I1[123]
           + (twooo2zn)*I4[83];
H_Cl->Val[124] += U[0][0]*I0[124] + U[4][0]*I1[124]
           + (twooo2zn)*I4[84];
H_Cl->Val[125] += U[0][0]*I0[125] + U[4][0]*I1[125]
           + (twooo2zn)*I4[85];
H_Cl->Val[126] += U[0][0]*I0[126] + U[4][0]*I1[126]
           + (oo2zn)*I4[86];
H_Cl->Val[127] += U[0][0]*I0[127] + U[4][0]*I1[127]
           + (oo2zn)*I4[87];
H_Cl->Val[128] += U[0][0]*I0[128] + U[4][0]*I1[128]
           + (oo2zn)*I4[88];
H_Cl->Val[129] += U[0][0]*I0[129] + U[4][0]*I1[129]
           + (oo2zn)*I4[89];
H_Cl->Val[130] += U[0][0]*I0[130] + U[4][0]*I1[130];
H_Cl->Val[131] += U[0][0]*I0[131] + U[4][0]*I1[131];
H_Cl->Val[132] += U[0][0]*I0[132] + U[4][0]*I1[132];
H_Cl->Val[133] += U[0][0]*I0[133] + U[4][0]*I1[133];
H_Cl->Val[134] += U[0][0]*I0[134] + U[4][0]*I1[134];
H_Cl->Val[135] += U[0][0]*I0[135] + U[4][0]*I1[135]
           + (fouroo2zn)*I4[90];
H_Cl->Val[136] += U[0][0]*I0[136] + U[4][0]*I1[136]
           + (threeoo2zn)*I4[91];
H_Cl->Val[137] += U[0][0]*I0[137] + U[4][0]*I1[137]
           + (threeoo2zn)*I4[92];
H_Cl->Val[138] += U[0][0]*I0[138] + U[4][0]*I1[138]
           + (twooo2zn)*I4[93];
H_Cl->Val[139] += U[0][0]*I0[139] + U[4][0]*I1[139]
           + (twooo2zn)*I4[94];
H_Cl->Val[140] += U[0][0]*I0[140] + U[4][0]*I1[140]
           + (twooo2zn)*I4[95];
H_Cl->Val[141] += U[0][0]*I0[141] + U[4][0]*I1[141]
           + (oo2zn)*I4[96];
H_Cl->Val[142] += U[0][0]*I0[142] + U[4][0]*I1[142]
           + (oo2zn)*I4[97];
H_Cl->Val[143] += U[0][0]*I0[143] + U[4][0]*I1[143]
           + (oo2zn)*I4[98];
H_Cl->Val[144] += U[0][0]*I0[144] + U[4][0]*I1[144]
           + (oo2zn)*I4[99];
H_Cl->Val[145] += U[0][0]*I0[145] + U[4][0]*I1[145];
H_Cl->Val[146] += U[0][0]*I0[146] + U[4][0]*I1[146];
H_Cl->Val[147] += U[0][0]*I0[147] + U[4][0]*I1[147];
H_Cl->Val[148] += U[0][0]*I0[148] + U[4][0]*I1[148];
H_Cl->Val[149] += U[0][0]*I0[149] + U[4][0]*I1[149];
H_Cl->Val[150] += U[0][1]*I0[90] + U[4][1]*I1[90]
           + (threeoo2z)*(I2[45] - (poz)*I3[45]);
H_Cl->Val[151] += U[0][1]*I0[91] + U[4][1]*I1[91]
           + (threeoo2z)*(I2[46] - (poz)*I3[46])
           + (oo2zn)*I4[60];
H_Cl->Val[152] += U[0][1]*I0[92] + U[4][1]*I1[92]
           + (threeoo2z)*(I2[47] - (poz)*I3[47]);
H_Cl->Val[153] += U[0][1]*I0[93] + U[4][1]*I1[93]
           + (threeoo2z)*(I2[48] - (poz)*I3[48])
           + (twooo2zn)*I4[61];
H_Cl->Val[154] += U[0][1]*I0[94] + U[4][1]*I1[94]
           + (threeoo2z)*(I2[49] - (poz)*I3[49])
           + (oo2zn)*I4[62];
H_Cl->Val[155] += U[0][1]*I0[95] + U[4][1]*I1[95]
           + (threeoo2z)*(I2[50] - (poz)*I3[50]);
H_Cl->Val[156] += U[0][1]*I0[96] + U[4][1]*I1[96]
           + (threeoo2z)*(I2[51] - (poz)*I3[51])
           + (threeoo2zn)*I4[63];
H_Cl->Val[157] += U[0][1]*I0[97] + U[4][1]*I1[97]
           + (threeoo2z)*(I2[52] - (poz)*I3[52])
           + (twooo2zn)*I4[64];
H_Cl->Val[158] += U[0][1]*I0[98] + U[4][1]*I1[98]
           + (threeoo2z)*(I2[53] - (poz)*I3[53])
           + (oo2zn)*I4[65];
H_Cl->Val[159] += U[0][1]*I0[99] + U[4][1]*I1[99]
           + (threeoo2z)*(I2[54] - (poz)*I3[54]);
H_Cl->Val[160] += U[0][1]*I0[100] + U[4][1]*I1[100]
           + (threeoo2z)*(I2[55] - (poz)*I3[55])
           + (fouroo2zn)*I4[66];
H_Cl->Val[161] += U[0][1]*I0[101] + U[4][1]*I1[101]
           + (threeoo2z)*(I2[56] - (poz)*I3[56])
           + (threeoo2zn)*I4[67];
H_Cl->Val[162] += U[0][1]*I0[102] + U[4][1]*I1[102]
           + (threeoo2z)*(I2[57] - (poz)*I3[57])
           + (twooo2zn)*I4[68];
H_Cl->Val[163] += U[0][1]*I0[103] + U[4][1]*I1[103]
           + (threeoo2z)*(I2[58] - (poz)*I3[58])
           + (oo2zn)*I4[69];
H_Cl->Val[164] += U[0][1]*I0[104] + U[4][1]*I1[104]
           + (threeoo2z)*(I2[59] - (poz)*I3[59]);
H_Cl->Val[165] += U[0][1]*I0[105] + U[4][1]*I1[105]
           + (twooo2z)*(I2[60] - (poz)*I3[60]);
H_Cl->Val[166] += U[0][1]*I0[106] + U[4][1]*I1[106]
           + (twooo2z)*(I2[61] - (poz)*I3[61])
           + (oo2zn)*I4[70];
H_Cl->Val[167] += U[0][1]*I0[107] + U[4][1]*I1[107]
           + (twooo2z)*(I2[62] - (poz)*I3[62]);
H_Cl->Val[168] += U[0][1]*I0[108] + U[4][1]*I1[108]
           + (twooo2z)*(I2[63] - (poz)*I3[63])
           + (twooo2zn)*I4[71];
H_Cl->Val[169] += U[0][1]*I0[109] + U[4][1]*I1[109]
           + (twooo2z)*(I2[64] - (poz)*I3[64])
           + (oo2zn)*I4[72];
H_Cl->Val[170] += U[0][1]*I0[110] + U[4][1]*I1[110]
           + (twooo2z)*(I2[65] - (poz)*I3[65]);
H_Cl->Val[171] += U[0][1]*I0[111] + U[4][1]*I1[111]
           + (twooo2z)*(I2[66] - (poz)*I3[66])
           + (threeoo2zn)*I4[73];
H_Cl->Val[172] += U[0][1]*I0[112] + U[4][1]*I1[112]
           + (twooo2z)*(I2[67] - (poz)*I3[67])
           + (twooo2zn)*I4[74];
H_Cl->Val[173] += U[0][1]*I0[113] + U[4][1]*I1[113]
           + (twooo2z)*(I2[68] - (poz)*I3[68])
           + (oo2zn)*I4[75];
H_Cl->Val[174] += U[0][1]*I0[114] + U[4][1]*I1[114]
           + (twooo2z)*(I2[69] - (poz)*I3[69]);
H_Cl->Val[175] += U[0][1]*I0[115] + U[4][1]*I1[115]
           + (twooo2z)*(I2[70] - (poz)*I3[70])
           + (fouroo2zn)*I4[76];
H_Cl->Val[176] += U[0][1]*I0[116] + U[4][1]*I1[116]
           + (twooo2z)*(I2[71] - (poz)*I3[71])
           + (threeoo2zn)*I4[77];
H_Cl->Val[177] += U[0][1]*I0[117] + U[4][1]*I1[117]
           + (twooo2z)*(I2[72] - (poz)*I3[72])
           + (twooo2zn)*I4[78];
H_Cl->Val[178] += U[0][1]*I0[118] + U[4][1]*I1[118]
           + (twooo2z)*(I2[73] - (poz)*I3[73])
           + (oo2zn)*I4[79];
H_Cl->Val[179] += U[0][1]*I0[119] + U[4][1]*I1[119]
           + (twooo2z)*(I2[74] - (poz)*I3[74]);
H_Cl->Val[180] += U[0][1]*I0[120] + U[4][1]*I1[120]
           + (oo2z)*(I2[75] - (poz)*I3[75]);
H_Cl->Val[181] += U[0][1]*I0[121] + U[4][1]*I1[121]
           + (oo2z)*(I2[76] - (poz)*I3[76])
           + (oo2zn)*I4[80];
H_Cl->Val[182] += U[0][1]*I0[122] + U[4][1]*I1[122]
           + (oo2z)*(I2[77] - (poz)*I3[77]);
H_Cl->Val[183] += U[0][1]*I0[123] + U[4][1]*I1[123]
           + (oo2z)*(I2[78] - (poz)*I3[78])
           + (twooo2zn)*I4[81];
H_Cl->Val[184] += U[0][1]*I0[124] + U[4][1]*I1[124]
           + (oo2z)*(I2[79] - (poz)*I3[79])
           + (oo2zn)*I4[82];
H_Cl->Val[185] += U[0][1]*I0[125] + U[4][1]*I1[125]
           + (oo2z)*(I2[80] - (poz)*I3[80]);
H_Cl->Val[186] += U[0][1]*I0[126] + U[4][1]*I1[126]
           + (oo2z)*(I2[81] - (poz)*I3[81])
           + (threeoo2zn)*I4[83];
H_Cl->Val[187] += U[0][1]*I0[127] + U[4][1]*I1[127]
           + (oo2z)*(I2[82] - (poz)*I3[82])
           + (twooo2zn)*I4[84];
H_Cl->Val[188] += U[0][1]*I0[128] + U[4][1]*I1[128]
           + (oo2z)*(I2[83] - (poz)*I3[83])
           + (oo2zn)*I4[85];
H_Cl->Val[189] += U[0][1]*I0[129] + U[4][1]*I1[129]
           + (oo2z)*(I2[84] - (poz)*I3[84]);
H_Cl->Val[190] += U[0][1]*I0[130] + U[4][1]*I1[130]
           + (oo2z)*(I2[85] - (poz)*I3[85])
           + (fouroo2zn)*I4[86];
H_Cl->Val[191] += U[0][1]*I0[131] + U[4][1]*I1[131]
           + (oo2z)*(I2[86] - (poz)*I3[86])
           + (threeoo2zn)*I4[87];
H_Cl->Val[192] += U[0][1]*I0[132] + U[4][1]*I1[132]
           + (oo2z)*(I2[87] - (poz)*I3[87])
           + (twooo2zn)*I4[88];
H_Cl->Val[193] += U[0][1]*I0[133] + U[4][1]*I1[133]
           + (oo2z)*(I2[88] - (poz)*I3[88])
           + (oo2zn)*I4[89];
H_Cl->Val[194] += U[0][1]*I0[134] + U[4][1]*I1[134]
           + (oo2z)*(I2[89] - (poz)*I3[89]);
H_Cl->Val[195] += U[0][1]*I0[135] + U[4][1]*I1[135];
H_Cl->Val[196] += U[0][1]*I0[136] + U[4][1]*I1[136]
           + (oo2zn)*I4[90];
H_Cl->Val[197] += U[0][1]*I0[137] + U[4][1]*I1[137];
H_Cl->Val[198] += U[0][1]*I0[138] + U[4][1]*I1[138]
           + (twooo2zn)*I4[91];
H_Cl->Val[199] += U[0][1]*I0[139] + U[4][1]*I1[139]
           + (oo2zn)*I4[92];
H_Cl->Val[200] += U[0][1]*I0[140] + U[4][1]*I1[140];
H_Cl->Val[201] += U[0][1]*I0[141] + U[4][1]*I1[141]
           + (threeoo2zn)*I4[93];
H_Cl->Val[202] += U[0][1]*I0[142] + U[4][1]*I1[142]
           + (twooo2zn)*I4[94];
H_Cl->Val[203] += U[0][1]*I0[143] + U[4][1]*I1[143]
           + (oo2zn)*I4[95];
H_Cl->Val[204] += U[0][1]*I0[144] + U[4][1]*I1[144];
H_Cl->Val[205] += U[0][1]*I0[145] + U[4][1]*I1[145]
           + (fouroo2zn)*I4[96];
H_Cl->Val[206] += U[0][1]*I0[146] + U[4][1]*I1[146]
           + (threeoo2zn)*I4[97];
H_Cl->Val[207] += U[0][1]*I0[147] + U[4][1]*I1[147]
           + (twooo2zn)*I4[98];
H_Cl->Val[208] += U[0][1]*I0[148] + U[4][1]*I1[148]
           + (oo2zn)*I4[99];
H_Cl->Val[209] += U[0][1]*I0[149] + U[4][1]*I1[149];
H_Cl->Val[210] += U[0][2]*I0[135] + U[4][2]*I1[135]
           + (threeoo2z)*(I2[75] - (poz)*I3[75]);
H_Cl->Val[211] += U[0][2]*I0[136] + U[4][2]*I1[136]
           + (threeoo2z)*(I2[76] - (poz)*I3[76]);
H_Cl->Val[212] += U[0][2]*I0[137] + U[4][2]*I1[137]
           + (threeoo2z)*(I2[77] - (poz)*I3[77])
           + (oo2zn)*I4[90];
H_Cl->Val[213] += U[0][2]*I0[138] + U[4][2]*I1[138]
           + (threeoo2z)*(I2[78] - (poz)*I3[78]);
H_Cl->Val[214] += U[0][2]*I0[139] + U[4][2]*I1[139]
           + (threeoo2z)*(I2[79] - (poz)*I3[79])
           + (oo2zn)*I4[91];
H_Cl->Val[215] += U[0][2]*I0[140] + U[4][2]*I1[140]
           + (threeoo2z)*(I2[80] - (poz)*I3[80])
           + (twooo2zn)*I4[92];
H_Cl->Val[216] += U[0][2]*I0[141] + U[4][2]*I1[141]
           + (threeoo2z)*(I2[81] - (poz)*I3[81]);
H_Cl->Val[217] += U[0][2]*I0[142] + U[4][2]*I1[142]
           + (threeoo2z)*(I2[82] - (poz)*I3[82])
           + (oo2zn)*I4[93];
H_Cl->Val[218] += U[0][2]*I0[143] + U[4][2]*I1[143]
           + (threeoo2z)*(I2[83] - (poz)*I3[83])
           + (twooo2zn)*I4[94];
H_Cl->Val[219] += U[0][2]*I0[144] + U[4][2]*I1[144]
           + (threeoo2z)*(I2[84] - (poz)*I3[84])
           + (threeoo2zn)*I4[95];
H_Cl->Val[220] += U[0][2]*I0[145] + U[4][2]*I1[145]
           + (threeoo2z)*(I2[85] - (poz)*I3[85]);
H_Cl->Val[221] += U[0][2]*I0[146] + U[4][2]*I1[146]
           + (threeoo2z)*(I2[86] - (poz)*I3[86])
           + (oo2zn)*I4[96];
H_Cl->Val[222] += U[0][2]*I0[147] + U[4][2]*I1[147]
           + (threeoo2z)*(I2[87] - (poz)*I3[87])
           + (twooo2zn)*I4[97];
H_Cl->Val[223] += U[0][2]*I0[148] + U[4][2]*I1[148]
           + (threeoo2z)*(I2[88] - (poz)*I3[88])
           + (threeoo2zn)*I4[98];
H_Cl->Val[224] += U[0][2]*I0[149] + U[4][2]*I1[149]
           + (threeoo2z)*(I2[89] - (poz)*I3[89])
           + (fouroo2zn)*I4[99];
}
