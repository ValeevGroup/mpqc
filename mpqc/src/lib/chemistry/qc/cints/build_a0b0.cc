
#include <chemistry/qc/cints/int2jf.h>

double *
TwoBodyIntJF::build_p000(iclass *V_Cl, int in) /* type = 2 */
{
  int m, mp1;
  if(V_done[in]) return V_Cl[in].Val;

  m = V_Cl[in].m;
  mp1 = m+1;
  V_Cl[in].Val[0] = U[0][0]*F[m] + U[4][0]*F[mp1];
  V_Cl[in].Val[1] = U[0][1]*F[m] + U[4][1]*F[mp1];
  V_Cl[in].Val[2] = U[0][2]*F[m] + U[4][2]*F[mp1];
  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_00p0(iclass *V_Cl, int in) /* type = 3 */
{
  int m, mp1;
  if(V_done[in]) return V_Cl[in].Val;

  m = V_Cl[in].m;
  mp1 = m+1;
  V_Cl[in].Val[0] = U[2][0]*F[m] + U[5][0]*F[mp1];
  V_Cl[in].Val[1] = U[2][1]*F[m] + U[5][1]*F[mp1];
  V_Cl[in].Val[2] = U[2][2]*F[m] + U[5][2]*F[mp1];
  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_d000(iclass *V_Cl, int in) /* type = 4 */
{
  double *I0, *I1;
  double K1;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_p000(V_Cl, V_Cl[in].operands[0]);
  I1 = build_p000(V_Cl, V_Cl[in].operands[1]);

  K1 = oo2z*(F[V_Cl[in].m]-poz*F[V_Cl[in].m+1]);

  V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0] + K1;
  V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1];
  V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2];
  V_Cl[in].Val[3] = U[0][1]*I0[1] + U[4][1]*I1[1] + K1;
  V_Cl[in].Val[4] = U[0][1]*I0[2] + U[4][1]*I1[2];
  V_Cl[in].Val[5] = U[0][2]*I0[2] + U[4][2]*I1[2] + K1;

  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_00d0(iclass *V_Cl, int in) /* type = 5 */
{
  double *I0, *I1;
  double K1;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_00p0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_00p0(V_Cl, V_Cl[in].operands[1]);

  K1 = oo2n*(F[V_Cl[in].m]-pon*F[V_Cl[in].m+1]);

  V_Cl[in].Val[0] = U[2][0]*I0[0] + U[5][0]*I1[0] + K1;
  V_Cl[in].Val[1] = U[2][0]*I0[1] + U[5][0]*I1[1];
  V_Cl[in].Val[2] = U[2][0]*I0[2] + U[5][0]*I1[2];
  V_Cl[in].Val[3] = U[2][1]*I0[1] + U[5][1]*I1[1] + K1;
  V_Cl[in].Val[4] = U[2][1]*I0[2] + U[5][1]*I1[2];
  V_Cl[in].Val[5] = U[2][2]*I0[2] + U[5][2]*I1[2] + K1;

  V_done[in] = 1;
  return V_Cl[in].Val;
}


double *
TwoBodyIntJF::build_p0p0(iclass *V_Cl, int in) /* type = 6 */
{
  double *I0, *I1;
  double K1;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_00p0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_00p0(V_Cl, V_Cl[in].operands[1]);
  K1 = oo2zn*F[V_Cl[in].m+1];

  V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0] + K1;
  V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1];
  V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2];
  V_Cl[in].Val[3] = U[0][1]*I0[0] + U[4][1]*I1[0];
  V_Cl[in].Val[4] = U[0][1]*I0[1] + U[4][1]*I1[1] + K1;
  V_Cl[in].Val[5] = U[0][1]*I0[2] + U[4][1]*I1[2];
  V_Cl[in].Val[6] = U[0][2]*I0[0] + U[4][2]*I1[0];
  V_Cl[in].Val[7] = U[0][2]*I0[1] + U[4][2]*I1[1];
  V_Cl[in].Val[8] = U[0][2]*I0[2] + U[4][2]*I1[2] + K1;

  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_f000(iclass *V_Cl, int in) /* type = 7 */
{
  double *I0, *I1, *I2, *I3;
  double K1 = 2*oo2z;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_d000(V_Cl, V_Cl[in].operands[0]);
  I1 = build_d000(V_Cl, V_Cl[in].operands[1]);
  I2 = build_p000(V_Cl, V_Cl[in].operands[2]);
  I3 = build_p000(V_Cl, V_Cl[in].operands[3]);

  V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0] + K1*(I2[0]-poz*I3[0]);
  V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1] + oo2z*(I2[1]-poz*I3[1]);
  V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2] + oo2z*(I2[2]-poz*I3[2]);
  V_Cl[in].Val[3] = U[0][0]*I0[3] + U[4][0]*I1[3];
  V_Cl[in].Val[4] = U[0][0]*I0[4] + U[4][0]*I1[4];
  V_Cl[in].Val[5] = U[0][0]*I0[5] + U[4][0]*I1[5];
  V_Cl[in].Val[6] = U[0][1]*I0[3] + U[4][1]*I1[3] + K1*(I2[1]-poz*I3[1]);
  V_Cl[in].Val[7] = U[0][1]*I0[4] + U[4][1]*I1[4] + oo2z*(I2[2]-poz*I3[2]);
  V_Cl[in].Val[8] = U[0][1]*I0[5] + U[4][1]*I1[5];
  V_Cl[in].Val[9] = U[0][2]*I0[5] + U[4][2]*I1[5] + K1*(I2[2]-poz*I3[2]);

  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_00f0(iclass *V_Cl, int in) /* type = 8 */
{
  double *I0, *I1, *I2, *I3;
  double K1 = 2*oo2n;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_00d0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_00d0(V_Cl, V_Cl[in].operands[1]);
  I2 = build_00p0(V_Cl, V_Cl[in].operands[2]);
  I3 = build_00p0(V_Cl, V_Cl[in].operands[3]);

  V_Cl[in].Val[0] = U[2][0]*I0[0] + U[5][0]*I1[0] + K1*(I2[0]-pon*I3[0]);
  V_Cl[in].Val[1] = U[2][0]*I0[1] + U[5][0]*I1[1] + oo2n*(I2[1]-pon*I3[1]);
  V_Cl[in].Val[2] = U[2][0]*I0[2] + U[5][0]*I1[2] + oo2n*(I2[2]-pon*I3[2]);
  V_Cl[in].Val[3] = U[2][0]*I0[3] + U[5][0]*I1[3];
  V_Cl[in].Val[4] = U[2][0]*I0[4] + U[5][0]*I1[4];
  V_Cl[in].Val[5] = U[2][0]*I0[5] + U[5][0]*I1[5];
  V_Cl[in].Val[6] = U[2][1]*I0[3] + U[5][1]*I1[3] + K1*(I2[1]-pon*I3[1]);
  V_Cl[in].Val[7] = U[2][1]*I0[4] + U[5][1]*I1[4] + oo2n*(I2[2]-pon*I3[2]);
  V_Cl[in].Val[8] = U[2][1]*I0[5] + U[5][1]*I1[5];
  V_Cl[in].Val[9] = U[2][2]*I0[5] + U[5][2]*I1[5] + K1*(I2[2]-pon*I3[2]);

  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_d0p0(iclass *V_Cl, int in) /* type = 9 */
{
  double *I0, *I1, *I2, *I3, *I4;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_p0p0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_p0p0(V_Cl, V_Cl[in].operands[1]);
  I2 = build_00p0(V_Cl, V_Cl[in].operands[2]);
  I3 = build_00p0(V_Cl, V_Cl[in].operands[3]);
  I4 = build_p000(V_Cl, V_Cl[in].operands[4]);

V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0] + oo2z*(I2[0] - poz*I3[0]) + oo2zn*I4[0];
V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1] + oo2z*(I2[1] - poz*I3[1]);
V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2] + oo2z*(I2[2] - poz*I3[2]);
V_Cl[in].Val[3] = U[0][0]*I0[3] + U[4][0]*I1[3] + oo2zn*I4[1];
V_Cl[in].Val[4] = U[0][0]*I0[4] + U[4][0]*I1[4];
V_Cl[in].Val[5] = U[0][0]*I0[5] + U[4][0]*I1[5];
V_Cl[in].Val[6] = U[0][0]*I0[6] + U[4][0]*I1[6] + oo2zn*I4[2];
V_Cl[in].Val[7] = U[0][0]*I0[7] + U[4][0]*I1[7];
V_Cl[in].Val[8] = U[0][0]*I0[8] + U[4][0]*I1[8];
V_Cl[in].Val[9] = U[0][1]*I0[3] + U[4][1]*I1[3] + oo2z*(I2[0] - poz*I3[0]);
V_Cl[in].Val[10] = U[0][1]*I0[4] + U[4][1]*I1[4] + oo2z*(I2[1] - poz*I3[1]) + oo2zn*I4[1];
V_Cl[in].Val[11] = U[0][1]*I0[5] + U[4][1]*I1[5] + oo2z*(I2[2] - poz*I3[2]);
V_Cl[in].Val[12] = U[0][1]*I0[6] + U[4][1]*I1[6];
V_Cl[in].Val[13] = U[0][1]*I0[7] + U[4][1]*I1[7] + oo2zn*I4[2];
V_Cl[in].Val[14] = U[0][1]*I0[8] + U[4][1]*I1[8];
V_Cl[in].Val[15] = U[0][2]*I0[6] + U[4][2]*I1[6] + oo2z*(I2[0] - poz*I3[0]);
V_Cl[in].Val[16] = U[0][2]*I0[7] + U[4][2]*I1[7] + oo2z*(I2[1] - poz*I3[1]);
V_Cl[in].Val[17] = U[0][2]*I0[8] + U[4][2]*I1[8] + oo2z*(I2[2] - poz*I3[2]) + oo2zn*I4[2];

  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_p0d0(iclass *V_Cl, int in) /* type = 10 */
{
  double *I0, *I1, *I4;
  double K1 = 2.0*oo2zn;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_00d0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_00d0(V_Cl, V_Cl[in].operands[1]);
  I4 = build_p000(V_Cl, V_Cl[in].operands[4]);

V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0] + K1*I4[0];
V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1] + oo2zn*I4[1];
V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2] + oo2zn*I4[2];
V_Cl[in].Val[3] = U[0][0]*I0[3] + U[4][0]*I1[3];
V_Cl[in].Val[4] = U[0][0]*I0[4] + U[4][0]*I1[4];
V_Cl[in].Val[5] = U[0][0]*I0[5] + U[4][0]*I1[5];
V_Cl[in].Val[6] = U[0][1]*I0[0] + U[4][1]*I1[0];
V_Cl[in].Val[7] = U[0][1]*I0[1] + U[4][1]*I1[1] + oo2zn*I4[0];
V_Cl[in].Val[8] = U[0][1]*I0[2] + U[4][1]*I1[2];
V_Cl[in].Val[9] = U[0][1]*I0[3] + U[4][1]*I1[3] + K1*I4[1];
V_Cl[in].Val[10] = U[0][1]*I0[4] + U[4][1]*I1[4] + oo2zn*I4[2];
V_Cl[in].Val[11] = U[0][1]*I0[5] + U[4][1]*I1[5];
V_Cl[in].Val[12] = U[0][2]*I0[0] + U[4][2]*I1[0];
V_Cl[in].Val[13] = U[0][2]*I0[1] + U[4][2]*I1[1];
V_Cl[in].Val[14] = U[0][2]*I0[2] + U[4][2]*I1[2] + oo2zn*I4[0];
V_Cl[in].Val[15] = U[0][2]*I0[3] + U[4][2]*I1[3];
V_Cl[in].Val[16] = U[0][2]*I0[4] + U[4][2]*I1[4] + oo2zn*I4[1];
V_Cl[in].Val[17] = U[0][2]*I0[5] + U[4][2]*I1[5] + K1*I4[2];


  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_f0p0(iclass *V_Cl, int in) /* type = 11 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double K1 = 2.0*oo2z;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_d0p0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_d0p0(V_Cl, V_Cl[in].operands[1]);
  I2 = build_p0p0(V_Cl, V_Cl[in].operands[2]);
  I3 = build_p0p0(V_Cl, V_Cl[in].operands[3]);
  I4 = build_d000(V_Cl, V_Cl[in].operands[4]);

V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0] + K1*(I2[0] - poz*I3[0]) +
         oo2zn*I4[0];
V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1] + K1*(I2[1] - poz*I3[1]);
V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2] + K1*(I2[2] - poz*I3[2]);
V_Cl[in].Val[3] = U[0][0]*I0[3] + U[4][0]*I1[3] + oo2z*(I2[3] - poz*I3[3]) +
         oo2zn*I4[1];
V_Cl[in].Val[4] = U[0][0]*I0[4] + U[4][0]*I1[4] + oo2z*(I2[4] - poz*I3[4]);
V_Cl[in].Val[5] = U[0][0]*I0[5] + U[4][0]*I1[5] + oo2z*(I2[5] - poz*I3[5]);
V_Cl[in].Val[6] = U[0][0]*I0[6] + U[4][0]*I1[6] + oo2z*(I2[6] - poz*I3[6]) +
         oo2zn*I4[2];
V_Cl[in].Val[7] = U[0][0]*I0[7] + U[4][0]*I1[7] + oo2z*(I2[7] - poz*I3[7]);
V_Cl[in].Val[8] = U[0][0]*I0[8] + U[4][0]*I1[8] + oo2z*(I2[8] - poz*I3[8]);
V_Cl[in].Val[9] = U[0][0]*I0[9] + U[4][0]*I1[9] + oo2zn*I4[3];
V_Cl[in].Val[10] = U[0][0]*I0[10] + U[4][0]*I1[10];
V_Cl[in].Val[11] = U[0][0]*I0[11] + U[4][0]*I1[11];
V_Cl[in].Val[12] = U[0][0]*I0[12] + U[4][0]*I1[12] + oo2zn*I4[4];
V_Cl[in].Val[13] = U[0][0]*I0[13] + U[4][0]*I1[13];
V_Cl[in].Val[14] = U[0][0]*I0[14] + U[4][0]*I1[14];
V_Cl[in].Val[15] = U[0][0]*I0[15] + U[4][0]*I1[15] + oo2zn*I4[5];
V_Cl[in].Val[16] = U[0][0]*I0[16] + U[4][0]*I1[16];
V_Cl[in].Val[17] = U[0][0]*I0[17] + U[4][0]*I1[17];
V_Cl[in].Val[18] = U[0][1]*I0[9] + U[4][1]*I1[9] + K1*(I2[3] - poz*I3[3]);
V_Cl[in].Val[19] = U[0][1]*I0[10] + U[4][1]*I1[10] + K1*(I2[4] - poz*I3[4]) +
         oo2zn*I4[3];
V_Cl[in].Val[20] = U[0][1]*I0[11] + U[4][1]*I1[11] + K1*(I2[5] - poz*I3[5]);
V_Cl[in].Val[21] = U[0][1]*I0[12] + U[4][1]*I1[12] + oo2z*(I2[6] - poz*I3[6]);
V_Cl[in].Val[22] = U[0][1]*I0[13] + U[4][1]*I1[13] + oo2z*(I2[7] - poz*I3[7]) +
         oo2zn*I4[4];
V_Cl[in].Val[23] = U[0][1]*I0[14] + U[4][1]*I1[14] + oo2z*(I2[8] - poz*I3[8]);
V_Cl[in].Val[24] = U[0][1]*I0[15] + U[4][1]*I1[15];
V_Cl[in].Val[25] = U[0][1]*I0[16] + U[4][1]*I1[16] + oo2zn*I4[5];
V_Cl[in].Val[26] = U[0][1]*I0[17] + U[4][1]*I1[17];
V_Cl[in].Val[27] = U[0][2]*I0[15] + U[4][2]*I1[15] + K1*(I2[6] - poz*I3[6]);
V_Cl[in].Val[28] = U[0][2]*I0[16] + U[4][2]*I1[16] + K1*(I2[7] - poz*I3[7]);
V_Cl[in].Val[29] = U[0][2]*I0[17] + U[4][2]*I1[17] + K1*(I2[8] - poz*I3[8]) +
         oo2zn*I4[5];

  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_p0f0(iclass *V_Cl, int in) /* type = 12 */
{
  double *I0, *I1, *I4;
  double K1 = 2.0*oo2zn;
  double K2 = 3.0*oo2zn;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_00f0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_00f0(V_Cl, V_Cl[in].operands[1]);
  I4 = build_00d0(V_Cl, V_Cl[in].operands[4]);

V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0] + K2*I4[0];
V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1] + K1*I4[1];
V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2] + K1*I4[2];
V_Cl[in].Val[3] = U[0][0]*I0[3] + U[4][0]*I1[3] + oo2zn*I4[3];
V_Cl[in].Val[4] = U[0][0]*I0[4] + U[4][0]*I1[4] + oo2zn*I4[4];
V_Cl[in].Val[5] = U[0][0]*I0[5] + U[4][0]*I1[5] + oo2zn*I4[5];
V_Cl[in].Val[6] = U[0][0]*I0[6] + U[4][0]*I1[6];
V_Cl[in].Val[7] = U[0][0]*I0[7] + U[4][0]*I1[7];
V_Cl[in].Val[8] = U[0][0]*I0[8] + U[4][0]*I1[8];
V_Cl[in].Val[9] = U[0][0]*I0[9] + U[4][0]*I1[9];
V_Cl[in].Val[10] = U[0][1]*I0[0] + U[4][1]*I1[0];
V_Cl[in].Val[11] = U[0][1]*I0[1] + U[4][1]*I1[1] + oo2zn*I4[0];
V_Cl[in].Val[12] = U[0][1]*I0[2] + U[4][1]*I1[2];
V_Cl[in].Val[13] = U[0][1]*I0[3] + U[4][1]*I1[3] + K1*I4[1];
V_Cl[in].Val[14] = U[0][1]*I0[4] + U[4][1]*I1[4] + oo2zn*I4[2];
V_Cl[in].Val[15] = U[0][1]*I0[5] + U[4][1]*I1[5];
V_Cl[in].Val[16] = U[0][1]*I0[6] + U[4][1]*I1[6] + K2*I4[3];
V_Cl[in].Val[17] = U[0][1]*I0[7] + U[4][1]*I1[7] + K1*I4[4];
V_Cl[in].Val[18] = U[0][1]*I0[8] + U[4][1]*I1[8] + oo2zn*I4[5];
V_Cl[in].Val[19] = U[0][1]*I0[9] + U[4][1]*I1[9];
V_Cl[in].Val[20] = U[0][2]*I0[0] + U[4][2]*I1[0];
V_Cl[in].Val[21] = U[0][2]*I0[1] + U[4][2]*I1[1];
V_Cl[in].Val[22] = U[0][2]*I0[2] + U[4][2]*I1[2] + oo2zn*I4[0];
V_Cl[in].Val[23] = U[0][2]*I0[3] + U[4][2]*I1[3];
V_Cl[in].Val[24] = U[0][2]*I0[4] + U[4][2]*I1[4] + oo2zn*I4[1];
V_Cl[in].Val[25] = U[0][2]*I0[5] + U[4][2]*I1[5] + K1*I4[2];
V_Cl[in].Val[26] = U[0][2]*I0[6] + U[4][2]*I1[6];
V_Cl[in].Val[27] = U[0][2]*I0[7] + U[4][2]*I1[7] + oo2zn*I4[3];
V_Cl[in].Val[28] = U[0][2]*I0[8] + U[4][2]*I1[8] + K1*I4[4];
V_Cl[in].Val[29] = U[0][2]*I0[9] + U[4][2]*I1[9] + K2*I4[5];

  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_d0d0(iclass *V_Cl, int in) /* type = 13 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double K1 = 2.0*oo2zn;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_p0d0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_p0d0(V_Cl, V_Cl[in].operands[1]);
  I2 = build_00d0(V_Cl, V_Cl[in].operands[2]);
  I3 = build_00d0(V_Cl, V_Cl[in].operands[3]);
  I4 = build_p0p0(V_Cl, V_Cl[in].operands[4]);

V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0] + oo2z*(I2[0] - poz*I3[0]) +
         K1*I4[0];
V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1] + oo2z*(I2[1] - poz*I3[1]) +
         oo2zn*I4[1];
V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2] + oo2z*(I2[2] - poz*I3[2]) +
         oo2zn*I4[2];
V_Cl[in].Val[3] = U[0][0]*I0[3] + U[4][0]*I1[3] + oo2z*(I2[3] - poz*I3[3]);
V_Cl[in].Val[4] = U[0][0]*I0[4] + U[4][0]*I1[4] + oo2z*(I2[4] - poz*I3[4]);
V_Cl[in].Val[5] = U[0][0]*I0[5] + U[4][0]*I1[5] + oo2z*(I2[5] - poz*I3[5]);
V_Cl[in].Val[6] = U[0][0]*I0[6] + U[4][0]*I1[6] + K1*I4[3];
V_Cl[in].Val[7] = U[0][0]*I0[7] + U[4][0]*I1[7] + oo2zn*I4[4];
V_Cl[in].Val[8] = U[0][0]*I0[8] + U[4][0]*I1[8] + oo2zn*I4[5];
V_Cl[in].Val[9] = U[0][0]*I0[9] + U[4][0]*I1[9];
V_Cl[in].Val[10] = U[0][0]*I0[10] + U[4][0]*I1[10];
V_Cl[in].Val[11] = U[0][0]*I0[11] + U[4][0]*I1[11];
V_Cl[in].Val[12] = U[0][0]*I0[12] + U[4][0]*I1[12] + K1*I4[6];
V_Cl[in].Val[13] = U[0][0]*I0[13] + U[4][0]*I1[13] + oo2zn*I4[7];
V_Cl[in].Val[14] = U[0][0]*I0[14] + U[4][0]*I1[14] + oo2zn*I4[8];
V_Cl[in].Val[15] = U[0][0]*I0[15] + U[4][0]*I1[15];
V_Cl[in].Val[16] = U[0][0]*I0[16] + U[4][0]*I1[16];
V_Cl[in].Val[17] = U[0][0]*I0[17] + U[4][0]*I1[17];
V_Cl[in].Val[18] = U[0][1]*I0[6] + U[4][1]*I1[6] + oo2z*(I2[0] - poz*I3[0]);
V_Cl[in].Val[19] = U[0][1]*I0[7] + U[4][1]*I1[7] + oo2z*(I2[1] - poz*I3[1]) +
         oo2zn*I4[3];
V_Cl[in].Val[20] = U[0][1]*I0[8] + U[4][1]*I1[8] + oo2z*(I2[2] - poz*I3[2]);
V_Cl[in].Val[21] = U[0][1]*I0[9] + U[4][1]*I1[9] + oo2z*(I2[3] - poz*I3[3]) +
         K1*I4[4];
V_Cl[in].Val[22] = U[0][1]*I0[10] + U[4][1]*I1[10] + oo2z*(I2[4] - poz*I3[4]) +
         oo2zn*I4[5];
V_Cl[in].Val[23] = U[0][1]*I0[11] + U[4][1]*I1[11] + oo2z*(I2[5] - poz*I3[5]);
V_Cl[in].Val[24] = U[0][1]*I0[12] + U[4][1]*I1[12];
V_Cl[in].Val[25] = U[0][1]*I0[13] + U[4][1]*I1[13] + oo2zn*I4[6];
V_Cl[in].Val[26] = U[0][1]*I0[14] + U[4][1]*I1[14];
V_Cl[in].Val[27] = U[0][1]*I0[15] + U[4][1]*I1[15] + K1*I4[7];
V_Cl[in].Val[28] = U[0][1]*I0[16] + U[4][1]*I1[16] + oo2zn*I4[8];
V_Cl[in].Val[29] = U[0][1]*I0[17] + U[4][1]*I1[17];
V_Cl[in].Val[30] = U[0][2]*I0[12] + U[4][2]*I1[12] + oo2z*(I2[0] - poz*I3[0]);
V_Cl[in].Val[31] = U[0][2]*I0[13] + U[4][2]*I1[13] + oo2z*(I2[1] - poz*I3[1]);
V_Cl[in].Val[32] = U[0][2]*I0[14] + U[4][2]*I1[14] + oo2z*(I2[2] - poz*I3[2]) +
         oo2zn*I4[6];
V_Cl[in].Val[33] = U[0][2]*I0[15] + U[4][2]*I1[15] + oo2z*(I2[3] - poz*I3[3]);
V_Cl[in].Val[34] = U[0][2]*I0[16] + U[4][2]*I1[16] + oo2z*(I2[4] - poz*I3[4]) +
         oo2zn*I4[7];
V_Cl[in].Val[35] = U[0][2]*I0[17] + U[4][2]*I1[17] + oo2z*(I2[5] - poz*I3[5]) +
         K1*I4[8];

  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_f0d0(iclass *V_Cl, int in) /* type = 14 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2z = 2.0*oo2z;
  double twooo2zn = 2.0*oo2zn;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_d0d0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_d0d0(V_Cl, V_Cl[in].operands[1]);
  I2 = build_p0p0(V_Cl, V_Cl[in].operands[2]);
  I3 = build_p0p0(V_Cl, V_Cl[in].operands[3]);
  I4 = build_d0p0(V_Cl, V_Cl[in].operands[4]);


V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0]
           + (twooo2z)*(I2[0] - (poz)*I3[0])
           + (twooo2zn)*I4[0];
V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1]
           + (twooo2z)*(I2[1] - (poz)*I3[1])
           + (oo2zn)*I4[1];
V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2]
           + (twooo2z)*(I2[2] - (poz)*I3[2])
           + (oo2zn)*I4[2];
V_Cl[in].Val[3] = U[0][0]*I0[3] + U[4][0]*I1[3]
           + (twooo2z)*(I2[3] - (poz)*I3[3]);
V_Cl[in].Val[4] = U[0][0]*I0[4] + U[4][0]*I1[4]
           + (twooo2z)*(I2[4] - (poz)*I3[4]);
V_Cl[in].Val[5] = U[0][0]*I0[5] + U[4][0]*I1[5]
           + (twooo2z)*(I2[5] - (poz)*I3[5]);
V_Cl[in].Val[6] = U[0][0]*I0[6] + U[4][0]*I1[6]
           + (oo2z)*(I2[6] - (poz)*I3[6])
           + (twooo2zn)*I4[3];
V_Cl[in].Val[7] = U[0][0]*I0[7] + U[4][0]*I1[7]
           + (oo2z)*(I2[7] - (poz)*I3[7])
           + (oo2zn)*I4[4];
V_Cl[in].Val[8] = U[0][0]*I0[8] + U[4][0]*I1[8]
           + (oo2z)*(I2[8] - (poz)*I3[8])
           + (oo2zn)*I4[5];
V_Cl[in].Val[9] = U[0][0]*I0[9] + U[4][0]*I1[9]
           + (oo2z)*(I2[9] - (poz)*I3[9]);
V_Cl[in].Val[10] = U[0][0]*I0[10] + U[4][0]*I1[10]
           + (oo2z)*(I2[10] - (poz)*I3[10]);
V_Cl[in].Val[11] = U[0][0]*I0[11] + U[4][0]*I1[11]
           + (oo2z)*(I2[11] - (poz)*I3[11]);
V_Cl[in].Val[12] = U[0][0]*I0[12] + U[4][0]*I1[12]
           + (oo2z)*(I2[12] - (poz)*I3[12])
           + (twooo2zn)*I4[6];
V_Cl[in].Val[13] = U[0][0]*I0[13] + U[4][0]*I1[13]
           + (oo2z)*(I2[13] - (poz)*I3[13])
           + (oo2zn)*I4[7];
V_Cl[in].Val[14] = U[0][0]*I0[14] + U[4][0]*I1[14]
           + (oo2z)*(I2[14] - (poz)*I3[14])
           + (oo2zn)*I4[8];
V_Cl[in].Val[15] = U[0][0]*I0[15] + U[4][0]*I1[15]
           + (oo2z)*(I2[15] - (poz)*I3[15]);
V_Cl[in].Val[16] = U[0][0]*I0[16] + U[4][0]*I1[16]
           + (oo2z)*(I2[16] - (poz)*I3[16]);
V_Cl[in].Val[17] = U[0][0]*I0[17] + U[4][0]*I1[17]
           + (oo2z)*(I2[17] - (poz)*I3[17]);
V_Cl[in].Val[18] = U[0][0]*I0[18] + U[4][0]*I1[18]
           + (twooo2zn)*I4[9];
V_Cl[in].Val[19] = U[0][0]*I0[19] + U[4][0]*I1[19]
           + (oo2zn)*I4[10];
V_Cl[in].Val[20] = U[0][0]*I0[20] + U[4][0]*I1[20]
           + (oo2zn)*I4[11];
V_Cl[in].Val[21] = U[0][0]*I0[21] + U[4][0]*I1[21];
V_Cl[in].Val[22] = U[0][0]*I0[22] + U[4][0]*I1[22];
V_Cl[in].Val[23] = U[0][0]*I0[23] + U[4][0]*I1[23];
V_Cl[in].Val[24] = U[0][0]*I0[24] + U[4][0]*I1[24]
           + (twooo2zn)*I4[12];
V_Cl[in].Val[25] = U[0][0]*I0[25] + U[4][0]*I1[25]
           + (oo2zn)*I4[13];
V_Cl[in].Val[26] = U[0][0]*I0[26] + U[4][0]*I1[26]
           + (oo2zn)*I4[14];
V_Cl[in].Val[27] = U[0][0]*I0[27] + U[4][0]*I1[27];
V_Cl[in].Val[28] = U[0][0]*I0[28] + U[4][0]*I1[28];
V_Cl[in].Val[29] = U[0][0]*I0[29] + U[4][0]*I1[29];
V_Cl[in].Val[30] = U[0][0]*I0[30] + U[4][0]*I1[30]
           + (twooo2zn)*I4[15];
V_Cl[in].Val[31] = U[0][0]*I0[31] + U[4][0]*I1[31]
           + (oo2zn)*I4[16];
V_Cl[in].Val[32] = U[0][0]*I0[32] + U[4][0]*I1[32]
           + (oo2zn)*I4[17];
V_Cl[in].Val[33] = U[0][0]*I0[33] + U[4][0]*I1[33];
V_Cl[in].Val[34] = U[0][0]*I0[34] + U[4][0]*I1[34];
V_Cl[in].Val[35] = U[0][0]*I0[35] + U[4][0]*I1[35];
V_Cl[in].Val[36] = U[0][1]*I0[18] + U[4][1]*I1[18]
           + (twooo2z)*(I2[6] - (poz)*I3[6]);
V_Cl[in].Val[37] = U[0][1]*I0[19] + U[4][1]*I1[19]
           + (twooo2z)*(I2[7] - (poz)*I3[7])
           + (oo2zn)*I4[9];
V_Cl[in].Val[38] = U[0][1]*I0[20] + U[4][1]*I1[20]
           + (twooo2z)*(I2[8] - (poz)*I3[8]);
V_Cl[in].Val[39] = U[0][1]*I0[21] + U[4][1]*I1[21]
           + (twooo2z)*(I2[9] - (poz)*I3[9])
           + (twooo2zn)*I4[10];
V_Cl[in].Val[40] = U[0][1]*I0[22] + U[4][1]*I1[22]
           + (twooo2z)*(I2[10] - (poz)*I3[10])
           + (oo2zn)*I4[11];
V_Cl[in].Val[41] = U[0][1]*I0[23] + U[4][1]*I1[23]
           + (twooo2z)*(I2[11] - (poz)*I3[11]);
V_Cl[in].Val[42] = U[0][1]*I0[24] + U[4][1]*I1[24]
           + (oo2z)*(I2[12] - (poz)*I3[12]);
V_Cl[in].Val[43] = U[0][1]*I0[25] + U[4][1]*I1[25]
           + (oo2z)*(I2[13] - (poz)*I3[13])
           + (oo2zn)*I4[12];
V_Cl[in].Val[44] = U[0][1]*I0[26] + U[4][1]*I1[26]
           + (oo2z)*(I2[14] - (poz)*I3[14]);
V_Cl[in].Val[45] = U[0][1]*I0[27] + U[4][1]*I1[27]
           + (oo2z)*(I2[15] - (poz)*I3[15])
           + (twooo2zn)*I4[13];
V_Cl[in].Val[46] = U[0][1]*I0[28] + U[4][1]*I1[28]
           + (oo2z)*(I2[16] - (poz)*I3[16])
           + (oo2zn)*I4[14];
V_Cl[in].Val[47] = U[0][1]*I0[29] + U[4][1]*I1[29]
           + (oo2z)*(I2[17] - (poz)*I3[17]);
V_Cl[in].Val[48] = U[0][1]*I0[30] + U[4][1]*I1[30];
V_Cl[in].Val[49] = U[0][1]*I0[31] + U[4][1]*I1[31]
           + (oo2zn)*I4[15];
V_Cl[in].Val[50] = U[0][1]*I0[32] + U[4][1]*I1[32];
V_Cl[in].Val[51] = U[0][1]*I0[33] + U[4][1]*I1[33]
           + (twooo2zn)*I4[16];
V_Cl[in].Val[52] = U[0][1]*I0[34] + U[4][1]*I1[34]
           + (oo2zn)*I4[17];
V_Cl[in].Val[53] = U[0][1]*I0[35] + U[4][1]*I1[35];
V_Cl[in].Val[54] = U[0][2]*I0[30] + U[4][2]*I1[30]
           + (twooo2z)*(I2[12] - (poz)*I3[12]);
V_Cl[in].Val[55] = U[0][2]*I0[31] + U[4][2]*I1[31]
           + (twooo2z)*(I2[13] - (poz)*I3[13]);
V_Cl[in].Val[56] = U[0][2]*I0[32] + U[4][2]*I1[32]
           + (twooo2z)*(I2[14] - (poz)*I3[14])
           + (oo2zn)*I4[15];
V_Cl[in].Val[57] = U[0][2]*I0[33] + U[4][2]*I1[33]
           + (twooo2z)*(I2[15] - (poz)*I3[15]);
V_Cl[in].Val[58] = U[0][2]*I0[34] + U[4][2]*I1[34]
           + (twooo2z)*(I2[16] - (poz)*I3[16])
           + (oo2zn)*I4[16];
V_Cl[in].Val[59] = U[0][2]*I0[35] + U[4][2]*I1[35]
           + (twooo2z)*(I2[17] - (poz)*I3[17])
           + (twooo2zn)*I4[17];
  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_d0f0(iclass *V_Cl, int in) /* type = 15 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2zn = 2.0*oo2zn;
  double threeoo2zn = 3.0*oo2zn;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_p0f0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_p0f0(V_Cl, V_Cl[in].operands[1]);
  I2 = build_00f0(V_Cl, V_Cl[in].operands[2]);
  I3 = build_00f0(V_Cl, V_Cl[in].operands[3]);
  I4 = build_p0d0(V_Cl, V_Cl[in].operands[4]);


V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0]
           + (oo2z)*(I2[0] - (poz)*I3[0])
           + (threeoo2zn)*I4[0];
V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1]
           + (oo2z)*(I2[1] - (poz)*I3[1])
           + (twooo2zn)*I4[1];
V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2]
           + (oo2z)*(I2[2] - (poz)*I3[2])
           + (twooo2zn)*I4[2];
V_Cl[in].Val[3] = U[0][0]*I0[3] + U[4][0]*I1[3]
           + (oo2z)*(I2[3] - (poz)*I3[3])
           + (oo2zn)*I4[3];
V_Cl[in].Val[4] = U[0][0]*I0[4] + U[4][0]*I1[4]
           + (oo2z)*(I2[4] - (poz)*I3[4])
           + (oo2zn)*I4[4];
V_Cl[in].Val[5] = U[0][0]*I0[5] + U[4][0]*I1[5]
           + (oo2z)*(I2[5] - (poz)*I3[5])
           + (oo2zn)*I4[5];
V_Cl[in].Val[6] = U[0][0]*I0[6] + U[4][0]*I1[6]
           + (oo2z)*(I2[6] - (poz)*I3[6]);
V_Cl[in].Val[7] = U[0][0]*I0[7] + U[4][0]*I1[7]
           + (oo2z)*(I2[7] - (poz)*I3[7]);
V_Cl[in].Val[8] = U[0][0]*I0[8] + U[4][0]*I1[8]
           + (oo2z)*(I2[8] - (poz)*I3[8]);
V_Cl[in].Val[9] = U[0][0]*I0[9] + U[4][0]*I1[9]
           + (oo2z)*(I2[9] - (poz)*I3[9]);
V_Cl[in].Val[10] = U[0][0]*I0[10] + U[4][0]*I1[10]
           + (threeoo2zn)*I4[6];
V_Cl[in].Val[11] = U[0][0]*I0[11] + U[4][0]*I1[11]
           + (twooo2zn)*I4[7];
V_Cl[in].Val[12] = U[0][0]*I0[12] + U[4][0]*I1[12]
           + (twooo2zn)*I4[8];
V_Cl[in].Val[13] = U[0][0]*I0[13] + U[4][0]*I1[13]
           + (oo2zn)*I4[9];
V_Cl[in].Val[14] = U[0][0]*I0[14] + U[4][0]*I1[14]
           + (oo2zn)*I4[10];
V_Cl[in].Val[15] = U[0][0]*I0[15] + U[4][0]*I1[15]
           + (oo2zn)*I4[11];
V_Cl[in].Val[16] = U[0][0]*I0[16] + U[4][0]*I1[16];
V_Cl[in].Val[17] = U[0][0]*I0[17] + U[4][0]*I1[17];
V_Cl[in].Val[18] = U[0][0]*I0[18] + U[4][0]*I1[18];
V_Cl[in].Val[19] = U[0][0]*I0[19] + U[4][0]*I1[19];
V_Cl[in].Val[20] = U[0][0]*I0[20] + U[4][0]*I1[20]
           + (threeoo2zn)*I4[12];
V_Cl[in].Val[21] = U[0][0]*I0[21] + U[4][0]*I1[21]
           + (twooo2zn)*I4[13];
V_Cl[in].Val[22] = U[0][0]*I0[22] + U[4][0]*I1[22]
           + (twooo2zn)*I4[14];
V_Cl[in].Val[23] = U[0][0]*I0[23] + U[4][0]*I1[23]
           + (oo2zn)*I4[15];
V_Cl[in].Val[24] = U[0][0]*I0[24] + U[4][0]*I1[24]
           + (oo2zn)*I4[16];
V_Cl[in].Val[25] = U[0][0]*I0[25] + U[4][0]*I1[25]
           + (oo2zn)*I4[17];
V_Cl[in].Val[26] = U[0][0]*I0[26] + U[4][0]*I1[26];
V_Cl[in].Val[27] = U[0][0]*I0[27] + U[4][0]*I1[27];
V_Cl[in].Val[28] = U[0][0]*I0[28] + U[4][0]*I1[28];
V_Cl[in].Val[29] = U[0][0]*I0[29] + U[4][0]*I1[29];
V_Cl[in].Val[30] = U[0][1]*I0[10] + U[4][1]*I1[10]
           + (oo2z)*(I2[0] - (poz)*I3[0]);
V_Cl[in].Val[31] = U[0][1]*I0[11] + U[4][1]*I1[11]
           + (oo2z)*(I2[1] - (poz)*I3[1])
           + (oo2zn)*I4[6];
V_Cl[in].Val[32] = U[0][1]*I0[12] + U[4][1]*I1[12]
           + (oo2z)*(I2[2] - (poz)*I3[2]);
V_Cl[in].Val[33] = U[0][1]*I0[13] + U[4][1]*I1[13]
           + (oo2z)*(I2[3] - (poz)*I3[3])
           + (twooo2zn)*I4[7];
V_Cl[in].Val[34] = U[0][1]*I0[14] + U[4][1]*I1[14]
           + (oo2z)*(I2[4] - (poz)*I3[4])
           + (oo2zn)*I4[8];
V_Cl[in].Val[35] = U[0][1]*I0[15] + U[4][1]*I1[15]
           + (oo2z)*(I2[5] - (poz)*I3[5]);
V_Cl[in].Val[36] = U[0][1]*I0[16] + U[4][1]*I1[16]
           + (oo2z)*(I2[6] - (poz)*I3[6])
           + (threeoo2zn)*I4[9];
V_Cl[in].Val[37] = U[0][1]*I0[17] + U[4][1]*I1[17]
           + (oo2z)*(I2[7] - (poz)*I3[7])
           + (twooo2zn)*I4[10];
V_Cl[in].Val[38] = U[0][1]*I0[18] + U[4][1]*I1[18]
           + (oo2z)*(I2[8] - (poz)*I3[8])
           + (oo2zn)*I4[11];
V_Cl[in].Val[39] = U[0][1]*I0[19] + U[4][1]*I1[19]
           + (oo2z)*(I2[9] - (poz)*I3[9]);
V_Cl[in].Val[40] = U[0][1]*I0[20] + U[4][1]*I1[20];
V_Cl[in].Val[41] = U[0][1]*I0[21] + U[4][1]*I1[21]
           + (oo2zn)*I4[12];
V_Cl[in].Val[42] = U[0][1]*I0[22] + U[4][1]*I1[22];
V_Cl[in].Val[43] = U[0][1]*I0[23] + U[4][1]*I1[23]
           + (twooo2zn)*I4[13];
V_Cl[in].Val[44] = U[0][1]*I0[24] + U[4][1]*I1[24]
           + (oo2zn)*I4[14];
V_Cl[in].Val[45] = U[0][1]*I0[25] + U[4][1]*I1[25];
V_Cl[in].Val[46] = U[0][1]*I0[26] + U[4][1]*I1[26]
           + (threeoo2zn)*I4[15];
V_Cl[in].Val[47] = U[0][1]*I0[27] + U[4][1]*I1[27]
           + (twooo2zn)*I4[16];
V_Cl[in].Val[48] = U[0][1]*I0[28] + U[4][1]*I1[28]
           + (oo2zn)*I4[17];
V_Cl[in].Val[49] = U[0][1]*I0[29] + U[4][1]*I1[29];
V_Cl[in].Val[50] = U[0][2]*I0[20] + U[4][2]*I1[20]
           + (oo2z)*(I2[0] - (poz)*I3[0]);
V_Cl[in].Val[51] = U[0][2]*I0[21] + U[4][2]*I1[21]
           + (oo2z)*(I2[1] - (poz)*I3[1]);
V_Cl[in].Val[52] = U[0][2]*I0[22] + U[4][2]*I1[22]
           + (oo2z)*(I2[2] - (poz)*I3[2])
           + (oo2zn)*I4[12];
V_Cl[in].Val[53] = U[0][2]*I0[23] + U[4][2]*I1[23]
           + (oo2z)*(I2[3] - (poz)*I3[3]);
V_Cl[in].Val[54] = U[0][2]*I0[24] + U[4][2]*I1[24]
           + (oo2z)*(I2[4] - (poz)*I3[4])
           + (oo2zn)*I4[13];
V_Cl[in].Val[55] = U[0][2]*I0[25] + U[4][2]*I1[25]
           + (oo2z)*(I2[5] - (poz)*I3[5])
           + (twooo2zn)*I4[14];
V_Cl[in].Val[56] = U[0][2]*I0[26] + U[4][2]*I1[26]
           + (oo2z)*(I2[6] - (poz)*I3[6]);
V_Cl[in].Val[57] = U[0][2]*I0[27] + U[4][2]*I1[27]
           + (oo2z)*(I2[7] - (poz)*I3[7])
           + (oo2zn)*I4[15];
V_Cl[in].Val[58] = U[0][2]*I0[28] + U[4][2]*I1[28]
           + (oo2z)*(I2[8] - (poz)*I3[8])
           + (twooo2zn)*I4[16];
V_Cl[in].Val[59] = U[0][2]*I0[29] + U[4][2]*I1[29]
           + (oo2z)*(I2[9] - (poz)*I3[9])
           + (threeoo2zn)*I4[17];
  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_f0f0(iclass *V_Cl, int in) /* type = 16 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2z = 2.0*oo2z;
  double twooo2zn = 2.0*oo2zn;
  double threeoo2zn = 3.0*oo2zn;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_d0f0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_d0f0(V_Cl, V_Cl[in].operands[1]);
  I2 = build_p0f0(V_Cl, V_Cl[in].operands[2]);
  I3 = build_p0f0(V_Cl, V_Cl[in].operands[3]);
  I4 = build_d0d0(V_Cl, V_Cl[in].operands[4]);


V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0]
           + (twooo2z)*(I2[0] - (poz)*I3[0])
           + (threeoo2zn)*I4[0];
V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1]
           + (twooo2z)*(I2[1] - (poz)*I3[1])
           + (twooo2zn)*I4[1];
V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2]
           + (twooo2z)*(I2[2] - (poz)*I3[2])
           + (twooo2zn)*I4[2];
V_Cl[in].Val[3] = U[0][0]*I0[3] + U[4][0]*I1[3]
           + (twooo2z)*(I2[3] - (poz)*I3[3])
           + (oo2zn)*I4[3];
V_Cl[in].Val[4] = U[0][0]*I0[4] + U[4][0]*I1[4]
           + (twooo2z)*(I2[4] - (poz)*I3[4])
           + (oo2zn)*I4[4];
V_Cl[in].Val[5] = U[0][0]*I0[5] + U[4][0]*I1[5]
           + (twooo2z)*(I2[5] - (poz)*I3[5])
           + (oo2zn)*I4[5];
V_Cl[in].Val[6] = U[0][0]*I0[6] + U[4][0]*I1[6]
           + (twooo2z)*(I2[6] - (poz)*I3[6]);
V_Cl[in].Val[7] = U[0][0]*I0[7] + U[4][0]*I1[7]
           + (twooo2z)*(I2[7] - (poz)*I3[7]);
V_Cl[in].Val[8] = U[0][0]*I0[8] + U[4][0]*I1[8]
           + (twooo2z)*(I2[8] - (poz)*I3[8]);
V_Cl[in].Val[9] = U[0][0]*I0[9] + U[4][0]*I1[9]
           + (twooo2z)*(I2[9] - (poz)*I3[9]);
V_Cl[in].Val[10] = U[0][0]*I0[10] + U[4][0]*I1[10]
           + (oo2z)*(I2[10] - (poz)*I3[10])
           + (threeoo2zn)*I4[6];
V_Cl[in].Val[11] = U[0][0]*I0[11] + U[4][0]*I1[11]
           + (oo2z)*(I2[11] - (poz)*I3[11])
           + (twooo2zn)*I4[7];
V_Cl[in].Val[12] = U[0][0]*I0[12] + U[4][0]*I1[12]
           + (oo2z)*(I2[12] - (poz)*I3[12])
           + (twooo2zn)*I4[8];
V_Cl[in].Val[13] = U[0][0]*I0[13] + U[4][0]*I1[13]
           + (oo2z)*(I2[13] - (poz)*I3[13])
           + (oo2zn)*I4[9];
V_Cl[in].Val[14] = U[0][0]*I0[14] + U[4][0]*I1[14]
           + (oo2z)*(I2[14] - (poz)*I3[14])
           + (oo2zn)*I4[10];
V_Cl[in].Val[15] = U[0][0]*I0[15] + U[4][0]*I1[15]
           + (oo2z)*(I2[15] - (poz)*I3[15])
           + (oo2zn)*I4[11];
V_Cl[in].Val[16] = U[0][0]*I0[16] + U[4][0]*I1[16]
           + (oo2z)*(I2[16] - (poz)*I3[16]);
V_Cl[in].Val[17] = U[0][0]*I0[17] + U[4][0]*I1[17]
           + (oo2z)*(I2[17] - (poz)*I3[17]);
V_Cl[in].Val[18] = U[0][0]*I0[18] + U[4][0]*I1[18]
           + (oo2z)*(I2[18] - (poz)*I3[18]);
V_Cl[in].Val[19] = U[0][0]*I0[19] + U[4][0]*I1[19]
           + (oo2z)*(I2[19] - (poz)*I3[19]);
V_Cl[in].Val[20] = U[0][0]*I0[20] + U[4][0]*I1[20]
           + (oo2z)*(I2[20] - (poz)*I3[20])
           + (threeoo2zn)*I4[12];
V_Cl[in].Val[21] = U[0][0]*I0[21] + U[4][0]*I1[21]
           + (oo2z)*(I2[21] - (poz)*I3[21])
           + (twooo2zn)*I4[13];
V_Cl[in].Val[22] = U[0][0]*I0[22] + U[4][0]*I1[22]
           + (oo2z)*(I2[22] - (poz)*I3[22])
           + (twooo2zn)*I4[14];
V_Cl[in].Val[23] = U[0][0]*I0[23] + U[4][0]*I1[23]
           + (oo2z)*(I2[23] - (poz)*I3[23])
           + (oo2zn)*I4[15];
V_Cl[in].Val[24] = U[0][0]*I0[24] + U[4][0]*I1[24]
           + (oo2z)*(I2[24] - (poz)*I3[24])
           + (oo2zn)*I4[16];
V_Cl[in].Val[25] = U[0][0]*I0[25] + U[4][0]*I1[25]
           + (oo2z)*(I2[25] - (poz)*I3[25])
           + (oo2zn)*I4[17];
V_Cl[in].Val[26] = U[0][0]*I0[26] + U[4][0]*I1[26]
           + (oo2z)*(I2[26] - (poz)*I3[26]);
V_Cl[in].Val[27] = U[0][0]*I0[27] + U[4][0]*I1[27]
           + (oo2z)*(I2[27] - (poz)*I3[27]);
V_Cl[in].Val[28] = U[0][0]*I0[28] + U[4][0]*I1[28]
           + (oo2z)*(I2[28] - (poz)*I3[28]);
V_Cl[in].Val[29] = U[0][0]*I0[29] + U[4][0]*I1[29]
           + (oo2z)*(I2[29] - (poz)*I3[29]);
V_Cl[in].Val[30] = U[0][0]*I0[30] + U[4][0]*I1[30]
           + (threeoo2zn)*I4[18];
V_Cl[in].Val[31] = U[0][0]*I0[31] + U[4][0]*I1[31]
           + (twooo2zn)*I4[19];
V_Cl[in].Val[32] = U[0][0]*I0[32] + U[4][0]*I1[32]
           + (twooo2zn)*I4[20];
V_Cl[in].Val[33] = U[0][0]*I0[33] + U[4][0]*I1[33]
           + (oo2zn)*I4[21];
V_Cl[in].Val[34] = U[0][0]*I0[34] + U[4][0]*I1[34]
           + (oo2zn)*I4[22];
V_Cl[in].Val[35] = U[0][0]*I0[35] + U[4][0]*I1[35]
           + (oo2zn)*I4[23];
V_Cl[in].Val[36] = U[0][0]*I0[36] + U[4][0]*I1[36];
V_Cl[in].Val[37] = U[0][0]*I0[37] + U[4][0]*I1[37];
V_Cl[in].Val[38] = U[0][0]*I0[38] + U[4][0]*I1[38];
V_Cl[in].Val[39] = U[0][0]*I0[39] + U[4][0]*I1[39];
V_Cl[in].Val[40] = U[0][0]*I0[40] + U[4][0]*I1[40]
           + (threeoo2zn)*I4[24];
V_Cl[in].Val[41] = U[0][0]*I0[41] + U[4][0]*I1[41]
           + (twooo2zn)*I4[25];
V_Cl[in].Val[42] = U[0][0]*I0[42] + U[4][0]*I1[42]
           + (twooo2zn)*I4[26];
V_Cl[in].Val[43] = U[0][0]*I0[43] + U[4][0]*I1[43]
           + (oo2zn)*I4[27];
V_Cl[in].Val[44] = U[0][0]*I0[44] + U[4][0]*I1[44]
           + (oo2zn)*I4[28];
V_Cl[in].Val[45] = U[0][0]*I0[45] + U[4][0]*I1[45]
           + (oo2zn)*I4[29];
V_Cl[in].Val[46] = U[0][0]*I0[46] + U[4][0]*I1[46];
V_Cl[in].Val[47] = U[0][0]*I0[47] + U[4][0]*I1[47];
V_Cl[in].Val[48] = U[0][0]*I0[48] + U[4][0]*I1[48];
V_Cl[in].Val[49] = U[0][0]*I0[49] + U[4][0]*I1[49];
V_Cl[in].Val[50] = U[0][0]*I0[50] + U[4][0]*I1[50]
           + (threeoo2zn)*I4[30];
V_Cl[in].Val[51] = U[0][0]*I0[51] + U[4][0]*I1[51]
           + (twooo2zn)*I4[31];
V_Cl[in].Val[52] = U[0][0]*I0[52] + U[4][0]*I1[52]
           + (twooo2zn)*I4[32];
V_Cl[in].Val[53] = U[0][0]*I0[53] + U[4][0]*I1[53]
           + (oo2zn)*I4[33];
V_Cl[in].Val[54] = U[0][0]*I0[54] + U[4][0]*I1[54]
           + (oo2zn)*I4[34];
V_Cl[in].Val[55] = U[0][0]*I0[55] + U[4][0]*I1[55]
           + (oo2zn)*I4[35];
V_Cl[in].Val[56] = U[0][0]*I0[56] + U[4][0]*I1[56];
V_Cl[in].Val[57] = U[0][0]*I0[57] + U[4][0]*I1[57];
V_Cl[in].Val[58] = U[0][0]*I0[58] + U[4][0]*I1[58];
V_Cl[in].Val[59] = U[0][0]*I0[59] + U[4][0]*I1[59];
V_Cl[in].Val[60] = U[0][1]*I0[30] + U[4][1]*I1[30]
           + (twooo2z)*(I2[10] - (poz)*I3[10]);
V_Cl[in].Val[61] = U[0][1]*I0[31] + U[4][1]*I1[31]
           + (twooo2z)*(I2[11] - (poz)*I3[11])
           + (oo2zn)*I4[18];
V_Cl[in].Val[62] = U[0][1]*I0[32] + U[4][1]*I1[32]
           + (twooo2z)*(I2[12] - (poz)*I3[12]);
V_Cl[in].Val[63] = U[0][1]*I0[33] + U[4][1]*I1[33]
           + (twooo2z)*(I2[13] - (poz)*I3[13])
           + (twooo2zn)*I4[19];
V_Cl[in].Val[64] = U[0][1]*I0[34] + U[4][1]*I1[34]
           + (twooo2z)*(I2[14] - (poz)*I3[14])
           + (oo2zn)*I4[20];
V_Cl[in].Val[65] = U[0][1]*I0[35] + U[4][1]*I1[35]
           + (twooo2z)*(I2[15] - (poz)*I3[15]);
V_Cl[in].Val[66] = U[0][1]*I0[36] + U[4][1]*I1[36]
           + (twooo2z)*(I2[16] - (poz)*I3[16])
           + (threeoo2zn)*I4[21];
V_Cl[in].Val[67] = U[0][1]*I0[37] + U[4][1]*I1[37]
           + (twooo2z)*(I2[17] - (poz)*I3[17])
           + (twooo2zn)*I4[22];
V_Cl[in].Val[68] = U[0][1]*I0[38] + U[4][1]*I1[38]
           + (twooo2z)*(I2[18] - (poz)*I3[18])
           + (oo2zn)*I4[23];
V_Cl[in].Val[69] = U[0][1]*I0[39] + U[4][1]*I1[39]
           + (twooo2z)*(I2[19] - (poz)*I3[19]);
V_Cl[in].Val[70] = U[0][1]*I0[40] + U[4][1]*I1[40]
           + (oo2z)*(I2[20] - (poz)*I3[20]);
V_Cl[in].Val[71] = U[0][1]*I0[41] + U[4][1]*I1[41]
           + (oo2z)*(I2[21] - (poz)*I3[21])
           + (oo2zn)*I4[24];
V_Cl[in].Val[72] = U[0][1]*I0[42] + U[4][1]*I1[42]
           + (oo2z)*(I2[22] - (poz)*I3[22]);
V_Cl[in].Val[73] = U[0][1]*I0[43] + U[4][1]*I1[43]
           + (oo2z)*(I2[23] - (poz)*I3[23])
           + (twooo2zn)*I4[25];
V_Cl[in].Val[74] = U[0][1]*I0[44] + U[4][1]*I1[44]
           + (oo2z)*(I2[24] - (poz)*I3[24])
           + (oo2zn)*I4[26];
V_Cl[in].Val[75] = U[0][1]*I0[45] + U[4][1]*I1[45]
           + (oo2z)*(I2[25] - (poz)*I3[25]);
V_Cl[in].Val[76] = U[0][1]*I0[46] + U[4][1]*I1[46]
           + (oo2z)*(I2[26] - (poz)*I3[26])
           + (threeoo2zn)*I4[27];
V_Cl[in].Val[77] = U[0][1]*I0[47] + U[4][1]*I1[47]
           + (oo2z)*(I2[27] - (poz)*I3[27])
           + (twooo2zn)*I4[28];
V_Cl[in].Val[78] = U[0][1]*I0[48] + U[4][1]*I1[48]
           + (oo2z)*(I2[28] - (poz)*I3[28])
           + (oo2zn)*I4[29];
V_Cl[in].Val[79] = U[0][1]*I0[49] + U[4][1]*I1[49]
           + (oo2z)*(I2[29] - (poz)*I3[29]);
V_Cl[in].Val[80] = U[0][1]*I0[50] + U[4][1]*I1[50];
V_Cl[in].Val[81] = U[0][1]*I0[51] + U[4][1]*I1[51]
           + (oo2zn)*I4[30];
V_Cl[in].Val[82] = U[0][1]*I0[52] + U[4][1]*I1[52];
V_Cl[in].Val[83] = U[0][1]*I0[53] + U[4][1]*I1[53]
           + (twooo2zn)*I4[31];
V_Cl[in].Val[84] = U[0][1]*I0[54] + U[4][1]*I1[54]
           + (oo2zn)*I4[32];
V_Cl[in].Val[85] = U[0][1]*I0[55] + U[4][1]*I1[55];
V_Cl[in].Val[86] = U[0][1]*I0[56] + U[4][1]*I1[56]
           + (threeoo2zn)*I4[33];
V_Cl[in].Val[87] = U[0][1]*I0[57] + U[4][1]*I1[57]
           + (twooo2zn)*I4[34];
V_Cl[in].Val[88] = U[0][1]*I0[58] + U[4][1]*I1[58]
           + (oo2zn)*I4[35];
V_Cl[in].Val[89] = U[0][1]*I0[59] + U[4][1]*I1[59];
V_Cl[in].Val[90] = U[0][2]*I0[50] + U[4][2]*I1[50]
           + (twooo2z)*(I2[20] - (poz)*I3[20]);
V_Cl[in].Val[91] = U[0][2]*I0[51] + U[4][2]*I1[51]
           + (twooo2z)*(I2[21] - (poz)*I3[21]);
V_Cl[in].Val[92] = U[0][2]*I0[52] + U[4][2]*I1[52]
           + (twooo2z)*(I2[22] - (poz)*I3[22])
           + (oo2zn)*I4[30];
V_Cl[in].Val[93] = U[0][2]*I0[53] + U[4][2]*I1[53]
           + (twooo2z)*(I2[23] - (poz)*I3[23]);
V_Cl[in].Val[94] = U[0][2]*I0[54] + U[4][2]*I1[54]
           + (twooo2z)*(I2[24] - (poz)*I3[24])
           + (oo2zn)*I4[31];
V_Cl[in].Val[95] = U[0][2]*I0[55] + U[4][2]*I1[55]
           + (twooo2z)*(I2[25] - (poz)*I3[25])
           + (twooo2zn)*I4[32];
V_Cl[in].Val[96] = U[0][2]*I0[56] + U[4][2]*I1[56]
           + (twooo2z)*(I2[26] - (poz)*I3[26]);
V_Cl[in].Val[97] = U[0][2]*I0[57] + U[4][2]*I1[57]
           + (twooo2z)*(I2[27] - (poz)*I3[27])
           + (oo2zn)*I4[33];
V_Cl[in].Val[98] = U[0][2]*I0[58] + U[4][2]*I1[58]
           + (twooo2z)*(I2[28] - (poz)*I3[28])
           + (twooo2zn)*I4[34];
V_Cl[in].Val[99] = U[0][2]*I0[59] + U[4][2]*I1[59]
           + (twooo2z)*(I2[29] - (poz)*I3[29])
           + (threeoo2zn)*I4[35];
  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_g000(iclass *V_Cl, int in) /* type = 17 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2z = 2.0*oo2z;
  double threeoo2z = 3.0*oo2z;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_f000(V_Cl, V_Cl[in].operands[0]);
  I1 = build_f000(V_Cl, V_Cl[in].operands[1]);
  I2 = build_d000(V_Cl, V_Cl[in].operands[2]);
  I3 = build_d000(V_Cl, V_Cl[in].operands[3]);


V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0]
           + (threeoo2z)*(I2[0] - (poz)*I3[0]);
V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1]
           + (twooo2z)*(I2[1] - (poz)*I3[1]);
V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2]
           + (twooo2z)*(I2[2] - (poz)*I3[2]);
V_Cl[in].Val[3] = U[0][0]*I0[3] + U[4][0]*I1[3]
           + (oo2z)*(I2[3] - (poz)*I3[3]);
V_Cl[in].Val[4] = U[0][0]*I0[4] + U[4][0]*I1[4]
           + (oo2z)*(I2[4] - (poz)*I3[4]);
V_Cl[in].Val[5] = U[0][0]*I0[5] + U[4][0]*I1[5]
           + (oo2z)*(I2[5] - (poz)*I3[5]);
V_Cl[in].Val[6] = U[0][0]*I0[6] + U[4][0]*I1[6];
V_Cl[in].Val[7] = U[0][0]*I0[7] + U[4][0]*I1[7];
V_Cl[in].Val[8] = U[0][0]*I0[8] + U[4][0]*I1[8];
V_Cl[in].Val[9] = U[0][0]*I0[9] + U[4][0]*I1[9];
V_Cl[in].Val[10] = U[0][1]*I0[6] + U[4][1]*I1[6]
           + (threeoo2z)*(I2[3] - (poz)*I3[3]);
V_Cl[in].Val[11] = U[0][1]*I0[7] + U[4][1]*I1[7]
           + (twooo2z)*(I2[4] - (poz)*I3[4]);
V_Cl[in].Val[12] = U[0][1]*I0[8] + U[4][1]*I1[8]
           + (oo2z)*(I2[5] - (poz)*I3[5]);
V_Cl[in].Val[13] = U[0][1]*I0[9] + U[4][1]*I1[9];
V_Cl[in].Val[14] = U[0][2]*I0[9] + U[4][2]*I1[9]
           + (threeoo2z)*(I2[5] - (poz)*I3[5]);
  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_00g0(iclass *V_Cl, int in) /* type = 18 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2n = 2.0*oo2n;
  double threeoo2n = 3.0*oo2n;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_00f0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_00f0(V_Cl, V_Cl[in].operands[1]);
  I2 = build_00d0(V_Cl, V_Cl[in].operands[2]);
  I3 = build_00d0(V_Cl, V_Cl[in].operands[3]);


V_Cl[in].Val[0] = U[2][0]*I0[0] + U[5][0]*I1[0]
           + (threeoo2n)*(I2[0] - (pon)*I3[0]);
V_Cl[in].Val[1] = U[2][0]*I0[1] + U[5][0]*I1[1]
           + (twooo2n)*(I2[1] - (pon)*I3[1]);
V_Cl[in].Val[2] = U[2][0]*I0[2] + U[5][0]*I1[2]
           + (twooo2n)*(I2[2] - (pon)*I3[2]);
V_Cl[in].Val[3] = U[2][0]*I0[3] + U[5][0]*I1[3]
           + (oo2n)*(I2[3] - (pon)*I3[3]);
V_Cl[in].Val[4] = U[2][0]*I0[4] + U[5][0]*I1[4]
           + (oo2n)*(I2[4] - (pon)*I3[4]);
V_Cl[in].Val[5] = U[2][0]*I0[5] + U[5][0]*I1[5]
           + (oo2n)*(I2[5] - (pon)*I3[5]);
V_Cl[in].Val[6] = U[2][0]*I0[6] + U[5][0]*I1[6];
V_Cl[in].Val[7] = U[2][0]*I0[7] + U[5][0]*I1[7];
V_Cl[in].Val[8] = U[2][0]*I0[8] + U[5][0]*I1[8];
V_Cl[in].Val[9] = U[2][0]*I0[9] + U[5][0]*I1[9];
V_Cl[in].Val[10] = U[2][1]*I0[6] + U[5][1]*I1[6]
           + (threeoo2n)*(I2[3] - (pon)*I3[3]);
V_Cl[in].Val[11] = U[2][1]*I0[7] + U[5][1]*I1[7]
           + (twooo2n)*(I2[4] - (pon)*I3[4]);
V_Cl[in].Val[12] = U[2][1]*I0[8] + U[5][1]*I1[8]
           + (oo2n)*(I2[5] - (pon)*I3[5]);
V_Cl[in].Val[13] = U[2][1]*I0[9] + U[5][1]*I1[9];
V_Cl[in].Val[14] = U[2][2]*I0[9] + U[5][2]*I1[9]
           + (threeoo2n)*(I2[5] - (pon)*I3[5]);
  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_g0p0(iclass *V_Cl, int in) /* type = 19 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2z = 2.0*oo2z;
  double threeoo2z = 3.0*oo2z;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_f0p0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_f0p0(V_Cl, V_Cl[in].operands[1]);
  I2 = build_d0p0(V_Cl, V_Cl[in].operands[2]);
  I3 = build_d0p0(V_Cl, V_Cl[in].operands[3]);
  I4 = build_f000(V_Cl, V_Cl[in].operands[4]);


V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0]
           + (threeoo2z)*(I2[0] - (poz)*I3[0])
           + (oo2zn)*I4[0];
V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1]
           + (threeoo2z)*(I2[1] - (poz)*I3[1]);
V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2]
           + (threeoo2z)*(I2[2] - (poz)*I3[2]);
V_Cl[in].Val[3] = U[0][0]*I0[3] + U[4][0]*I1[3]
           + (twooo2z)*(I2[3] - (poz)*I3[3])
           + (oo2zn)*I4[1];
V_Cl[in].Val[4] = U[0][0]*I0[4] + U[4][0]*I1[4]
           + (twooo2z)*(I2[4] - (poz)*I3[4]);
V_Cl[in].Val[5] = U[0][0]*I0[5] + U[4][0]*I1[5]
           + (twooo2z)*(I2[5] - (poz)*I3[5]);
V_Cl[in].Val[6] = U[0][0]*I0[6] + U[4][0]*I1[6]
           + (twooo2z)*(I2[6] - (poz)*I3[6])
           + (oo2zn)*I4[2];
V_Cl[in].Val[7] = U[0][0]*I0[7] + U[4][0]*I1[7]
           + (twooo2z)*(I2[7] - (poz)*I3[7]);
V_Cl[in].Val[8] = U[0][0]*I0[8] + U[4][0]*I1[8]
           + (twooo2z)*(I2[8] - (poz)*I3[8]);
V_Cl[in].Val[9] = U[0][0]*I0[9] + U[4][0]*I1[9]
           + (oo2z)*(I2[9] - (poz)*I3[9])
           + (oo2zn)*I4[3];
V_Cl[in].Val[10] = U[0][0]*I0[10] + U[4][0]*I1[10]
           + (oo2z)*(I2[10] - (poz)*I3[10]);
V_Cl[in].Val[11] = U[0][0]*I0[11] + U[4][0]*I1[11]
           + (oo2z)*(I2[11] - (poz)*I3[11]);
V_Cl[in].Val[12] = U[0][0]*I0[12] + U[4][0]*I1[12]
           + (oo2z)*(I2[12] - (poz)*I3[12])
           + (oo2zn)*I4[4];
V_Cl[in].Val[13] = U[0][0]*I0[13] + U[4][0]*I1[13]
           + (oo2z)*(I2[13] - (poz)*I3[13]);
V_Cl[in].Val[14] = U[0][0]*I0[14] + U[4][0]*I1[14]
           + (oo2z)*(I2[14] - (poz)*I3[14]);
V_Cl[in].Val[15] = U[0][0]*I0[15] + U[4][0]*I1[15]
           + (oo2z)*(I2[15] - (poz)*I3[15])
           + (oo2zn)*I4[5];
V_Cl[in].Val[16] = U[0][0]*I0[16] + U[4][0]*I1[16]
           + (oo2z)*(I2[16] - (poz)*I3[16]);
V_Cl[in].Val[17] = U[0][0]*I0[17] + U[4][0]*I1[17]
           + (oo2z)*(I2[17] - (poz)*I3[17]);
V_Cl[in].Val[18] = U[0][0]*I0[18] + U[4][0]*I1[18]
           + (oo2zn)*I4[6];
V_Cl[in].Val[19] = U[0][0]*I0[19] + U[4][0]*I1[19];
V_Cl[in].Val[20] = U[0][0]*I0[20] + U[4][0]*I1[20];
V_Cl[in].Val[21] = U[0][0]*I0[21] + U[4][0]*I1[21]
           + (oo2zn)*I4[7];
V_Cl[in].Val[22] = U[0][0]*I0[22] + U[4][0]*I1[22];
V_Cl[in].Val[23] = U[0][0]*I0[23] + U[4][0]*I1[23];
V_Cl[in].Val[24] = U[0][0]*I0[24] + U[4][0]*I1[24]
           + (oo2zn)*I4[8];
V_Cl[in].Val[25] = U[0][0]*I0[25] + U[4][0]*I1[25];
V_Cl[in].Val[26] = U[0][0]*I0[26] + U[4][0]*I1[26];
V_Cl[in].Val[27] = U[0][0]*I0[27] + U[4][0]*I1[27]
           + (oo2zn)*I4[9];
V_Cl[in].Val[28] = U[0][0]*I0[28] + U[4][0]*I1[28];
V_Cl[in].Val[29] = U[0][0]*I0[29] + U[4][0]*I1[29];
V_Cl[in].Val[30] = U[0][1]*I0[18] + U[4][1]*I1[18]
           + (threeoo2z)*(I2[9] - (poz)*I3[9]);
V_Cl[in].Val[31] = U[0][1]*I0[19] + U[4][1]*I1[19]
           + (threeoo2z)*(I2[10] - (poz)*I3[10])
           + (oo2zn)*I4[6];
V_Cl[in].Val[32] = U[0][1]*I0[20] + U[4][1]*I1[20]
           + (threeoo2z)*(I2[11] - (poz)*I3[11]);
V_Cl[in].Val[33] = U[0][1]*I0[21] + U[4][1]*I1[21]
           + (twooo2z)*(I2[12] - (poz)*I3[12]);
V_Cl[in].Val[34] = U[0][1]*I0[22] + U[4][1]*I1[22]
           + (twooo2z)*(I2[13] - (poz)*I3[13])
           + (oo2zn)*I4[7];
V_Cl[in].Val[35] = U[0][1]*I0[23] + U[4][1]*I1[23]
           + (twooo2z)*(I2[14] - (poz)*I3[14]);
V_Cl[in].Val[36] = U[0][1]*I0[24] + U[4][1]*I1[24]
           + (oo2z)*(I2[15] - (poz)*I3[15]);
V_Cl[in].Val[37] = U[0][1]*I0[25] + U[4][1]*I1[25]
           + (oo2z)*(I2[16] - (poz)*I3[16])
           + (oo2zn)*I4[8];
V_Cl[in].Val[38] = U[0][1]*I0[26] + U[4][1]*I1[26]
           + (oo2z)*(I2[17] - (poz)*I3[17]);
V_Cl[in].Val[39] = U[0][1]*I0[27] + U[4][1]*I1[27];
V_Cl[in].Val[40] = U[0][1]*I0[28] + U[4][1]*I1[28]
           + (oo2zn)*I4[9];
V_Cl[in].Val[41] = U[0][1]*I0[29] + U[4][1]*I1[29];
V_Cl[in].Val[42] = U[0][2]*I0[27] + U[4][2]*I1[27]
           + (threeoo2z)*(I2[15] - (poz)*I3[15]);
V_Cl[in].Val[43] = U[0][2]*I0[28] + U[4][2]*I1[28]
           + (threeoo2z)*(I2[16] - (poz)*I3[16]);
V_Cl[in].Val[44] = U[0][2]*I0[29] + U[4][2]*I1[29]
           + (threeoo2z)*(I2[17] - (poz)*I3[17])
           + (oo2zn)*I4[9];
  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_p0g0(iclass *V_Cl, int in) /* type = 20 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2zn = 2.0*oo2zn;
  double threeoo2zn = 3.0*oo2zn;
  double fouroo2zn = 4.0*oo2zn;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_00g0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_00g0(V_Cl, V_Cl[in].operands[1]);
  I4 = build_00f0(V_Cl, V_Cl[in].operands[4]);


V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0]
           + (fouroo2zn)*I4[0];
V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1]
           + (threeoo2zn)*I4[1];
V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2]
           + (threeoo2zn)*I4[2];
V_Cl[in].Val[3] = U[0][0]*I0[3] + U[4][0]*I1[3]
           + (twooo2zn)*I4[3];
V_Cl[in].Val[4] = U[0][0]*I0[4] + U[4][0]*I1[4]
           + (twooo2zn)*I4[4];
V_Cl[in].Val[5] = U[0][0]*I0[5] + U[4][0]*I1[5]
           + (twooo2zn)*I4[5];
V_Cl[in].Val[6] = U[0][0]*I0[6] + U[4][0]*I1[6]
           + (oo2zn)*I4[6];
V_Cl[in].Val[7] = U[0][0]*I0[7] + U[4][0]*I1[7]
           + (oo2zn)*I4[7];
V_Cl[in].Val[8] = U[0][0]*I0[8] + U[4][0]*I1[8]
           + (oo2zn)*I4[8];
V_Cl[in].Val[9] = U[0][0]*I0[9] + U[4][0]*I1[9]
           + (oo2zn)*I4[9];
V_Cl[in].Val[10] = U[0][0]*I0[10] + U[4][0]*I1[10];
V_Cl[in].Val[11] = U[0][0]*I0[11] + U[4][0]*I1[11];
V_Cl[in].Val[12] = U[0][0]*I0[12] + U[4][0]*I1[12];
V_Cl[in].Val[13] = U[0][0]*I0[13] + U[4][0]*I1[13];
V_Cl[in].Val[14] = U[0][0]*I0[14] + U[4][0]*I1[14];
V_Cl[in].Val[15] = U[0][1]*I0[0] + U[4][1]*I1[0];
V_Cl[in].Val[16] = U[0][1]*I0[1] + U[4][1]*I1[1]
           + (oo2zn)*I4[0];
V_Cl[in].Val[17] = U[0][1]*I0[2] + U[4][1]*I1[2];
V_Cl[in].Val[18] = U[0][1]*I0[3] + U[4][1]*I1[3]
           + (twooo2zn)*I4[1];
V_Cl[in].Val[19] = U[0][1]*I0[4] + U[4][1]*I1[4]
           + (oo2zn)*I4[2];
V_Cl[in].Val[20] = U[0][1]*I0[5] + U[4][1]*I1[5];
V_Cl[in].Val[21] = U[0][1]*I0[6] + U[4][1]*I1[6]
           + (threeoo2zn)*I4[3];
V_Cl[in].Val[22] = U[0][1]*I0[7] + U[4][1]*I1[7]
           + (twooo2zn)*I4[4];
V_Cl[in].Val[23] = U[0][1]*I0[8] + U[4][1]*I1[8]
           + (oo2zn)*I4[5];
V_Cl[in].Val[24] = U[0][1]*I0[9] + U[4][1]*I1[9];
V_Cl[in].Val[25] = U[0][1]*I0[10] + U[4][1]*I1[10]
           + (fouroo2zn)*I4[6];
V_Cl[in].Val[26] = U[0][1]*I0[11] + U[4][1]*I1[11]
           + (threeoo2zn)*I4[7];
V_Cl[in].Val[27] = U[0][1]*I0[12] + U[4][1]*I1[12]
           + (twooo2zn)*I4[8];
V_Cl[in].Val[28] = U[0][1]*I0[13] + U[4][1]*I1[13]
           + (oo2zn)*I4[9];
V_Cl[in].Val[29] = U[0][1]*I0[14] + U[4][1]*I1[14];
V_Cl[in].Val[30] = U[0][2]*I0[0] + U[4][2]*I1[0];
V_Cl[in].Val[31] = U[0][2]*I0[1] + U[4][2]*I1[1];
V_Cl[in].Val[32] = U[0][2]*I0[2] + U[4][2]*I1[2]
           + (oo2zn)*I4[0];
V_Cl[in].Val[33] = U[0][2]*I0[3] + U[4][2]*I1[3];
V_Cl[in].Val[34] = U[0][2]*I0[4] + U[4][2]*I1[4]
           + (oo2zn)*I4[1];
V_Cl[in].Val[35] = U[0][2]*I0[5] + U[4][2]*I1[5]
           + (twooo2zn)*I4[2];
V_Cl[in].Val[36] = U[0][2]*I0[6] + U[4][2]*I1[6];
V_Cl[in].Val[37] = U[0][2]*I0[7] + U[4][2]*I1[7]
           + (oo2zn)*I4[3];
V_Cl[in].Val[38] = U[0][2]*I0[8] + U[4][2]*I1[8]
           + (twooo2zn)*I4[4];
V_Cl[in].Val[39] = U[0][2]*I0[9] + U[4][2]*I1[9]
           + (threeoo2zn)*I4[5];
V_Cl[in].Val[40] = U[0][2]*I0[10] + U[4][2]*I1[10];
V_Cl[in].Val[41] = U[0][2]*I0[11] + U[4][2]*I1[11]
           + (oo2zn)*I4[6];
V_Cl[in].Val[42] = U[0][2]*I0[12] + U[4][2]*I1[12]
           + (twooo2zn)*I4[7];
V_Cl[in].Val[43] = U[0][2]*I0[13] + U[4][2]*I1[13]
           + (threeoo2zn)*I4[8];
V_Cl[in].Val[44] = U[0][2]*I0[14] + U[4][2]*I1[14]
           + (fouroo2zn)*I4[9];
  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_d0g0(iclass *V_Cl, int in) /* type = 21 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2zn = 2.0*oo2zn;
  double threeoo2zn = 3.0*oo2zn;
  double fouroo2zn = 4.0*oo2zn;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_p0g0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_p0g0(V_Cl, V_Cl[in].operands[1]);
  I2 = build_00g0(V_Cl, V_Cl[in].operands[2]);
  I3 = build_00g0(V_Cl, V_Cl[in].operands[3]);
  I4 = build_p0f0(V_Cl, V_Cl[in].operands[4]);


V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0]
           + (oo2z)*(I2[0] - (poz)*I3[0])
           + (fouroo2zn)*I4[0];
V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1]
           + (oo2z)*(I2[1] - (poz)*I3[1])
           + (threeoo2zn)*I4[1];
V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2]
           + (oo2z)*(I2[2] - (poz)*I3[2])
           + (threeoo2zn)*I4[2];
V_Cl[in].Val[3] = U[0][0]*I0[3] + U[4][0]*I1[3]
           + (oo2z)*(I2[3] - (poz)*I3[3])
           + (twooo2zn)*I4[3];
V_Cl[in].Val[4] = U[0][0]*I0[4] + U[4][0]*I1[4]
           + (oo2z)*(I2[4] - (poz)*I3[4])
           + (twooo2zn)*I4[4];
V_Cl[in].Val[5] = U[0][0]*I0[5] + U[4][0]*I1[5]
           + (oo2z)*(I2[5] - (poz)*I3[5])
           + (twooo2zn)*I4[5];
V_Cl[in].Val[6] = U[0][0]*I0[6] + U[4][0]*I1[6]
           + (oo2z)*(I2[6] - (poz)*I3[6])
           + (oo2zn)*I4[6];
V_Cl[in].Val[7] = U[0][0]*I0[7] + U[4][0]*I1[7]
           + (oo2z)*(I2[7] - (poz)*I3[7])
           + (oo2zn)*I4[7];
V_Cl[in].Val[8] = U[0][0]*I0[8] + U[4][0]*I1[8]
           + (oo2z)*(I2[8] - (poz)*I3[8])
           + (oo2zn)*I4[8];
V_Cl[in].Val[9] = U[0][0]*I0[9] + U[4][0]*I1[9]
           + (oo2z)*(I2[9] - (poz)*I3[9])
           + (oo2zn)*I4[9];
V_Cl[in].Val[10] = U[0][0]*I0[10] + U[4][0]*I1[10]
           + (oo2z)*(I2[10] - (poz)*I3[10]);
V_Cl[in].Val[11] = U[0][0]*I0[11] + U[4][0]*I1[11]
           + (oo2z)*(I2[11] - (poz)*I3[11]);
V_Cl[in].Val[12] = U[0][0]*I0[12] + U[4][0]*I1[12]
           + (oo2z)*(I2[12] - (poz)*I3[12]);
V_Cl[in].Val[13] = U[0][0]*I0[13] + U[4][0]*I1[13]
           + (oo2z)*(I2[13] - (poz)*I3[13]);
V_Cl[in].Val[14] = U[0][0]*I0[14] + U[4][0]*I1[14]
           + (oo2z)*(I2[14] - (poz)*I3[14]);
V_Cl[in].Val[15] = U[0][0]*I0[15] + U[4][0]*I1[15]
           + (fouroo2zn)*I4[10];
V_Cl[in].Val[16] = U[0][0]*I0[16] + U[4][0]*I1[16]
           + (threeoo2zn)*I4[11];
V_Cl[in].Val[17] = U[0][0]*I0[17] + U[4][0]*I1[17]
           + (threeoo2zn)*I4[12];
V_Cl[in].Val[18] = U[0][0]*I0[18] + U[4][0]*I1[18]
           + (twooo2zn)*I4[13];
V_Cl[in].Val[19] = U[0][0]*I0[19] + U[4][0]*I1[19]
           + (twooo2zn)*I4[14];
V_Cl[in].Val[20] = U[0][0]*I0[20] + U[4][0]*I1[20]
           + (twooo2zn)*I4[15];
V_Cl[in].Val[21] = U[0][0]*I0[21] + U[4][0]*I1[21]
           + (oo2zn)*I4[16];
V_Cl[in].Val[22] = U[0][0]*I0[22] + U[4][0]*I1[22]
           + (oo2zn)*I4[17];
V_Cl[in].Val[23] = U[0][0]*I0[23] + U[4][0]*I1[23]
           + (oo2zn)*I4[18];
V_Cl[in].Val[24] = U[0][0]*I0[24] + U[4][0]*I1[24]
           + (oo2zn)*I4[19];
V_Cl[in].Val[25] = U[0][0]*I0[25] + U[4][0]*I1[25];
V_Cl[in].Val[26] = U[0][0]*I0[26] + U[4][0]*I1[26];
V_Cl[in].Val[27] = U[0][0]*I0[27] + U[4][0]*I1[27];
V_Cl[in].Val[28] = U[0][0]*I0[28] + U[4][0]*I1[28];
V_Cl[in].Val[29] = U[0][0]*I0[29] + U[4][0]*I1[29];
V_Cl[in].Val[30] = U[0][0]*I0[30] + U[4][0]*I1[30]
           + (fouroo2zn)*I4[20];
V_Cl[in].Val[31] = U[0][0]*I0[31] + U[4][0]*I1[31]
           + (threeoo2zn)*I4[21];
V_Cl[in].Val[32] = U[0][0]*I0[32] + U[4][0]*I1[32]
           + (threeoo2zn)*I4[22];
V_Cl[in].Val[33] = U[0][0]*I0[33] + U[4][0]*I1[33]
           + (twooo2zn)*I4[23];
V_Cl[in].Val[34] = U[0][0]*I0[34] + U[4][0]*I1[34]
           + (twooo2zn)*I4[24];
V_Cl[in].Val[35] = U[0][0]*I0[35] + U[4][0]*I1[35]
           + (twooo2zn)*I4[25];
V_Cl[in].Val[36] = U[0][0]*I0[36] + U[4][0]*I1[36]
           + (oo2zn)*I4[26];
V_Cl[in].Val[37] = U[0][0]*I0[37] + U[4][0]*I1[37]
           + (oo2zn)*I4[27];
V_Cl[in].Val[38] = U[0][0]*I0[38] + U[4][0]*I1[38]
           + (oo2zn)*I4[28];
V_Cl[in].Val[39] = U[0][0]*I0[39] + U[4][0]*I1[39]
           + (oo2zn)*I4[29];
V_Cl[in].Val[40] = U[0][0]*I0[40] + U[4][0]*I1[40];
V_Cl[in].Val[41] = U[0][0]*I0[41] + U[4][0]*I1[41];
V_Cl[in].Val[42] = U[0][0]*I0[42] + U[4][0]*I1[42];
V_Cl[in].Val[43] = U[0][0]*I0[43] + U[4][0]*I1[43];
V_Cl[in].Val[44] = U[0][0]*I0[44] + U[4][0]*I1[44];
V_Cl[in].Val[45] = U[0][1]*I0[15] + U[4][1]*I1[15]
           + (oo2z)*(I2[0] - (poz)*I3[0]);
V_Cl[in].Val[46] = U[0][1]*I0[16] + U[4][1]*I1[16]
           + (oo2z)*(I2[1] - (poz)*I3[1])
           + (oo2zn)*I4[10];
V_Cl[in].Val[47] = U[0][1]*I0[17] + U[4][1]*I1[17]
           + (oo2z)*(I2[2] - (poz)*I3[2]);
V_Cl[in].Val[48] = U[0][1]*I0[18] + U[4][1]*I1[18]
           + (oo2z)*(I2[3] - (poz)*I3[3])
           + (twooo2zn)*I4[11];
V_Cl[in].Val[49] = U[0][1]*I0[19] + U[4][1]*I1[19]
           + (oo2z)*(I2[4] - (poz)*I3[4])
           + (oo2zn)*I4[12];
V_Cl[in].Val[50] = U[0][1]*I0[20] + U[4][1]*I1[20]
           + (oo2z)*(I2[5] - (poz)*I3[5]);
V_Cl[in].Val[51] = U[0][1]*I0[21] + U[4][1]*I1[21]
           + (oo2z)*(I2[6] - (poz)*I3[6])
           + (threeoo2zn)*I4[13];
V_Cl[in].Val[52] = U[0][1]*I0[22] + U[4][1]*I1[22]
           + (oo2z)*(I2[7] - (poz)*I3[7])
           + (twooo2zn)*I4[14];
V_Cl[in].Val[53] = U[0][1]*I0[23] + U[4][1]*I1[23]
           + (oo2z)*(I2[8] - (poz)*I3[8])
           + (oo2zn)*I4[15];
V_Cl[in].Val[54] = U[0][1]*I0[24] + U[4][1]*I1[24]
           + (oo2z)*(I2[9] - (poz)*I3[9]);
V_Cl[in].Val[55] = U[0][1]*I0[25] + U[4][1]*I1[25]
           + (oo2z)*(I2[10] - (poz)*I3[10])
           + (fouroo2zn)*I4[16];
V_Cl[in].Val[56] = U[0][1]*I0[26] + U[4][1]*I1[26]
           + (oo2z)*(I2[11] - (poz)*I3[11])
           + (threeoo2zn)*I4[17];
V_Cl[in].Val[57] = U[0][1]*I0[27] + U[4][1]*I1[27]
           + (oo2z)*(I2[12] - (poz)*I3[12])
           + (twooo2zn)*I4[18];
V_Cl[in].Val[58] = U[0][1]*I0[28] + U[4][1]*I1[28]
           + (oo2z)*(I2[13] - (poz)*I3[13])
           + (oo2zn)*I4[19];
V_Cl[in].Val[59] = U[0][1]*I0[29] + U[4][1]*I1[29]
           + (oo2z)*(I2[14] - (poz)*I3[14]);
V_Cl[in].Val[60] = U[0][1]*I0[30] + U[4][1]*I1[30];
V_Cl[in].Val[61] = U[0][1]*I0[31] + U[4][1]*I1[31]
           + (oo2zn)*I4[20];
V_Cl[in].Val[62] = U[0][1]*I0[32] + U[4][1]*I1[32];
V_Cl[in].Val[63] = U[0][1]*I0[33] + U[4][1]*I1[33]
           + (twooo2zn)*I4[21];
V_Cl[in].Val[64] = U[0][1]*I0[34] + U[4][1]*I1[34]
           + (oo2zn)*I4[22];
V_Cl[in].Val[65] = U[0][1]*I0[35] + U[4][1]*I1[35];
V_Cl[in].Val[66] = U[0][1]*I0[36] + U[4][1]*I1[36]
           + (threeoo2zn)*I4[23];
V_Cl[in].Val[67] = U[0][1]*I0[37] + U[4][1]*I1[37]
           + (twooo2zn)*I4[24];
V_Cl[in].Val[68] = U[0][1]*I0[38] + U[4][1]*I1[38]
           + (oo2zn)*I4[25];
V_Cl[in].Val[69] = U[0][1]*I0[39] + U[4][1]*I1[39];
V_Cl[in].Val[70] = U[0][1]*I0[40] + U[4][1]*I1[40]
           + (fouroo2zn)*I4[26];
V_Cl[in].Val[71] = U[0][1]*I0[41] + U[4][1]*I1[41]
           + (threeoo2zn)*I4[27];
V_Cl[in].Val[72] = U[0][1]*I0[42] + U[4][1]*I1[42]
           + (twooo2zn)*I4[28];
V_Cl[in].Val[73] = U[0][1]*I0[43] + U[4][1]*I1[43]
           + (oo2zn)*I4[29];
V_Cl[in].Val[74] = U[0][1]*I0[44] + U[4][1]*I1[44];
V_Cl[in].Val[75] = U[0][2]*I0[30] + U[4][2]*I1[30]
           + (oo2z)*(I2[0] - (poz)*I3[0]);
V_Cl[in].Val[76] = U[0][2]*I0[31] + U[4][2]*I1[31]
           + (oo2z)*(I2[1] - (poz)*I3[1]);
V_Cl[in].Val[77] = U[0][2]*I0[32] + U[4][2]*I1[32]
           + (oo2z)*(I2[2] - (poz)*I3[2])
           + (oo2zn)*I4[20];
V_Cl[in].Val[78] = U[0][2]*I0[33] + U[4][2]*I1[33]
           + (oo2z)*(I2[3] - (poz)*I3[3]);
V_Cl[in].Val[79] = U[0][2]*I0[34] + U[4][2]*I1[34]
           + (oo2z)*(I2[4] - (poz)*I3[4])
           + (oo2zn)*I4[21];
V_Cl[in].Val[80] = U[0][2]*I0[35] + U[4][2]*I1[35]
           + (oo2z)*(I2[5] - (poz)*I3[5])
           + (twooo2zn)*I4[22];
V_Cl[in].Val[81] = U[0][2]*I0[36] + U[4][2]*I1[36]
           + (oo2z)*(I2[6] - (poz)*I3[6]);
V_Cl[in].Val[82] = U[0][2]*I0[37] + U[4][2]*I1[37]
           + (oo2z)*(I2[7] - (poz)*I3[7])
           + (oo2zn)*I4[23];
V_Cl[in].Val[83] = U[0][2]*I0[38] + U[4][2]*I1[38]
           + (oo2z)*(I2[8] - (poz)*I3[8])
           + (twooo2zn)*I4[24];
V_Cl[in].Val[84] = U[0][2]*I0[39] + U[4][2]*I1[39]
           + (oo2z)*(I2[9] - (poz)*I3[9])
           + (threeoo2zn)*I4[25];
V_Cl[in].Val[85] = U[0][2]*I0[40] + U[4][2]*I1[40]
           + (oo2z)*(I2[10] - (poz)*I3[10]);
V_Cl[in].Val[86] = U[0][2]*I0[41] + U[4][2]*I1[41]
           + (oo2z)*(I2[11] - (poz)*I3[11])
           + (oo2zn)*I4[26];
V_Cl[in].Val[87] = U[0][2]*I0[42] + U[4][2]*I1[42]
           + (oo2z)*(I2[12] - (poz)*I3[12])
           + (twooo2zn)*I4[27];
V_Cl[in].Val[88] = U[0][2]*I0[43] + U[4][2]*I1[43]
           + (oo2z)*(I2[13] - (poz)*I3[13])
           + (threeoo2zn)*I4[28];
V_Cl[in].Val[89] = U[0][2]*I0[44] + U[4][2]*I1[44]
           + (oo2z)*(I2[14] - (poz)*I3[14])
           + (fouroo2zn)*I4[29];
  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_g0d0(iclass *V_Cl, int in) /* type = 22 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2z = 2.0*oo2z;
  double threeoo2z = 3.0*oo2z;
  double twooo2zn = 2.0*oo2zn;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_f0d0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_f0d0(V_Cl, V_Cl[in].operands[1]);
  I2 = build_d0d0(V_Cl, V_Cl[in].operands[2]);
  I3 = build_d0d0(V_Cl, V_Cl[in].operands[3]);
  I4 = build_f0p0(V_Cl, V_Cl[in].operands[4]);


V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0]
           + (threeoo2z)*(I2[0] - (poz)*I3[0])
           + (twooo2zn)*I4[0];
V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1]
           + (threeoo2z)*(I2[1] - (poz)*I3[1])
           + (oo2zn)*I4[1];
V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2]
           + (threeoo2z)*(I2[2] - (poz)*I3[2])
           + (oo2zn)*I4[2];
V_Cl[in].Val[3] = U[0][0]*I0[3] + U[4][0]*I1[3]
           + (threeoo2z)*(I2[3] - (poz)*I3[3]);
V_Cl[in].Val[4] = U[0][0]*I0[4] + U[4][0]*I1[4]
           + (threeoo2z)*(I2[4] - (poz)*I3[4]);
V_Cl[in].Val[5] = U[0][0]*I0[5] + U[4][0]*I1[5]
           + (threeoo2z)*(I2[5] - (poz)*I3[5]);
V_Cl[in].Val[6] = U[0][0]*I0[6] + U[4][0]*I1[6]
           + (twooo2z)*(I2[6] - (poz)*I3[6])
           + (twooo2zn)*I4[3];
V_Cl[in].Val[7] = U[0][0]*I0[7] + U[4][0]*I1[7]
           + (twooo2z)*(I2[7] - (poz)*I3[7])
           + (oo2zn)*I4[4];
V_Cl[in].Val[8] = U[0][0]*I0[8] + U[4][0]*I1[8]
           + (twooo2z)*(I2[8] - (poz)*I3[8])
           + (oo2zn)*I4[5];
V_Cl[in].Val[9] = U[0][0]*I0[9] + U[4][0]*I1[9]
           + (twooo2z)*(I2[9] - (poz)*I3[9]);
V_Cl[in].Val[10] = U[0][0]*I0[10] + U[4][0]*I1[10]
           + (twooo2z)*(I2[10] - (poz)*I3[10]);
V_Cl[in].Val[11] = U[0][0]*I0[11] + U[4][0]*I1[11]
           + (twooo2z)*(I2[11] - (poz)*I3[11]);
V_Cl[in].Val[12] = U[0][0]*I0[12] + U[4][0]*I1[12]
           + (twooo2z)*(I2[12] - (poz)*I3[12])
           + (twooo2zn)*I4[6];
V_Cl[in].Val[13] = U[0][0]*I0[13] + U[4][0]*I1[13]
           + (twooo2z)*(I2[13] - (poz)*I3[13])
           + (oo2zn)*I4[7];
V_Cl[in].Val[14] = U[0][0]*I0[14] + U[4][0]*I1[14]
           + (twooo2z)*(I2[14] - (poz)*I3[14])
           + (oo2zn)*I4[8];
V_Cl[in].Val[15] = U[0][0]*I0[15] + U[4][0]*I1[15]
           + (twooo2z)*(I2[15] - (poz)*I3[15]);
V_Cl[in].Val[16] = U[0][0]*I0[16] + U[4][0]*I1[16]
           + (twooo2z)*(I2[16] - (poz)*I3[16]);
V_Cl[in].Val[17] = U[0][0]*I0[17] + U[4][0]*I1[17]
           + (twooo2z)*(I2[17] - (poz)*I3[17]);
V_Cl[in].Val[18] = U[0][0]*I0[18] + U[4][0]*I1[18]
           + (oo2z)*(I2[18] - (poz)*I3[18])
           + (twooo2zn)*I4[9];
V_Cl[in].Val[19] = U[0][0]*I0[19] + U[4][0]*I1[19]
           + (oo2z)*(I2[19] - (poz)*I3[19])
           + (oo2zn)*I4[10];
V_Cl[in].Val[20] = U[0][0]*I0[20] + U[4][0]*I1[20]
           + (oo2z)*(I2[20] - (poz)*I3[20])
           + (oo2zn)*I4[11];
V_Cl[in].Val[21] = U[0][0]*I0[21] + U[4][0]*I1[21]
           + (oo2z)*(I2[21] - (poz)*I3[21]);
V_Cl[in].Val[22] = U[0][0]*I0[22] + U[4][0]*I1[22]
           + (oo2z)*(I2[22] - (poz)*I3[22]);
V_Cl[in].Val[23] = U[0][0]*I0[23] + U[4][0]*I1[23]
           + (oo2z)*(I2[23] - (poz)*I3[23]);
V_Cl[in].Val[24] = U[0][0]*I0[24] + U[4][0]*I1[24]
           + (oo2z)*(I2[24] - (poz)*I3[24])
           + (twooo2zn)*I4[12];
V_Cl[in].Val[25] = U[0][0]*I0[25] + U[4][0]*I1[25]
           + (oo2z)*(I2[25] - (poz)*I3[25])
           + (oo2zn)*I4[13];
V_Cl[in].Val[26] = U[0][0]*I0[26] + U[4][0]*I1[26]
           + (oo2z)*(I2[26] - (poz)*I3[26])
           + (oo2zn)*I4[14];
V_Cl[in].Val[27] = U[0][0]*I0[27] + U[4][0]*I1[27]
           + (oo2z)*(I2[27] - (poz)*I3[27]);
V_Cl[in].Val[28] = U[0][0]*I0[28] + U[4][0]*I1[28]
           + (oo2z)*(I2[28] - (poz)*I3[28]);
V_Cl[in].Val[29] = U[0][0]*I0[29] + U[4][0]*I1[29]
           + (oo2z)*(I2[29] - (poz)*I3[29]);
V_Cl[in].Val[30] = U[0][0]*I0[30] + U[4][0]*I1[30]
           + (oo2z)*(I2[30] - (poz)*I3[30])
           + (twooo2zn)*I4[15];
V_Cl[in].Val[31] = U[0][0]*I0[31] + U[4][0]*I1[31]
           + (oo2z)*(I2[31] - (poz)*I3[31])
           + (oo2zn)*I4[16];
V_Cl[in].Val[32] = U[0][0]*I0[32] + U[4][0]*I1[32]
           + (oo2z)*(I2[32] - (poz)*I3[32])
           + (oo2zn)*I4[17];
V_Cl[in].Val[33] = U[0][0]*I0[33] + U[4][0]*I1[33]
           + (oo2z)*(I2[33] - (poz)*I3[33]);
V_Cl[in].Val[34] = U[0][0]*I0[34] + U[4][0]*I1[34]
           + (oo2z)*(I2[34] - (poz)*I3[34]);
V_Cl[in].Val[35] = U[0][0]*I0[35] + U[4][0]*I1[35]
           + (oo2z)*(I2[35] - (poz)*I3[35]);
V_Cl[in].Val[36] = U[0][0]*I0[36] + U[4][0]*I1[36]
           + (twooo2zn)*I4[18];
V_Cl[in].Val[37] = U[0][0]*I0[37] + U[4][0]*I1[37]
           + (oo2zn)*I4[19];
V_Cl[in].Val[38] = U[0][0]*I0[38] + U[4][0]*I1[38]
           + (oo2zn)*I4[20];
V_Cl[in].Val[39] = U[0][0]*I0[39] + U[4][0]*I1[39];
V_Cl[in].Val[40] = U[0][0]*I0[40] + U[4][0]*I1[40];
V_Cl[in].Val[41] = U[0][0]*I0[41] + U[4][0]*I1[41];
V_Cl[in].Val[42] = U[0][0]*I0[42] + U[4][0]*I1[42]
           + (twooo2zn)*I4[21];
V_Cl[in].Val[43] = U[0][0]*I0[43] + U[4][0]*I1[43]
           + (oo2zn)*I4[22];
V_Cl[in].Val[44] = U[0][0]*I0[44] + U[4][0]*I1[44]
           + (oo2zn)*I4[23];
V_Cl[in].Val[45] = U[0][0]*I0[45] + U[4][0]*I1[45];
V_Cl[in].Val[46] = U[0][0]*I0[46] + U[4][0]*I1[46];
V_Cl[in].Val[47] = U[0][0]*I0[47] + U[4][0]*I1[47];
V_Cl[in].Val[48] = U[0][0]*I0[48] + U[4][0]*I1[48]
           + (twooo2zn)*I4[24];
V_Cl[in].Val[49] = U[0][0]*I0[49] + U[4][0]*I1[49]
           + (oo2zn)*I4[25];
V_Cl[in].Val[50] = U[0][0]*I0[50] + U[4][0]*I1[50]
           + (oo2zn)*I4[26];
V_Cl[in].Val[51] = U[0][0]*I0[51] + U[4][0]*I1[51];
V_Cl[in].Val[52] = U[0][0]*I0[52] + U[4][0]*I1[52];
V_Cl[in].Val[53] = U[0][0]*I0[53] + U[4][0]*I1[53];
V_Cl[in].Val[54] = U[0][0]*I0[54] + U[4][0]*I1[54]
           + (twooo2zn)*I4[27];
V_Cl[in].Val[55] = U[0][0]*I0[55] + U[4][0]*I1[55]
           + (oo2zn)*I4[28];
V_Cl[in].Val[56] = U[0][0]*I0[56] + U[4][0]*I1[56]
           + (oo2zn)*I4[29];
V_Cl[in].Val[57] = U[0][0]*I0[57] + U[4][0]*I1[57];
V_Cl[in].Val[58] = U[0][0]*I0[58] + U[4][0]*I1[58];
V_Cl[in].Val[59] = U[0][0]*I0[59] + U[4][0]*I1[59];
V_Cl[in].Val[60] = U[0][1]*I0[36] + U[4][1]*I1[36]
           + (threeoo2z)*(I2[18] - (poz)*I3[18]);
V_Cl[in].Val[61] = U[0][1]*I0[37] + U[4][1]*I1[37]
           + (threeoo2z)*(I2[19] - (poz)*I3[19])
           + (oo2zn)*I4[18];
V_Cl[in].Val[62] = U[0][1]*I0[38] + U[4][1]*I1[38]
           + (threeoo2z)*(I2[20] - (poz)*I3[20]);
V_Cl[in].Val[63] = U[0][1]*I0[39] + U[4][1]*I1[39]
           + (threeoo2z)*(I2[21] - (poz)*I3[21])
           + (twooo2zn)*I4[19];
V_Cl[in].Val[64] = U[0][1]*I0[40] + U[4][1]*I1[40]
           + (threeoo2z)*(I2[22] - (poz)*I3[22])
           + (oo2zn)*I4[20];
V_Cl[in].Val[65] = U[0][1]*I0[41] + U[4][1]*I1[41]
           + (threeoo2z)*(I2[23] - (poz)*I3[23]);
V_Cl[in].Val[66] = U[0][1]*I0[42] + U[4][1]*I1[42]
           + (twooo2z)*(I2[24] - (poz)*I3[24]);
V_Cl[in].Val[67] = U[0][1]*I0[43] + U[4][1]*I1[43]
           + (twooo2z)*(I2[25] - (poz)*I3[25])
           + (oo2zn)*I4[21];
V_Cl[in].Val[68] = U[0][1]*I0[44] + U[4][1]*I1[44]
           + (twooo2z)*(I2[26] - (poz)*I3[26]);
V_Cl[in].Val[69] = U[0][1]*I0[45] + U[4][1]*I1[45]
           + (twooo2z)*(I2[27] - (poz)*I3[27])
           + (twooo2zn)*I4[22];
V_Cl[in].Val[70] = U[0][1]*I0[46] + U[4][1]*I1[46]
           + (twooo2z)*(I2[28] - (poz)*I3[28])
           + (oo2zn)*I4[23];
V_Cl[in].Val[71] = U[0][1]*I0[47] + U[4][1]*I1[47]
           + (twooo2z)*(I2[29] - (poz)*I3[29]);
V_Cl[in].Val[72] = U[0][1]*I0[48] + U[4][1]*I1[48]
           + (oo2z)*(I2[30] - (poz)*I3[30]);
V_Cl[in].Val[73] = U[0][1]*I0[49] + U[4][1]*I1[49]
           + (oo2z)*(I2[31] - (poz)*I3[31])
           + (oo2zn)*I4[24];
V_Cl[in].Val[74] = U[0][1]*I0[50] + U[4][1]*I1[50]
           + (oo2z)*(I2[32] - (poz)*I3[32]);
V_Cl[in].Val[75] = U[0][1]*I0[51] + U[4][1]*I1[51]
           + (oo2z)*(I2[33] - (poz)*I3[33])
           + (twooo2zn)*I4[25];
V_Cl[in].Val[76] = U[0][1]*I0[52] + U[4][1]*I1[52]
           + (oo2z)*(I2[34] - (poz)*I3[34])
           + (oo2zn)*I4[26];
V_Cl[in].Val[77] = U[0][1]*I0[53] + U[4][1]*I1[53]
           + (oo2z)*(I2[35] - (poz)*I3[35]);
V_Cl[in].Val[78] = U[0][1]*I0[54] + U[4][1]*I1[54];
V_Cl[in].Val[79] = U[0][1]*I0[55] + U[4][1]*I1[55]
           + (oo2zn)*I4[27];
V_Cl[in].Val[80] = U[0][1]*I0[56] + U[4][1]*I1[56];
V_Cl[in].Val[81] = U[0][1]*I0[57] + U[4][1]*I1[57]
           + (twooo2zn)*I4[28];
V_Cl[in].Val[82] = U[0][1]*I0[58] + U[4][1]*I1[58]
           + (oo2zn)*I4[29];
V_Cl[in].Val[83] = U[0][1]*I0[59] + U[4][1]*I1[59];
V_Cl[in].Val[84] = U[0][2]*I0[54] + U[4][2]*I1[54]
           + (threeoo2z)*(I2[30] - (poz)*I3[30]);
V_Cl[in].Val[85] = U[0][2]*I0[55] + U[4][2]*I1[55]
           + (threeoo2z)*(I2[31] - (poz)*I3[31]);
V_Cl[in].Val[86] = U[0][2]*I0[56] + U[4][2]*I1[56]
           + (threeoo2z)*(I2[32] - (poz)*I3[32])
           + (oo2zn)*I4[27];
V_Cl[in].Val[87] = U[0][2]*I0[57] + U[4][2]*I1[57]
           + (threeoo2z)*(I2[33] - (poz)*I3[33]);
V_Cl[in].Val[88] = U[0][2]*I0[58] + U[4][2]*I1[58]
           + (threeoo2z)*(I2[34] - (poz)*I3[34])
           + (oo2zn)*I4[28];
V_Cl[in].Val[89] = U[0][2]*I0[59] + U[4][2]*I1[59]
           + (threeoo2z)*(I2[35] - (poz)*I3[35])
           + (twooo2zn)*I4[29];
  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_g0f0(iclass *V_Cl, int in) /* type = 23 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2z = 2.0*oo2z;
  double threeoo2z = 3.0*oo2z;
  double twooo2zn = 2.0*oo2zn;
  double threeoo2zn = 3.0*oo2zn;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_f0f0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_f0f0(V_Cl, V_Cl[in].operands[1]);
  I2 = build_d0f0(V_Cl, V_Cl[in].operands[2]);
  I3 = build_d0f0(V_Cl, V_Cl[in].operands[3]);
  I4 = build_f0d0(V_Cl, V_Cl[in].operands[4]);


V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0]
           + (threeoo2z)*(I2[0] - (poz)*I3[0])
           + (threeoo2zn)*I4[0];
V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1]
           + (threeoo2z)*(I2[1] - (poz)*I3[1])
           + (twooo2zn)*I4[1];
V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2]
           + (threeoo2z)*(I2[2] - (poz)*I3[2])
           + (twooo2zn)*I4[2];
V_Cl[in].Val[3] = U[0][0]*I0[3] + U[4][0]*I1[3]
           + (threeoo2z)*(I2[3] - (poz)*I3[3])
           + (oo2zn)*I4[3];
V_Cl[in].Val[4] = U[0][0]*I0[4] + U[4][0]*I1[4]
           + (threeoo2z)*(I2[4] - (poz)*I3[4])
           + (oo2zn)*I4[4];
V_Cl[in].Val[5] = U[0][0]*I0[5] + U[4][0]*I1[5]
           + (threeoo2z)*(I2[5] - (poz)*I3[5])
           + (oo2zn)*I4[5];
V_Cl[in].Val[6] = U[0][0]*I0[6] + U[4][0]*I1[6]
           + (threeoo2z)*(I2[6] - (poz)*I3[6]);
V_Cl[in].Val[7] = U[0][0]*I0[7] + U[4][0]*I1[7]
           + (threeoo2z)*(I2[7] - (poz)*I3[7]);
V_Cl[in].Val[8] = U[0][0]*I0[8] + U[4][0]*I1[8]
           + (threeoo2z)*(I2[8] - (poz)*I3[8]);
V_Cl[in].Val[9] = U[0][0]*I0[9] + U[4][0]*I1[9]
           + (threeoo2z)*(I2[9] - (poz)*I3[9]);
V_Cl[in].Val[10] = U[0][0]*I0[10] + U[4][0]*I1[10]
           + (twooo2z)*(I2[10] - (poz)*I3[10])
           + (threeoo2zn)*I4[6];
V_Cl[in].Val[11] = U[0][0]*I0[11] + U[4][0]*I1[11]
           + (twooo2z)*(I2[11] - (poz)*I3[11])
           + (twooo2zn)*I4[7];
V_Cl[in].Val[12] = U[0][0]*I0[12] + U[4][0]*I1[12]
           + (twooo2z)*(I2[12] - (poz)*I3[12])
           + (twooo2zn)*I4[8];
V_Cl[in].Val[13] = U[0][0]*I0[13] + U[4][0]*I1[13]
           + (twooo2z)*(I2[13] - (poz)*I3[13])
           + (oo2zn)*I4[9];
V_Cl[in].Val[14] = U[0][0]*I0[14] + U[4][0]*I1[14]
           + (twooo2z)*(I2[14] - (poz)*I3[14])
           + (oo2zn)*I4[10];
V_Cl[in].Val[15] = U[0][0]*I0[15] + U[4][0]*I1[15]
           + (twooo2z)*(I2[15] - (poz)*I3[15])
           + (oo2zn)*I4[11];
V_Cl[in].Val[16] = U[0][0]*I0[16] + U[4][0]*I1[16]
           + (twooo2z)*(I2[16] - (poz)*I3[16]);
V_Cl[in].Val[17] = U[0][0]*I0[17] + U[4][0]*I1[17]
           + (twooo2z)*(I2[17] - (poz)*I3[17]);
V_Cl[in].Val[18] = U[0][0]*I0[18] + U[4][0]*I1[18]
           + (twooo2z)*(I2[18] - (poz)*I3[18]);
V_Cl[in].Val[19] = U[0][0]*I0[19] + U[4][0]*I1[19]
           + (twooo2z)*(I2[19] - (poz)*I3[19]);
V_Cl[in].Val[20] = U[0][0]*I0[20] + U[4][0]*I1[20]
           + (twooo2z)*(I2[20] - (poz)*I3[20])
           + (threeoo2zn)*I4[12];
V_Cl[in].Val[21] = U[0][0]*I0[21] + U[4][0]*I1[21]
           + (twooo2z)*(I2[21] - (poz)*I3[21])
           + (twooo2zn)*I4[13];
V_Cl[in].Val[22] = U[0][0]*I0[22] + U[4][0]*I1[22]
           + (twooo2z)*(I2[22] - (poz)*I3[22])
           + (twooo2zn)*I4[14];
V_Cl[in].Val[23] = U[0][0]*I0[23] + U[4][0]*I1[23]
           + (twooo2z)*(I2[23] - (poz)*I3[23])
           + (oo2zn)*I4[15];
V_Cl[in].Val[24] = U[0][0]*I0[24] + U[4][0]*I1[24]
           + (twooo2z)*(I2[24] - (poz)*I3[24])
           + (oo2zn)*I4[16];
V_Cl[in].Val[25] = U[0][0]*I0[25] + U[4][0]*I1[25]
           + (twooo2z)*(I2[25] - (poz)*I3[25])
           + (oo2zn)*I4[17];
V_Cl[in].Val[26] = U[0][0]*I0[26] + U[4][0]*I1[26]
           + (twooo2z)*(I2[26] - (poz)*I3[26]);
V_Cl[in].Val[27] = U[0][0]*I0[27] + U[4][0]*I1[27]
           + (twooo2z)*(I2[27] - (poz)*I3[27]);
V_Cl[in].Val[28] = U[0][0]*I0[28] + U[4][0]*I1[28]
           + (twooo2z)*(I2[28] - (poz)*I3[28]);
V_Cl[in].Val[29] = U[0][0]*I0[29] + U[4][0]*I1[29]
           + (twooo2z)*(I2[29] - (poz)*I3[29]);
V_Cl[in].Val[30] = U[0][0]*I0[30] + U[4][0]*I1[30]
           + (oo2z)*(I2[30] - (poz)*I3[30])
           + (threeoo2zn)*I4[18];
V_Cl[in].Val[31] = U[0][0]*I0[31] + U[4][0]*I1[31]
           + (oo2z)*(I2[31] - (poz)*I3[31])
           + (twooo2zn)*I4[19];
V_Cl[in].Val[32] = U[0][0]*I0[32] + U[4][0]*I1[32]
           + (oo2z)*(I2[32] - (poz)*I3[32])
           + (twooo2zn)*I4[20];
V_Cl[in].Val[33] = U[0][0]*I0[33] + U[4][0]*I1[33]
           + (oo2z)*(I2[33] - (poz)*I3[33])
           + (oo2zn)*I4[21];
V_Cl[in].Val[34] = U[0][0]*I0[34] + U[4][0]*I1[34]
           + (oo2z)*(I2[34] - (poz)*I3[34])
           + (oo2zn)*I4[22];
V_Cl[in].Val[35] = U[0][0]*I0[35] + U[4][0]*I1[35]
           + (oo2z)*(I2[35] - (poz)*I3[35])
           + (oo2zn)*I4[23];
V_Cl[in].Val[36] = U[0][0]*I0[36] + U[4][0]*I1[36]
           + (oo2z)*(I2[36] - (poz)*I3[36]);
V_Cl[in].Val[37] = U[0][0]*I0[37] + U[4][0]*I1[37]
           + (oo2z)*(I2[37] - (poz)*I3[37]);
V_Cl[in].Val[38] = U[0][0]*I0[38] + U[4][0]*I1[38]
           + (oo2z)*(I2[38] - (poz)*I3[38]);
V_Cl[in].Val[39] = U[0][0]*I0[39] + U[4][0]*I1[39]
           + (oo2z)*(I2[39] - (poz)*I3[39]);
V_Cl[in].Val[40] = U[0][0]*I0[40] + U[4][0]*I1[40]
           + (oo2z)*(I2[40] - (poz)*I3[40])
           + (threeoo2zn)*I4[24];
V_Cl[in].Val[41] = U[0][0]*I0[41] + U[4][0]*I1[41]
           + (oo2z)*(I2[41] - (poz)*I3[41])
           + (twooo2zn)*I4[25];
V_Cl[in].Val[42] = U[0][0]*I0[42] + U[4][0]*I1[42]
           + (oo2z)*(I2[42] - (poz)*I3[42])
           + (twooo2zn)*I4[26];
V_Cl[in].Val[43] = U[0][0]*I0[43] + U[4][0]*I1[43]
           + (oo2z)*(I2[43] - (poz)*I3[43])
           + (oo2zn)*I4[27];
V_Cl[in].Val[44] = U[0][0]*I0[44] + U[4][0]*I1[44]
           + (oo2z)*(I2[44] - (poz)*I3[44])
           + (oo2zn)*I4[28];
V_Cl[in].Val[45] = U[0][0]*I0[45] + U[4][0]*I1[45]
           + (oo2z)*(I2[45] - (poz)*I3[45])
           + (oo2zn)*I4[29];
V_Cl[in].Val[46] = U[0][0]*I0[46] + U[4][0]*I1[46]
           + (oo2z)*(I2[46] - (poz)*I3[46]);
V_Cl[in].Val[47] = U[0][0]*I0[47] + U[4][0]*I1[47]
           + (oo2z)*(I2[47] - (poz)*I3[47]);
V_Cl[in].Val[48] = U[0][0]*I0[48] + U[4][0]*I1[48]
           + (oo2z)*(I2[48] - (poz)*I3[48]);
V_Cl[in].Val[49] = U[0][0]*I0[49] + U[4][0]*I1[49]
           + (oo2z)*(I2[49] - (poz)*I3[49]);
V_Cl[in].Val[50] = U[0][0]*I0[50] + U[4][0]*I1[50]
           + (oo2z)*(I2[50] - (poz)*I3[50])
           + (threeoo2zn)*I4[30];
V_Cl[in].Val[51] = U[0][0]*I0[51] + U[4][0]*I1[51]
           + (oo2z)*(I2[51] - (poz)*I3[51])
           + (twooo2zn)*I4[31];
V_Cl[in].Val[52] = U[0][0]*I0[52] + U[4][0]*I1[52]
           + (oo2z)*(I2[52] - (poz)*I3[52])
           + (twooo2zn)*I4[32];
V_Cl[in].Val[53] = U[0][0]*I0[53] + U[4][0]*I1[53]
           + (oo2z)*(I2[53] - (poz)*I3[53])
           + (oo2zn)*I4[33];
V_Cl[in].Val[54] = U[0][0]*I0[54] + U[4][0]*I1[54]
           + (oo2z)*(I2[54] - (poz)*I3[54])
           + (oo2zn)*I4[34];
V_Cl[in].Val[55] = U[0][0]*I0[55] + U[4][0]*I1[55]
           + (oo2z)*(I2[55] - (poz)*I3[55])
           + (oo2zn)*I4[35];
V_Cl[in].Val[56] = U[0][0]*I0[56] + U[4][0]*I1[56]
           + (oo2z)*(I2[56] - (poz)*I3[56]);
V_Cl[in].Val[57] = U[0][0]*I0[57] + U[4][0]*I1[57]
           + (oo2z)*(I2[57] - (poz)*I3[57]);
V_Cl[in].Val[58] = U[0][0]*I0[58] + U[4][0]*I1[58]
           + (oo2z)*(I2[58] - (poz)*I3[58]);
V_Cl[in].Val[59] = U[0][0]*I0[59] + U[4][0]*I1[59]
           + (oo2z)*(I2[59] - (poz)*I3[59]);
V_Cl[in].Val[60] = U[0][0]*I0[60] + U[4][0]*I1[60]
           + (threeoo2zn)*I4[36];
V_Cl[in].Val[61] = U[0][0]*I0[61] + U[4][0]*I1[61]
           + (twooo2zn)*I4[37];
V_Cl[in].Val[62] = U[0][0]*I0[62] + U[4][0]*I1[62]
           + (twooo2zn)*I4[38];
V_Cl[in].Val[63] = U[0][0]*I0[63] + U[4][0]*I1[63]
           + (oo2zn)*I4[39];
V_Cl[in].Val[64] = U[0][0]*I0[64] + U[4][0]*I1[64]
           + (oo2zn)*I4[40];
V_Cl[in].Val[65] = U[0][0]*I0[65] + U[4][0]*I1[65]
           + (oo2zn)*I4[41];
V_Cl[in].Val[66] = U[0][0]*I0[66] + U[4][0]*I1[66];
V_Cl[in].Val[67] = U[0][0]*I0[67] + U[4][0]*I1[67];
V_Cl[in].Val[68] = U[0][0]*I0[68] + U[4][0]*I1[68];
V_Cl[in].Val[69] = U[0][0]*I0[69] + U[4][0]*I1[69];
V_Cl[in].Val[70] = U[0][0]*I0[70] + U[4][0]*I1[70]
           + (threeoo2zn)*I4[42];
V_Cl[in].Val[71] = U[0][0]*I0[71] + U[4][0]*I1[71]
           + (twooo2zn)*I4[43];
V_Cl[in].Val[72] = U[0][0]*I0[72] + U[4][0]*I1[72]
           + (twooo2zn)*I4[44];
V_Cl[in].Val[73] = U[0][0]*I0[73] + U[4][0]*I1[73]
           + (oo2zn)*I4[45];
V_Cl[in].Val[74] = U[0][0]*I0[74] + U[4][0]*I1[74]
           + (oo2zn)*I4[46];
V_Cl[in].Val[75] = U[0][0]*I0[75] + U[4][0]*I1[75]
           + (oo2zn)*I4[47];
V_Cl[in].Val[76] = U[0][0]*I0[76] + U[4][0]*I1[76];
V_Cl[in].Val[77] = U[0][0]*I0[77] + U[4][0]*I1[77];
V_Cl[in].Val[78] = U[0][0]*I0[78] + U[4][0]*I1[78];
V_Cl[in].Val[79] = U[0][0]*I0[79] + U[4][0]*I1[79];
V_Cl[in].Val[80] = U[0][0]*I0[80] + U[4][0]*I1[80]
           + (threeoo2zn)*I4[48];
V_Cl[in].Val[81] = U[0][0]*I0[81] + U[4][0]*I1[81]
           + (twooo2zn)*I4[49];
V_Cl[in].Val[82] = U[0][0]*I0[82] + U[4][0]*I1[82]
           + (twooo2zn)*I4[50];
V_Cl[in].Val[83] = U[0][0]*I0[83] + U[4][0]*I1[83]
           + (oo2zn)*I4[51];
V_Cl[in].Val[84] = U[0][0]*I0[84] + U[4][0]*I1[84]
           + (oo2zn)*I4[52];
V_Cl[in].Val[85] = U[0][0]*I0[85] + U[4][0]*I1[85]
           + (oo2zn)*I4[53];
V_Cl[in].Val[86] = U[0][0]*I0[86] + U[4][0]*I1[86];
V_Cl[in].Val[87] = U[0][0]*I0[87] + U[4][0]*I1[87];
V_Cl[in].Val[88] = U[0][0]*I0[88] + U[4][0]*I1[88];
V_Cl[in].Val[89] = U[0][0]*I0[89] + U[4][0]*I1[89];
V_Cl[in].Val[90] = U[0][0]*I0[90] + U[4][0]*I1[90]
           + (threeoo2zn)*I4[54];
V_Cl[in].Val[91] = U[0][0]*I0[91] + U[4][0]*I1[91]
           + (twooo2zn)*I4[55];
V_Cl[in].Val[92] = U[0][0]*I0[92] + U[4][0]*I1[92]
           + (twooo2zn)*I4[56];
V_Cl[in].Val[93] = U[0][0]*I0[93] + U[4][0]*I1[93]
           + (oo2zn)*I4[57];
V_Cl[in].Val[94] = U[0][0]*I0[94] + U[4][0]*I1[94]
           + (oo2zn)*I4[58];
V_Cl[in].Val[95] = U[0][0]*I0[95] + U[4][0]*I1[95]
           + (oo2zn)*I4[59];
V_Cl[in].Val[96] = U[0][0]*I0[96] + U[4][0]*I1[96];
V_Cl[in].Val[97] = U[0][0]*I0[97] + U[4][0]*I1[97];
V_Cl[in].Val[98] = U[0][0]*I0[98] + U[4][0]*I1[98];
V_Cl[in].Val[99] = U[0][0]*I0[99] + U[4][0]*I1[99];
V_Cl[in].Val[100] = U[0][1]*I0[60] + U[4][1]*I1[60]
           + (threeoo2z)*(I2[30] - (poz)*I3[30]);
V_Cl[in].Val[101] = U[0][1]*I0[61] + U[4][1]*I1[61]
           + (threeoo2z)*(I2[31] - (poz)*I3[31])
           + (oo2zn)*I4[36];
V_Cl[in].Val[102] = U[0][1]*I0[62] + U[4][1]*I1[62]
           + (threeoo2z)*(I2[32] - (poz)*I3[32]);
V_Cl[in].Val[103] = U[0][1]*I0[63] + U[4][1]*I1[63]
           + (threeoo2z)*(I2[33] - (poz)*I3[33])
           + (twooo2zn)*I4[37];
V_Cl[in].Val[104] = U[0][1]*I0[64] + U[4][1]*I1[64]
           + (threeoo2z)*(I2[34] - (poz)*I3[34])
           + (oo2zn)*I4[38];
V_Cl[in].Val[105] = U[0][1]*I0[65] + U[4][1]*I1[65]
           + (threeoo2z)*(I2[35] - (poz)*I3[35]);
V_Cl[in].Val[106] = U[0][1]*I0[66] + U[4][1]*I1[66]
           + (threeoo2z)*(I2[36] - (poz)*I3[36])
           + (threeoo2zn)*I4[39];
V_Cl[in].Val[107] = U[0][1]*I0[67] + U[4][1]*I1[67]
           + (threeoo2z)*(I2[37] - (poz)*I3[37])
           + (twooo2zn)*I4[40];
V_Cl[in].Val[108] = U[0][1]*I0[68] + U[4][1]*I1[68]
           + (threeoo2z)*(I2[38] - (poz)*I3[38])
           + (oo2zn)*I4[41];
V_Cl[in].Val[109] = U[0][1]*I0[69] + U[4][1]*I1[69]
           + (threeoo2z)*(I2[39] - (poz)*I3[39]);
V_Cl[in].Val[110] = U[0][1]*I0[70] + U[4][1]*I1[70]
           + (twooo2z)*(I2[40] - (poz)*I3[40]);
V_Cl[in].Val[111] = U[0][1]*I0[71] + U[4][1]*I1[71]
           + (twooo2z)*(I2[41] - (poz)*I3[41])
           + (oo2zn)*I4[42];
V_Cl[in].Val[112] = U[0][1]*I0[72] + U[4][1]*I1[72]
           + (twooo2z)*(I2[42] - (poz)*I3[42]);
V_Cl[in].Val[113] = U[0][1]*I0[73] + U[4][1]*I1[73]
           + (twooo2z)*(I2[43] - (poz)*I3[43])
           + (twooo2zn)*I4[43];
V_Cl[in].Val[114] = U[0][1]*I0[74] + U[4][1]*I1[74]
           + (twooo2z)*(I2[44] - (poz)*I3[44])
           + (oo2zn)*I4[44];
V_Cl[in].Val[115] = U[0][1]*I0[75] + U[4][1]*I1[75]
           + (twooo2z)*(I2[45] - (poz)*I3[45]);
V_Cl[in].Val[116] = U[0][1]*I0[76] + U[4][1]*I1[76]
           + (twooo2z)*(I2[46] - (poz)*I3[46])
           + (threeoo2zn)*I4[45];
V_Cl[in].Val[117] = U[0][1]*I0[77] + U[4][1]*I1[77]
           + (twooo2z)*(I2[47] - (poz)*I3[47])
           + (twooo2zn)*I4[46];
V_Cl[in].Val[118] = U[0][1]*I0[78] + U[4][1]*I1[78]
           + (twooo2z)*(I2[48] - (poz)*I3[48])
           + (oo2zn)*I4[47];
V_Cl[in].Val[119] = U[0][1]*I0[79] + U[4][1]*I1[79]
           + (twooo2z)*(I2[49] - (poz)*I3[49]);
V_Cl[in].Val[120] = U[0][1]*I0[80] + U[4][1]*I1[80]
           + (oo2z)*(I2[50] - (poz)*I3[50]);
V_Cl[in].Val[121] = U[0][1]*I0[81] + U[4][1]*I1[81]
           + (oo2z)*(I2[51] - (poz)*I3[51])
           + (oo2zn)*I4[48];
V_Cl[in].Val[122] = U[0][1]*I0[82] + U[4][1]*I1[82]
           + (oo2z)*(I2[52] - (poz)*I3[52]);
V_Cl[in].Val[123] = U[0][1]*I0[83] + U[4][1]*I1[83]
           + (oo2z)*(I2[53] - (poz)*I3[53])
           + (twooo2zn)*I4[49];
V_Cl[in].Val[124] = U[0][1]*I0[84] + U[4][1]*I1[84]
           + (oo2z)*(I2[54] - (poz)*I3[54])
           + (oo2zn)*I4[50];
V_Cl[in].Val[125] = U[0][1]*I0[85] + U[4][1]*I1[85]
           + (oo2z)*(I2[55] - (poz)*I3[55]);
V_Cl[in].Val[126] = U[0][1]*I0[86] + U[4][1]*I1[86]
           + (oo2z)*(I2[56] - (poz)*I3[56])
           + (threeoo2zn)*I4[51];
V_Cl[in].Val[127] = U[0][1]*I0[87] + U[4][1]*I1[87]
           + (oo2z)*(I2[57] - (poz)*I3[57])
           + (twooo2zn)*I4[52];
V_Cl[in].Val[128] = U[0][1]*I0[88] + U[4][1]*I1[88]
           + (oo2z)*(I2[58] - (poz)*I3[58])
           + (oo2zn)*I4[53];
V_Cl[in].Val[129] = U[0][1]*I0[89] + U[4][1]*I1[89]
           + (oo2z)*(I2[59] - (poz)*I3[59]);
V_Cl[in].Val[130] = U[0][1]*I0[90] + U[4][1]*I1[90];
V_Cl[in].Val[131] = U[0][1]*I0[91] + U[4][1]*I1[91]
           + (oo2zn)*I4[54];
V_Cl[in].Val[132] = U[0][1]*I0[92] + U[4][1]*I1[92];
V_Cl[in].Val[133] = U[0][1]*I0[93] + U[4][1]*I1[93]
           + (twooo2zn)*I4[55];
V_Cl[in].Val[134] = U[0][1]*I0[94] + U[4][1]*I1[94]
           + (oo2zn)*I4[56];
V_Cl[in].Val[135] = U[0][1]*I0[95] + U[4][1]*I1[95];
V_Cl[in].Val[136] = U[0][1]*I0[96] + U[4][1]*I1[96]
           + (threeoo2zn)*I4[57];
V_Cl[in].Val[137] = U[0][1]*I0[97] + U[4][1]*I1[97]
           + (twooo2zn)*I4[58];
V_Cl[in].Val[138] = U[0][1]*I0[98] + U[4][1]*I1[98]
           + (oo2zn)*I4[59];
V_Cl[in].Val[139] = U[0][1]*I0[99] + U[4][1]*I1[99];
V_Cl[in].Val[140] = U[0][2]*I0[90] + U[4][2]*I1[90]
           + (threeoo2z)*(I2[50] - (poz)*I3[50]);
V_Cl[in].Val[141] = U[0][2]*I0[91] + U[4][2]*I1[91]
           + (threeoo2z)*(I2[51] - (poz)*I3[51]);
V_Cl[in].Val[142] = U[0][2]*I0[92] + U[4][2]*I1[92]
           + (threeoo2z)*(I2[52] - (poz)*I3[52])
           + (oo2zn)*I4[54];
V_Cl[in].Val[143] = U[0][2]*I0[93] + U[4][2]*I1[93]
           + (threeoo2z)*(I2[53] - (poz)*I3[53]);
V_Cl[in].Val[144] = U[0][2]*I0[94] + U[4][2]*I1[94]
           + (threeoo2z)*(I2[54] - (poz)*I3[54])
           + (oo2zn)*I4[55];
V_Cl[in].Val[145] = U[0][2]*I0[95] + U[4][2]*I1[95]
           + (threeoo2z)*(I2[55] - (poz)*I3[55])
           + (twooo2zn)*I4[56];
V_Cl[in].Val[146] = U[0][2]*I0[96] + U[4][2]*I1[96]
           + (threeoo2z)*(I2[56] - (poz)*I3[56]);
V_Cl[in].Val[147] = U[0][2]*I0[97] + U[4][2]*I1[97]
           + (threeoo2z)*(I2[57] - (poz)*I3[57])
           + (oo2zn)*I4[57];
V_Cl[in].Val[148] = U[0][2]*I0[98] + U[4][2]*I1[98]
           + (threeoo2z)*(I2[58] - (poz)*I3[58])
           + (twooo2zn)*I4[58];
V_Cl[in].Val[149] = U[0][2]*I0[99] + U[4][2]*I1[99]
           + (threeoo2z)*(I2[59] - (poz)*I3[59])
           + (threeoo2zn)*I4[59];
  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_f0g0(iclass *V_Cl, int in) /* type = 24 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2z = 2.0*oo2z;
  double twooo2zn = 2.0*oo2zn;
  double threeoo2zn = 3.0*oo2zn;
  double fouroo2zn = 4.0*oo2zn;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_d0g0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_d0g0(V_Cl, V_Cl[in].operands[1]);
  I2 = build_p0g0(V_Cl, V_Cl[in].operands[2]);
  I3 = build_p0g0(V_Cl, V_Cl[in].operands[3]);
  I4 = build_d0f0(V_Cl, V_Cl[in].operands[4]);


V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0]
           + (twooo2z)*(I2[0] - (poz)*I3[0])
           + (fouroo2zn)*I4[0];
V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1]
           + (twooo2z)*(I2[1] - (poz)*I3[1])
           + (threeoo2zn)*I4[1];
V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2]
           + (twooo2z)*(I2[2] - (poz)*I3[2])
           + (threeoo2zn)*I4[2];
V_Cl[in].Val[3] = U[0][0]*I0[3] + U[4][0]*I1[3]
           + (twooo2z)*(I2[3] - (poz)*I3[3])
           + (twooo2zn)*I4[3];
V_Cl[in].Val[4] = U[0][0]*I0[4] + U[4][0]*I1[4]
           + (twooo2z)*(I2[4] - (poz)*I3[4])
           + (twooo2zn)*I4[4];
V_Cl[in].Val[5] = U[0][0]*I0[5] + U[4][0]*I1[5]
           + (twooo2z)*(I2[5] - (poz)*I3[5])
           + (twooo2zn)*I4[5];
V_Cl[in].Val[6] = U[0][0]*I0[6] + U[4][0]*I1[6]
           + (twooo2z)*(I2[6] - (poz)*I3[6])
           + (oo2zn)*I4[6];
V_Cl[in].Val[7] = U[0][0]*I0[7] + U[4][0]*I1[7]
           + (twooo2z)*(I2[7] - (poz)*I3[7])
           + (oo2zn)*I4[7];
V_Cl[in].Val[8] = U[0][0]*I0[8] + U[4][0]*I1[8]
           + (twooo2z)*(I2[8] - (poz)*I3[8])
           + (oo2zn)*I4[8];
V_Cl[in].Val[9] = U[0][0]*I0[9] + U[4][0]*I1[9]
           + (twooo2z)*(I2[9] - (poz)*I3[9])
           + (oo2zn)*I4[9];
V_Cl[in].Val[10] = U[0][0]*I0[10] + U[4][0]*I1[10]
           + (twooo2z)*(I2[10] - (poz)*I3[10]);
V_Cl[in].Val[11] = U[0][0]*I0[11] + U[4][0]*I1[11]
           + (twooo2z)*(I2[11] - (poz)*I3[11]);
V_Cl[in].Val[12] = U[0][0]*I0[12] + U[4][0]*I1[12]
           + (twooo2z)*(I2[12] - (poz)*I3[12]);
V_Cl[in].Val[13] = U[0][0]*I0[13] + U[4][0]*I1[13]
           + (twooo2z)*(I2[13] - (poz)*I3[13]);
V_Cl[in].Val[14] = U[0][0]*I0[14] + U[4][0]*I1[14]
           + (twooo2z)*(I2[14] - (poz)*I3[14]);
V_Cl[in].Val[15] = U[0][0]*I0[15] + U[4][0]*I1[15]
           + (oo2z)*(I2[15] - (poz)*I3[15])
           + (fouroo2zn)*I4[10];
V_Cl[in].Val[16] = U[0][0]*I0[16] + U[4][0]*I1[16]
           + (oo2z)*(I2[16] - (poz)*I3[16])
           + (threeoo2zn)*I4[11];
V_Cl[in].Val[17] = U[0][0]*I0[17] + U[4][0]*I1[17]
           + (oo2z)*(I2[17] - (poz)*I3[17])
           + (threeoo2zn)*I4[12];
V_Cl[in].Val[18] = U[0][0]*I0[18] + U[4][0]*I1[18]
           + (oo2z)*(I2[18] - (poz)*I3[18])
           + (twooo2zn)*I4[13];
V_Cl[in].Val[19] = U[0][0]*I0[19] + U[4][0]*I1[19]
           + (oo2z)*(I2[19] - (poz)*I3[19])
           + (twooo2zn)*I4[14];
V_Cl[in].Val[20] = U[0][0]*I0[20] + U[4][0]*I1[20]
           + (oo2z)*(I2[20] - (poz)*I3[20])
           + (twooo2zn)*I4[15];
V_Cl[in].Val[21] = U[0][0]*I0[21] + U[4][0]*I1[21]
           + (oo2z)*(I2[21] - (poz)*I3[21])
           + (oo2zn)*I4[16];
V_Cl[in].Val[22] = U[0][0]*I0[22] + U[4][0]*I1[22]
           + (oo2z)*(I2[22] - (poz)*I3[22])
           + (oo2zn)*I4[17];
V_Cl[in].Val[23] = U[0][0]*I0[23] + U[4][0]*I1[23]
           + (oo2z)*(I2[23] - (poz)*I3[23])
           + (oo2zn)*I4[18];
V_Cl[in].Val[24] = U[0][0]*I0[24] + U[4][0]*I1[24]
           + (oo2z)*(I2[24] - (poz)*I3[24])
           + (oo2zn)*I4[19];
V_Cl[in].Val[25] = U[0][0]*I0[25] + U[4][0]*I1[25]
           + (oo2z)*(I2[25] - (poz)*I3[25]);
V_Cl[in].Val[26] = U[0][0]*I0[26] + U[4][0]*I1[26]
           + (oo2z)*(I2[26] - (poz)*I3[26]);
V_Cl[in].Val[27] = U[0][0]*I0[27] + U[4][0]*I1[27]
           + (oo2z)*(I2[27] - (poz)*I3[27]);
V_Cl[in].Val[28] = U[0][0]*I0[28] + U[4][0]*I1[28]
           + (oo2z)*(I2[28] - (poz)*I3[28]);
V_Cl[in].Val[29] = U[0][0]*I0[29] + U[4][0]*I1[29]
           + (oo2z)*(I2[29] - (poz)*I3[29]);
V_Cl[in].Val[30] = U[0][0]*I0[30] + U[4][0]*I1[30]
           + (oo2z)*(I2[30] - (poz)*I3[30])
           + (fouroo2zn)*I4[20];
V_Cl[in].Val[31] = U[0][0]*I0[31] + U[4][0]*I1[31]
           + (oo2z)*(I2[31] - (poz)*I3[31])
           + (threeoo2zn)*I4[21];
V_Cl[in].Val[32] = U[0][0]*I0[32] + U[4][0]*I1[32]
           + (oo2z)*(I2[32] - (poz)*I3[32])
           + (threeoo2zn)*I4[22];
V_Cl[in].Val[33] = U[0][0]*I0[33] + U[4][0]*I1[33]
           + (oo2z)*(I2[33] - (poz)*I3[33])
           + (twooo2zn)*I4[23];
V_Cl[in].Val[34] = U[0][0]*I0[34] + U[4][0]*I1[34]
           + (oo2z)*(I2[34] - (poz)*I3[34])
           + (twooo2zn)*I4[24];
V_Cl[in].Val[35] = U[0][0]*I0[35] + U[4][0]*I1[35]
           + (oo2z)*(I2[35] - (poz)*I3[35])
           + (twooo2zn)*I4[25];
V_Cl[in].Val[36] = U[0][0]*I0[36] + U[4][0]*I1[36]
           + (oo2z)*(I2[36] - (poz)*I3[36])
           + (oo2zn)*I4[26];
V_Cl[in].Val[37] = U[0][0]*I0[37] + U[4][0]*I1[37]
           + (oo2z)*(I2[37] - (poz)*I3[37])
           + (oo2zn)*I4[27];
V_Cl[in].Val[38] = U[0][0]*I0[38] + U[4][0]*I1[38]
           + (oo2z)*(I2[38] - (poz)*I3[38])
           + (oo2zn)*I4[28];
V_Cl[in].Val[39] = U[0][0]*I0[39] + U[4][0]*I1[39]
           + (oo2z)*(I2[39] - (poz)*I3[39])
           + (oo2zn)*I4[29];
V_Cl[in].Val[40] = U[0][0]*I0[40] + U[4][0]*I1[40]
           + (oo2z)*(I2[40] - (poz)*I3[40]);
V_Cl[in].Val[41] = U[0][0]*I0[41] + U[4][0]*I1[41]
           + (oo2z)*(I2[41] - (poz)*I3[41]);
V_Cl[in].Val[42] = U[0][0]*I0[42] + U[4][0]*I1[42]
           + (oo2z)*(I2[42] - (poz)*I3[42]);
V_Cl[in].Val[43] = U[0][0]*I0[43] + U[4][0]*I1[43]
           + (oo2z)*(I2[43] - (poz)*I3[43]);
V_Cl[in].Val[44] = U[0][0]*I0[44] + U[4][0]*I1[44]
           + (oo2z)*(I2[44] - (poz)*I3[44]);
V_Cl[in].Val[45] = U[0][0]*I0[45] + U[4][0]*I1[45]
           + (fouroo2zn)*I4[30];
V_Cl[in].Val[46] = U[0][0]*I0[46] + U[4][0]*I1[46]
           + (threeoo2zn)*I4[31];
V_Cl[in].Val[47] = U[0][0]*I0[47] + U[4][0]*I1[47]
           + (threeoo2zn)*I4[32];
V_Cl[in].Val[48] = U[0][0]*I0[48] + U[4][0]*I1[48]
           + (twooo2zn)*I4[33];
V_Cl[in].Val[49] = U[0][0]*I0[49] + U[4][0]*I1[49]
           + (twooo2zn)*I4[34];
V_Cl[in].Val[50] = U[0][0]*I0[50] + U[4][0]*I1[50]
           + (twooo2zn)*I4[35];
V_Cl[in].Val[51] = U[0][0]*I0[51] + U[4][0]*I1[51]
           + (oo2zn)*I4[36];
V_Cl[in].Val[52] = U[0][0]*I0[52] + U[4][0]*I1[52]
           + (oo2zn)*I4[37];
V_Cl[in].Val[53] = U[0][0]*I0[53] + U[4][0]*I1[53]
           + (oo2zn)*I4[38];
V_Cl[in].Val[54] = U[0][0]*I0[54] + U[4][0]*I1[54]
           + (oo2zn)*I4[39];
V_Cl[in].Val[55] = U[0][0]*I0[55] + U[4][0]*I1[55];
V_Cl[in].Val[56] = U[0][0]*I0[56] + U[4][0]*I1[56];
V_Cl[in].Val[57] = U[0][0]*I0[57] + U[4][0]*I1[57];
V_Cl[in].Val[58] = U[0][0]*I0[58] + U[4][0]*I1[58];
V_Cl[in].Val[59] = U[0][0]*I0[59] + U[4][0]*I1[59];
V_Cl[in].Val[60] = U[0][0]*I0[60] + U[4][0]*I1[60]
           + (fouroo2zn)*I4[40];
V_Cl[in].Val[61] = U[0][0]*I0[61] + U[4][0]*I1[61]
           + (threeoo2zn)*I4[41];
V_Cl[in].Val[62] = U[0][0]*I0[62] + U[4][0]*I1[62]
           + (threeoo2zn)*I4[42];
V_Cl[in].Val[63] = U[0][0]*I0[63] + U[4][0]*I1[63]
           + (twooo2zn)*I4[43];
V_Cl[in].Val[64] = U[0][0]*I0[64] + U[4][0]*I1[64]
           + (twooo2zn)*I4[44];
V_Cl[in].Val[65] = U[0][0]*I0[65] + U[4][0]*I1[65]
           + (twooo2zn)*I4[45];
V_Cl[in].Val[66] = U[0][0]*I0[66] + U[4][0]*I1[66]
           + (oo2zn)*I4[46];
V_Cl[in].Val[67] = U[0][0]*I0[67] + U[4][0]*I1[67]
           + (oo2zn)*I4[47];
V_Cl[in].Val[68] = U[0][0]*I0[68] + U[4][0]*I1[68]
           + (oo2zn)*I4[48];
V_Cl[in].Val[69] = U[0][0]*I0[69] + U[4][0]*I1[69]
           + (oo2zn)*I4[49];
V_Cl[in].Val[70] = U[0][0]*I0[70] + U[4][0]*I1[70];
V_Cl[in].Val[71] = U[0][0]*I0[71] + U[4][0]*I1[71];
V_Cl[in].Val[72] = U[0][0]*I0[72] + U[4][0]*I1[72];
V_Cl[in].Val[73] = U[0][0]*I0[73] + U[4][0]*I1[73];
V_Cl[in].Val[74] = U[0][0]*I0[74] + U[4][0]*I1[74];
V_Cl[in].Val[75] = U[0][0]*I0[75] + U[4][0]*I1[75]
           + (fouroo2zn)*I4[50];
V_Cl[in].Val[76] = U[0][0]*I0[76] + U[4][0]*I1[76]
           + (threeoo2zn)*I4[51];
V_Cl[in].Val[77] = U[0][0]*I0[77] + U[4][0]*I1[77]
           + (threeoo2zn)*I4[52];
V_Cl[in].Val[78] = U[0][0]*I0[78] + U[4][0]*I1[78]
           + (twooo2zn)*I4[53];
V_Cl[in].Val[79] = U[0][0]*I0[79] + U[4][0]*I1[79]
           + (twooo2zn)*I4[54];
V_Cl[in].Val[80] = U[0][0]*I0[80] + U[4][0]*I1[80]
           + (twooo2zn)*I4[55];
V_Cl[in].Val[81] = U[0][0]*I0[81] + U[4][0]*I1[81]
           + (oo2zn)*I4[56];
V_Cl[in].Val[82] = U[0][0]*I0[82] + U[4][0]*I1[82]
           + (oo2zn)*I4[57];
V_Cl[in].Val[83] = U[0][0]*I0[83] + U[4][0]*I1[83]
           + (oo2zn)*I4[58];
V_Cl[in].Val[84] = U[0][0]*I0[84] + U[4][0]*I1[84]
           + (oo2zn)*I4[59];
V_Cl[in].Val[85] = U[0][0]*I0[85] + U[4][0]*I1[85];
V_Cl[in].Val[86] = U[0][0]*I0[86] + U[4][0]*I1[86];
V_Cl[in].Val[87] = U[0][0]*I0[87] + U[4][0]*I1[87];
V_Cl[in].Val[88] = U[0][0]*I0[88] + U[4][0]*I1[88];
V_Cl[in].Val[89] = U[0][0]*I0[89] + U[4][0]*I1[89];
V_Cl[in].Val[90] = U[0][1]*I0[45] + U[4][1]*I1[45]
           + (twooo2z)*(I2[15] - (poz)*I3[15]);
V_Cl[in].Val[91] = U[0][1]*I0[46] + U[4][1]*I1[46]
           + (twooo2z)*(I2[16] - (poz)*I3[16])
           + (oo2zn)*I4[30];
V_Cl[in].Val[92] = U[0][1]*I0[47] + U[4][1]*I1[47]
           + (twooo2z)*(I2[17] - (poz)*I3[17]);
V_Cl[in].Val[93] = U[0][1]*I0[48] + U[4][1]*I1[48]
           + (twooo2z)*(I2[18] - (poz)*I3[18])
           + (twooo2zn)*I4[31];
V_Cl[in].Val[94] = U[0][1]*I0[49] + U[4][1]*I1[49]
           + (twooo2z)*(I2[19] - (poz)*I3[19])
           + (oo2zn)*I4[32];
V_Cl[in].Val[95] = U[0][1]*I0[50] + U[4][1]*I1[50]
           + (twooo2z)*(I2[20] - (poz)*I3[20]);
V_Cl[in].Val[96] = U[0][1]*I0[51] + U[4][1]*I1[51]
           + (twooo2z)*(I2[21] - (poz)*I3[21])
           + (threeoo2zn)*I4[33];
V_Cl[in].Val[97] = U[0][1]*I0[52] + U[4][1]*I1[52]
           + (twooo2z)*(I2[22] - (poz)*I3[22])
           + (twooo2zn)*I4[34];
V_Cl[in].Val[98] = U[0][1]*I0[53] + U[4][1]*I1[53]
           + (twooo2z)*(I2[23] - (poz)*I3[23])
           + (oo2zn)*I4[35];
V_Cl[in].Val[99] = U[0][1]*I0[54] + U[4][1]*I1[54]
           + (twooo2z)*(I2[24] - (poz)*I3[24]);
V_Cl[in].Val[100] = U[0][1]*I0[55] + U[4][1]*I1[55]
           + (twooo2z)*(I2[25] - (poz)*I3[25])
           + (fouroo2zn)*I4[36];
V_Cl[in].Val[101] = U[0][1]*I0[56] + U[4][1]*I1[56]
           + (twooo2z)*(I2[26] - (poz)*I3[26])
           + (threeoo2zn)*I4[37];
V_Cl[in].Val[102] = U[0][1]*I0[57] + U[4][1]*I1[57]
           + (twooo2z)*(I2[27] - (poz)*I3[27])
           + (twooo2zn)*I4[38];
V_Cl[in].Val[103] = U[0][1]*I0[58] + U[4][1]*I1[58]
           + (twooo2z)*(I2[28] - (poz)*I3[28])
           + (oo2zn)*I4[39];
V_Cl[in].Val[104] = U[0][1]*I0[59] + U[4][1]*I1[59]
           + (twooo2z)*(I2[29] - (poz)*I3[29]);
V_Cl[in].Val[105] = U[0][1]*I0[60] + U[4][1]*I1[60]
           + (oo2z)*(I2[30] - (poz)*I3[30]);
V_Cl[in].Val[106] = U[0][1]*I0[61] + U[4][1]*I1[61]
           + (oo2z)*(I2[31] - (poz)*I3[31])
           + (oo2zn)*I4[40];
V_Cl[in].Val[107] = U[0][1]*I0[62] + U[4][1]*I1[62]
           + (oo2z)*(I2[32] - (poz)*I3[32]);
V_Cl[in].Val[108] = U[0][1]*I0[63] + U[4][1]*I1[63]
           + (oo2z)*(I2[33] - (poz)*I3[33])
           + (twooo2zn)*I4[41];
V_Cl[in].Val[109] = U[0][1]*I0[64] + U[4][1]*I1[64]
           + (oo2z)*(I2[34] - (poz)*I3[34])
           + (oo2zn)*I4[42];
V_Cl[in].Val[110] = U[0][1]*I0[65] + U[4][1]*I1[65]
           + (oo2z)*(I2[35] - (poz)*I3[35]);
V_Cl[in].Val[111] = U[0][1]*I0[66] + U[4][1]*I1[66]
           + (oo2z)*(I2[36] - (poz)*I3[36])
           + (threeoo2zn)*I4[43];
V_Cl[in].Val[112] = U[0][1]*I0[67] + U[4][1]*I1[67]
           + (oo2z)*(I2[37] - (poz)*I3[37])
           + (twooo2zn)*I4[44];
V_Cl[in].Val[113] = U[0][1]*I0[68] + U[4][1]*I1[68]
           + (oo2z)*(I2[38] - (poz)*I3[38])
           + (oo2zn)*I4[45];
V_Cl[in].Val[114] = U[0][1]*I0[69] + U[4][1]*I1[69]
           + (oo2z)*(I2[39] - (poz)*I3[39]);
V_Cl[in].Val[115] = U[0][1]*I0[70] + U[4][1]*I1[70]
           + (oo2z)*(I2[40] - (poz)*I3[40])
           + (fouroo2zn)*I4[46];
V_Cl[in].Val[116] = U[0][1]*I0[71] + U[4][1]*I1[71]
           + (oo2z)*(I2[41] - (poz)*I3[41])
           + (threeoo2zn)*I4[47];
V_Cl[in].Val[117] = U[0][1]*I0[72] + U[4][1]*I1[72]
           + (oo2z)*(I2[42] - (poz)*I3[42])
           + (twooo2zn)*I4[48];
V_Cl[in].Val[118] = U[0][1]*I0[73] + U[4][1]*I1[73]
           + (oo2z)*(I2[43] - (poz)*I3[43])
           + (oo2zn)*I4[49];
V_Cl[in].Val[119] = U[0][1]*I0[74] + U[4][1]*I1[74]
           + (oo2z)*(I2[44] - (poz)*I3[44]);
V_Cl[in].Val[120] = U[0][1]*I0[75] + U[4][1]*I1[75];
V_Cl[in].Val[121] = U[0][1]*I0[76] + U[4][1]*I1[76]
           + (oo2zn)*I4[50];
V_Cl[in].Val[122] = U[0][1]*I0[77] + U[4][1]*I1[77];
V_Cl[in].Val[123] = U[0][1]*I0[78] + U[4][1]*I1[78]
           + (twooo2zn)*I4[51];
V_Cl[in].Val[124] = U[0][1]*I0[79] + U[4][1]*I1[79]
           + (oo2zn)*I4[52];
V_Cl[in].Val[125] = U[0][1]*I0[80] + U[4][1]*I1[80];
V_Cl[in].Val[126] = U[0][1]*I0[81] + U[4][1]*I1[81]
           + (threeoo2zn)*I4[53];
V_Cl[in].Val[127] = U[0][1]*I0[82] + U[4][1]*I1[82]
           + (twooo2zn)*I4[54];
V_Cl[in].Val[128] = U[0][1]*I0[83] + U[4][1]*I1[83]
           + (oo2zn)*I4[55];
V_Cl[in].Val[129] = U[0][1]*I0[84] + U[4][1]*I1[84];
V_Cl[in].Val[130] = U[0][1]*I0[85] + U[4][1]*I1[85]
           + (fouroo2zn)*I4[56];
V_Cl[in].Val[131] = U[0][1]*I0[86] + U[4][1]*I1[86]
           + (threeoo2zn)*I4[57];
V_Cl[in].Val[132] = U[0][1]*I0[87] + U[4][1]*I1[87]
           + (twooo2zn)*I4[58];
V_Cl[in].Val[133] = U[0][1]*I0[88] + U[4][1]*I1[88]
           + (oo2zn)*I4[59];
V_Cl[in].Val[134] = U[0][1]*I0[89] + U[4][1]*I1[89];
V_Cl[in].Val[135] = U[0][2]*I0[75] + U[4][2]*I1[75]
           + (twooo2z)*(I2[30] - (poz)*I3[30]);
V_Cl[in].Val[136] = U[0][2]*I0[76] + U[4][2]*I1[76]
           + (twooo2z)*(I2[31] - (poz)*I3[31]);
V_Cl[in].Val[137] = U[0][2]*I0[77] + U[4][2]*I1[77]
           + (twooo2z)*(I2[32] - (poz)*I3[32])
           + (oo2zn)*I4[50];
V_Cl[in].Val[138] = U[0][2]*I0[78] + U[4][2]*I1[78]
           + (twooo2z)*(I2[33] - (poz)*I3[33]);
V_Cl[in].Val[139] = U[0][2]*I0[79] + U[4][2]*I1[79]
           + (twooo2z)*(I2[34] - (poz)*I3[34])
           + (oo2zn)*I4[51];
V_Cl[in].Val[140] = U[0][2]*I0[80] + U[4][2]*I1[80]
           + (twooo2z)*(I2[35] - (poz)*I3[35])
           + (twooo2zn)*I4[52];
V_Cl[in].Val[141] = U[0][2]*I0[81] + U[4][2]*I1[81]
           + (twooo2z)*(I2[36] - (poz)*I3[36]);
V_Cl[in].Val[142] = U[0][2]*I0[82] + U[4][2]*I1[82]
           + (twooo2z)*(I2[37] - (poz)*I3[37])
           + (oo2zn)*I4[53];
V_Cl[in].Val[143] = U[0][2]*I0[83] + U[4][2]*I1[83]
           + (twooo2z)*(I2[38] - (poz)*I3[38])
           + (twooo2zn)*I4[54];
V_Cl[in].Val[144] = U[0][2]*I0[84] + U[4][2]*I1[84]
           + (twooo2z)*(I2[39] - (poz)*I3[39])
           + (threeoo2zn)*I4[55];
V_Cl[in].Val[145] = U[0][2]*I0[85] + U[4][2]*I1[85]
           + (twooo2z)*(I2[40] - (poz)*I3[40]);
V_Cl[in].Val[146] = U[0][2]*I0[86] + U[4][2]*I1[86]
           + (twooo2z)*(I2[41] - (poz)*I3[41])
           + (oo2zn)*I4[56];
V_Cl[in].Val[147] = U[0][2]*I0[87] + U[4][2]*I1[87]
           + (twooo2z)*(I2[42] - (poz)*I3[42])
           + (twooo2zn)*I4[57];
V_Cl[in].Val[148] = U[0][2]*I0[88] + U[4][2]*I1[88]
           + (twooo2z)*(I2[43] - (poz)*I3[43])
           + (threeoo2zn)*I4[58];
V_Cl[in].Val[149] = U[0][2]*I0[89] + U[4][2]*I1[89]
           + (twooo2z)*(I2[44] - (poz)*I3[44])
           + (fouroo2zn)*I4[59];
  V_done[in] = 1;
  return V_Cl[in].Val;
}

double *
TwoBodyIntJF::build_g0g0(iclass *V_Cl, int in) /* type = 25 */
{
  double *I0, *I1, *I2, *I3, *I4;
  double twooo2z = 2.0*oo2z;
  double threeoo2z = 3.0*oo2z;
  double twooo2zn = 2.0*oo2zn;
  double threeoo2zn = 3.0*oo2zn;
  double fouroo2zn = 4.0*oo2zn;

  if(V_done[in]) return V_Cl[in].Val;

  I0 = build_f0g0(V_Cl, V_Cl[in].operands[0]);
  I1 = build_f0g0(V_Cl, V_Cl[in].operands[1]);
  I2 = build_d0g0(V_Cl, V_Cl[in].operands[2]);
  I3 = build_d0g0(V_Cl, V_Cl[in].operands[3]);
  I4 = build_f0f0(V_Cl, V_Cl[in].operands[4]);


V_Cl[in].Val[0] = U[0][0]*I0[0] + U[4][0]*I1[0]
           + (threeoo2z)*(I2[0] - (poz)*I3[0])
           + (fouroo2zn)*I4[0];
V_Cl[in].Val[1] = U[0][0]*I0[1] + U[4][0]*I1[1]
           + (threeoo2z)*(I2[1] - (poz)*I3[1])
           + (threeoo2zn)*I4[1];
V_Cl[in].Val[2] = U[0][0]*I0[2] + U[4][0]*I1[2]
           + (threeoo2z)*(I2[2] - (poz)*I3[2])
           + (threeoo2zn)*I4[2];
V_Cl[in].Val[3] = U[0][0]*I0[3] + U[4][0]*I1[3]
           + (threeoo2z)*(I2[3] - (poz)*I3[3])
           + (twooo2zn)*I4[3];
V_Cl[in].Val[4] = U[0][0]*I0[4] + U[4][0]*I1[4]
           + (threeoo2z)*(I2[4] - (poz)*I3[4])
           + (twooo2zn)*I4[4];
V_Cl[in].Val[5] = U[0][0]*I0[5] + U[4][0]*I1[5]
           + (threeoo2z)*(I2[5] - (poz)*I3[5])
           + (twooo2zn)*I4[5];
V_Cl[in].Val[6] = U[0][0]*I0[6] + U[4][0]*I1[6]
           + (threeoo2z)*(I2[6] - (poz)*I3[6])
           + (oo2zn)*I4[6];
V_Cl[in].Val[7] = U[0][0]*I0[7] + U[4][0]*I1[7]
           + (threeoo2z)*(I2[7] - (poz)*I3[7])
           + (oo2zn)*I4[7];
V_Cl[in].Val[8] = U[0][0]*I0[8] + U[4][0]*I1[8]
           + (threeoo2z)*(I2[8] - (poz)*I3[8])
           + (oo2zn)*I4[8];
V_Cl[in].Val[9] = U[0][0]*I0[9] + U[4][0]*I1[9]
           + (threeoo2z)*(I2[9] - (poz)*I3[9])
           + (oo2zn)*I4[9];
V_Cl[in].Val[10] = U[0][0]*I0[10] + U[4][0]*I1[10]
           + (threeoo2z)*(I2[10] - (poz)*I3[10]);
V_Cl[in].Val[11] = U[0][0]*I0[11] + U[4][0]*I1[11]
           + (threeoo2z)*(I2[11] - (poz)*I3[11]);
V_Cl[in].Val[12] = U[0][0]*I0[12] + U[4][0]*I1[12]
           + (threeoo2z)*(I2[12] - (poz)*I3[12]);
V_Cl[in].Val[13] = U[0][0]*I0[13] + U[4][0]*I1[13]
           + (threeoo2z)*(I2[13] - (poz)*I3[13]);
V_Cl[in].Val[14] = U[0][0]*I0[14] + U[4][0]*I1[14]
           + (threeoo2z)*(I2[14] - (poz)*I3[14]);
V_Cl[in].Val[15] = U[0][0]*I0[15] + U[4][0]*I1[15]
           + (twooo2z)*(I2[15] - (poz)*I3[15])
           + (fouroo2zn)*I4[10];
V_Cl[in].Val[16] = U[0][0]*I0[16] + U[4][0]*I1[16]
           + (twooo2z)*(I2[16] - (poz)*I3[16])
           + (threeoo2zn)*I4[11];
V_Cl[in].Val[17] = U[0][0]*I0[17] + U[4][0]*I1[17]
           + (twooo2z)*(I2[17] - (poz)*I3[17])
           + (threeoo2zn)*I4[12];
V_Cl[in].Val[18] = U[0][0]*I0[18] + U[4][0]*I1[18]
           + (twooo2z)*(I2[18] - (poz)*I3[18])
           + (twooo2zn)*I4[13];
V_Cl[in].Val[19] = U[0][0]*I0[19] + U[4][0]*I1[19]
           + (twooo2z)*(I2[19] - (poz)*I3[19])
           + (twooo2zn)*I4[14];
V_Cl[in].Val[20] = U[0][0]*I0[20] + U[4][0]*I1[20]
           + (twooo2z)*(I2[20] - (poz)*I3[20])
           + (twooo2zn)*I4[15];
V_Cl[in].Val[21] = U[0][0]*I0[21] + U[4][0]*I1[21]
           + (twooo2z)*(I2[21] - (poz)*I3[21])
           + (oo2zn)*I4[16];
V_Cl[in].Val[22] = U[0][0]*I0[22] + U[4][0]*I1[22]
           + (twooo2z)*(I2[22] - (poz)*I3[22])
           + (oo2zn)*I4[17];
V_Cl[in].Val[23] = U[0][0]*I0[23] + U[4][0]*I1[23]
           + (twooo2z)*(I2[23] - (poz)*I3[23])
           + (oo2zn)*I4[18];
V_Cl[in].Val[24] = U[0][0]*I0[24] + U[4][0]*I1[24]
           + (twooo2z)*(I2[24] - (poz)*I3[24])
           + (oo2zn)*I4[19];
V_Cl[in].Val[25] = U[0][0]*I0[25] + U[4][0]*I1[25]
           + (twooo2z)*(I2[25] - (poz)*I3[25]);
V_Cl[in].Val[26] = U[0][0]*I0[26] + U[4][0]*I1[26]
           + (twooo2z)*(I2[26] - (poz)*I3[26]);
V_Cl[in].Val[27] = U[0][0]*I0[27] + U[4][0]*I1[27]
           + (twooo2z)*(I2[27] - (poz)*I3[27]);
V_Cl[in].Val[28] = U[0][0]*I0[28] + U[4][0]*I1[28]
           + (twooo2z)*(I2[28] - (poz)*I3[28]);
V_Cl[in].Val[29] = U[0][0]*I0[29] + U[4][0]*I1[29]
           + (twooo2z)*(I2[29] - (poz)*I3[29]);
V_Cl[in].Val[30] = U[0][0]*I0[30] + U[4][0]*I1[30]
           + (twooo2z)*(I2[30] - (poz)*I3[30])
           + (fouroo2zn)*I4[20];
V_Cl[in].Val[31] = U[0][0]*I0[31] + U[4][0]*I1[31]
           + (twooo2z)*(I2[31] - (poz)*I3[31])
           + (threeoo2zn)*I4[21];
V_Cl[in].Val[32] = U[0][0]*I0[32] + U[4][0]*I1[32]
           + (twooo2z)*(I2[32] - (poz)*I3[32])
           + (threeoo2zn)*I4[22];
V_Cl[in].Val[33] = U[0][0]*I0[33] + U[4][0]*I1[33]
           + (twooo2z)*(I2[33] - (poz)*I3[33])
           + (twooo2zn)*I4[23];
V_Cl[in].Val[34] = U[0][0]*I0[34] + U[4][0]*I1[34]
           + (twooo2z)*(I2[34] - (poz)*I3[34])
           + (twooo2zn)*I4[24];
V_Cl[in].Val[35] = U[0][0]*I0[35] + U[4][0]*I1[35]
           + (twooo2z)*(I2[35] - (poz)*I3[35])
           + (twooo2zn)*I4[25];
V_Cl[in].Val[36] = U[0][0]*I0[36] + U[4][0]*I1[36]
           + (twooo2z)*(I2[36] - (poz)*I3[36])
           + (oo2zn)*I4[26];
V_Cl[in].Val[37] = U[0][0]*I0[37] + U[4][0]*I1[37]
           + (twooo2z)*(I2[37] - (poz)*I3[37])
           + (oo2zn)*I4[27];
V_Cl[in].Val[38] = U[0][0]*I0[38] + U[4][0]*I1[38]
           + (twooo2z)*(I2[38] - (poz)*I3[38])
           + (oo2zn)*I4[28];
V_Cl[in].Val[39] = U[0][0]*I0[39] + U[4][0]*I1[39]
           + (twooo2z)*(I2[39] - (poz)*I3[39])
           + (oo2zn)*I4[29];
V_Cl[in].Val[40] = U[0][0]*I0[40] + U[4][0]*I1[40]
           + (twooo2z)*(I2[40] - (poz)*I3[40]);
V_Cl[in].Val[41] = U[0][0]*I0[41] + U[4][0]*I1[41]
           + (twooo2z)*(I2[41] - (poz)*I3[41]);
V_Cl[in].Val[42] = U[0][0]*I0[42] + U[4][0]*I1[42]
           + (twooo2z)*(I2[42] - (poz)*I3[42]);
V_Cl[in].Val[43] = U[0][0]*I0[43] + U[4][0]*I1[43]
           + (twooo2z)*(I2[43] - (poz)*I3[43]);
V_Cl[in].Val[44] = U[0][0]*I0[44] + U[4][0]*I1[44]
           + (twooo2z)*(I2[44] - (poz)*I3[44]);
V_Cl[in].Val[45] = U[0][0]*I0[45] + U[4][0]*I1[45]
           + (oo2z)*(I2[45] - (poz)*I3[45])
           + (fouroo2zn)*I4[30];
V_Cl[in].Val[46] = U[0][0]*I0[46] + U[4][0]*I1[46]
           + (oo2z)*(I2[46] - (poz)*I3[46])
           + (threeoo2zn)*I4[31];
V_Cl[in].Val[47] = U[0][0]*I0[47] + U[4][0]*I1[47]
           + (oo2z)*(I2[47] - (poz)*I3[47])
           + (threeoo2zn)*I4[32];
V_Cl[in].Val[48] = U[0][0]*I0[48] + U[4][0]*I1[48]
           + (oo2z)*(I2[48] - (poz)*I3[48])
           + (twooo2zn)*I4[33];
V_Cl[in].Val[49] = U[0][0]*I0[49] + U[4][0]*I1[49]
           + (oo2z)*(I2[49] - (poz)*I3[49])
           + (twooo2zn)*I4[34];
V_Cl[in].Val[50] = U[0][0]*I0[50] + U[4][0]*I1[50]
           + (oo2z)*(I2[50] - (poz)*I3[50])
           + (twooo2zn)*I4[35];
V_Cl[in].Val[51] = U[0][0]*I0[51] + U[4][0]*I1[51]
           + (oo2z)*(I2[51] - (poz)*I3[51])
           + (oo2zn)*I4[36];
V_Cl[in].Val[52] = U[0][0]*I0[52] + U[4][0]*I1[52]
           + (oo2z)*(I2[52] - (poz)*I3[52])
           + (oo2zn)*I4[37];
V_Cl[in].Val[53] = U[0][0]*I0[53] + U[4][0]*I1[53]
           + (oo2z)*(I2[53] - (poz)*I3[53])
           + (oo2zn)*I4[38];
V_Cl[in].Val[54] = U[0][0]*I0[54] + U[4][0]*I1[54]
           + (oo2z)*(I2[54] - (poz)*I3[54])
           + (oo2zn)*I4[39];
V_Cl[in].Val[55] = U[0][0]*I0[55] + U[4][0]*I1[55]
           + (oo2z)*(I2[55] - (poz)*I3[55]);
V_Cl[in].Val[56] = U[0][0]*I0[56] + U[4][0]*I1[56]
           + (oo2z)*(I2[56] - (poz)*I3[56]);
V_Cl[in].Val[57] = U[0][0]*I0[57] + U[4][0]*I1[57]
           + (oo2z)*(I2[57] - (poz)*I3[57]);
V_Cl[in].Val[58] = U[0][0]*I0[58] + U[4][0]*I1[58]
           + (oo2z)*(I2[58] - (poz)*I3[58]);
V_Cl[in].Val[59] = U[0][0]*I0[59] + U[4][0]*I1[59]
           + (oo2z)*(I2[59] - (poz)*I3[59]);
V_Cl[in].Val[60] = U[0][0]*I0[60] + U[4][0]*I1[60]
           + (oo2z)*(I2[60] - (poz)*I3[60])
           + (fouroo2zn)*I4[40];
V_Cl[in].Val[61] = U[0][0]*I0[61] + U[4][0]*I1[61]
           + (oo2z)*(I2[61] - (poz)*I3[61])
           + (threeoo2zn)*I4[41];
V_Cl[in].Val[62] = U[0][0]*I0[62] + U[4][0]*I1[62]
           + (oo2z)*(I2[62] - (poz)*I3[62])
           + (threeoo2zn)*I4[42];
V_Cl[in].Val[63] = U[0][0]*I0[63] + U[4][0]*I1[63]
           + (oo2z)*(I2[63] - (poz)*I3[63])
           + (twooo2zn)*I4[43];
V_Cl[in].Val[64] = U[0][0]*I0[64] + U[4][0]*I1[64]
           + (oo2z)*(I2[64] - (poz)*I3[64])
           + (twooo2zn)*I4[44];
V_Cl[in].Val[65] = U[0][0]*I0[65] + U[4][0]*I1[65]
           + (oo2z)*(I2[65] - (poz)*I3[65])
           + (twooo2zn)*I4[45];
V_Cl[in].Val[66] = U[0][0]*I0[66] + U[4][0]*I1[66]
           + (oo2z)*(I2[66] - (poz)*I3[66])
           + (oo2zn)*I4[46];
V_Cl[in].Val[67] = U[0][0]*I0[67] + U[4][0]*I1[67]
           + (oo2z)*(I2[67] - (poz)*I3[67])
           + (oo2zn)*I4[47];
V_Cl[in].Val[68] = U[0][0]*I0[68] + U[4][0]*I1[68]
           + (oo2z)*(I2[68] - (poz)*I3[68])
           + (oo2zn)*I4[48];
V_Cl[in].Val[69] = U[0][0]*I0[69] + U[4][0]*I1[69]
           + (oo2z)*(I2[69] - (poz)*I3[69])
           + (oo2zn)*I4[49];
V_Cl[in].Val[70] = U[0][0]*I0[70] + U[4][0]*I1[70]
           + (oo2z)*(I2[70] - (poz)*I3[70]);
V_Cl[in].Val[71] = U[0][0]*I0[71] + U[4][0]*I1[71]
           + (oo2z)*(I2[71] - (poz)*I3[71]);
V_Cl[in].Val[72] = U[0][0]*I0[72] + U[4][0]*I1[72]
           + (oo2z)*(I2[72] - (poz)*I3[72]);
V_Cl[in].Val[73] = U[0][0]*I0[73] + U[4][0]*I1[73]
           + (oo2z)*(I2[73] - (poz)*I3[73]);
V_Cl[in].Val[74] = U[0][0]*I0[74] + U[4][0]*I1[74]
           + (oo2z)*(I2[74] - (poz)*I3[74]);
V_Cl[in].Val[75] = U[0][0]*I0[75] + U[4][0]*I1[75]
           + (oo2z)*(I2[75] - (poz)*I3[75])
           + (fouroo2zn)*I4[50];
V_Cl[in].Val[76] = U[0][0]*I0[76] + U[4][0]*I1[76]
           + (oo2z)*(I2[76] - (poz)*I3[76])
           + (threeoo2zn)*I4[51];
V_Cl[in].Val[77] = U[0][0]*I0[77] + U[4][0]*I1[77]
           + (oo2z)*(I2[77] - (poz)*I3[77])
           + (threeoo2zn)*I4[52];
V_Cl[in].Val[78] = U[0][0]*I0[78] + U[4][0]*I1[78]
           + (oo2z)*(I2[78] - (poz)*I3[78])
           + (twooo2zn)*I4[53];
V_Cl[in].Val[79] = U[0][0]*I0[79] + U[4][0]*I1[79]
           + (oo2z)*(I2[79] - (poz)*I3[79])
           + (twooo2zn)*I4[54];
V_Cl[in].Val[80] = U[0][0]*I0[80] + U[4][0]*I1[80]
           + (oo2z)*(I2[80] - (poz)*I3[80])
           + (twooo2zn)*I4[55];
V_Cl[in].Val[81] = U[0][0]*I0[81] + U[4][0]*I1[81]
           + (oo2z)*(I2[81] - (poz)*I3[81])
           + (oo2zn)*I4[56];
V_Cl[in].Val[82] = U[0][0]*I0[82] + U[4][0]*I1[82]
           + (oo2z)*(I2[82] - (poz)*I3[82])
           + (oo2zn)*I4[57];
V_Cl[in].Val[83] = U[0][0]*I0[83] + U[4][0]*I1[83]
           + (oo2z)*(I2[83] - (poz)*I3[83])
           + (oo2zn)*I4[58];
V_Cl[in].Val[84] = U[0][0]*I0[84] + U[4][0]*I1[84]
           + (oo2z)*(I2[84] - (poz)*I3[84])
           + (oo2zn)*I4[59];
V_Cl[in].Val[85] = U[0][0]*I0[85] + U[4][0]*I1[85]
           + (oo2z)*(I2[85] - (poz)*I3[85]);
V_Cl[in].Val[86] = U[0][0]*I0[86] + U[4][0]*I1[86]
           + (oo2z)*(I2[86] - (poz)*I3[86]);
V_Cl[in].Val[87] = U[0][0]*I0[87] + U[4][0]*I1[87]
           + (oo2z)*(I2[87] - (poz)*I3[87]);
V_Cl[in].Val[88] = U[0][0]*I0[88] + U[4][0]*I1[88]
           + (oo2z)*(I2[88] - (poz)*I3[88]);
V_Cl[in].Val[89] = U[0][0]*I0[89] + U[4][0]*I1[89]
           + (oo2z)*(I2[89] - (poz)*I3[89]);
V_Cl[in].Val[90] = U[0][0]*I0[90] + U[4][0]*I1[90]
           + (fouroo2zn)*I4[60];
V_Cl[in].Val[91] = U[0][0]*I0[91] + U[4][0]*I1[91]
           + (threeoo2zn)*I4[61];
V_Cl[in].Val[92] = U[0][0]*I0[92] + U[4][0]*I1[92]
           + (threeoo2zn)*I4[62];
V_Cl[in].Val[93] = U[0][0]*I0[93] + U[4][0]*I1[93]
           + (twooo2zn)*I4[63];
V_Cl[in].Val[94] = U[0][0]*I0[94] + U[4][0]*I1[94]
           + (twooo2zn)*I4[64];
V_Cl[in].Val[95] = U[0][0]*I0[95] + U[4][0]*I1[95]
           + (twooo2zn)*I4[65];
V_Cl[in].Val[96] = U[0][0]*I0[96] + U[4][0]*I1[96]
           + (oo2zn)*I4[66];
V_Cl[in].Val[97] = U[0][0]*I0[97] + U[4][0]*I1[97]
           + (oo2zn)*I4[67];
V_Cl[in].Val[98] = U[0][0]*I0[98] + U[4][0]*I1[98]
           + (oo2zn)*I4[68];
V_Cl[in].Val[99] = U[0][0]*I0[99] + U[4][0]*I1[99]
           + (oo2zn)*I4[69];
V_Cl[in].Val[100] = U[0][0]*I0[100] + U[4][0]*I1[100];
V_Cl[in].Val[101] = U[0][0]*I0[101] + U[4][0]*I1[101];
V_Cl[in].Val[102] = U[0][0]*I0[102] + U[4][0]*I1[102];
V_Cl[in].Val[103] = U[0][0]*I0[103] + U[4][0]*I1[103];
V_Cl[in].Val[104] = U[0][0]*I0[104] + U[4][0]*I1[104];
V_Cl[in].Val[105] = U[0][0]*I0[105] + U[4][0]*I1[105]
           + (fouroo2zn)*I4[70];
V_Cl[in].Val[106] = U[0][0]*I0[106] + U[4][0]*I1[106]
           + (threeoo2zn)*I4[71];
V_Cl[in].Val[107] = U[0][0]*I0[107] + U[4][0]*I1[107]
           + (threeoo2zn)*I4[72];
V_Cl[in].Val[108] = U[0][0]*I0[108] + U[4][0]*I1[108]
           + (twooo2zn)*I4[73];
V_Cl[in].Val[109] = U[0][0]*I0[109] + U[4][0]*I1[109]
           + (twooo2zn)*I4[74];
V_Cl[in].Val[110] = U[0][0]*I0[110] + U[4][0]*I1[110]
           + (twooo2zn)*I4[75];
V_Cl[in].Val[111] = U[0][0]*I0[111] + U[4][0]*I1[111]
           + (oo2zn)*I4[76];
V_Cl[in].Val[112] = U[0][0]*I0[112] + U[4][0]*I1[112]
           + (oo2zn)*I4[77];
V_Cl[in].Val[113] = U[0][0]*I0[113] + U[4][0]*I1[113]
           + (oo2zn)*I4[78];
V_Cl[in].Val[114] = U[0][0]*I0[114] + U[4][0]*I1[114]
           + (oo2zn)*I4[79];
V_Cl[in].Val[115] = U[0][0]*I0[115] + U[4][0]*I1[115];
V_Cl[in].Val[116] = U[0][0]*I0[116] + U[4][0]*I1[116];
V_Cl[in].Val[117] = U[0][0]*I0[117] + U[4][0]*I1[117];
V_Cl[in].Val[118] = U[0][0]*I0[118] + U[4][0]*I1[118];
V_Cl[in].Val[119] = U[0][0]*I0[119] + U[4][0]*I1[119];
V_Cl[in].Val[120] = U[0][0]*I0[120] + U[4][0]*I1[120]
           + (fouroo2zn)*I4[80];
V_Cl[in].Val[121] = U[0][0]*I0[121] + U[4][0]*I1[121]
           + (threeoo2zn)*I4[81];
V_Cl[in].Val[122] = U[0][0]*I0[122] + U[4][0]*I1[122]
           + (threeoo2zn)*I4[82];
V_Cl[in].Val[123] = U[0][0]*I0[123] + U[4][0]*I1[123]
           + (twooo2zn)*I4[83];
V_Cl[in].Val[124] = U[0][0]*I0[124] + U[4][0]*I1[124]
           + (twooo2zn)*I4[84];
V_Cl[in].Val[125] = U[0][0]*I0[125] + U[4][0]*I1[125]
           + (twooo2zn)*I4[85];
V_Cl[in].Val[126] = U[0][0]*I0[126] + U[4][0]*I1[126]
           + (oo2zn)*I4[86];
V_Cl[in].Val[127] = U[0][0]*I0[127] + U[4][0]*I1[127]
           + (oo2zn)*I4[87];
V_Cl[in].Val[128] = U[0][0]*I0[128] + U[4][0]*I1[128]
           + (oo2zn)*I4[88];
V_Cl[in].Val[129] = U[0][0]*I0[129] + U[4][0]*I1[129]
           + (oo2zn)*I4[89];
V_Cl[in].Val[130] = U[0][0]*I0[130] + U[4][0]*I1[130];
V_Cl[in].Val[131] = U[0][0]*I0[131] + U[4][0]*I1[131];
V_Cl[in].Val[132] = U[0][0]*I0[132] + U[4][0]*I1[132];
V_Cl[in].Val[133] = U[0][0]*I0[133] + U[4][0]*I1[133];
V_Cl[in].Val[134] = U[0][0]*I0[134] + U[4][0]*I1[134];
V_Cl[in].Val[135] = U[0][0]*I0[135] + U[4][0]*I1[135]
           + (fouroo2zn)*I4[90];
V_Cl[in].Val[136] = U[0][0]*I0[136] + U[4][0]*I1[136]
           + (threeoo2zn)*I4[91];
V_Cl[in].Val[137] = U[0][0]*I0[137] + U[4][0]*I1[137]
           + (threeoo2zn)*I4[92];
V_Cl[in].Val[138] = U[0][0]*I0[138] + U[4][0]*I1[138]
           + (twooo2zn)*I4[93];
V_Cl[in].Val[139] = U[0][0]*I0[139] + U[4][0]*I1[139]
           + (twooo2zn)*I4[94];
V_Cl[in].Val[140] = U[0][0]*I0[140] + U[4][0]*I1[140]
           + (twooo2zn)*I4[95];
V_Cl[in].Val[141] = U[0][0]*I0[141] + U[4][0]*I1[141]
           + (oo2zn)*I4[96];
V_Cl[in].Val[142] = U[0][0]*I0[142] + U[4][0]*I1[142]
           + (oo2zn)*I4[97];
V_Cl[in].Val[143] = U[0][0]*I0[143] + U[4][0]*I1[143]
           + (oo2zn)*I4[98];
V_Cl[in].Val[144] = U[0][0]*I0[144] + U[4][0]*I1[144]
           + (oo2zn)*I4[99];
V_Cl[in].Val[145] = U[0][0]*I0[145] + U[4][0]*I1[145];
V_Cl[in].Val[146] = U[0][0]*I0[146] + U[4][0]*I1[146];
V_Cl[in].Val[147] = U[0][0]*I0[147] + U[4][0]*I1[147];
V_Cl[in].Val[148] = U[0][0]*I0[148] + U[4][0]*I1[148];
V_Cl[in].Val[149] = U[0][0]*I0[149] + U[4][0]*I1[149];
V_Cl[in].Val[150] = U[0][1]*I0[90] + U[4][1]*I1[90]
           + (threeoo2z)*(I2[45] - (poz)*I3[45]);
V_Cl[in].Val[151] = U[0][1]*I0[91] + U[4][1]*I1[91]
           + (threeoo2z)*(I2[46] - (poz)*I3[46])
           + (oo2zn)*I4[60];
V_Cl[in].Val[152] = U[0][1]*I0[92] + U[4][1]*I1[92]
           + (threeoo2z)*(I2[47] - (poz)*I3[47]);
V_Cl[in].Val[153] = U[0][1]*I0[93] + U[4][1]*I1[93]
           + (threeoo2z)*(I2[48] - (poz)*I3[48])
           + (twooo2zn)*I4[61];
V_Cl[in].Val[154] = U[0][1]*I0[94] + U[4][1]*I1[94]
           + (threeoo2z)*(I2[49] - (poz)*I3[49])
           + (oo2zn)*I4[62];
V_Cl[in].Val[155] = U[0][1]*I0[95] + U[4][1]*I1[95]
           + (threeoo2z)*(I2[50] - (poz)*I3[50]);
V_Cl[in].Val[156] = U[0][1]*I0[96] + U[4][1]*I1[96]
           + (threeoo2z)*(I2[51] - (poz)*I3[51])
           + (threeoo2zn)*I4[63];
V_Cl[in].Val[157] = U[0][1]*I0[97] + U[4][1]*I1[97]
           + (threeoo2z)*(I2[52] - (poz)*I3[52])
           + (twooo2zn)*I4[64];
V_Cl[in].Val[158] = U[0][1]*I0[98] + U[4][1]*I1[98]
           + (threeoo2z)*(I2[53] - (poz)*I3[53])
           + (oo2zn)*I4[65];
V_Cl[in].Val[159] = U[0][1]*I0[99] + U[4][1]*I1[99]
           + (threeoo2z)*(I2[54] - (poz)*I3[54]);
V_Cl[in].Val[160] = U[0][1]*I0[100] + U[4][1]*I1[100]
           + (threeoo2z)*(I2[55] - (poz)*I3[55])
           + (fouroo2zn)*I4[66];
V_Cl[in].Val[161] = U[0][1]*I0[101] + U[4][1]*I1[101]
           + (threeoo2z)*(I2[56] - (poz)*I3[56])
           + (threeoo2zn)*I4[67];
V_Cl[in].Val[162] = U[0][1]*I0[102] + U[4][1]*I1[102]
           + (threeoo2z)*(I2[57] - (poz)*I3[57])
           + (twooo2zn)*I4[68];
V_Cl[in].Val[163] = U[0][1]*I0[103] + U[4][1]*I1[103]
           + (threeoo2z)*(I2[58] - (poz)*I3[58])
           + (oo2zn)*I4[69];
V_Cl[in].Val[164] = U[0][1]*I0[104] + U[4][1]*I1[104]
           + (threeoo2z)*(I2[59] - (poz)*I3[59]);
V_Cl[in].Val[165] = U[0][1]*I0[105] + U[4][1]*I1[105]
           + (twooo2z)*(I2[60] - (poz)*I3[60]);
V_Cl[in].Val[166] = U[0][1]*I0[106] + U[4][1]*I1[106]
           + (twooo2z)*(I2[61] - (poz)*I3[61])
           + (oo2zn)*I4[70];
V_Cl[in].Val[167] = U[0][1]*I0[107] + U[4][1]*I1[107]
           + (twooo2z)*(I2[62] - (poz)*I3[62]);
V_Cl[in].Val[168] = U[0][1]*I0[108] + U[4][1]*I1[108]
           + (twooo2z)*(I2[63] - (poz)*I3[63])
           + (twooo2zn)*I4[71];
V_Cl[in].Val[169] = U[0][1]*I0[109] + U[4][1]*I1[109]
           + (twooo2z)*(I2[64] - (poz)*I3[64])
           + (oo2zn)*I4[72];
V_Cl[in].Val[170] = U[0][1]*I0[110] + U[4][1]*I1[110]
           + (twooo2z)*(I2[65] - (poz)*I3[65]);
V_Cl[in].Val[171] = U[0][1]*I0[111] + U[4][1]*I1[111]
           + (twooo2z)*(I2[66] - (poz)*I3[66])
           + (threeoo2zn)*I4[73];
V_Cl[in].Val[172] = U[0][1]*I0[112] + U[4][1]*I1[112]
           + (twooo2z)*(I2[67] - (poz)*I3[67])
           + (twooo2zn)*I4[74];
V_Cl[in].Val[173] = U[0][1]*I0[113] + U[4][1]*I1[113]
           + (twooo2z)*(I2[68] - (poz)*I3[68])
           + (oo2zn)*I4[75];
V_Cl[in].Val[174] = U[0][1]*I0[114] + U[4][1]*I1[114]
           + (twooo2z)*(I2[69] - (poz)*I3[69]);
V_Cl[in].Val[175] = U[0][1]*I0[115] + U[4][1]*I1[115]
           + (twooo2z)*(I2[70] - (poz)*I3[70])
           + (fouroo2zn)*I4[76];
V_Cl[in].Val[176] = U[0][1]*I0[116] + U[4][1]*I1[116]
           + (twooo2z)*(I2[71] - (poz)*I3[71])
           + (threeoo2zn)*I4[77];
V_Cl[in].Val[177] = U[0][1]*I0[117] + U[4][1]*I1[117]
           + (twooo2z)*(I2[72] - (poz)*I3[72])
           + (twooo2zn)*I4[78];
V_Cl[in].Val[178] = U[0][1]*I0[118] + U[4][1]*I1[118]
           + (twooo2z)*(I2[73] - (poz)*I3[73])
           + (oo2zn)*I4[79];
V_Cl[in].Val[179] = U[0][1]*I0[119] + U[4][1]*I1[119]
           + (twooo2z)*(I2[74] - (poz)*I3[74]);
V_Cl[in].Val[180] = U[0][1]*I0[120] + U[4][1]*I1[120]
           + (oo2z)*(I2[75] - (poz)*I3[75]);
V_Cl[in].Val[181] = U[0][1]*I0[121] + U[4][1]*I1[121]
           + (oo2z)*(I2[76] - (poz)*I3[76])
           + (oo2zn)*I4[80];
V_Cl[in].Val[182] = U[0][1]*I0[122] + U[4][1]*I1[122]
           + (oo2z)*(I2[77] - (poz)*I3[77]);
V_Cl[in].Val[183] = U[0][1]*I0[123] + U[4][1]*I1[123]
           + (oo2z)*(I2[78] - (poz)*I3[78])
           + (twooo2zn)*I4[81];
V_Cl[in].Val[184] = U[0][1]*I0[124] + U[4][1]*I1[124]
           + (oo2z)*(I2[79] - (poz)*I3[79])
           + (oo2zn)*I4[82];
V_Cl[in].Val[185] = U[0][1]*I0[125] + U[4][1]*I1[125]
           + (oo2z)*(I2[80] - (poz)*I3[80]);
V_Cl[in].Val[186] = U[0][1]*I0[126] + U[4][1]*I1[126]
           + (oo2z)*(I2[81] - (poz)*I3[81])
           + (threeoo2zn)*I4[83];
V_Cl[in].Val[187] = U[0][1]*I0[127] + U[4][1]*I1[127]
           + (oo2z)*(I2[82] - (poz)*I3[82])
           + (twooo2zn)*I4[84];
V_Cl[in].Val[188] = U[0][1]*I0[128] + U[4][1]*I1[128]
           + (oo2z)*(I2[83] - (poz)*I3[83])
           + (oo2zn)*I4[85];
V_Cl[in].Val[189] = U[0][1]*I0[129] + U[4][1]*I1[129]
           + (oo2z)*(I2[84] - (poz)*I3[84]);
V_Cl[in].Val[190] = U[0][1]*I0[130] + U[4][1]*I1[130]
           + (oo2z)*(I2[85] - (poz)*I3[85])
           + (fouroo2zn)*I4[86];
V_Cl[in].Val[191] = U[0][1]*I0[131] + U[4][1]*I1[131]
           + (oo2z)*(I2[86] - (poz)*I3[86])
           + (threeoo2zn)*I4[87];
V_Cl[in].Val[192] = U[0][1]*I0[132] + U[4][1]*I1[132]
           + (oo2z)*(I2[87] - (poz)*I3[87])
           + (twooo2zn)*I4[88];
V_Cl[in].Val[193] = U[0][1]*I0[133] + U[4][1]*I1[133]
           + (oo2z)*(I2[88] - (poz)*I3[88])
           + (oo2zn)*I4[89];
V_Cl[in].Val[194] = U[0][1]*I0[134] + U[4][1]*I1[134]
           + (oo2z)*(I2[89] - (poz)*I3[89]);
V_Cl[in].Val[195] = U[0][1]*I0[135] + U[4][1]*I1[135];
V_Cl[in].Val[196] = U[0][1]*I0[136] + U[4][1]*I1[136]
           + (oo2zn)*I4[90];
V_Cl[in].Val[197] = U[0][1]*I0[137] + U[4][1]*I1[137];
V_Cl[in].Val[198] = U[0][1]*I0[138] + U[4][1]*I1[138]
           + (twooo2zn)*I4[91];
V_Cl[in].Val[199] = U[0][1]*I0[139] + U[4][1]*I1[139]
           + (oo2zn)*I4[92];
V_Cl[in].Val[200] = U[0][1]*I0[140] + U[4][1]*I1[140];
V_Cl[in].Val[201] = U[0][1]*I0[141] + U[4][1]*I1[141]
           + (threeoo2zn)*I4[93];
V_Cl[in].Val[202] = U[0][1]*I0[142] + U[4][1]*I1[142]
           + (twooo2zn)*I4[94];
V_Cl[in].Val[203] = U[0][1]*I0[143] + U[4][1]*I1[143]
           + (oo2zn)*I4[95];
V_Cl[in].Val[204] = U[0][1]*I0[144] + U[4][1]*I1[144];
V_Cl[in].Val[205] = U[0][1]*I0[145] + U[4][1]*I1[145]
           + (fouroo2zn)*I4[96];
V_Cl[in].Val[206] = U[0][1]*I0[146] + U[4][1]*I1[146]
           + (threeoo2zn)*I4[97];
V_Cl[in].Val[207] = U[0][1]*I0[147] + U[4][1]*I1[147]
           + (twooo2zn)*I4[98];
V_Cl[in].Val[208] = U[0][1]*I0[148] + U[4][1]*I1[148]
           + (oo2zn)*I4[99];
V_Cl[in].Val[209] = U[0][1]*I0[149] + U[4][1]*I1[149];
V_Cl[in].Val[210] = U[0][2]*I0[135] + U[4][2]*I1[135]
           + (threeoo2z)*(I2[75] - (poz)*I3[75]);
V_Cl[in].Val[211] = U[0][2]*I0[136] + U[4][2]*I1[136]
           + (threeoo2z)*(I2[76] - (poz)*I3[76]);
V_Cl[in].Val[212] = U[0][2]*I0[137] + U[4][2]*I1[137]
           + (threeoo2z)*(I2[77] - (poz)*I3[77])
           + (oo2zn)*I4[90];
V_Cl[in].Val[213] = U[0][2]*I0[138] + U[4][2]*I1[138]
           + (threeoo2z)*(I2[78] - (poz)*I3[78]);
V_Cl[in].Val[214] = U[0][2]*I0[139] + U[4][2]*I1[139]
           + (threeoo2z)*(I2[79] - (poz)*I3[79])
           + (oo2zn)*I4[91];
V_Cl[in].Val[215] = U[0][2]*I0[140] + U[4][2]*I1[140]
           + (threeoo2z)*(I2[80] - (poz)*I3[80])
           + (twooo2zn)*I4[92];
V_Cl[in].Val[216] = U[0][2]*I0[141] + U[4][2]*I1[141]
           + (threeoo2z)*(I2[81] - (poz)*I3[81]);
V_Cl[in].Val[217] = U[0][2]*I0[142] + U[4][2]*I1[142]
           + (threeoo2z)*(I2[82] - (poz)*I3[82])
           + (oo2zn)*I4[93];
V_Cl[in].Val[218] = U[0][2]*I0[143] + U[4][2]*I1[143]
           + (threeoo2z)*(I2[83] - (poz)*I3[83])
           + (twooo2zn)*I4[94];
V_Cl[in].Val[219] = U[0][2]*I0[144] + U[4][2]*I1[144]
           + (threeoo2z)*(I2[84] - (poz)*I3[84])
           + (threeoo2zn)*I4[95];
V_Cl[in].Val[220] = U[0][2]*I0[145] + U[4][2]*I1[145]
           + (threeoo2z)*(I2[85] - (poz)*I3[85]);
V_Cl[in].Val[221] = U[0][2]*I0[146] + U[4][2]*I1[146]
           + (threeoo2z)*(I2[86] - (poz)*I3[86])
           + (oo2zn)*I4[96];
V_Cl[in].Val[222] = U[0][2]*I0[147] + U[4][2]*I1[147]
           + (threeoo2z)*(I2[87] - (poz)*I3[87])
           + (twooo2zn)*I4[97];
V_Cl[in].Val[223] = U[0][2]*I0[148] + U[4][2]*I1[148]
           + (threeoo2z)*(I2[88] - (poz)*I3[88])
           + (threeoo2zn)*I4[98];
V_Cl[in].Val[224] = U[0][2]*I0[149] + U[4][2]*I1[149]
           + (threeoo2z)*(I2[89] - (poz)*I3[89])
           + (fouroo2zn)*I4[99];
  V_done[in] = 1;
  return V_Cl[in].Val;
}
