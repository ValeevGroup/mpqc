
#include <chemistry/qc/cints/int2jf.h>

static inline int
is_eq(int in[5], struct iclass *cl)
{
  if (in[0] != cl->am[0]) return 0;
  if (in[1] != cl->am[1]) return 0;
  if (in[2] != cl->am[2]) return 0;
  if (in[3] != cl->am[3]) return 0;
  if (in[4] != cl->m) return 0;

  return 1;
}

int
TwoBodyIntJF::List_HRR(iclass *C_stack, int in[5], int &last, double*& dp)
{
  int i, n=last;
  int OK, test;
  int am[4], out1[5], out2[5];

  for (i=0; i < n; i++)
    if (is_eq(in, &C_stack[i]))
      return i;

  // add it to stack and build the operands
  last++;
  int size = 1;

  for (i=0; i < 4; i++) {
    C_stack[n].am[i] = in[i];
    size *= ioff(in[i]+1);
  }

  C_stack[n].Val = dp;
  C_stack[n].done = 0;
  C_stack[n].m = 0;

  memset(C_stack[n].Val, 0, size*sizeof(double));
  dp += size;

  // decide how to decrement
  if (in[1]>0) {
    out1[0] = in[0]+1;
    out1[1] = in[1]-1;
    out1[2] = in[2];
    out1[3] = in[3];
    out1[4] = in[4];
    out2[0] = in[0];
    out2[1] = in[1]-1;
    out2[2] = in[2];
    out2[3] = in[3];
    out2[4] = in[4];
    C_stack[n].type = 0; /* for HRR build */

  } else if (in[3]>0) {
    out1[0] = in[0];
    out1[1] = in[1];
    out1[2] = in[2]+1;
    out1[3] = in[3]-1;
    out1[4] = in[4];
    out2[0] = in[0];
    out2[1] = in[1];
    out2[2] = in[2];
    out2[3] = in[3]-1;
    out2[4] = in[4];
    C_stack[n].type = 0; /* for HRR build */

  } else {
    C_stack[n].type = 1; /* need VRR here - will get specific later */
  }

  if (!(C_stack[n].type)) {
    C_stack[n].operands[0] = List_HRR(C_stack, out1, last, dp);
    C_stack[n].operands[1] = List_HRR(C_stack, out2, last, dp);
  }

  return n;
}

void
TwoBodyIntJF::Top_VRR(int n, iclass *H, iclass *V, int& last, double*& dp)
{
  int i, a, m, tot_am;
  int size;
  int am[4], out0[5], out1[5], out2[5], out3[5], out4[5];
  int foo;

  for (i=0; i < 4; i++) {
    H[n].operands[i] = 0;
    am[i] = H[n].am[i];
    out0[i] = am[i];
    out1[i] = am[i];
    out2[i] = am[i];
    out3[i] = am[i];
    out4[i] = am[i];
  }

  if (i < 4000) foo = 4;

  // this is the spiffy control thing
  if(am[0]==1&&am[2]==0)      H[n].type = 2;  /* (p000) */
  else if(am[0]==0&&am[2]==1) H[n].type = 3;  /* (00p0) */
  else if(am[0]==2&&am[2]==0) H[n].type = 4;  /* (d000) */
  else if(am[0]==0&&am[2]==2) H[n].type = 5;  /* (00d0) */
  else if(am[0]==1&&am[2]==1) H[n].type = 6;  /* (p0p0) */
  else if(am[0]==3&&am[2]==0) H[n].type = 7;  /* (f000) */
  else if(am[0]==0&&am[2]==3) H[n].type = 8;  /* (00f0) */
  else if(am[0]==2&&am[2]==1) H[n].type = 9;  /* (d0p0) */
  else if(am[0]==1&&am[2]==2) H[n].type = 10; /* (p0d0) */
  else if(am[0]==3&&am[2]==1) H[n].type = 11; /* (f0p0) */
  else if(am[0]==1&&am[2]==3) H[n].type = 12; /* (p0f0) */
  else if(am[0]==2&&am[2]==2) H[n].type = 13; /* (d0d0) */
  else if(am[0]==3&&am[2]==2) H[n].type = 14; /* (f0d0) */
  else if(am[0]==2&&am[2]==3) H[n].type = 15; /* (d0f0) */
  else if(am[0]==3&&am[2]==3) H[n].type = 16; /* (f0f0) */
  else if(am[0]==4&&am[2]==0) H[n].type = 17; /* (g000) */
  else if(am[0]==0&&am[2]==4) H[n].type = 18; /* (00g0) */
  else if(am[0]==4&&am[2]==1) H[n].type = 19; /* (g0p0) */
  else if(am[0]==1&&am[2]==4) H[n].type = 20; /* (p0g0) */
  else if(am[0]==2&&am[2]==4) H[n].type = 21; /* (d0g0) */
  else if(am[0]==4&&am[2]==2) H[n].type = 22; /* (g0d0) */
  else if(am[0]==4&&am[2]==3) H[n].type = 23; /* (g0f0) */
  else if(am[0]==3&&am[2]==4) H[n].type = 24; /* (f0g0) */
  else if(am[0]==4&&am[2]==4) H[n].type = 25; /* (g0g0) */
  else H[n].type = 1; /* something else */

  H[n].operands[4] = 0;
  out0[foo] = H[n].m;
  out1[foo] = H[n].m;
  out2[foo] = H[n].m;
  out3[foo] = H[n].m;
  out4[foo] = H[n].m;
 
  // fill in the operands
  a=0;
  if (am[a]==0) a=2;

  if (am[a]>0) {
    out0[a] = am[a]-1;
    out1[a] = am[a]-1;
    out1[foo] = H[n].m+1;
    H[n].operands[0] = List_VRR(out0, V, last, dp);
    H[n].operands[1] = List_VRR(out1, V, last, dp);
    if (am[a]>1) {
      out2[a] = am[a]-2;
      out3[a] = am[a]-2;
      out3[foo] = H[n].m+1;
      H[n].operands[2] = List_VRR(out2, V, last, dp);
      H[n].operands[3] = List_VRR(out3, V, last, dp);
    }
    if (am[a^2]>0) {
      out4[a] = am[a]-1;
      out4[a^2] = am[a^2]-1;
      out4[foo] = H[n].m+1;
      H[n].operands[4] = List_VRR(out4, V, last, dp);
    }
  }
}

int
TwoBodyIntJF::List_VRR(int in[5], iclass *V_stack, int& last, double*& dp)
{
  int i, a, m, tot_am;
  int n = last;
  int foo, size;
  int doit = 0;
  int out0[5], out1[5], out2[5], out3[5], out4[5];

  for (i=1 ; i < last; i++)
    if (is_eq(in, &V_stack[i]))
      return i;

  // should handle if it's a leaf
  if (in[0]+in[2] < 0)
    exit(0);

  if (in[0]+in[2] == 0)
    return 0;

  // add it to stack and build the operands
  last++;
  size = 1;
  for (i=0; i < 4; i++) {
    V_stack[n].am[i] = in[i];
    size *= ioff(in[i]+1);
  }

  V_stack[n].Val = dp;
  V_stack[n].done = 0;
  V_stack[n].m = in[4];

  memset(V_stack[n].Val, 0, size*sizeof(double));
  dp += size;

  // this is the spiffy control thing
  if(in[0]==1&&in[2]==0)      V_stack[n].type = 2;  /* (p000) */
  else if(in[0]==0&&in[2]==1) V_stack[n].type = 3;  /* (00p0) */
  else if(in[0]==2&&in[2]==0) V_stack[n].type = 4;  /* (d000) */
  else if(in[0]==0&&in[2]==2) V_stack[n].type = 5;  /* (00d0) */
  else if(in[0]==1&&in[2]==1) V_stack[n].type = 6;  /* (p0p0) */
  else if(in[0]==3&&in[2]==0) V_stack[n].type = 7;  /* (f000) */
  else if(in[0]==0&&in[2]==3) V_stack[n].type = 8;  /* (00f0) */
  else if(in[0]==2&&in[2]==1) V_stack[n].type = 9;  /* (d0p0) */
  else if(in[0]==1&&in[2]==2) V_stack[n].type = 10; /* (p0d0) */
  else if(in[0]==3&&in[2]==1) V_stack[n].type = 11; /* (f0p0) */
  else if(in[0]==1&&in[2]==3) V_stack[n].type = 12; /* (p0f0) */
  else if(in[0]==2&&in[2]==2) V_stack[n].type = 13; /* (d0d0) */
  else if(in[0]==3&&in[2]==2) V_stack[n].type = 14; /* (f0d0) */
  else if(in[0]==2&&in[2]==3) V_stack[n].type = 15; /* (d0f0) */
  else if(in[0]==3&&in[2]==3) V_stack[n].type = 16; /* (f0f0) */
  else if(in[0]==4&&in[2]==0) V_stack[n].type = 17; /* (g000) */
  else if(in[0]==0&&in[2]==4) V_stack[n].type = 18; /* (00g0) */
  else if(in[0]==4&&in[2]==1) V_stack[n].type = 19; /* (g0p0) */
  else if(in[0]==1&&in[2]==4) V_stack[n].type = 20; /* (p0g0) */
  else if(in[0]==2&&in[2]==4) V_stack[n].type = 21; /* (d0g0) */
  else if(in[0]==4&&in[2]==2) V_stack[n].type = 22; /* (g0d0) */
  else if(in[0]==4&&in[2]==3) V_stack[n].type = 23; /* (g0f0) */
  else if(in[0]==3&&in[2]==4) V_stack[n].type = 24; /* (f0g0) */
  else if(in[0]==4&&in[2]==4) V_stack[n].type = 25; /* (g0g0) */
  else V_stack[n].type = 1; /* something else */

  a=0;
  if (in[0]==0) a=2;
  if (a<4000) foo = 4;

  for (i=0; i < 5; i++) {
    V_stack[n].operands[i] = 0;
    out0[i] = in[i];
    out1[i] = in[i];
    out2[i] = in[i];
    out3[i] = in[i];
    out4[i] = in[i];
  }

  if (in[a]>0) {
    out0[a] = in[a]-1;
    out1[a] = in[a]-1;
    out1[foo] = in[foo]+1;
    V_stack[n].operands[0] = List_VRR(out0, V_stack, last, dp);
    V_stack[n].operands[1] = List_VRR(out1, V_stack, last, dp);
    if (in[a] > 1) {
      out2[a] = in[a]-2;
      out3[a] = in[a]-2;
      out3[foo] = in[foo]+1;
      V_stack[n].operands[2] = List_VRR(out2, V_stack, last, dp);
      V_stack[n].operands[3] = List_VRR(out3, V_stack, last, dp);
    }
    if (in[a^2] > 0) {
      out4[a] = in[a]-1;
      out4[a^2] = in[a^2]-1;
      out4[foo] = in[foo]+1;
      V_stack[n].operands[4] = List_VRR(out4, V_stack, last, dp);
    }
  }

  return n;
}
