
#include <chemistry/qc/cints/int2jf.h>

static int io[] = {1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153};

void
TwoBodyIntJF::Init_VRR(iclass *H_Cl, iclass *V_Cl, int in)
{
  if (H_Cl[in].am[0]+H_Cl[in].am[1]+H_Cl[in].am[2]+H_Cl[in].am[3]==0) {
    H_Cl[in].Val[0] += F[H_Cl[in].m];
    return;
  }

  // these will be the top level recursive calls
  switch(H_Cl[in].type) {
    case(1):
      Top_VRR_build(H_Cl, V_Cl, in);
      break;
    case(2):
      top_build_p000(&H_Cl[in], V_Cl);
      break;
    case(3):
      top_build_00p0(&H_Cl[in], V_Cl);
      break;
    case(4):
      top_build_d000(&H_Cl[in], V_Cl);
      break;
    case(5):
      top_build_00d0(&H_Cl[in], V_Cl);
      break;
    case(6):
      top_build_p0p0(&H_Cl[in], V_Cl);
      break;
    case(7):
      top_build_f000(&H_Cl[in], V_Cl);
      break;
    case(8):
      top_build_00f0(&H_Cl[in], V_Cl);
      break;
    case(9):
      top_build_d0p0(&H_Cl[in], V_Cl);
      break;
    case(10):
      top_build_p0d0(&H_Cl[in], V_Cl);
      break;
    case(11):
      top_build_f0p0(&H_Cl[in], V_Cl);
      break;
    case(12):
      top_build_p0f0(&H_Cl[in], V_Cl);
      break;
    case(13):
      top_build_d0d0(&H_Cl[in], V_Cl);
      break;
    case(14):
      top_build_f0d0(&H_Cl[in], V_Cl);
      break;
    case(15):
      top_build_d0f0(&H_Cl[in], V_Cl);
      break;
    case(16):
      top_build_f0f0(&H_Cl[in], V_Cl);
      break;
    case(17):
      top_build_g000(&H_Cl[in], V_Cl);
      break;
    case(18):
      top_build_00g0(&H_Cl[in], V_Cl);
      break;
    case(19):
      top_build_g0p0(&H_Cl[in], V_Cl);
      break;
    case(20):
      top_build_p0g0(&H_Cl[in], V_Cl);
      break;
    case(21):
      top_build_d0g0(&H_Cl[in], V_Cl);
      break;
    case(22):
      top_build_g0d0(&H_Cl[in], V_Cl);
      break;
    case(23):
      top_build_g0f0(&H_Cl[in], V_Cl);
      break;
    case(24):
      top_build_f0g0(&H_Cl[in], V_Cl);
      break;
    case(25):
      top_build_g0g0(&H_Cl[in], V_Cl);
      break;
    default:
      fprintf(stderr, "something is very wrong, and"
              " it is concentrated in Top_VRR_build()\n");
      abort();
  }
}

void
TwoBodyIntJF::Top_VRR_build(iclass *H_Cl, iclass *V_Cl, int in)
{
  double *hclinv = H_Cl[in].Val;
  double *hclinvi = hclinv;
  int L0[3];
  int L2[3];
  int k;
  double K1, K2;
  int i1, j1, k1;
  int i2, j2, k2;
  int ii, jj, kk, ll;
  int a=0;
  int b=0;
  int c=0;
  int d=0;
  int i, j;
  double Ub[3];
  double Uk[3];
  int t1, t2, t3, t4;
  int T1, T2, T3, T4;
  int T2max;
  int am[4], ampass[4];
  double *I[5]; /* these for getting children easier */

  if (H_Cl[in].am[0]+H_Cl[in].am[1]+H_Cl[in].am[2]+H_Cl[in].am[3]==0) {
    H_Cl[in].Val[0] += F[H_Cl[in].m];
    return;
  }

  /* these will be the recursive calls */  
  for (i=0;i<5;i++) {
    if(H_Cl[in].operands[i]){
      switch(V_Cl[H_Cl[in].operands[i]].type) {
        case(1):
          I[i] = VRR_Build(V_Cl, H_Cl[in].operands[i]);
          break;
        case(2):
          I[i] = build_p000(V_Cl, H_Cl[in].operands[i]);
          break;
        case(3):
          I[i] = build_00p0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(4):
          I[i] = build_d000(V_Cl, H_Cl[in].operands[i]);
          break;
        case(5):
          I[i] = build_00d0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(6):
          I[i] = build_p0p0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(7):
          I[i] = build_f000(V_Cl, H_Cl[in].operands[i]);
          break;
        case(8):
          I[i] = build_00f0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(9):
          I[i] = build_d0p0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(10):
          I[i] = build_p0d0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(11):
          I[i] = build_f0p0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(12):
          I[i] = build_p0f0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(13):
          I[i] = build_d0d0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(14):
          I[i] = build_f0d0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(15):
          I[i] = build_d0f0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(16):
          I[i] = build_f0f0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(17):
          I[i] = build_g000(V_Cl, H_Cl[in].operands[i]);
          break;
        case(18):
          I[i] = build_00g0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(19):
          I[i] = build_g0p0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(20):
          I[i] = build_p0g0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(21):
          I[i] = build_d0g0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(22):
          I[i] = build_g0d0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(23):
          I[i] = build_g0f0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(24):
          I[i] = build_f0g0(V_Cl, H_Cl[in].operands[i]);
          break;
        case(25):
          I[i] = build_g0g0(V_Cl, H_Cl[in].operands[i]);
          break;
        default:
          fprintf(stderr, "something is very wrong, and"
                  " it is concentrated in Top_VRR_build()\n");
          abort();
        }
      }
    else {
       if(i==0||i==2) I[i] = &F[H_Cl[in].m];
       if(i==1||i==3||i==4) I[i] = &F[H_Cl[in].m+1];
       }
    }

  /* zero L - this /is/ necessary... */
  bzero(L0,3*sizeof(int));
  bzero(L2,3*sizeof(int));

  /* set up am vector */
  for(i=0;i<4;i++){
    ampass[i] = H_Cl[in].am[i];
    am[i] = H_Cl[in].am[i];
    }

  /* b for picking which to decrement */
  b=0;
  if(H_Cl[in].am[0]==0) b=2;

  if(b==2){
    T2max = (am[2])*(am[2]+1)/2;
    T1 = T2 = T3 = T4 = 0;

/* set up constants based on which center b refers to */
    Ub[0] = U[2][0];
    Ub[1] = U[2][1];
    Ub[2] = U[2][2];
    Uk[0] = U[5][0];
    Uk[1] = U[5][1];
    Uk[2] = U[5][2];

    K1 = oo2n;
    K2 = pon;

/* build L */
    t1 = 0;
    for(ii = 0; ii <= am[0]; ii++){
      L0[0] = am[0] - ii;
      for(jj = 0; jj <= ii; jj++){
        L0[1] = ii - jj;
        L0[2] = jj ;

        for(kk = 0; kk <= am[2]; kk++){
          L2[0] = am[2] - kk;
          for(ll = 0; ll <= kk; ll++){
            L2[1] = kk - ll;
            L2[2] = ll ;

          /* d to find xyz */
            if(L2[2]) d=2;
            if(L2[1]) d=1;
            if(L2[0]) d=0;

            if(T2==T2max){
              L2[d] = L2[d] - 1;
              ampass[2] = ampass[2] - 1;
              a = c = 0;
              if(ampass[0]){
                i = ampass[0]-L0[0];
                a = i + ioff(i) - L0[1];
                }
              if(ampass[2]){
                i = ampass[2]-L2[0];
                c = i + ioff(i) - L2[1];
                }
              T2 = a*io[ampass[2]]+c;
              L2[d] = L2[d] - 1;
              ampass[2] = ampass[b] - 1;
              a = c = 0;
              if(ampass[0]){
                i = ampass[0]-L0[0];
                a = i + ioff(i) - L0[1];
               }
              if(ampass[2]){
                i = ampass[2]-L2[0];
                c = i + ioff(i) - L2[1];
                }
              T3 = a*io[ampass[2]]+c;
              L2[d] = L2[d] + 1;
              ampass[2] = ampass[2] + 1;
              L0[0] = L0[0] - 1;
              ampass[0] = ampass[0] - 1;
              a = c = 0;
              if(ampass[0]){
                i = ampass[0]-L0[0];
                a = i + ioff(i) - L0[1];
                }
              if(ampass[2]){
                i = ampass[2]-L2[0];
                c = i + ioff(i) - L2[1];
                }
              T4 = a*io[ampass[2]]+c;
              L0[0] = L0[0] + 1;
              ampass[0] = ampass[0] + 1;
              L2[d] = L2[d] + 1;
              ampass[2] = ampass[2] + 1;
              }
  
            if(fabs(Ub[d])> 1.0e-15)
              *hclinvi += Ub[d]*I[0][T2];
            if(fabs(Uk[d])> 1.0e-15)
              *hclinvi += Uk[d]*I[1][T2];
  
            if(L2[d]>1){
              *hclinvi += K1*(L2[d]-1)*
                          (I[2][T3] - K2*I[3][T3]);
              T3++;
              }
            if(L0[d]){
              *hclinvi += (L0[d])*oo2zn*I[4][T4];
              T4++;
              }

            T2++;
            hclinvi++;
            }
          }
        }
      }
    return;
    }
  else if(b==0){
/* for function-less indexing */
    T2max = (am[0]*(am[0]+1)*(am[2]+1)*(am[2]+2))/4;
    T1 = T2 = T3 = T4 = 0;

/* set up constants based on which center b refers to */
    k=4;
    K1=oo2z;
    K2=poz;
    Ub[0] = U[0][0];
    Ub[1] = U[0][1];
    Ub[2] = U[0][2];
    Uk[0] = U[4][0];
    Uk[1] = U[4][1];
    Uk[2] = U[4][2];

/* build L */
    t1 = 0;
    for(ii = 0; ii <= am[0]; ii++){
      L0[0] = am[0] - ii;
      for(jj = 0; jj <= ii; jj++){
        L0[1] = ii - jj;
        L0[2] = jj ;
  
      /* d to find xyz */
        if(L0[0]) d=0;
        else if(L0[1]) d=1;
        else if(L0[2]) d=2;

        for(kk = 0; kk <= am[2]; kk++){
          L2[0] = am[2] - kk;
          for(ll = 0; ll <= kk; ll++){
            L2[1] = kk - ll;
            L2[2] = ll ;

            if(T2==T2max){
              L0[d] = L0[d] - 1;
              ampass[0] = ampass[0] - 1;
              a = c = 0;
              if(ampass[0]){
                i = ampass[0]-L0[0];
                a = i + ioff(i) - L0[1];
                }
              if(ampass[2]){
                i = ampass[2]-L2[0];
                c = i + ioff(i) - L2[1];
                }
              T2 = a*io[ampass[2]]+c;
              L0[d] = L0[d] - 1;
              ampass[0] = ampass[0] - 1;
              a = c = 0;
              if(ampass[0]){
                i = ampass[0]-L0[0];
                a = i + ioff(i) - L0[1];
                }
              if(ampass[2]){
                i = ampass[2]-L2[0];
                c = i + ioff(i) - L2[1];
                }
              T3 = a*io[ampass[2]]+c;
              L0[d] = L0[d] + 1;
              ampass[0] = ampass[0] + 1;
              L2[0] = L2[0] - 1;
              ampass[2] = ampass[2] - 1;
              a = c = 0;
              if(ampass[0]){
                i = ampass[0]-L0[0];
                a = i + ioff(i) - L0[1];
                }
              if(ampass[2]){
                i = ampass[2]-L2[0];
                c = i + ioff(i) - L2[1];
                }
              T4 = a*io[ampass[2]]+c;
              L2[0] = L2[0] + 1;
              ampass[2] = ampass[2] + 1;
              L0[d] = L0[d] + 1;
              ampass[0] = ampass[0] + 1;
              }

            if(fabs(Ub[d])> 1.0e-15)
              *hclinvi += Ub[d]*I[0][T2];
            if(fabs(Uk[d])> 1.0e-15)
              *hclinvi += Uk[d]*I[1][T2];

            if(L0[d]>1){
              *hclinvi += K1*(L0[d]-1)*
                        (I[2][T3] - K2*I[3][T3]);
              T3++;
              }
            if(L2[d]){
              *hclinvi += (L2[d])*oo2zn*I[4][T4];
              T4++;
              }
          
            T2++;
            hclinvi++;
            }
          }
        }
      }
    }
}

double *
TwoBodyIntJF::VRR_Build(iclass *V_Cl, int in)
{
  struct am_str L;
  int i1, j1, k1;
  int i2, j2, k2;
  int ii, jj, kk, ll;
  int k, index;
  double K1, K2;
  int a=0;
  int b=0;
  int c=0;
  int d=0;
  int i, j;
  int count = 0;
  double alpha[4];
  double t;
  int t1, t2, t3;
  int T1, T2, T3, T4;
  int T2max;
  int am[4], ampass[4];
  double *I[5]; /* these for getting children easier */

  if(V_done[in]) return V_Cl[in].Val;

/* print out what is needed so I can make special routines for
   those cases most commonly called */
  /*printf("    VRR_Build called with am = %d %d %d %d\n",
         V_Cl[in].am[0], V_Cl[in].am[1], V_Cl[in].am[2], V_Cl[in].am[3]); */

/* these will be the recursive calls */  
  for(i=0;i<5;i++){
    if(V_Cl[in].operands[i]){
      switch(V_Cl[V_Cl[in].operands[i]].type) {
        case(1):
          I[i] = VRR_Build(V_Cl, V_Cl[in].operands[i]);
          break;
        case(2):
          I[i] = build_p000(V_Cl, V_Cl[in].operands[i]);
          break;
        case(3):
          I[i] = build_00p0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(4):
          I[i] = build_d000(V_Cl, V_Cl[in].operands[i]);
          break;
        case(5):
          I[i] = build_00d0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(6):
          I[i] = build_p0p0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(7):
          I[i] = build_f000(V_Cl, V_Cl[in].operands[i]);
          break;
        case(8):
          I[i] = build_00f0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(9):
          I[i] = build_d0p0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(10):
          I[i] = build_p0d0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(11):
          I[i] = build_f0p0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(12):
          I[i] = build_p0f0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(13):
          I[i] = build_d0d0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(14):
          I[i] = build_f0d0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(15):
          I[i] = build_d0f0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(16):
          I[i] = build_f0f0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(17):
          I[i] = build_g000(V_Cl, V_Cl[in].operands[i]);
          break;
        case(18):
          I[i] = build_00g0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(19):
          I[i] = build_g0p0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(20):
          I[i] = build_p0g0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(21):
          I[i] = build_d0g0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(22):
          I[i] = build_g0d0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(23):
          I[i] = build_g0f0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(24):
          I[i] = build_f0g0(V_Cl, V_Cl[in].operands[i]);
          break;
        case(25):
          I[i] = build_g0g0(V_Cl, V_Cl[in].operands[i]);
          break;
        default:
          fprintf(stderr, "something is very wrong, and"
                  " it is concentrated in VRR_Build()\n");
          abort();
        }
      }
    else {
       if(i==0||i==2) I[i] = &F[V_Cl[in].m];
       if(i==1||i==3||i==4) I[i] = &F[V_Cl[in].m+1];
       }
    }

  /* zero L - this /is/ necessary... */
  bzero(L.index,12*sizeof(int));

  /* set up am vector */
  for(i=0;i<4;i++){
    ampass[i] = V_Cl[in].am[i];
    am[i] = V_Cl[in].am[i];
    }

  /* b for picking which to decrement */
  b=0;
  if(V_Cl[in].am[0]==0) b=2;

  /* for function-less indexing */
  if(b==0) T2max = (am[0]*(am[0]+1)*(am[2]+1)*(am[2]+2))/4;
  if(b==2) T2max = (am[2])*(am[2]+1)/2;
  T1 = T2 = T3 = T4 = 0;

/* set up constants based on which center b refers to */
  k=5;
  K1 = oo2n;
  K2 = pon;
  if(b<2) {k=4;K1=oo2z;K2=poz;}

/* build L */
  for(ii = 0; ii <= am[0]; ii++){
    L.index[0][0] = am[0] - ii;
    for(jj = 0; jj <= ii; jj++){
      L.index[0][1] = ii - jj;
      L.index[0][2] = jj ;

      for(kk = 0; kk <= am[2]; kk++){
        L.index[2][0] = am[2] - kk;
        for(ll = 0; ll <= kk; ll++){
          L.index[2][1] = kk - ll;
          L.index[2][2] = ll ;

          /* have L defining element of V_Cl[in] */
          V_Cl[in].Val[T1] = 0;
          
          /* d to find xyz */
          if(L.index[b][2]) d=2;
          if(L.index[b][1]) d=1;
          if(L.index[b][0]) d=0;

          if(T2==T2max){
            L.index[b][d] = L.index[b][d] - 1;
            ampass[b] = ampass[b] - 1;
            a = c = 0;
            if(ampass[0]){
              i = ampass[0]-L.index[0][0];
              a = i + ioff(i) - L.index[0][1];
              }
            if(ampass[2]){
              i = ampass[2]-L.index[2][0];
              c = i + ioff(i) - L.index[2][1];
              }
            T2 = a*io[ampass[2]]+c;
            L.index[b][d] = L.index[b][d] - 1;
            ampass[b] = ampass[b] - 1;
            a = c = 0;
            if(ampass[0]){
              i = ampass[0]-L.index[0][0];
              a = i + ioff(i) - L.index[0][1];
              }
            if(ampass[2]){
              i = ampass[2]-L.index[2][0];
              c = i + ioff(i) - L.index[2][1];
              }
            T3 = a*io[ampass[2]]+c;
            L.index[b][d] = L.index[b][d] + 1;
            ampass[b] = ampass[b] + 1;
            L.index[b^2][0] = L.index[b^2][0] - 1;
            ampass[b^2] = ampass[b^2] - 1;
            a = c = 0;
            if(ampass[0]){
              i = ampass[0]-L.index[0][0];
              a = i + ioff(i) - L.index[0][1];
              }
            if(ampass[2]){
              i = ampass[2]-L.index[2][0];
              c = i + ioff(i) - L.index[2][1];
              }
            T4 = a*io[ampass[2]]+c;
            L.index[b^2][0] = L.index[b^2][0] + 1;
            ampass[b^2] = ampass[b^2] + 1;
            L.index[b][d] = L.index[b][d] + 1;
            ampass[b] = ampass[b] + 1;
            }

          /* build that one from pieces already calculated */
          if(U[b][d]> 1.0e-15 || U[b][d]< -1.0e-15)
            V_Cl[in].Val[T1] += U[b][d]*I[0][T2];
          if(U[k][d]> 1.0e-15 || U[k][d] < -1.0e-15)
            V_Cl[in].Val[T1] += U[k][d]*I[1][T2];

          if(L.index[b][d]>1){
            V_Cl[in].Val[T1] += K1*(L.index[b][d]-1)*(I[2][T3] - K2*I[3][T3]);
            T3++;
            }
          if(L.index[b^2][d]){
            V_Cl[in].Val[T1] += (L.index[b^2][d])*oo2zn*I[4][T4];
            T4++;
            }

          T2++;
          T1++;
          }
        }
      }
    }
  V_done[in] = 1;
  return V_Cl[in].Val;
}


int
TwoBodyIntJF::Fill_data(iclass *H_Cl, int in)
{
  int index;
  int i, j, k, l, m, n, o, p;
  double K1, K2;
  int ii, jj, kk, ll;
  int am[4];

  for (i=0;i<4;i++){
    am[i] = H_Cl[in].am[i];
  }
  
  index = 0;
  ii = -1;
  for(i = 0; i <= am[0]; i++){
    for(j = 0; j <= i; j++){
      ii++;

      jj = -1;
      for(k = 0; k <= am[1]; k++){
        for(l = 0; l <= k; l++){
          jj++;

          kk = -1;
          for(m = 0; m <= am[2]; m++){
            for(n = 0; n <= m; n++){
              kk++;

              ll = -1;
              for(o = 0; o <= am[3]; o++){
                for(p = 0; p <= o; p++){
                  ll++;

                  if(fabs(H_Cl[in].Val[index])> 1.0e-15) {
                    printf("%d %d %d %d %lf\n",
                           ii,jj,kk,ll,H_Cl[in].Val[index]);
                  }
                  index++;
                  /* end of all the index generation loops */
                }
              }
            }
          }
        }
      }
    }
  }
   return index;
}
  
double *
TwoBodyIntJF::HRR_build(iclass *H_Cl, int in, double AB[3], double CD[3])
{
  if (H_Cl[in].type || H_done[in])
    return H_Cl[in].Val;

  if (H_Cl[in].am[1])
    return HRR_build_on1(H_Cl, in, AB, CD);

  return HRR_build_on3(H_Cl, in, AB, CD);
}

double *
TwoBodyIntJF::HRR_build_on3(iclass *H_Cl, int in, double AB[3], double CD[3])
{
  double *hclinv = H_Cl[in].Val;
  double *hclinvi=hclinv;
  int index;
  register char i, j, k, l, m, n, o, p;
  int t1, t2;
  int am[4];
  double *I0, *I1; /* these for getting children easier */
  int c11, c12, c13a, c13b, c14;
  int ioa2, ioa3a, ioa3b, ioa4;
  int h1, h2;
  double C[3];
  int xyz;

  I0 = HRR_build(H_Cl, H_Cl[in].operands[0], AB, CD);
  I1 = HRR_build(H_Cl, H_Cl[in].operands[1], AB, CD);

  H_done[in] = 1;
  memcpy(am,H_Cl[in].am,sizeof(int)*4);

  if (!am[3])
    exit(0);

  C[0] = CD[0];
  C[1] = CD[1];
  C[2] = CD[2];

  ioa2 = io[am[1]];
  ioa3a = io[am[2]+1];
  ioa3b = io[am[2]];
  ioa4 = io[am[3]-1];

  index = 0;
  for (i = 0; i <= am[0]; i++) {
    for (j = 0; j <= i; j++) {
      c11 = io[i] - i + j - 1;
      t1 =  c11*ioa2;

      for (k = 0; k <= am[1]; k++) {
        for (l = 0; l <= k; l++) {
          c12 = io[k] - k + l - 1;
          t2 = t1+c12;

          for (m = 0; m <= am[2]; m++) {
            for (n = 0; n <= m; n++) {

              for (o = 0; o <= am[3]; o++) {
                for (p = 0; p <= o; p++) {

                  /* assuming a==3 */

                  c13b = io[m] - m + n - 1;  /* not incremented */
                  if(am[3]-o){ /* x build */
                    c13a = c13b;             /* incremented x */
                    c14 = io[o] - o + p - 1;  /* decremented x */
                    xyz = 0;
                    }
                  else if(o-p){   /* y build */
                    c13a = io[m+1] - m + n - 2;
                    c14 = io[o-1] - o + p;
                    xyz = 1;
                    }
                  else {          /* z build */
                    c13a = io[m+1] - m + n - 1;
                    c14 = io[o-1] - o + p - 1;
                    xyz = 2;
                    }

                  h1 = (t2*ioa3a+c13a)*ioa4+c14;
                  h2 = (t2*ioa3b+c13b)*ioa4+c14;

                  *hclinvi = I0[h1] + C[xyz] * I1[h2];
                  hclinvi++;
                }
              }
            }
          }
        }
      }
    }
  }

  return hclinv;
}


double *
TwoBodyIntJF::HRR_build_on1(iclass *H_Cl, int in, double AB[3], double CD[3])
{
  double *hclinv = H_Cl[in].Val;
  double *hclinvi=hclinv;
  int index;
  register char i, j, k, l, m, n, o, p;
  int t1, t2, t3, t4;
  int am[4];
  double *I0, *I1; /* these for getting children easier */
  double C[3];
  int c11a, c11b, c12, c13, c14;
  int ioa2, ioa3, ioa4;
  int h1, h2;
  int xyz;


  I0 = HRR_build(H_Cl, H_Cl[in].operands[0], AB, CD);
  I1 = HRR_build(H_Cl, H_Cl[in].operands[1], AB, CD);

  H_done[in] = 1;
  memcpy(am,H_Cl[in].am,sizeof(int)*4);

  if (!am[1])
    exit(0);

  C[0] = AB[0];
  C[1] = AB[1];
  C[2] = AB[2];

  ioa2 = io[am[1]-1];
  ioa3 = io[am[2]];
  ioa4 = io[am[3]];


  index = 0;
  for (i = 0; i <= am[0]; i++) {
    for (j = 0; j <= i; j++) {

      for (k = 0; k <= am[1]; k++) {
        for (l = 0; l <= k; l++) {

          /* assuming a==1 */

          c11b = io[i] - i + j - 1;  /* not incremented */
          if(am[1]-k){ /* x build */
            c11a = c11b;             /* incremented x */
            c12 = io[k] - k + l - 1;  /* decremented x */
            xyz = 0;
            }
          else if(k-l){   /* y build */
            c11a = io[i+1] - i + j - 2;
            c12 = io[k-1] - k + l;
            xyz = 1;
            }
          else {          /* z build */
            c11a = io[i+1] - i + j - 1;
            c12 = io[k-1] - k + l - 1;
            xyz = 2;
            }

          t1 = (c11a*ioa2+c12)*ioa3;
          t2 = (c11b*ioa2+c12)*ioa3;

          for (m = 0; m <= am[2]; m++) {
            for (n = 0; n <= m; n++) {
              c13 = io[m] - m + n - 1;
              t3 = (t1+c13)*ioa4;
              t4 = (t2+c13)*ioa4;

              for (o = 0; o <= am[3]; o++) {
                for (p = 0; p <= o; p++) {
                  c14 = io[o] - o + p - 1;

                  h1 = t3+c14;
                  h2 = t4+c14;

                  *hclinvi = I0[h1] + C[xyz] * I1[h2];

                  hclinvi++;
                }
              }
            }
          }
        }
      }
    }
  }

  return hclinv;
}
