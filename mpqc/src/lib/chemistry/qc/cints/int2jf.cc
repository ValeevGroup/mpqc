
#include <math.h>

#include <chemistry/qc/cints/int2jf.h>


TwoBodyIntJF::TwoBodyIntJF(const RefGaussianBasisSet& gbs) :
  gbs_(gbs),
  fjt(20)
{
  Classes = new iclass[MAXCLASS];
  V_done = new char[MAXCLASS];
  H_done = new char[MAXCLASS];

  memset(H_done,0,MAXCLASS);
  memset(V_done,0,MAXCLASS);

  // use Size to count how much room needs to be allocated to the stack
  int sz = 0;
  int i;
  for (i=0; i < MAXAM*2; i++)
    for (int j=0; j < MAXAM*2; j++)
      sz += ioff(i+1)*ioff(j+1)*MAXAM*8;

  DP = new double[sz];

}

TwoBodyIntJF::~TwoBodyIntJF()
{
  if (Classes) {
    delete[] Classes;
    Classes=0;
  }

  if (V_done) {
    delete[] V_done;
    V_done=0;
  }

  if (H_done) {
    delete[] H_done;
    H_done=0;
  }

  if (DP) {
    delete[] DP;
    DP=0;
  }
}

static void
swap(shell_stuff*& a, shell_stuff*& b)
{
  shell_stuff *t = a;
  a = b;
  b = t;
}

void
TwoBodyIntJF::compute_shell(int sii, int sjj, int skk, int sll, double *buf)
{
  const double F0[20] = {1.0,  1.0/3.0,  1.0/5.0,  1.0/7.0,  1.0/9.0,
                  1.0/11.0, 1.0/13.0, 1.0/15.0, 1.0/17.0, 1.0/19.0,
                  1.0/21.0, 1.0/23.0, 1.0/25.0, 1.0/27.0, 1.0/29.0,
                  1.0/31.0, 1.0/33.0, 1.0/35.0, 1.0/37.0, 1.0/39.0};

  GaussianBasisSet& gbs = *gbs_.pointer();
  Molecule& mol = *gbs.molecule().pointer();

  GaussianShell& gsi = gbs(sii);
  GaussianShell& gsj = gbs(sjj);
  GaussianShell& gsk = gbs(skk);
  GaussianShell& gsl = gbs(sll);
  
  for (int ci=0; ci < gsi.ncontraction(); ci++) {
    shell_stuff shli(gbs,sii,ci);
    
    for (int cj=0; cj < gsj.ncontraction(); cj++) {
      shell_stuff shlj(gbs,sjj,cj);

      if ((shli.center==shlj.center) && (cj > ci))
        break;

      for (int ck=0; ck < gsk.ncontraction(); ck++) {
        shell_stuff shlk(gbs,skk,ck);
    
        if ((shli.center==shlk.center) && (ck > ci))
          break;

        for (int cl=0; cl < gsl.ncontraction(); cl++) {
          shell_stuff shll(gbs,sll,cl);
    
          if (((shli.center==shlk.center) && (shlj.center==shll.center) &&
               (cl > cj)) ||
              (shlk.center==shll.center) && (cl > ck))
            break;

          // need to decide if we even need to calculate this one... odd am=no?
          int total_am = shli.am + shlj.am + shlk.am + shll.am;

          if ((total_am%2) && (shli.center==shlj.center) &&
              (shlj.center==shlk.center) && (shlk.center==shll.center))
            continue;

          // place in "descending" angular mom-
          // my simple way of optimizing PHG recursion (VRR)
          shell_stuff *shpi = &shli;
          shell_stuff *shpj = &shlj;
          shell_stuff *shpk = &shlk;
          shell_stuff *shpl = &shll;
          
          if (shpi->am < shpj->am)
            swap(shpi,shpj);

          if (shpk->am < shpl->am)
            swap(shpk,shpl);

          if (shpi->am < shpk->am){
            swap(shpi,shpk);
            swap(shpj,shpl);
          }

          int ni = shpi->gs.nfunction(ci);
          int nj = shpj->gs.nfunction(cj);
          int nk = shpk->gs.nfunction(ck);
          int nl = shpl->gs.nfunction(cl);
          int len = ni*nj*nk*nl;
          
          double ab[3], cd[3];
          
          ab[0] = shpi->ac[0] - shpj->ac[0];
          ab[1] = shpi->ac[1] - shpj->ac[1];
          ab[2] = shpi->ac[2] - shpj->ac[2];
          cd[0] = shpk->ac[0] - shpl->ac[0];
          cd[1] = shpk->ac[1] - shpl->ac[1];
          cd[2] = shpk->ac[2] - shpl->ac[2];
  
          double ab2 = ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2];
          double cd2 = cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2];

          int orig_am[5];
          orig_am[0] = shpi->am;
          orig_am[1] = shpj->am;
          orig_am[2] = shpk->am;
          orig_am[3] = shpl->am;
    
          // usr HRR to get all classes needed to build (ab|cd) - 
          // i.e. all (e0|f0) classes.  push onto stack of classes.

          double *dp_use = DP;
          int n_hrr = 0;
          iclass *H_Cl = Classes;

          // recursively list (and allocate room for) HRR generated
          // intermediates
          List_HRR(H_Cl, orig_am, n_hrr, dp_use);

          int first_vrr = 1;
          int last_vrr = 1;

          // set V_Cl to after end of hrr generated intermediates
          iclass *V_Cl = &(Classes[n_hrr]);

          // zero flags for HRR generated intermediates
          memset(H_done,0,n_hrr*sizeof(char));

          // now loop over all elements in Classes (e0|f0) -
          // apply vrr to push all those classes onto stack of classes
          for (int class_it = 0; class_it < n_hrr; class_it++)
            if (H_Cl[class_it].type != 0)
              Top_VRR(class_it, H_Cl, V_Cl, last_vrr, dp_use);

          // contract by primitives out here
          for (int pi = 0; pi < shpi->gs.nprimitive(); pi++) {
            inorm = shpi->coef(pi);
            expi = shpi->gs.exponent(pi);
            
            for (int pj = 0; pj < shpj->gs.nprimitive(); pj++) {
              jnorm = shpj->coef(pj);
              expj = shpj->gs.exponent(pj);
              zeta = expi+expj;
              oo2z = 1.0/zeta;

              double Sovlp1 =
                inorm*jnorm*exp(-expi*expj*ab2*oo2z)*pow(M_PI*oo2z,1.5);

              double Px = (shpi->ac[0] * expi + shpj->ac[0] * expj) * oo2z;
              double Py = (shpi->ac[1] * expi + shpj->ac[1] * expj) * oo2z;
              double Pz = (shpi->ac[2] * expi + shpj->ac[2] * expj) * oo2z;

              oo2z *= 0.5;
              
              U[0][0] = Px - shpi->ac[0];
              U[0][1] = Py - shpi->ac[1];
              U[0][2] = Pz - shpi->ac[2];
              U[1][0] = Px - shpj->ac[0];
              U[1][1] = Py - shpj->ac[1];
              U[1][2] = Pz - shpj->ac[2];

              for (int pk = 0; pk < shpk->gs.nprimitive(); pk++) {
                knorm = shpk->coef(pk);
                expk = shpk->gs.exponent(pk);

                for (int pl = 0; pl < shpl->gs.nprimitive(); pl++) {
                  lnorm = shpl->coef(pl);
                  expl = shpl->gs.exponent(pl);
                  eta = expk + expl;
                  oo2n = 1.0/eta;

                  double Sovlp2 =
                    knorm*lnorm*exp(-expk*expl*cd2*oo2n)*pow(M_PI*oo2n,1.5);

                  double Qx = (shpk->ac[0] * expk + shpl->ac[0] * expl) * oo2n;
                  double Qy = (shpk->ac[1] * expk + shpl->ac[1] * expl) * oo2n;
                  double Qz = (shpk->ac[2] * expk + shpl->ac[2] * expl) * oo2n;

                  oo2n *= 0.5;
                  
                  U[2][0] = Qx - shpk->ac[0];
                  U[2][1] = Qy - shpk->ac[1];
                  U[2][2] = Qz - shpk->ac[2];
                  U[3][0] = Qx - shpl->ac[0];
                  U[3][1] = Qy - shpl->ac[1];
                  U[3][2] = Qz - shpl->ac[2];

                  double PQx = Px - Qx;
                  double PQy = Py - Qy;
                  double PQz = Pz - Qz;

                  double Wx = (Px*zeta + Qx*eta)/(zeta+eta);
                  double Wy = (Py*zeta + Qy*eta)/(zeta+eta);
                  double Wz = (Pz*zeta + Qz*eta)/(zeta+eta);

                  U[4][0] = Wx - Px;
                  U[4][1] = Wy - Py;
                  U[4][2] = Wz - Pz;
                  U[5][0] = Wx - Qx;
                  U[5][1] = Wy - Qy;
                  U[5][2] = Wz - Qz;

                  oo2zn = 0.5/(zeta+eta);
                  rho = zeta*eta / (zeta+eta);
                  poz = rho/zeta;
                  pon = rho/eta;
         
                  double coef1 = 2.0 * sqrt(rho/M_PI) * Sovlp1 * Sovlp2;

                  double pq2 = PQx*PQx + PQy*PQy + PQz*PQz;
                  
                  if (fabs(pq2) < 1.0e-15) {
                    for (int i=0; i <= total_am; i++)
                      F[i] = F0[i]*coef1;
                  } else {
                    fjt.fjt(total_am, rho*pq2);
                    for (int i=0; i <= total_am; i++)
                      F[i] = fjt.table_value(i) * coef1;
                  }
                      
                  // zero flags for vrr generated intermediates
                  memset(V_done, 0, last_vrr*sizeof(char));

                  /* call function to begin build by VRR of top classes */
                  for (int class_it = 0; class_it < n_hrr; class_it++)
                    if (H_Cl[class_it].type != 0)
                      Init_VRR(H_Cl, V_Cl, class_it);
                }
              }
            }
          } // end getting unique primitive set

          // call a function to accumulate all the HRR generated classes
          double *data = HRR_build(Classes, 0, ab, cd);
          //int num = Fill_data(Classes, 0);

#if 0
          int index=0;
          for (int fi = 0; fi < ioff(Classes[0].am[0]+1); fi++) {
            int bfi = fi+shpi->func0;
            
            for (int fj = 0; fj < ioff(Classes[0].am[1]+1); fj++) {
              int bfj = fj+shpj->func0;

              for (int fk = 0; fk < ioff(Classes[0].am[2]+1); fk++) {
                int bfk = fk+shpk->func0;

                for (int fl = 0; fl < ioff(Classes[0].am[3]+1); fl++) {
                  int bfl = fl+shpl->func0;
                  double v = Classes[0].Val[index];

                  if (fabs(v) > 1.0e-12)
                    printf("%5d %5d %5d %5d %20.15f\n",bfi,bfj,bfk,bfl,v);

                  index++;
                }
              }
            }
          }
#endif
        }
      }
    }
  }
}
