
#include <chemistry/qc/cints/cints.h>
#include <chemistry/qc/cints/int2jf.h>

void
TwoBodyIntJF::init()
{
  fjt.reset(20);

  V_done = new char[MAXCLASS];
  memset(V_done,0,sizeof(char)*MAXCLASS);
  
  H_done = new char[MAXCLASS];
  memset(H_done,0,sizeof(char)*MAXCLASS);

  Classes = new iclass[MAXCLASS];
  
  int sz=ioff(MAXAM);
  sz *= sz;
  sz *= sz;
  
  tot_data = new base_eri[sz];
  memset(tot_data,0,sizeof(base_eri)*sz);
  
  sz=0;
  for (int i=0; i < MAXAM*2; i++)
    for (int j=0; j < MAXAM*2; j++)
      sz += ioff(i+1)*ioff(j+1)*MAXAM*8;
  
  dp_use = new double[sz];
  
}

TwoBodyIntJF::TwoBodyIntJF(const RefGaussianBasisSet& gbs) :
  gbs_(gbs),
  V_done(0),
  H_done(0),
  dp_use(0),
  Classes(0),
  tot_data(0)
{
  init();
}

TwoBodyIntJF::~TwoBodyIntJF()
{
  if (V_done) {
    delete[] V_done;
    V_done = 0;
  }
  if (H_done) {
    delete[] H_done;
    H_done = 0;
  }
  if (dp_use) {
    delete[] dp_use;
    dp_use = 0;
  }
  if (Classes) {
    delete[] Classes;
    Classes = 0;
  }
  if (tot_data) {
    delete[] tot_data;
    tot_data = 0;
  }
}

static void
swap(double& i, double& j)
{
  double t=j;
  j=i;
  i=t;
}

static void
swap(shell*& i, shell*& j)
{
  shell* t=j;
  j=i;
  i=t;
}

static void
swap(int& i, int& j)
{
  int t=j;
  j=i;
  i=t;
}

void
TwoBodyIntJF::compute_shell(int sii, int sjj, int skk, int sll, double *buf)
{
  int i,j,k,l;
  
  GaussianBasisSet& gbs = *gbs_.pointer();
  
  GaussianShell* gsi = &gbs(sii);
  GaussianShell* gsj = &gbs(sjj);
  GaussianShell* gsk = &gbs(skk);
  GaussianShell* gsl = &gbs(sll);

  int ati = gbs.shell_to_center(sii);
  int atj = gbs.shell_to_center(sjj);
  int atk = gbs.shell_to_center(skk);
  int atl = gbs.shell_to_center(sll);
  
  int ati_e_atj = (ati==atj);
  int atj_e_atk = (atj==atk);
  int atk_e_atl = (atk==atl);
  
  for (int ci=0; ci < gbs(sii).ncontraction(); ci++) {
    int ami = gsi->am(ci);
    shell shli(ati, sii, ci, gbs);
    
    for (int cj=0; cj < gbs(sjj).ncontraction(); cj++) {
      int amj = gsj->am(cj);
      int amij = ami+amj;
      shell shlj(atj, sjj, cj, gbs);

      for (int ck=0; ck < gbs(skk).ncontraction(); ck++) {
        int amk = gsk->am(ck);
        int amijk = amij + amk;
        shell shlk(atk, skk, ck, gbs);

        for (int cl=0; cl < gbs(sll).ncontraction(); cl++) {
          int aml = gsl->am(cl);
          int amijkl = amijk + aml;
          shell shll(atl, sll, cl, gbs);

          if (amijkl%2 && ati_e_atj && atj_e_atk && atk_e_atl)
            continue;
          
          shell *lsi=&shli, *lsj=&shlj, *lsk=&shlk, *lsl=&shll;
          
          // place in "descending" angular mom-
          // my simple way of optimizing PHG recursion (VRR)
          if (lsi->am < lsj->am)
            swap(lsi,lsj);

          if (lsk->am < lsl->am)
            swap(lsk,lsl);

          if (lsi->am < lsk->am) {
            swap(lsi,lsk);
            swap(lsj,lsl);
          }

          double AB[3], CD[3];
          AB[0] = lsi->p[0] - lsj->p[0];
          AB[1] = lsi->p[1] - lsj->p[1];
          AB[2] = lsi->p[2] - lsj->p[2];
          CD[0] = lsk->p[0] - lsl->p[0];
          CD[1] = lsk->p[1] - lsl->p[1];
          CD[2] = lsk->p[2] - lsl->p[2];
  
          double AB2 = AB[0]*AB[0]+AB[1]*AB[1]+AB[2]*AB[2];
          double CD2 = CD[0]*CD[0]+CD[1]*CD[1]+CD[2]*CD[2];

          int orig_am[4];
          orig_am[0] = lsi->am();
          orig_am[1] = lsj->am();
          orig_am[2] = lsk->am();
          orig_am[3] = lsl->am();
                      
          // usr HRR to get all classes needed to build (ab|cd) - 
          // i.e. all (e0|f0) classes.  push onto stack of classes.

          double *DP = dp_use;
          int n_hrr = 0;
          iclass *H_Cl = Classes;

          // recursively list (and allocate room for) HRR generated
          // intermediates
          List_HRR(H_Cl, orig_am, n_hrr, DP);

          int first_vrr = 1;
          int last_vrr = 1;

          // set V_Cl to after end of hrr generated intermediates
          iclass *V_Cl = Classes + n_hrr;

          // zero flags for HRR generated intermediates
          memset(H_done,0,n_hrr*sizeof(char));

          // now loop over all elements in Classes (e0|f0) -
          // apply vrr to push all those classes onto stack of classes
          int class_it;
          for (class_it = 0; class_it < n_hrr; class_it++)
            if (H_Cl[class_it].type != 0)
              Top_VRR(class_it, H_Cl, V_Cl, last_vrr, DP);

          // contract by primitives out here
          for (int pi = 0; pi < lsi->nprim(); pi++) {
            for (int pj = 0; pj < lsj->nprim(); pj++) {
              for (int pk = 0; pk < lsk->nprim(); pk++) {
                for (int pl = 0; pl < lsl->nprim(); pl++) {
         
                  /* zero flags for vrr generated intermediates */
                  memset(V_done,0,last_vrr*sizeof(char));

                  //pii = pi + shells[si].fprim-1;
                  //pjj = pj + shells[sj].fprim-1;
                  //pkk = pk + shells[sk].fprim-1;
                  //pll = pl + shells[sl].fprim-1;

                  //Shell_init(AB2, CD2, ericent, shells,
                  //           si, sj, sk, sl, pii, pjj, pkk, pll, cgtos);

                  /* call function to begin build by VRR of top classes */
                  for (class_it = 0; class_it < n_hrr; class_it++)
                    if (H_Cl[class_it].type != 0)
                      Init_VRR(H_Cl, V_Cl, class_it);

                }
              }
            }
          } // end getting unique primitive set

#if 0
          /*  call a function to accumulate all the HRR generated classes */
          data = HRR_build(Classes, 0, &AB, &CD);
          num = Fill_data(Classes, 0, tot_data, ioffset, joffset, 
                          koffset, loffset, intfile);

          /* write out shell info */
          /*fprintf(eriout, "shell %2d %2d %2d %2d\n", si, sj, sk, sl);*/
          sz = sizeof(struct tebuf);
#if 0
          for(n=0; n<num; n++){
            if(fabs(tot_data[n].val)>ZERO){
              fprintf(eriout, "%2d %2d %2d %2d %20.13lf\n",
                      tot_data[n].i+ioffset, 
                      tot_data[n].j+joffset, 
                      tot_data[n].k+koffset, 
                      tot_data[n].l+loffset,
                      tot_data[n].val);
            }
            tot_data[n].val = 0.0;
          }
          fflush(eriout);
#endif
#endif

        } /* end if we need to do this shell block */
      }
    }
  }
}


