
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

  int dum = ioff(MAXAM);
  dum *= dum;
  dum *= dum;
     
  tot_data = new base_eri[dum];
  memset(tot_data,0,sizeof(base_eri)*dum);
 
  // use Size to count how much room needs to be allocated to the stack
  int sz = 0;
  int i;
  for (i=0; i < MAXAM*2; i++)
    for (int j=0; j < MAXAM*2; j++)
      sz += ioff(i+1)*ioff(j+1)*MAXAM*8;

  DP = (double *) malloc(sz*sizeof(double));

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

  if (tot_data) {
    delete[] tot_data;
    tot_data=0;
  }

  if (DP) {
    delete[] DP;
    DP=0;
  }
}

void
TwoBodyIntJF::compute_shell(int sii, int sjj, int skk, int sll, double *buf)
{

  double *dp_use = DP;

  int ioffset, joffset, koffset, loffset;
  int slmax ;
  struct am_str L ;
  struct iclass *Classes;
  struct iclass *H_Cl;
  struct iclass *V_Cl;
  double AB2, CD2;
  struct coordinates P, PA, PB, AB;
  struct coordinates Q, QC, QD, CD;
  int count ;
  int dum;
  struct base_eri *tot_data; /* accum. for contracted fn's */
  double *data;
  int n;
  int num;  /* number of base_eri's returned by shell_eri */
  double inorm, jnorm, knorm, lnorm ;
  struct coordinates ericent[4];
  int total_am;
  int n_hrr = 0;
  int orig_am[4];
  int first_vrr, last_vrr;
  int class_it;

  
  // need to decide if we even need to calculate this one... odd am=no?
  total_am = shells[sii].am+shells[sjj].am+shells[skk].am+shells[sll].am;
  if (!(total_am%2)||
      (shells[sii].center!=shells[sjj].center)||
      (shells[sjj].center!=shells[skk].center)||
      (shells[skk].center!=shells[sll].center)){

    /* si, sj, sk, sl refer to shell numbers here */
    /* place in "descending" angular mom-
       my simple way of optimizing PHG recursion (VRR) */
    si = sii; sj = sjj; sk = skk; sl = sll;
    if (shells[si].am < shells[sj].am){
      dum = si;
      si = sj;
      sj = dum;
    }
    if(shells[sk].am < shells[sl].am){
      dum = sk;
      sk = sl;
      sl = dum;
    }
    if(shells[si].am < shells[sk].am){
      dum = si;
      si = sk;
      sk = dum;
      dum = sj;
      sj = sl;
      sl = dum;
    }
    /* for numbering the AO's */
    ioffset = joffset = koffset = loffset = 0;
    for(i=0;i<si;i++) ioffset += ioff[shells[i].am];
    for(i=0;i<sj;i++) joffset += ioff[shells[i].am];
    for(i=0;i<sk;i++) koffset += ioff[shells[i].am];
    for(i=0;i<sl;i++) loffset += ioff[shells[i].am];


    ericent[0].x = centers[shells[si].center-1].x;
    ericent[0].y = centers[shells[si].center-1].y;
    ericent[0].z = centers[shells[si].center-1].z;
    ericent[0].Z_nuc = centers[shells[si].center-1].Z_nuc;
    ericent[1].x = centers[shells[sj].center-1].x;
    ericent[1].y = centers[shells[sj].center-1].y;
    ericent[1].z = centers[shells[sj].center-1].z;
    ericent[1].Z_nuc = centers[shells[sj].center-1].Z_nuc;
    ericent[2].x = centers[shells[sk].center-1].x;
    ericent[2].y = centers[shells[sk].center-1].y;
    ericent[2].z = centers[shells[sk].center-1].z;
    ericent[2].Z_nuc = centers[shells[sk].center-1].Z_nuc;
    ericent[3].x = centers[shells[sl].center-1].x;
    ericent[3].y = centers[shells[sl].center-1].y;
    ericent[3].z = centers[shells[sl].center-1].z;
    ericent[3].Z_nuc = centers[shells[sl].center-1].Z_nuc;

    AB.x = ericent[0].x-ericent[1].x;
    AB.y = ericent[0].y-ericent[1].y;
    AB.z = ericent[0].z-ericent[1].z;
    CD.x = ericent[2].x-ericent[3].x;
    CD.y = ericent[2].y-ericent[3].y;
    CD.z = ericent[2].z-ericent[3].z;
  
    AB2 = AB.x*AB.x+AB.y*AB.y+AB.z*AB.z;
    CD2 = CD.x*CD.x+CD.y*CD.y+CD.z*CD.z;

    orig_am[0] = shells[si].am-1;
    orig_am[1] = shells[sj].am-1;
    orig_am[2] = shells[sk].am-1;
    orig_am[3] = shells[sl].am-1;
    
    /* usr HRR to get all classes needed to build (ab|cd) - 
       i.e. all (e0|f0) classes.  push onto stack of classes. */

    /*printf("\n\nBegining new Shell Quartet:\n");*/
    dp_use = DP;
    n_hrr = 0;
    H_Cl = Classes;

    /* recursively list (and allocate room for) HRR generated intermediates */
    dum = List_HRR(H_Cl, orig_am, &n_hrr, &dp_use);

    first_vrr = 1;
    last_vrr = 1;

    /* set V_Cl to after end of hrr generated intermediates */
    V_Cl = &(Classes[n_hrr]);

    /* zero flags for HRR generated intermediates */
    bzero(H_done,n_hrr*sizeof(char));

    /* now loop over all elements in Classes (e0|f0) -
       apply vrr to push all those classes onto stack of classes */
    for(class_it = 0; class_it < n_hrr; class_it++){
      if(H_Cl[class_it].type != 0){
        Top_VRR(class_it, H_Cl, V_Cl, &last_vrr, &dp_use);
      }
    }

    /* contract by primitives out here */
    for (pi = 0; pi < shells[si].n_prims; pi++){
      for (pj = 0; pj < shells[sj].n_prims; pj++){
        for (pk = 0; pk < shells[sk].n_prims; pk++){
          for (pl = 0; pl < shells[sl].n_prims; pl++){
         
            /* zero flags for vrr generated intermediates */
            bzero(V_done,last_vrr*sizeof(char));

            pii = pi + shells[si].fprim-1;
            pjj = pj + shells[sj].fprim-1;
            pkk = pk + shells[sk].fprim-1;
            pll = pl + shells[sl].fprim-1;

            Shell_init(AB2, CD2, ericent, shells,
                       si, sj, sk, sl, pii, pjj, pkk, pll, cgtos);

            /* call function to begin build by VRR of top classes */
            for(class_it = 0; class_it < n_hrr; class_it++){
              if(H_Cl[class_it].type != 0){
                Init_VRR(H_Cl, V_Cl, class_it);
              }
            }
 
          }
        }
      }
    } /* end getting unique primitive set */

    /*  call a function to accumulate all the HRR generated classes */
    data = HRR_build(Classes, 0, &AB, &CD);
    num = Fill_data(Classes, 0, tot_data, ioffset, joffset, 
                    koffset, loffset, intfile);

    /* write out shell info */
    /*fprintf(eriout, "shell %2d %2d %2d %2d\n", si, sj, sk, sl);*/
    sz = sizeof(struct tebuf);

  } /* end if we need to do this shell block */
}
