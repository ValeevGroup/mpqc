

/* $Log$
 * Revision 1.1  1993/12/29 12:53:13  etseidl
 * Initial revision
 *
 * Revision 1.2  1992/09/11  12:14:43  seidl
 * add nrot and ntrans to char table
 *
 * Revision 1.1.1.1  1992/03/17  17:10:51  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:10:49  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/21  11:31:53  seidl
 * use chemistry/qc/intv2/int_libv2
 *
 * Revision 1.1  1991/12/20  15:45:23  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/17  21:55:18  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/03  16:50:10  etseidl
 * Initial revision
 *
 * Revision 1.2  1991/12/02  19:53:33  seidl
 * *** empty log message ***
 *
 * Revision 1.1  1991/11/22  18:28:13  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>

#ifndef IOFF
#define IOFF(a,b) (a>b)?a*(a+1)/2+b:b*(b+1)/2+a
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "symm.h"
#include "symm_mac.h"
#include "symmallc.h"
#include "symmfree.h"
#include "symerr.gbl"

#include "mkrpd.gbl"
#include "mkrpd.lcl"

#if 1
static double _pure_d[6][5] = {
  -0.28867513459481288225, -0.86602540378443864676, 0.0, 0.0, 0.0,
   0.0,                     0.0,                    0.0, 0.0, 1.0,
   0.57735026918962576450,  0.0,                    0.0, 0.0, 0.0,
   0.0,                     0.0,                    1.0, 0.0, 0.0,
   0.0,                     0.0,                    0.0, 1.0, 0.0,
  -0.28867513459481288225,  0.86602540378443864676, 0.0, 0.0, 0.0};
#else
static double _pure_d[6][5] = {
  -0.40824829046386301636, -0.70710678118654752440, 0.0, 0.0, 0.0,
   0.0,                     0.0,                    0.0, 0.0, 1.0,
   0.81649658092772603273,  0.0,                    0.0, 0.0, 0.0,
   0.0,                     0.0,                    1.0, 0.0, 0.0,
   0.0,                     0.0,                    0.0, 1.0, 0.0,
  -0.40824829046386301636,  0.70710678118654752440, 0.0, 0.0, 0.0};
#endif

static double *_pd[6] = 
  {_pure_d[0], _pure_d[1], _pure_d[2], _pure_d[3], _pure_d[4], _pure_d[5]};

GLOBAL_FUNCTION int
make_rp_d(si, n, theta, pg, f_exist, g_exist, _outfile)
sym_struct_t *si;
int n;
double theta;
enum pgroups pg;
int f_exist;
int g_exist;
FILE *_outfile;
{
  int i,j,k,l,m;
  int i0,j0;
  int jk,jkl,jklm;
  int errcod;
  char *label;
  double ct,st,c2t,s2t;
  double *rot=(double *) malloc(sizeof(double)*si->g);
  double *trans=(double *) malloc(sizeof(double)*si->g);
  double ineg1;
  double **rp, **rd, **rf, **rg;
  double **rpd, **rpf, **rpg;
  double **rpn;
#if 1
  double xynorm = sqrt(3.0);
  double xynormi = 1.0/sqrt(3.0);
#else
  double xynorm = sqrt(2.0);
  double xynormi = 1.0/sqrt(2.0);
#endif

  check_alloc(rot,"rot in make_rp_d");
  bzero(rot,sizeof(double)*si->g);

  check_alloc(trans,"trans in make_rp_d");
  bzero(trans,sizeof(double)*si->g);

/* This function will, given the point group and order, generate the
 * tranformation matrices for p and higher functions (up to f right now)
 * and will also generate a character table for the point group.
 * This character table is in the order that symmetry operations are
 * generated, not in Cotton order. If this is a problem, tough.
 */

/* please note that 0=y, 1=z, and 2=x because of the way integrals are
 * calculated.  for d shells order is 0=yy, 1=yz, 2=zz, 3=xy, 4=xz, 5=xx.
 * f will be 0=yyy, 1=yyz, 2=yzz, 3=zzz, 4=xyy, 5=xyz, 6=xzz, 7=xxy, 8=xxz,
 * 9=xxx. And so on.  
 *
 * This is a hard-wired hack job. It would be nice to one day make this 
 * completely * general (ie  any angular momentum, any axes, etc.).
 *
 * Also note, I am forming the transpose of the R matrices in Dupuis&King's
 * paper so that I can perform a more efficient matrix multiply in 
 * sym_sym_mat().
 */
  
  switch(pg) {
  case _PG_C1: /* no symmetry */
    si->Rp[0][0][0] = si->Rp[0][1][1] = si->Rp[0][2][2] = 1.0;

    errcod = allocbn_character(&si->ct.gamma[0],"g degen label nrot ntrans",
                                                          si->g,1,"A",3,3);
    if(errcod != 0) return rperr(_outfile,0);

    si->ct.gamma[0].rep[0]=1.0;
    break;

  case _PG_CI: /* equivalent to S2 about the z axis */
    si->Rp[0][0][0] = si->Rp[0][1][1] = si->Rp[0][2][2] = 1.0;
    si->Rp[1][0][0] = si->Rp[1][1][1] = si->Rp[1][2][2] = -1.0;

    errcod = allocbn_character(&si->ct.gamma[0],"g degen label nrot ntrans",
                                                        si->g,1,"Ag",3,0);
    if(errcod != 0) return rperr(_outfile,0);
    errcod = allocbn_character(&si->ct.gamma[1],"g degen label nrot ntrans",
                                                        si->g,1,"Au",0,3);
    if(errcod != 0) return rperr(_outfile,1);

    si->ct.gamma[0].rep[0]=  1.0;
    si->ct.gamma[0].rep[1]=  1.0;
    si->ct.gamma[1].rep[0]=  1.0;
    si->ct.gamma[1].rep[1]= -1.0;
    break;

  case _PG_CS: /* reflection through the xy plane */
    si->Rp[0][0][0] = si->Rp[0][1][1] = si->Rp[0][2][2] = 1.0;
    si->Rp[1][0][0] = si->Rp[1][2][2] = 1.0;
    si->Rp[1][1][1] = -1.0;

    errcod = allocbn_character(&si->ct.gamma[0],"g degen label nrot ntrans",
                                                         si->g,1,"A\'",1,2);
    if(errcod != 0) return rperr(_outfile,0);
    errcod = allocbn_character(&si->ct.gamma[1],"g degen label nrot ntrans",
                                                         si->g,1,"A\"",2,1);
    if(errcod != 0) return rperr(_outfile,1);

    si->ct.gamma[0].rep[0]=  1.0;
    si->ct.gamma[0].rep[1]=  1.0;
    si->ct.gamma[1].rep[0]=  1.0;
    si->ct.gamma[1].rep[1]= -1.0;
    break;

  case _PG_CN: /* clockwise rotation about z axis by theta*i radians */
    for(i=0; i < n ; i++) {
      ct = cos((double)theta*i);
      st = sin((double)theta*i);
      rp = si->Rp[i];

      rp[0][0] = rp[2][2] = ct;
      rp[0][2] = st;
      rp[2][0] = -st;
      rp[1][1] = 1.0;

      rot[i] = trans[i] = rp[0][0]+rp[1][1]+rp[2][2];
      }

    errcod = allocbn_character(&si->ct.gamma[0],"g degen label",si->g,1,"A");
    if(errcod != 0) return rperr(_outfile,0);
    for(i=0; i < si->g ; i++) si->ct.gamma[0].rep[i] = 1.0;

    if(n%2) {
      for(i=1; i < si->nirrep ; i++) {
        label = (char *) malloc(4);
        if(n==3) sprintf(label,"E");
        else sprintf(label,"E%d",i);
        errcod = 
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);

        si->ct.gamma[i].rep[0] = 2.0;
        for(j=1; j < si->g ; j++) 
          si->ct.gamma[i].rep[j] = 2.0*cos((double)theta*j*i);
        }
      }
    else {
      errcod = allocbn_character(&si->ct.gamma[1],"g degen label",si->g,1,"B");
      if(errcod != 0) return rperr(_outfile,1);
      for(i=0; i < si->g ; i++) si->ct.gamma[1].rep[i] = pow(-1.0,(double)i);

      for(i=2; i < si->nirrep ; i++) {
        label = (char *) malloc(4);
        if(n==4) sprintf(label,"E");
        else sprintf(label,"E%d",i-1);
        errcod = 
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);
        for(j=0; j < si->g ; j++) 
          si->ct.gamma[i].rep[j] = 2.0*cos((double)theta*j*(i-1));
        }
      }

    break;
  case _PG_CNV: /* clockwise rotation about z axis by theta*i radians, then
             *  reflect through xz plane */
    for(i=0; i < n ; i++) {
      ct = cos((double)theta*i);
      st = sin((double)theta*i);
      rp = si->Rp[i];

      rp[0][0] = rp[2][2] = ct;
      rp[0][2] = st;
      rp[2][0] = -st;
      rp[1][1] = 1.0;

      rot[i] = trans[i] = rp[0][0]+rp[1][1]+rp[2][2];
      }

    for(i=0,j=n; i < n ; i++,j++) {
      ct = cos((double)theta*i);
      st = sin((double)theta*i);
      rp = si->Rp[j];

      rp[0][0] = -ct;
      rp[2][2] = ct;
      rp[2][0] = rp[0][2] = st;
      rp[1][1] = 1.0;

      trans[j] = rp[0][0]+rp[1][1]+rp[2][2];
      rot[j] = -trans[j];
      }

    errcod = allocbn_character(&si->ct.gamma[0],"g degen label",si->g,1,"A1");
    if(errcod != 0) return rperr(_outfile,0);
    for(i=0; i < si->g ; i++) si->ct.gamma[0].rep[i] = 1.0;

    errcod = allocbn_character(&si->ct.gamma[1],"g degen label",si->g,1,"A2");
    if(errcod != 0) return rperr(_outfile,1);
    for(i=0; i < si->g/2 ; i++) si->ct.gamma[1].rep[i] = 1.0;
    for(i=si->g/2; i < si->g ; i++) si->ct.gamma[1].rep[i] = -1.0;

    if(n%2) {
      for(i=2; i < si->nirrep ; i++) {
        label = (char *) malloc(4);
        if(n==3) sprintf(label,"E");
        else sprintf(label,"E%d",i-1);
        errcod = 
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);

        si->ct.gamma[i].rep[0] = 2.0;
        for(j=1; j < si->g/2 ; j++) 
          si->ct.gamma[i].rep[j] = 2.0*cos((double)theta*j*(i-1));
        for(j=si->g/2; j < si->g ; j++) si->ct.gamma[i].rep[j] = 0.0;
        }
      }
    else {
      errcod = allocbn_character(&si->ct.gamma[2],"g degen label",si->g,1,"B1");
      if(errcod != 0) return rperr(_outfile,2);
      for(i=0; i < si->g/2 ; i++) si->ct.gamma[2].rep[i] = pow(-1.0,(double)i);
      for(i=si->g/2; i < si->g ; i++) si->ct.gamma[2].rep[i] = 
                                                           pow(-1.0,(double)i);
      errcod = allocbn_character(&si->ct.gamma[3],"g degen label",si->g,1,"B2");
      if(errcod != 0) return rperr(_outfile,3);
      for(i=0; i < si->g/2 ; i++) si->ct.gamma[3].rep[i] = pow(-1.0,(double)i);
      for(i=si->g/2; i < si->g ; i++) si->ct.gamma[3].rep[i] = 
                                                   pow(-1.0,(double)(i+1));

      for(i=4; i < si->nirrep; i++) {
        label = (char *) malloc(4);
        if(n==4) sprintf(label,"E");
        else sprintf(label,"E%d",i-3);
        errcod = 
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);
        for(j=0; j < si->g/2 ; j++) 
          si->ct.gamma[i].rep[j] = 2.0*cos((double)theta*j*(i-3));
        for(j=si->g/2; j < si->g ; j++) si->ct.gamma[i].rep[j] = 0.0;
        }
      }

    break;

  case _PG_CNH: /* clockwise rotation about z axis by theta*i radians, 
             * as well as rotation-reflection about same axis */
    for(i=0; i < n ; i++) {
      ct = cos((double)theta*i);
      st = sin((double)theta*i);
      rp = si->Rp[i]; rpn = si->Rp[i+n];

      rp[0][0] = rp[2][2] = rpn[0][0] = rpn[2][2] = ct;
      rp[0][2] = rpn[0][2] = st;
      rp[2][0] = rpn[2][0] = -st;
      rp[1][1] = 1.0; rpn[1][1] = -1.0;

      rot[i] = trans[i] = rp[0][0]+rp[1][1]+rp[2][2];
      trans[i+n] = rpn[0][0]+rpn[1][1]+rpn[2][2];
      rot[i+n] = -trans[i+n];
      }

    if(n%2) {
      errcod = 
        allocbn_character(&si->ct.gamma[0],"g degen label",si->g,1,"A\'");
      if(errcod != 0) return rperr(_outfile,0);
      for(i=0; i < si->g ; i++) si->ct.gamma[0].rep[i] = 1.0;

      for(i=1; i < si->nirrep/2 ; i++) {
        label = (char *) malloc(4);
        if(n==3) sprintf(label,"E\'");
        else sprintf(label,"E%d\'",i);
        errcod = 
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);

        si->ct.gamma[i].rep[0] = 2.0;
        for(j=1; j < si->g/2 ; j++) 
                        si->ct.gamma[i].rep[j] = 2.0*cos((double)theta*j*i);
        for(j=si->g/2; j < si->g ; j++) si->ct.gamma[i].rep[j] = 
                                           si->ct.gamma[i].rep[j-(si->g/2)];
        }

      i=si->nirrep/2;
      errcod = 
        allocbn_character(&si->ct.gamma[i],"g degen label",si->g,1,"A\"");
      if(errcod != 0) return rperr(_outfile,i);
      for(j=0; j < si->g/2 ; j++) si->ct.gamma[i].rep[j] = 1.0;
      for(j=si->g/2; j < si->g ; j++) si->ct.gamma[i].rep[j] = -1.0;

      for(i=i+1; i < si->nirrep ; i++) {
        i0 = i-(si->nirrep/2);
        label = (char *) malloc(4);
        if(n==3) sprintf(label,"E\"");
        else sprintf(label,"E%d\"",i0);
        errcod = 
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);
        for(j=0; j < si->g/2 ; j++) si->ct.gamma[i].rep[j] = 
                                                       si->ct.gamma[i0].rep[j];
        for(j=si->g/2; j < si->g ; j++) si->ct.gamma[i].rep[j] = 
                                                      -si->ct.gamma[i0].rep[j];
        }
      }
    else {
      double ineg = pow(-1.0,(double)n/2);

      errcod = 
        allocbn_character(&si->ct.gamma[0],"g degen label",si->g,1,"Ag");
      if(errcod != 0) return rperr(_outfile,0);
      for(i=0; i < si->g ; i++) si->ct.gamma[0].rep[i] = 1.0;
      errcod = 
        allocbn_character(&si->ct.gamma[1],"g degen label",si->g,1,"Bg");
      if(errcod != 0) return rperr(_outfile,1);
      for(i=0; i < si->g/2 ; i++) si->ct.gamma[1].rep[i] = 
                                                        pow(-1.0,(double) i);
      for(i=si->g/2; i < si->g ; i++) si->ct.gamma[1].rep[i] = 
                                                   ineg*pow(-1.0,(double) i);

      for(i=2; i < si->nirrep/2 ; i++) {
        ineg=pow(-1.0,(double) (i-1));
        label = (char *) malloc(4);
        if(n==4) sprintf(label,"Eg");
        else sprintf(label,"E%dg",i-1);
        errcod = 
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);

        si->ct.gamma[i].rep[0] = 2.0;
        for(j=1; j < si->g/2 ; j++) 
          si->ct.gamma[i].rep[j] = 2.0*cos((double)theta*j*(i-1));
        for(j=si->g/2; j < si->g ; j++) si->ct.gamma[i].rep[j] = 
                                        ineg*si->ct.gamma[i].rep[j-(si->g/2)];
        }

      i=si->nirrep/2;
      errcod = 
        allocbn_character(&si->ct.gamma[i],"g degen label",si->g,1,"Au");
      if(errcod != 0) return rperr(_outfile,i);
      for(j=0; j < si->g/2 ; j++) si->ct.gamma[i].rep[j] = 1.0;
      for(j=si->g/2; j < si->g ; j++) si->ct.gamma[i].rep[j] = -1.0;
      i++;
      errcod = 
        allocbn_character(&si->ct.gamma[i],"g degen label",si->g,1,"Bu");
      if(errcod != 0) return rperr(_outfile,i);
      for(j=0; j < si->g/2 ; j++) si->ct.gamma[i].rep[j] = 
                                                        si->ct.gamma[1].rep[j];
      for(j=si->g/2; j < si->g ; j++) si->ct.gamma[i].rep[j] =
                                                       -si->ct.gamma[1].rep[j];

      for(i=i+1; i <= si->nirrep ; i++) {
        i0=i-(si->nirrep/2);
        label = (char *) malloc(4);
        if(n==4) sprintf(label,"Eu");
        else sprintf(label,"E%du",i0-1);
        errcod = 
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);

        for(j=0; j < si->g/2 ; j++) si->ct.gamma[i].rep[j] = 
                                                     si->ct.gamma[i0].rep[j];
        for(j=si->g/2; j < si->g ; j++) si->ct.gamma[i].rep[j] = 
                                                    -si->ct.gamma[i0].rep[j];
        }
      }
    break;

  case _PG_SN: /* clockwise rotation-reflection by theta*i radians about 
            * z axis */
    for(i=0; i < n ; i++) {
      ct = cos((double)theta*i);
      st = sin((double)theta*i);
      rp = si->Rp[i];
      ineg1 = pow(-1.0,(double) i);

      rp[0][0] = rp[2][2] = ct;
      rp[0][2] = st;
      rp[2][0] = -st;
      rp[1][1] = ineg1;

      trans[i] = rp[0][0]+rp[1][1]+rp[2][2];
      rot[i] = (i%2) ? -trans[i] : trans[i];
      }

    if((n/2)%2) {
      errcod = allocbn_character(&si->ct.gamma[0],"g degen label",si->g,1,"Ag");
      if(errcod != 0) return rperr(_outfile,0);
      for(i=0; i < si->g ; i++) si->ct.gamma[0].rep[i] = 1.0;

      for(i=1; i < si->nirrep/2 ; i++) {
        label = (char *) malloc(4);
        if(n==6) sprintf(label,"Eg");
        else sprintf(label,"E%dg",i);
        errcod = 
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);

        si->ct.gamma[i].rep[0] = 2.0;
        for(j=1; j < si->g ; j++) 
          si->ct.gamma[i].rep[j] = 2.0*cos((double)2.0*theta*j*i);
        }

      i=si->nirrep/2;
      errcod = allocbn_character(&si->ct.gamma[i],"g degen label",si->g,1,"Au");
      if(errcod != 0) return rperr(_outfile,i);
      for(j=0; j < si->g ; j++) si->ct.gamma[i].rep[j] = pow(-1.0,(double)j);

      for(i=i+1; i < si->nirrep ; i++) {
        i0=i-si->nirrep/2;
        label = (char *) malloc(4);
        if(n==6) sprintf(label,"Eu");
        else sprintf(label,"E%du",i);
        errcod = 
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);

        for(j=0; j < si->g ; j++) si->ct.gamma[i].rep[j] = 
                                pow(-1.0,(double) j)*si->ct.gamma[i0].rep[j];
        }
      }
    else {
      errcod = allocbn_character(&si->ct.gamma[0],"g degen label",si->g,1,"A");
      if(errcod != 0) return rperr(_outfile,0);
      for(i=0; i < si->g ; i++) si->ct.gamma[0].rep[i] = 1.0;

      errcod = allocbn_character(&si->ct.gamma[1],"g degen label",si->g,1,"B");
      if(errcod != 0) return rperr(_outfile,1);
      for(i=0; i < si->g ; i++) si->ct.gamma[1].rep[i] = pow(-1.0,(double)i);

      for(i=2; i < si->nirrep ; i++) {
        label = (char *) malloc(4);
        if(n==4) sprintf(label,"E");
        else sprintf(label,"E%d",i-1);
        errcod = 
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);
        for(j=0; j < si->g ; j++) 
          si->ct.gamma[i].rep[j] = 2.0*cos((double)theta*j*(i-1));
        }
      }
    break;

  case _PG_DN: /* clockwise rotation about z axis, followed by C2 about 
            * x axis */
    for(i=0; i < n ; i++) {
      ct = cos((double)theta*i);
      st = sin((double)theta*i);
      rp = si->Rp[i];

      rp[0][0] = rp[2][2] = ct;
      rp[0][2] = st;
      rp[2][0] = -st;
      rp[1][1] = 1.0;

      rot[i] = trans[i] = rp[0][0]+rp[1][1]+rp[2][2];
      }

    for(i=0,j=n; i < n ; i++,j++) {
      ct = cos((double)theta*i);
      st = sin((double)theta*i);
      rp = si->Rp[j];

      rp[0][0] = -ct;
      rp[2][2] = ct;
      rp[0][2] = rp[2][0] = -st;
      rp[1][1] = -1.0;

      rot[j] = trans[j] = rp[0][0]+rp[1][1]+rp[2][2];
      }

    if(n==2) {
      errcod = allocbn_character(&si->ct.gamma[0],"g degen label",si->g,1,"A");
      if(errcod != 0) return rperr(_outfile,0);
      for(i=0; i < si->g ; i++) si->ct.gamma[0].rep[i] = 1.0;

      errcod = allocbn_character(&si->ct.gamma[1],"g degen label",si->g,1,"B1");
      if(errcod != 0) return rperr(_outfile,1);
      si->ct.gamma[1].rep[0] = si->ct.gamma[1].rep[1] = 1.0;
      si->ct.gamma[1].rep[2] = si->ct.gamma[1].rep[3] = -1.0;

      errcod = allocbn_character(&si->ct.gamma[2],"g degen label",si->g,1,"B2");
      if(errcod != 0) return rperr(_outfile,2);
      si->ct.gamma[2].rep[0] = si->ct.gamma[2].rep[3] = 1.0;
      si->ct.gamma[2].rep[1] = si->ct.gamma[2].rep[2] = -1.0;

      errcod = allocbn_character(&si->ct.gamma[3],"g degen label",si->g,1,"B3");
      if(errcod != 0) return rperr(_outfile,3);
      si->ct.gamma[3].rep[0] = si->ct.gamma[3].rep[2] = 1.0;
      si->ct.gamma[3].rep[1] = si->ct.gamma[3].rep[3] = -1.0;
      }
    else if(n%2) {
      errcod = allocbn_character(&si->ct.gamma[0],"g degen label",si->g,1,"A1");
      if(errcod != 0) return rperr(_outfile,0);
      for(i=0; i < si->g ; i++) si->ct.gamma[0].rep[i] = 1.0;

      errcod = allocbn_character(&si->ct.gamma[1],"g degen label",si->g,1,"A2");
      if(errcod != 0) return rperr(_outfile,1);
      for(i=0; i < si->g/2 ; i++) si->ct.gamma[1].rep[i] = 1.0;
      for(i=si->g/2; i < si->g ; i++) si->ct.gamma[1].rep[i] = -1.0;

      for(i=2; i < si->nirrep ; i++) {
        label = (char *) malloc(4);
        if(n==3) sprintf(label,"E");
        else sprintf(label,"E%d",i-1);
        errcod = 
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);

        si->ct.gamma[i].rep[0] = 2.0;
        for(j=1; j < si->g/2 ; j++) 
          si->ct.gamma[i].rep[j] = 2.0*cos((double)theta*j*(i-1));
        for(j=si->g/2; j < si->g ; j++) si->ct.gamma[i].rep[j] = 0.0;
        }
      }
    else {
      errcod = allocbn_character(&si->ct.gamma[0],"g degen label",si->g,1,"A1");
      if(errcod != 0) return rperr(_outfile,0);
      for(i=0; i < si->g ; i++) si->ct.gamma[0].rep[i] = 1.0;

      errcod = allocbn_character(&si->ct.gamma[1],"g degen label",si->g,1,"A2");
      if(errcod != 0) return rperr(_outfile,1);
      for(i=0; i < si->g/2 ; i++) si->ct.gamma[1].rep[i] = 1.0;
      for(i=si->g/2; i < si->g ; i++) si->ct.gamma[1].rep[i] = -1.0;

      errcod = allocbn_character(&si->ct.gamma[2],"g degen label",si->g,1,"B1");
      if(errcod != 0) return rperr(_outfile,2);
      for(i=0; i < si->g/2 ; i++) si->ct.gamma[2].rep[i] = pow(-1.0,(double)i);
      for(i=si->g/2; i < si->g ; i++) si->ct.gamma[2].rep[i] = 
                                                           pow(-1.0,(double)i);

      errcod = allocbn_character(&si->ct.gamma[3],"g degen label",si->g,1,"B2");
      if(errcod != 0) return rperr(_outfile,3);
      for(i=0; i < si->g/2 ; i++) si->ct.gamma[3].rep[i] = pow(-1.0,(double)i);
      for(i=si->g/2; i < si->g ; i++) si->ct.gamma[3].rep[i] = 
                                                   pow(-1.0,(double)(i+1));

      for(i=4; i < si->nirrep ; i++) {
        label = (char *) malloc(4);
        if(n==4) sprintf(label,"E");
        else sprintf(label,"E%d",i-3);
        errcod = 
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);
        for(j=0; j < si->g/2 ; j++) 
          si->ct.gamma[i].rep[j] = 2.0*cos((double)theta*j*(i-3));
        for(j=si->g/2; j < si->g ; j++) si->ct.gamma[i].rep[j] = 0.0;
        }
      }
    break;

  case _PG_DND: /* rotation reflection about z axis by theta/2 radians, followed
             * by c2 about x axis, then reflection through yz plane */
    for(i=0; i < 2*n ; i++) {
      ct = cos((double)theta*i*0.5);
      st = sin((double)theta*i*0.5);
      rp = si->Rp[i];
      ineg1 = pow(-1.0,(double) i);

      rp[0][0] = rp[2][2] = ct;
      rp[0][2] = st;
      rp[2][0] = -st;
      rp[1][1] = ineg1;

      trans[i] = rp[0][0]+rp[1][1]+rp[2][2];
      rot[i] = (i%2) ? -trans[i] : trans[i];
      }

    for(i=0,j=2*n; i < n ; i++,j++) {
      ct = cos((double)theta*i);
      st = sin((double)theta*i);
      rp = si->Rp[j];

      rp[0][0] = -ct;
      rp[2][2] = ct;
      rp[0][2] = rp[2][0] = -st;
      rp[1][1] = -1.0;

      rot[j] = trans[j] = rp[0][0]+rp[1][1]+rp[2][2];
      }

    for(i=0,j=3*n; i < n ; i++,j++) {
      ct = cos(((double) i+0.5)*theta);
      st = sin(((double) i+0.5)*theta);
      rp = si->Rp[j];

      rp[0][0] = -ct;
      rp[2][2] = ct;
      rp[0][2] = rp[2][0] = -st;
      rp[1][1] = 1.0;

      trans[j] = rp[0][0]+rp[1][1]+rp[2][2];
      rot[j] = -trans[j];
      }

    if(n%2) {
      errcod = 
        allocbn_character(&si->ct.gamma[0],"g degen label",si->g,1,"A1g");
      if(errcod != 0) return rperr(_outfile,0);
      for(i=0; i < si->g ; i++) si->ct.gamma[0].rep[i] = 1.0;

      errcod = 
        allocbn_character(&si->ct.gamma[1],"g degen label",si->g,1,"A2g");
      if(errcod != 0) return rperr(_outfile,1);
      for(i=0; i < si->g/2 ; i++) si->ct.gamma[1].rep[i] = 1.0;
      for(i=si->g/2; i < si->g ; i++) si->ct.gamma[1].rep[i] = -1.0;

      for(i=2; i < si->nirrep/2 ; i++) {
        label = (char *) malloc(4);
        if(n==3) sprintf(label,"Eg");
        else sprintf(label,"E%dg",i-1);
        errcod =
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);

        si->ct.gamma[i].rep[0] = 2.0;
        for(j=1; j < si->g/2 ; j++)
          si->ct.gamma[i].rep[j] = 2.0*cos((double)2.0*theta*j*(i-1));
        for(j=si->g/2; j < si->g ; j++) si->ct.gamma[i].rep[j] = 0.0;
        }

      i=si->nirrep/2;
      errcod = 
        allocbn_character(&si->ct.gamma[i],"g degen label",si->g,1,"A1u");
      if(errcod != 0) return rperr(_outfile,i);
      i++;
      errcod = 
        allocbn_character(&si->ct.gamma[i],"g degen label",si->g,1,"A2u");
      if(errcod != 0) return rperr(_outfile,i);
      i++;

      for(; i < si->nirrep ; i++) {
        i0=i-si->nirrep/2;
        label = (char *) malloc(4);
        if(n==3) sprintf(label,"Eu");
        else sprintf(label,"E%du",i0-1);
        errcod =
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);
        }

      for(i=si->nirrep/2,i0=0; i < si->nirrep; i++,i0++) {
        for(j=0; j < si->g/2 ; j++) si->ct.gamma[i].rep[j] =
                                pow(-1.0,(double) j)*si->ct.gamma[i0].rep[j];
        for(j=si->g/2; j < 3*si->g/4 ; j++) si->ct.gamma[i].rep[j] =
                                                     si->ct.gamma[i0].rep[j];
        for(j=3*si->g/4; j < si->g ; j++) si->ct.gamma[i].rep[j] =
                                                    -si->ct.gamma[i0].rep[j];
        }
      }
    else {
      errcod = 
        allocbn_character(&si->ct.gamma[0],"g degen label",si->g,1,"A1");
      if(errcod != 0) return rperr(_outfile,0);
      for(i=0; i < si->g ; i++) si->ct.gamma[0].rep[i] = 1.0;

      errcod = 
        allocbn_character(&si->ct.gamma[1],"g degen label",si->g,1,"A2");
      if(errcod != 0) return rperr(_outfile,1);
      for(i=0; i < si->g/2 ; i++) si->ct.gamma[1].rep[i] = 1.0;
      for(i=si->g/2; i < si->g ; i++) si->ct.gamma[1].rep[i] = -1.0;

      errcod = 
        allocbn_character(&si->ct.gamma[2],"g degen label",si->g,1,"B1");
      if(errcod != 0) return rperr(_outfile,2);
      for(i=0; i < si->g/2 ; i++) si->ct.gamma[2].rep[i] = 
                                                        pow(-1.0,(double) i);
      for(i=si->g/2; i < 3*si->g/4 ; i++) si->ct.gamma[2].rep[i] =  1.0;
      for(i=3*si->g/4; i < si->g ; i++) si->ct.gamma[2].rep[i] = -1.0;
                                                   
      errcod = 
        allocbn_character(&si->ct.gamma[3],"g degen label",si->g,1,"B2");
      if(errcod != 0) return rperr(_outfile,3);
      for(i=0; i < si->g/2 ; i++) si->ct.gamma[3].rep[i] = 
                                                        pow(-1.0,(double) i);
      for(i=si->g/2; i < 3*si->g/4 ; i++) si->ct.gamma[3].rep[i] = -1.0;
      for(i=3*si->g/4; i < si->g ; i++) si->ct.gamma[3].rep[i] = 1.0;

      for(i=4; i < si->nirrep ; i++) {
        label = (char *) malloc(4);
        if(n==2) sprintf(label,"E");
        else sprintf(label,"E%d",i-3);
        errcod =
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);

        si->ct.gamma[i].rep[0] = 2.0;
        for(j=1; j < si->g/2 ; j++)
          si->ct.gamma[i].rep[j] = 2.0*cos((double)theta*j*(i-3)*0.5);
        for(j=si->g/2; j < si->g ; j++) si->ct.gamma[i].rep[j] = 0.0;
        }
      }
                                                   
    break;

  case _PG_DNH: /* clockwise rotation and rotation-reflection about z axis,
             * followed by c2 about x axis and then reflection 
             * through xz */
    for(i=0; i < n ; i++) {
      ct = cos((double)theta*i);
      st = sin((double)theta*i);
      rp = si->Rp[i]; rpn = si->Rp[i+n];

      rp[0][0] = rp[2][2] = rpn[0][0] = rpn[2][2] = ct;
      rp[0][2] = rpn[0][2] = st;
      rp[2][0] = rpn[2][0] = -st;
      rp[1][1] = 1.0; rpn[1][1] = -1.0;

      rot[i] = trans[i] = rp[0][0]+rp[1][1]+rp[2][2];
      trans[i+n] = rpn[0][0]+rpn[1][1]+rpn[2][2];
      rot[i+n] = -trans[i+n];
      }

    for(i=0,j=2*n; i < n ; i++,j++) {
      ct = cos((double)theta*i);
      st = sin((double)theta*i);
      rp = si->Rp[j];

      rp[0][0] = -ct;
      rp[2][2] = ct;
      rp[0][2] = rp[2][0] = -st;
      rp[1][1] = -1.0;

      rot[j] = trans[j] = rp[0][0]+rp[1][1]+rp[2][2];
      }

    for(i=0,j=3*n; i < n ; i++,j++) {
      ct = cos((double)theta*i);
      st = sin((double)theta*i);
      rp = si->Rp[j];

      rp[0][0] = -ct;
      rp[2][2] = ct;
      rp[2][0] = rp[0][2] = st;
      rp[1][1] = 1.0;

      trans[j] = rp[0][0]+rp[1][1]+rp[2][2];
      rot[j] = -trans[j];
      }

    if(n==2) {
      errcod = allocbn_character(&si->ct.gamma[0],"g degen label",si->g,1,"Ag");
      if(errcod != 0) return rperr(_outfile,0);
      for(i=0; i < 8 ; i++) si->ct.gamma[0].rep[i] = 1.0;

      errcod = 
        allocbn_character(&si->ct.gamma[1],"g degen label",si->g,1,"B1g");
      if(errcod != 0) return rperr(_outfile,1);
      for(i=0; i < 4 ; i++) si->ct.gamma[1].rep[i] = 1.0;
      for(i=4; i < 8 ; i++) si->ct.gamma[1].rep[i] = -1.0;

      errcod = 
        allocbn_character(&si->ct.gamma[2],"g degen label",si->g,1,"B2g");
      if(errcod != 0) return rperr(_outfile,2);
      si->ct.gamma[2].rep[0] = si->ct.gamma[2].rep[3] =
       si->ct.gamma[2].rep[5] = si->ct.gamma[2].rep[6] = 1.0;
      si->ct.gamma[2].rep[1] = si->ct.gamma[2].rep[2] =
       si->ct.gamma[2].rep[4] = si->ct.gamma[2].rep[7] = -1.0;

      errcod = 
        allocbn_character(&si->ct.gamma[3],"g degen label",si->g,1,"B3g");
      if(errcod != 0) return rperr(_outfile,3);
      si->ct.gamma[3].rep[0] = si->ct.gamma[3].rep[3] =
       si->ct.gamma[3].rep[4] = si->ct.gamma[3].rep[7] = 1.0;
      si->ct.gamma[3].rep[1] = si->ct.gamma[3].rep[2] =
       si->ct.gamma[3].rep[5] = si->ct.gamma[3].rep[6] = -1.0;

      errcod = allocbn_character(&si->ct.gamma[4],"g degen label",si->g,1,"Au");
      if(errcod != 0) return rperr(_outfile,4);
      for(i=0; i < 2 ; i++) si->ct.gamma[4].rep[i] = 1.0;
      for(i=2; i < 4 ; i++) si->ct.gamma[4].rep[i] = -1.0;
      for(i=4; i < 6 ; i++) si->ct.gamma[4].rep[i] = 1.0;
      for(i=6; i < 8 ; i++) si->ct.gamma[4].rep[i] = -1.0;

      errcod = 
        allocbn_character(&si->ct.gamma[5],"g degen label",si->g,1,"B1u");
      if(errcod != 0) return rperr(_outfile,5);
      for(i=0; i < 2 ; i++) si->ct.gamma[5].rep[i] = 1.0;
      for(i=2; i < 6 ; i++) si->ct.gamma[5].rep[i] = -1.0;
      for(i=6; i < 8 ; i++) si->ct.gamma[5].rep[i] = 1.0;

      errcod = 
        allocbn_character(&si->ct.gamma[6],"g degen label",si->g,1,"B2u");
      if(errcod != 0) return rperr(_outfile,6);
      for(i=0; i < 4 ; i++) si->ct.gamma[6].rep[i] =  pow(-1.0,(double)i);
      for(i=4; i < 8 ; i++) si->ct.gamma[6].rep[i] = -pow(-1.0,(double)i);

      errcod = 
        allocbn_character(&si->ct.gamma[7],"g degen label",si->g,1,"B3u");
      if(errcod != 0) return rperr(_outfile,7);
      for(i=0; i < 8 ; i++) si->ct.gamma[7].rep[i] = pow(-1.0,(double)i);
      }
    else if(n%2) {
      errcod = 
        allocbn_character(&si->ct.gamma[0],"g degen label",si->g,1,"A1\'");
      if(errcod != 0) return rperr(_outfile,0);
      for(i=0; i < si->g ; i++) si->ct.gamma[0].rep[i] = 1.0;

      errcod = 
        allocbn_character(&si->ct.gamma[1],"g degen label",si->g,1,"A2\'");
      if(errcod != 0) return rperr(_outfile,1);
      for(i=0; i < si->g/2 ; i++) si->ct.gamma[1].rep[i] = 1.0;
      for(i=si->g/2,i0=0; i < si->g ; i++,i0++) si->ct.gamma[1].rep[i] = -1.0;

      for(i=2; i < si->nirrep/2 ; i++) {
        label = (char *) malloc(4);
        if(n==3) sprintf(label,"E\'");
        else sprintf(label,"E%d\'",i-1);
        errcod =
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);

        si->ct.gamma[i].rep[0] = 2.0;
        for(j=1; j < si->g/4 ; j++)
          si->ct.gamma[i].rep[j] = 2.0*cos((double)theta*j*(i-1));
        for(j=si->g/4; j < si->g/2 ; j++) si->ct.gamma[i].rep[j] = 
                                                2.0*cos((double)theta*j*(i-1));
        for(j=si->g/2,j0=0; j < si->g ; j++,j0++) si->ct.gamma[i].rep[j] = 0.0;
        }

      i=si->nirrep/2;
      errcod = 
        allocbn_character(&si->ct.gamma[i],"g degen label",si->g,1,"A1\"");
      if(errcod != 0) return rperr(_outfile,i);
      errcod = 
        allocbn_character(&si->ct.gamma[i+1],"g degen label",si->g,1,"A2\"");
      if(errcod != 0) return rperr(_outfile,i+1);

      for(i=i+2; i < si->nirrep ; i++) {
        i0=i-(si->nirrep/2);
        label = (char *) malloc(4);
        if(n==3) sprintf(label,"E\"");
        else sprintf(label,"E%d\"",i0-1);
        errcod =
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);
        }

      for(i=si->nirrep/2,i0=0; i < si->nirrep ; i0++,i++) {
        for(j=0; j < si->g/4 ; j++) si->ct.gamma[i].rep[j] = 
                                                       si->ct.gamma[i0].rep[j];
        for(j=si->g/4; j < si->g/2 ; j++) si->ct.gamma[i].rep[j] =
                                                      -si->ct.gamma[i0].rep[j];
        for(j=si->g/2; j < 3*si->g/4 ; j++) si->ct.gamma[i].rep[j] =
                                                       si->ct.gamma[i0].rep[j];
        for(j=3*si->g/4; j < si->g ; j++) si->ct.gamma[i].rep[j] =
                                                      -si->ct.gamma[i0].rep[j];
        }

      }
    else {
      errcod = 
        allocbn_character(&si->ct.gamma[0],"g degen label",si->g,1,"A1g");
      if(errcod != 0) return rperr(_outfile,0);
      for(i=0; i < si->g ; i++) si->ct.gamma[0].rep[i] = 1.0;

      errcod = 
        allocbn_character(&si->ct.gamma[1],"g degen label",si->g,1,"A2g");
      if(errcod != 0) return rperr(_outfile,1);
      for(i=0; i < si->g/2 ; i++) si->ct.gamma[1].rep[i] = 1.0;
      for(i=si->g/2; i < si->g ; i++) si->ct.gamma[1].rep[i] = -1.0;

      errcod = 
        allocbn_character(&si->ct.gamma[2],"g degen label",si->g,1,"B1g");
      if(errcod != 0) return rperr(_outfile,2);
      for(i=0; i < si->g ; i++) si->ct.gamma[2].rep[i] = pow(-1.0,(double)i);

      errcod =
        allocbn_character(&si->ct.gamma[3],"g degen label",si->g,1,"B2g");
      if(errcod != 0) return rperr(_outfile,3);
      for(i=0; i < si->g/2 ; i++) si->ct.gamma[3].rep[i] = pow(-1.0,(double)i);
      for(i=si->g/2; i < si->g ; i++) si->ct.gamma[3].rep[i] = 
                                                       pow(-1.0,(double)(i+1));

      for(i=4; i < si->nirrep/2 ; i++) {
        label = (char *) malloc(4);
        if(n==4) sprintf(label,"Eg");
        else sprintf(label,"E%dg",i-3);
        errcod =
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);
        for(j=0; j < si->g/4 ; j++)
          si->ct.gamma[i].rep[j] = 2.0*cos((double)theta*j*(i-3));
        for(j=si->g/4; j < si->g/2 ; j++) si->ct.gamma[i].rep[j] = 
                      pow(-1.0,(double)(i+1))*2.0*cos((double)theta*j*(i-3));
        for(j=si->g/2; j < si->g ; j++) si->ct.gamma[i].rep[j] = 0.0;
        }

      i=si->nirrep/2;
      errcod =
        allocbn_character(&si->ct.gamma[i],"g degen label",si->g,1,"A1u");
      if(errcod != 0) return rperr(_outfile,i);
      i++;
      errcod =
        allocbn_character(&si->ct.gamma[i],"g degen label",si->g,1,"A2u");
      if(errcod != 0) return rperr(_outfile,i);
      i++;
      errcod =
        allocbn_character(&si->ct.gamma[i],"g degen label",si->g,1,"B1u");
      if(errcod != 0) return rperr(_outfile,i);
      i++;
      errcod =
        allocbn_character(&si->ct.gamma[i],"g degen label",si->g,1,"B2u");
      if(errcod != 0) return rperr(_outfile,i);
      i++;
      for(; i < si->nirrep ; i++) {
        i0=i-(si->nirrep/2);
        label = (char *) malloc(4);
        if(n==4) sprintf(label,"Eu");
        else sprintf(label,"E%du",i0-3);
        errcod =
          allocbn_character(&si->ct.gamma[i],"g degen label",si->g,2,label);
        if(errcod != 0) return rperr(_outfile,i);
        }

      for(i=si->nirrep/2,i0=0; i < si->nirrep ; i0++,i++) {
        for(j=0; j < si->g/4 ; j++) si->ct.gamma[i].rep[j] =
                                                       si->ct.gamma[i0].rep[j];
        for(j=si->g/4; j < si->g/2 ; j++) si->ct.gamma[i].rep[j] =
                                                      -si->ct.gamma[i0].rep[j];
        for(j=si->g/2; j < 3*si->g/4 ; j++) si->ct.gamma[i].rep[j] =
                                                       si->ct.gamma[i0].rep[j];
        for(j=3*si->g/4; j < si->g ; j++) si->ct.gamma[i].rep[j] =
                                                      -si->ct.gamma[i0].rep[j];
        }
      }
    break;
  default:
    return(-1);
    }

/* ok, we have the reducible representation of the rotations and translations,
 * now let's use projection operators to find out how many rotations and
 * translations there are in each irrep
 */

  if(pg != _PG_C1 && pg != _PG_CI && pg != _PG_CS ) {
    for(i=0; i < si->ct.nirrep; i++) {
      double nr=0; double nt=0;
      for(j=0; j < si->ct.gamma[i].g; j++) {
        nr += rot[j]*si->ct.gamma[i].rep[j];
        nt += trans[j]*si->ct.gamma[i].rep[j];
        }
      if(pg==_PG_CN || pg==_PG_CNH || pg==_PG_SN) {
        si->ct.gamma[i].nrot =
             (nr+0.5)/(si->ct.gamma[i].g*si->ct.gamma[i].degen);
        si->ct.gamma[i].ntrans =
             (nt+0.5)/(si->ct.gamma[i].g*si->ct.gamma[i].degen);
        }
      else {
        si->ct.gamma[i].nrot = (nr+0.5)/si->ct.gamma[i].g;
        si->ct.gamma[i].ntrans = (nt+0.5)/si->ct.gamma[i].g;
        }
      }
    }

/* form the Rd matrices from the Rp matrices */

  for(i=0; i < si->g ; i++) {
    rp = si->Rp[i];
    rd = si->Rd[i];

    for(l=0; l < 6 ; l++) {
      switch (l) {
      case yy:
        j=y; k=y; jk = yy;
        break;
      case yz:
        j=y; k=z; jk = yz;
        break;
      case zz:
        j=z; k=z; jk = zz;
        break;
      case xy:
        j=x; k=y; jk = xy;
        break;
      case xz:
        j=x; k=z; jk = xz;
        break;
      case xx:
        j=x; k=x; jk = xx;
        break;
        }

      rd[jk][xx] = rp[j][x]*rp[k][x];
      rd[jk][yy] = rp[j][y]*rp[k][y];
      rd[jk][zz] = rp[j][z]*rp[k][z];
      rd[jk][xy] = rp[j][x]*rp[k][y]+rp[j][y]*rp[k][x];
      rd[jk][xz] = rp[j][x]*rp[k][z]+rp[j][z]*rp[k][x];
      rd[jk][yz] = rp[j][y]*rp[k][z]+rp[j][z]*rp[k][y];
      }

  /* and normalize a bit */
    rd[xy][xx] *= xynorm;
    rd[xy][yy] *= xynorm;
    rd[xy][zz] *= xynorm;
    rd[xz][xx] *= xynorm;
    rd[xz][yy] *= xynorm;
    rd[xz][zz] *= xynorm;
    rd[yz][xx] *= xynorm;
    rd[yz][yy] *= xynorm;
    rd[yz][zz] *= xynorm;

    rd[xx][xy] *= xynormi;
    rd[yy][xy] *= xynormi;
    rd[zz][xy] *= xynormi;
    rd[xx][xz] *= xynormi;
    rd[yy][xz] *= xynormi;
    rd[zz][xz] *= xynormi;
    rd[xx][yz] *= xynormi;
    rd[yy][yz] *= xynormi;
    rd[zz][yz] *= xynormi;
    }

 /* and form pure d trans matrices */

  for(i=0; i < si->g ; i++) {
#if 1
    math_double_mxm(si->Rd[i],1,_pd,0,si->Rdp[i],0,6,6,5,0);
#else
    rd = si->Rd[i];
    rpd = si->Rdp[i];

    rpd[xx][pzz]=-rd[zz][zz]/sqrt(6.0);
    rpd[yy][pzz]=-rd[zz][zz]/sqrt(6.0);
    rpd[zz][pzz]=2.0*rd[zz][zz]/sqrt(6.0);
    
    rpd[xx][xmy]=rd[xy][xy]/sqrt(2.0);
    rpd[yy][xmy]=-rd[xy][xy]/sqrt(2.0);
    rpd[xy][xmy]=-rd[yy][xy]*sqrt(3.0);

    rpd[xx][pxy]=rd[xy][xx]*sqrt(2.0)/sqrt(3.0);
    rpd[yy][pxy]=rd[xy][yy]*sqrt(2.0)/sqrt(3.0);
    rpd[xy][pxy]=rd[xy][xy];

    rpd[xz][pxz]=rd[xz][xz];
    rpd[yz][pyz]=rd[yz][yz];
    rpd[xz][pyz]=rd[yz][xz];
    rpd[yz][pxz]=rd[xz][yz];
#endif
    }
  

  if(f_exist) {
    double sq_3 = sqrt(3.0);
    double sq_15 = sqrt(15.0);
    double sq_1_3 = sqrt(1.0/3.0);
    double sq_1_15 = sqrt(1.0/15.0);
    double sq_3_15 = sqrt(3.0/15.0);
    double sq_15_3 = sqrt(15.0/3.0);

    for(i=0; i < si->g ; i++) {
      rp = si->Rp[i];
      rd = si->Rd[i];
      rf = si->Rf[i];

      for(m=0; m < 10 ; m++) {
        switch (m) {
        case yyy:
          j=y; k=y; l=y; jkl = yyy;
          break;
        case yyz:
          j=y; k=y; l=z; jkl = yyz;
          break;
        case yzz:
          j=y; k=z; l=z; jkl = yzz;
          break;
        case zzz:
          j=z; k=z; l=z; jkl = zzz;
          break;
        case xyy:
          j=x; k=y; l=y; jkl = xyy;
          break;
        case xyz:
          j=x; k=y; l=z; jkl = xyz;
          break;
        case xzz:
          j=x; k=z; l=z; jkl = xzz;
          break;
        case xxy:
          j=x; k=x; l=y; jkl = xxy;
          break;
        case xxz:
          j=x; k=x; l=z; jkl = xxz;
          break;
        case xxx:
          j=x; k=x; l=x; jkl = xxx;
          break; 
          }

        rf[jkl][xxx] = rp[j][x]*rp[k][x]*rp[l][x];
        rf[jkl][yyy] = rp[j][y]*rp[k][y]*rp[l][y];
        rf[jkl][zzz] = rp[j][z]*rp[k][z]*rp[l][z];

        rf[jkl][xxy] = rp[j][x]*rp[k][x]*rp[l][y]+
                       rp[j][x]*rp[k][y]*rp[l][x]+
                       rp[j][y]*rp[k][x]*rp[l][x];
        rf[jkl][xxz] = rp[j][x]*rp[k][x]*rp[l][z]+
                       rp[j][x]*rp[k][z]*rp[l][x]+
                       rp[j][z]*rp[k][x]*rp[l][x];
        rf[jkl][xyy] = rp[j][x]*rp[k][y]*rp[l][y]+
                       rp[j][y]*rp[k][x]*rp[l][y]+
                       rp[j][y]*rp[k][y]*rp[l][x];
        rf[jkl][yyz] = rp[j][z]*rp[k][y]*rp[l][y]+
                       rp[j][y]*rp[k][z]*rp[l][y]+
                       rp[j][y]*rp[k][y]*rp[l][z];
        rf[jkl][xzz] = rp[j][x]*rp[k][z]*rp[l][z]+
                       rp[j][z]*rp[k][x]*rp[l][z]+
                       rp[j][z]*rp[k][z]*rp[l][x];
        rf[jkl][yzz] = rp[j][y]*rp[k][z]*rp[l][z]+
                       rp[j][z]*rp[k][y]*rp[l][z]+
                       rp[j][z]*rp[k][z]*rp[l][y];
        rf[jkl][xyz] = rp[j][x]*rp[k][y]*rp[l][z]+
                       rp[j][z]*rp[k][x]*rp[l][y]+
                       rp[j][z]*rp[k][y]*rp[l][x]+
                       rp[j][y]*rp[k][x]*rp[l][z]+
                       rp[j][y]*rp[k][z]*rp[l][x]+
                       rp[j][x]*rp[k][z]*rp[l][y];
        }
      rf[xxx][xxy] *= sq_3_15;
      rf[xxx][xxz] *= sq_3_15;
      rf[xxx][xyy] *= sq_3_15;
      rf[xxx][yyz] *= sq_3_15;
      rf[xxx][xzz] *= sq_3_15;
      rf[xxx][yzz] *= sq_3_15;
      rf[xxx][xyz] *= sq_1_15;
      rf[yyy][xxy] *= sq_3_15;
      rf[yyy][xxz] *= sq_3_15;
      rf[yyy][xyy] *= sq_3_15;
      rf[yyy][yyz] *= sq_3_15;
      rf[yyy][xzz] *= sq_3_15;
      rf[yyy][yzz] *= sq_3_15;
      rf[yyy][xyz] *= sq_1_15;
      rf[zzz][xxy] *= sq_3_15;
      rf[zzz][xxz] *= sq_3_15;
      rf[zzz][xyy] *= sq_3_15;
      rf[zzz][yyz] *= sq_3_15;
      rf[zzz][xzz] *= sq_3_15;
      rf[zzz][yzz] *= sq_3_15;
      rf[zzz][xyz] *= sq_1_15;

      rf[xxy][xxx] *= sq_15_3;
      rf[xxz][xxx] *= sq_15_3;
      rf[xyy][xxx] *= sq_15_3;
      rf[yyz][xxx] *= sq_15_3;
      rf[xzz][xxx] *= sq_15_3;
      rf[yzz][xxx] *= sq_15_3;
      rf[xxy][yyy] *= sq_15_3;
      rf[xxz][yyy] *= sq_15_3;
      rf[xyy][yyy] *= sq_15_3;
      rf[yyz][yyy] *= sq_15_3;
      rf[xzz][yyy] *= sq_15_3;
      rf[yzz][yyy] *= sq_15_3;
      rf[xxy][zzz] *= sq_15_3;
      rf[xxz][zzz] *= sq_15_3;
      rf[xyy][zzz] *= sq_15_3;
      rf[yyz][zzz] *= sq_15_3;
      rf[xzz][zzz] *= sq_15_3;
      rf[yzz][zzz] *= sq_15_3;
      rf[xxy][xyz] *= sq_1_3;
      rf[xxz][xyz] *= sq_1_3;
      rf[xyy][xyz] *= sq_1_3;
      rf[yyz][xyz] *= sq_1_3;
      rf[xzz][xyz] *= sq_1_3;
      rf[yzz][xyz] *= sq_1_3;

      rf[xyz][xxx] *= sq_15;
      rf[xyz][yyy] *= sq_15;
      rf[xyz][zzz] *= sq_15;
      rf[xyz][xxy] *= sq_3;
      rf[xyz][xxz] *= sq_3;
      rf[xyz][xyy] *= sq_3;
      rf[xyz][yyz] *= sq_3;
      rf[xyz][xzz] *= sq_3;
      rf[xyz][yzz] *= sq_3;
      }
    }

  free(rot);
  free(trans);

  return(0);
  }

LOCAL_FUNCTION int
rperr(_outfile,i)
FILE *_outfile;
int i;
{
  fprintf(_outfile,"sym_init_centers: make_rp_d:\n");
  fprintf(_outfile,"could not allocate memory for gamma[%d]\n",i);
  fflush(_outfile);
  return(-1);
  }
