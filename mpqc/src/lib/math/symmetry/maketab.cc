
/* maketab.cc --
 *
 *      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
 *      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
 *      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
 *      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
 *      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
 *      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
 *
 *  Author:
 *      E. T. Seidl
 *      Bldg. 12A, Rm. 2033
 *      Computer Systems Laboratory
 *      Division of Computer Research and Technology
 *      National Institutes of Health
 *      Bethesda, Maryland 20892
 *      Internet: seidl@alw.nih.gov
 *      June, 1993
 */

#include <math.h>
#include <string.h>

#include <math/symmetry/pointgrp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * This function will generate a character table for the point group.
 * This character table is in the order that symmetry operations are
 * generated, not in Cotton order. If this is a problem, tough.
 * Also generate the transformation matrices.
 */

int CharacterTable::make_table()
{
  int i,j,i0;
  char label[4];

  if (!g) return 0;

  gamma_ = new IrreducibleRepresentation[nirrep_];

  for (i=0; i < nirrep_; i++)
    gamma_[i].init();

  symop = new SymmetryOperation[g];

  //for (i=0; i < g; i++) {
  //  symop[i].zero();
  //  }

 // this array forms a reducible representation for rotations about x,y,z
  double *rot = new double[g];
  memset(rot,'\0',sizeof(double)*g);

 // this array forms a reducible representation for translations along x,y,z
  double *trans = new double[g];
  memset(trans,'\0',sizeof(double)*g);

 // the angle to rotate about the principal axis
  double theta = (nt) ? 2.0*M_PI/nt : 2.0*M_PI;

  switch (pg) {
  case C1: // no symmetry
    { IrreducibleRepresentation ir(1,1,"A"); gamma_[0] = ir; }

    gamma_[0].nrot_ = 3; gamma_[0].ntrans_ = 3;
    gamma_[0].rep[0] = 1;

    symop[0](0,0) = symop[0](1,1) = symop[0](2,2) = 1;
    break;

  case CI:  // equivalent to S2 about the z axis
    { IrreducibleRepresentation ir(2,1,"Ag"); gamma_[0] = ir; }
    gamma_[0].rep[0] = 1; gamma_[0].rep[1] = 1; gamma_[0].nrot_=3;

    symop[0](0,0) = symop[0](1,1) = symop[0](2,2) = 1;

    { IrreducibleRepresentation ir(2,1,"Au"); gamma_[1] = ir; }
    gamma_[1].rep[0] = 1; gamma_[1].rep[1] = -1; gamma_[1].ntrans_=3;

    symop[1](0,0) = symop[1](1,1) = symop[1](2,2) = -1;

    break;

  case CS: // reflection through the xy plane
    { IrreducibleRepresentation ir(2,1,"A'"); gamma_[0] = ir; }
    gamma_[0].rep[0] = 1; gamma_[0].rep[1] = 1;
    gamma_[0].nrot_=1; gamma_[0].ntrans_=2;

    symop[0](0,0) = symop[0](1,1) = symop[0](2,2) = 1;

    { IrreducibleRepresentation ir(2,1,"A\""); gamma_[1] = ir; }
    gamma_[1].rep[0] = 1; gamma_[1].rep[1] = -1;
    gamma_[1].nrot_=2; gamma_[1].ntrans_=1;

    symop[1](0,0) = symop[1](1,1) = 1;
    symop[1](2,2) = -1;

    break;

  case CN: // clockwise rotation about z axis by theta*i radians
    { IrreducibleRepresentation ir(g,1,"A"); gamma_[0] = ir; }
    for (i=0; i < g; i++) gamma_[0].rep[i] = 1;

    if (nt%2) {
      for (i=1; i < nirrep_; i++) {
        if (nt==3) sprintf(label,"E");
        else sprintf(label,"E%d",i);

	IrreducibleRepresentation ir(g,2,label);

        ir.rep[0] = 2.0;
        for (j=1; j < g; j++) ir.rep[j] = 2.0*cos((double)theta*j*i);

        gamma_[i] = ir;
        }
      }
    else {
      { IrreducibleRepresentation ir(g,1,"B"); gamma_[1] = ir; }
      for (i=0; i < g; i++) gamma_[1].rep[i] = pow(-1.0,(double)i);

      for (i=2; i < nirrep_; i++) {
        if (nt==4) sprintf(label,"E");
        else sprintf(label,"E%d",i-1);

	IrreducibleRepresentation ir(g,2,label);

        for (j=0; j < g; j++) ir.rep[j] = 2.0*cos((double)theta*j*(i-1));

        gamma_[i] = ir;
        }
      }

    for (i=0; i < nt ; i++) {
      symop[i][0][0] = symop[i][1][1] = cos((double)theta*i);
      symop[i][0][1] = sin((double)theta*i);
      symop[i][1][0] = -sin((double)theta*i);
      symop[i][2][2] = 1.0;

      rot[i] = trans[i] = symop[i].trace();
      }

    break;

  case CNV: // clockwise rotation about z axis by theta*i radians, then
            // reflect through the xz plane
    { IrreducibleRepresentation ir(g,1,"A1"); gamma_[0] = ir; }
    for (i=0; i < g; i++) gamma_[0].rep[i] = 1;

    { IrreducibleRepresentation ir(g,1,"A2"); gamma_[1] = ir; }
    for (i=0; i < g/2; i++) gamma_[1].rep[i] = 1;
    for (; i < g; i++) gamma_[1].rep[i] = -1;

    if (nt%2) { // cnv with odd nt has only a1, a2, and e's
      for (i=2; i < nirrep_; i++) {
        char lab[4];
        if (nt==3) sprintf(lab,"E");
        else sprintf(lab,"E%d",i-1);

        IrreducibleRepresentation ir(g,2,lab);

        ir.rep[0] = 2;
        for (j=0; j < g/2; j++)
          ir.rep[j] = 2.0*cos((double)theta*j*(i-1));
        for (; j < g; j++)
          ir.rep[j] = 0;

        gamma_[i] = ir;
        }
      }
    else {  // ok, we've got a1, a2, b1, b2, and e's
      { IrreducibleRepresentation ir(g,1,"B1"); gamma_[2] = ir; }
      for (i=0; i < g; i++) gamma_[2].rep[i] = pow(-1.0,(double)i);

      { IrreducibleRepresentation ir(g,1,"B2"); gamma_[3] = ir; }
      for (i=0; i < g/2 ; i++) gamma_[3].rep[i] = pow(-1.0,(double)i);
      for (; i < g ; i++) gamma_[3].rep[i] = pow(-1.0,(double)(i+1));

      for (i=4; i < nirrep_; i++) {
        char lab[4];
        if (nt==4) sprintf(lab,"E");
        else sprintf(lab,"E%d",i-3);

        IrreducibleRepresentation ir(g,2,lab);

        ir.rep[0] = 2;
        for (j=0; j < g/2; j++)
          ir.rep[j] = 2.0*cos((double)theta*j*(i-3));
        for (; j < g; j++)
          ir.rep[j] = 0;

        gamma_[i] = ir;
        }
      }

    for (i=0; i < nt ; i++) {
      symop[i][0][0] = symop[i][1][1] = cos((double)theta*i);
      symop[i][0][1] = sin((double)theta*i);
      symop[i][1][0] = -sin((double)theta*i);
      symop[i][2][2] = 1.0;

      rot[i] = trans[i] = symop[i].trace();
      }

    for (i=0; i < nt ; i++) {
      symop[i+nt][0][0] = cos((double)theta*i);
      symop[i+nt][1][1] = -cos((double)theta*i);
      symop[i+nt][1][0] = symop[i+nt][0][1] = sin((double)theta*i);
      symop[i+nt][2][2] = 1.0;

      trans[i+nt] = symop[i+nt].trace();
      rot[i+nt] = -trans[i+nt];
      }

    break;

  case CNH: // lockwise rotation about z axis by theta*i radians,
            // as well as rotation-reflection about same axis

    if (nt%2) {
      { IrreducibleRepresentation ir(g,1,"A'"); gamma_[0] = ir; }
      for (i=0; i < g ; i++) gamma_[0].rep[i] = 1;

      for (i=1; i < nirrep_/2 ; i++) {
        if (nt==3) sprintf(label,"E'");
        else sprintf(label,"E%d'",i);

        IrreducibleRepresentation ir(g,2,label);

        ir.rep[0] = 2.0;

        for (j=1; j < g/2; j++) ir.rep[j] = 2.0*cos((double)theta*j*i);
        for (; j < g; j++) ir.rep[j] = ir.rep[j-(g/2)];

        gamma_[i] = ir;
        }

      i=nirrep_/2;
      { IrreducibleRepresentation ir(g,1,"A\""); gamma_[i] = ir; }
      for (j=0; j < g/2; j++) gamma_[i].rep[j] = 1.0;
      for ( ; j < g; j++) gamma_[i].rep[j] = -1.0;

      for (i=i+1; i < nirrep_ ; i++) {
        i0 = i-(nirrep_/2);
        if (nt==3) sprintf(label,"E\"");
        else sprintf(label,"E%d\"",i0);

        IrreducibleRepresentation ir(g,2,label);

        for (j=0; j < g/2; j++) ir.rep[j] = gamma_[i0].rep[j];
        for (j=g/2; j < g ; j++) ir.rep[j] = -gamma_[i0].rep[j];

        gamma_[i] = ir;
        }
      }
    else {
      double ineg = pow(-1.0,(double)nt/2);

      { IrreducibleRepresentation ir(g,1,"Ag"); gamma_[0] = ir; }
      for (i=0; i < g; i++) gamma_[0].rep[i] = 1.0;

      { IrreducibleRepresentation ir(g,1,"Bg"); gamma_[1] = ir; }
      for (i=0; i < g/2; i++) gamma_[1].rep[i] = pow(-1.0,(double) i);
      for (; i < g; i++) gamma_[1].rep[i] = ineg*pow(-1.0,(double) i);

      for (i=2; i < nirrep_/2; i++) {
        ineg=pow(-1.0,(double) (i-1));

        if (nt==4) sprintf(label,"Eg");
        else sprintf(label,"E%dg",i-1);

        IrreducibleRepresentation ir(g,2,label);

        ir.rep[0] = 2.0;
        for (j=1; j < g/2; j++) ir.rep[j] = 2.0*cos((double)theta*j*(i-1));
        for (; j < g; j++) ir.rep[j] = ineg*ir.rep[j-(g/2)];

        gamma_[i] = ir;
        }

      i=nirrep_/2;
      { IrreducibleRepresentation ir(g,1,"Au"); gamma_[i] = ir; }
      for (j=0; j < g/2; j++) gamma_[i].rep[j] = 1.0;
      for (j=g/2; j < g; j++) gamma_[i].rep[j] = -1.0;

      i++;
      { IrreducibleRepresentation ir(g,1,"Bu"); gamma_[i] = ir; }
      for (j=0; j < g/2; j++) gamma_[i].rep[j] = gamma_[1].rep[j];
      for (j=g/2; j < g; j++) gamma_[i].rep[j] = -gamma_[1].rep[j];

      for (i=i+1; i < nirrep_ ; i++) {
        i0=i-(nirrep_/2);

        if (nt==4) sprintf(label,"Eu");
        else sprintf(label,"E%du",i0-1);

        IrreducibleRepresentation ir(g,2,label);

        for (j=0; j < g/2; j++) ir.rep[j] =  gamma_[i0].rep[j];
        for (j=g/2; j < g; j++) ir.rep[j] = -gamma_[i0].rep[j];

        gamma_[i] = ir;
        }
      }

    for (i=0; i < nt ; i++) {
      symop[i][0][0] = symop[i][1][1] =
      symop[i+nt][0][0] = symop[i+nt][1][1] = cos((double)theta*i);
      symop[i][0][1] = symop[i+nt][0][1] = sin((double)theta*i);
      symop[i][1][0] = symop[i+nt][1][0] = -sin((double)theta*i);
      symop[i][2][2] = 1.0;
      symop[i+nt][2][2] = -1.0;

      rot[i] = trans[i] = symop[i].trace();
      trans[i+nt] = symop[i+nt].trace();
      rot[i+nt] = -trans[i+nt];
      }

    break;

  case SN: // clockwise rotation-reflection by theta*i radians about z axis

    if ((nt/2)%2) {
      { IrreducibleRepresentation ir(g,1,"Ag"); gamma_[0] = ir; }
      for (i=0; i < g ; i++) gamma_[0].rep[i] = 1.0;

      for (i=1; i < nirrep_/2 ; i++) {
        if (nt==6) sprintf(label,"Eg");
        else sprintf(label,"E%dg",i);

        IrreducibleRepresentation ir(g,2,label);

        ir.rep[0] = 2.0;
        for (j=1; j < g ; j++) ir.rep[j] = 2.0*cos((double)2.0*theta*j*i);

        gamma_[i] = ir;
        }

      { IrreducibleRepresentation ir(g,1,"Au"); gamma_[i] = ir; }
      for (j=0; j < g ; j++) gamma_[i].rep[j] = pow(-1.0,(double)j);

      for (i=i+1; i < nirrep_ ; i++) {
        i0=i-nirrep_/2;

        if (nt==6) sprintf(label,"Eu");
        else sprintf(label,"E%du",i);

        IrreducibleRepresentation ir(g,2,label);

        for (j=0; j < g; j++)
          ir.rep[j] = pow(-1.0,(double) j)*gamma_[i0].rep[j];

        gamma_[i] = ir;
        }
      }
    else {
      { IrreducibleRepresentation ir(g,1,"A"); gamma_[0] = ir; }
      for (i=0; i < g; i++) gamma_[0].rep[i] = 1.0;

      { IrreducibleRepresentation ir(g,1,"B"); gamma_[1] = ir; }
      for (i=0; i < g; i++) gamma_[1].rep[i] = pow(-1.0,(double)i);

      for (i=2; i < nirrep_ ; i++) {
        if (nt==4) sprintf(label,"E");
        else sprintf(label,"E%d",i-1);

        IrreducibleRepresentation ir(g,2,label);

        for (j=0; j < g; j++) ir.rep[j] = 2.0*cos((double)theta*j*(i-1));

        gamma_[i] = ir;
        }
      }

    for (i=0; i < nt ; i++) {
      symop[i][0][0] = symop[i][1][1] = cos((double)theta*i);
      symop[i][0][1] = sin((double)theta*i);
      symop[i][1][0] = -sin((double)theta*i);
      symop[i][2][2] = pow(-1.0,(double) i);

      trans[i] = symop[i].trace();
      rot[i] = (i%2) ? -trans[i] : trans[i];
      }

    break;

  case DN: // clockwise rotation about z axis, followed by C2 about
           // x axis

    if (nt==2) { // D2 is a special case
      { IrreducibleRepresentation ir(g,1,"A"); gamma_[0] = ir; }
      for (i=0; i < g; i++) gamma_[0].rep[i] = 1.0;

      { IrreducibleRepresentation ir(g,1,"B1"); gamma_[1] = ir; }
      gamma_[1].rep[0] = gamma_[1].rep[1] = 1.0;
      gamma_[1].rep[2] = gamma_[1].rep[3] = -1.0;

      { IrreducibleRepresentation ir(g,1,"B2"); gamma_[2] = ir; }
      gamma_[2].rep[0] = gamma_[2].rep[3] = 1.0;
      gamma_[2].rep[1] = gamma_[2].rep[2] = -1.0;

      { IrreducibleRepresentation ir(g,1,"B3"); gamma_[3] = ir; }
      gamma_[3].rep[0] = gamma_[3].rep[2] = 1.0;
      gamma_[3].rep[1] = gamma_[3].rep[3] = -1.0;
      }
    else if (nt%2) {
      { IrreducibleRepresentation ir(g,1,"A1"); gamma_[0] = ir; }
      for (i=0; i < g; i++) gamma_[0].rep[i] = 1.0;

      { IrreducibleRepresentation ir(g,1,"A2"); gamma_[1] = ir; }
      for (i=0; i < g/2; i++) gamma_[1].rep[i] = 1.0;
      for (i=g/2; i < g; i++) gamma_[1].rep[i] = -1.0;

      for (i=2; i < nirrep_ ; i++) {
        if (nt==3) sprintf(label,"E");
        else sprintf(label,"E%d",i-1);

        IrreducibleRepresentation ir(g,2,label);

        ir.rep[0] = 2.0;
        for (j=1; j < g/2; j++) ir.rep[j] = 2.0*cos((double)theta*j*(i-1));
        for (j=g/2; j < g; j++) ir.rep[j] = 0.0;

        gamma_[i] = ir;
        }
      }
    else {
      { IrreducibleRepresentation ir(g,1,"A1"); gamma_[0] = ir; }
      for (i=0; i < g; i++) gamma_[0].rep[i] = 1.0;

      { IrreducibleRepresentation ir(g,1,"A2"); gamma_[1] = ir; }
      for (i=0; i < g/2; i++) gamma_[1].rep[i] = 1.0;
      for (i=g/2; i < g; i++) gamma_[1].rep[i] = -1.0;

      { IrreducibleRepresentation ir(g,1,"B1"); gamma_[2] = ir; }
      for (i=0; i < g; i++) gamma_[2].rep[i] = pow(-1.0,(double)i);

      { IrreducibleRepresentation ir(g,1,"B2"); gamma_[3] = ir; }
      for (i=0; i < g/2; i++) gamma_[3].rep[i] = pow(-1.0,(double)i);
      for (i=g/2; i < g; i++) gamma_[3].rep[i] = pow(-1.0,(double)(i+1));

      for (i=4; i < nirrep_; i++) {
        if (nt==4) sprintf(label,"E");
        else sprintf(label,"E%d",i-3);

        IrreducibleRepresentation ir(g,2,label);

        for (j=0; j < g/2; j++) ir.rep[j] = 2.0*cos((double)theta*j*(i-3));
        for (j=g/2; j < g; j++) ir.rep[j] = 0.0;

        gamma_[i] = ir;
        }
      }

    for (i=0; i < nt ; i++) {
      symop[i][0][0] = symop[i][1][1] = cos((double) theta*i);
      symop[i][0][1] = sin((double) theta*i);
      symop[i][1][0] = -sin((double) theta*i);
      symop[i][2][2] = 1.0;

      rot[i] = trans[i] = symop[i].trace();
      }

    for (i=0,j=nt; i < nt ; i++,j++) {
      symop[j][0][0] = cos((double)theta*i);
      symop[j][1][1] = -cos((double)theta*i);
      symop[j][1][0] = symop[j][0][1] = -sin((double)theta*i);
      symop[j][2][2] = -1.0;

      rot[j] = trans[j] = symop[j].trace();
      }

    break;

  case DND: // rotation reflection about z axis by theta/2 radians, followed
            // by c2 about x axis, then reflection through yz plane
    if (nt%2) {
      { IrreducibleRepresentation ir(g,1,"A1g"); gamma_[0] = ir; }
      for (i=0; i < g; i++) gamma_[0].rep[i] = 1.0;

      { IrreducibleRepresentation ir(g,1,"A2g"); gamma_[1] = ir; }
      for (i=0; i < g/2; i++) gamma_[1].rep[i] = 1.0;
      for (i=g/2; i < g; i++) gamma_[1].rep[i] = -1.0;

      for (i=2; i < nirrep_/2 ; i++) {
        if (nt==3) sprintf(label,"Eg");
        else sprintf(label,"E%dg",i-1);

        IrreducibleRepresentation ir(g,2,label);

        ir.rep[0] = 2.0;
        for (j=1; j < g/2; j++) ir.rep[j] = 2.0*cos((double)2.0*theta*j*(i-1));
        for (j=g/2; j < g; j++) ir.rep[j] = 0.0;

        gamma_[i] = ir;
        }

      { IrreducibleRepresentation ir(g,1,"A1u"); gamma_[i] = ir; i++; }
      { IrreducibleRepresentation ir(g,1,"A2u"); gamma_[i] = ir; i++; }

      for (; i < nirrep_; i++) {
        i0=i-nirrep_/2;
        
        if (nt==3) sprintf(label,"Eu");
        else sprintf(label,"E%du",i0-1);

        IrreducibleRepresentation ir(g,2,label);
        gamma_[i] = ir;
        }

      for (i=nirrep_/2,i0=0; i < nirrep_; i++,i0++) {
        for (j=0; j < g/2; j++) gamma_[i].rep[j] =
                              pow(-1.0,(double) j)*gamma_[i0].rep[j];
        for (; j < 3*g/4; j++)  gamma_[i].rep[j] =  gamma_[i0].rep[j];
        for (; j < g    ; j++)  gamma_[i].rep[j] = -gamma_[i0].rep[j];
        }
      }
    else {
      { IrreducibleRepresentation ir(g,1,"A1"); gamma_[0] = ir; }
      for (i=0; i < g; i++) gamma_[0].rep[i] = 1.0;

      { IrreducibleRepresentation ir(g,1,"A2"); gamma_[1] = ir; }
      for (i=0; i < g/2; i++) gamma_[1].rep[i] = 1.0;
      for (i=g/2; i < g; i++) gamma_[1].rep[i] = -1.0;

      { IrreducibleRepresentation ir(g,1,"B1"); gamma_[2] = ir; }
      for (i=0; i < g/2; i++) gamma_[2].rep[i] = pow(-1.0,(double) i);
      for ( ; i < 3*g/4; i++) gamma_[2].rep[i] =  1.0;
      for ( ; i < g    ; i++) gamma_[2].rep[i] = -1.0;

      { IrreducibleRepresentation ir(g,1,"B2"); gamma_[3] = ir; }
      for (i=0; i < g/2; i++) gamma_[3].rep[i] = pow(-1.0,(double) i);
      for ( ; i < 3*g/4; i++) gamma_[3].rep[i] = -1.0;
      for ( ; i < g    ; i++) gamma_[3].rep[i] = 1.0;

      for (i=4; i < nirrep_; i++) {
        if (nt==2) sprintf(label,"E");
        else sprintf(label,"E%d",i-3);

        IrreducibleRepresentation ir(g,2,label);

        ir.rep[0] = 2.0;
        for (j=1; j < g/2; j++) ir.rep[j] = 2.0*cos((double)theta*j*(i-3)*0.5);
        for (j=g/2; j < g; j++) ir.rep[j] = 0.0;

        gamma_[i] = ir;
        }
      }

    for (i=0; i < 2*nt ; i++) {
      symop[i][0][0] = symop[i][1][1] = cos((double)theta*i*0.5);
      symop[i][0][1] = sin((double)theta*i*0.5);
      symop[i][1][0] = -sin((double)theta*i*0.5);
      symop[i][2][2] = pow(-1.0,(double) i);

      trans[i] = symop[i].trace();
      rot[i] = (i%2) ? -trans[i] : trans[i];
      }

    for (i=0,j=2*nt; i < nt ; i++,j++) {
      symop[j][0][0] = cos((double)theta*i);
      symop[j][1][1] = -cos((double)theta*i);
      symop[j][1][0] = symop[j][0][1] = -sin((double)theta*i);
      symop[j][2][2] = -1.0;

      rot[j] = trans[j] = symop[j].trace();
      }

    for (i=0,j=3*nt; i < nt ; i++,j++) {
      symop[j][0][0] = cos((double)theta*i+theta*0.5);
      symop[j][1][1] = -cos((double)theta*i+theta*0.5);
      symop[j][1][0] = symop[j][0][1] = -sin((double)theta*i+theta*0.5);
      symop[j][2][2] = 1.0;

      trans[j] = symop[j].trace();
      rot[j] = -trans[j];
      }

    break;

  case DNH: // clockwise rotation and rotation-reflection about z axis,
            // followed by c2 about x axis and then reflection
            // through xz 

    if (nt==2) { // d2h is a special case
      { IrreducibleRepresentation ir(g,1,"Ag"); gamma_[0] = ir; }
      for (i=0; i < 8; i++) gamma_[0].rep[i] = 1.0;

      { IrreducibleRepresentation ir(g,1,"B1g"); gamma_[1] = ir; }
      for (i=0; i < 4; i++) gamma_[1].rep[i] = 1.0;
      for (i=4; i < 8; i++) gamma_[1].rep[i] = -1.0;

      { IrreducibleRepresentation ir(g,1,"B2g"); gamma_[2] = ir; }
      gamma_[2].rep[0] = gamma_[2].rep[3] =
                         gamma_[2].rep[5] = gamma_[2].rep[6] = 1.0;
      gamma_[2].rep[1] = gamma_[2].rep[2] =
                         gamma_[2].rep[4] = gamma_[2].rep[7] = -1.0;

      { IrreducibleRepresentation ir(g,1,"B3g"); gamma_[3] = ir; }
      gamma_[3].rep[0] = gamma_[3].rep[3] =
                         gamma_[3].rep[4] = gamma_[3].rep[7] = 1.0;
      gamma_[3].rep[1] = gamma_[3].rep[2] =
                         gamma_[3].rep[5] = gamma_[3].rep[6] = -1.0;

      { IrreducibleRepresentation ir(g,1,"Au"); gamma_[4] = ir; }
      for (i=0; i < 2; i++) gamma_[4].rep[i] = 1.0;
      for (i=2; i < 4; i++) gamma_[4].rep[i] = -1.0;
      for (i=4; i < 6; i++) gamma_[4].rep[i] = 1.0;
      for (i=6; i < 8; i++) gamma_[4].rep[i] = -1.0;

      { IrreducibleRepresentation ir(g,1,"B1u"); gamma_[5] = ir; }
      for (i=0; i < 2; i++) gamma_[5].rep[i] = 1.0;
      for (i=2; i < 6; i++) gamma_[5].rep[i] = -1.0;
      for (i=6; i < 8; i++) gamma_[5].rep[i] = 1.0;

      { IrreducibleRepresentation ir(g,1,"B2u"); gamma_[6] = ir; }
      for (i=0; i < 4; i++) gamma_[6].rep[i] =  pow(-1.0,(double)i);
      for (i=4; i < 8; i++) gamma_[6].rep[i] = -pow(-1.0,(double)i);

      { IrreducibleRepresentation ir(g,1,"B3u"); gamma_[7] = ir; }
      for (i=0; i < 8; i++) gamma_[7].rep[i] = pow(-1.0,(double)i);

      }

    else if (nt%2) {
      { IrreducibleRepresentation ir(g,1,"A1'"); gamma_[0] = ir; }
      for (i=0; i < g; i++) gamma_[0].rep[i] = 1.0;

      { IrreducibleRepresentation ir(g,1,"A2'"); gamma_[1] = ir; }
      for (i=0; i < g/2; i++) gamma_[1].rep[i] = 1.0;
      for (i=g/2; i < g; i++) gamma_[1].rep[i] = -1.0;

      for (i=2; i < nirrep_/2 ; i++) {
        if (nt==3) sprintf(label,"E'");
        else sprintf(label,"E%d'",i-1);

        IrreducibleRepresentation ir(g,2,label);

        ir.rep[0] = 2.0;
        for (j=1; j < g/4; j++) ir.rep[j] = 2.0*cos((double)theta*j*(i-1));
        for (   ; j < g/2; j++) ir.rep[j] = 2.0*cos((double)theta*j*(i-1));
        for (   ; j < g  ; j++) ir.rep[j] = 0.0;

        gamma_[i] = ir;
        }

      { IrreducibleRepresentation ir(g,1,"A1\""); gamma_[i] = ir; i++; }
      { IrreducibleRepresentation ir(g,1,"A2\""); gamma_[i] = ir; i++; }

      for ( ; i < nirrep_; i++) {
        i0=i-(nirrep_/2);

        if (nt==3) sprintf(label,"E\"");
        else sprintf(label,"E%d\"",i0-1);

        IrreducibleRepresentation ir(g,2,label);

        gamma_[i] = ir;
        }

      for (i=nirrep_/2,i0=0; i < nirrep_; i0++,i++) {
        for (j=0; j < g/4; j++) gamma_[i].rep[j] =  gamma_[i0].rep[j];
        for (   ; j < g/2; j++) gamma_[i].rep[j] = -gamma_[i0].rep[j];
        for ( ; j < 3*g/4; j++) gamma_[i].rep[j] =  gamma_[i0].rep[j];
        for ( ; j < g    ; j++) gamma_[i].rep[j] = -gamma_[i0].rep[j];
        }
      }

    else {
      { IrreducibleRepresentation ir(g,1,"A1g"); gamma_[0] = ir; }
      for (i=0; i < g; i++) gamma_[0].rep[i] = 1.0;

      { IrreducibleRepresentation ir(g,1,"A2g"); gamma_[1] = ir; }
      for (i=0; i < g/2; i++) gamma_[1].rep[i] = 1.0;
      for (i=g/2; i < g; i++) gamma_[1].rep[i] = -1.0;

      { IrreducibleRepresentation ir(g,1,"B1g"); gamma_[2] = ir; }
      for (i=0; i < g; i++) gamma_[2].rep[i] = pow(-1.0,(double)i);

      { IrreducibleRepresentation ir(g,1,"B2g"); gamma_[3] = ir; }
      for (i=0; i < g/2; i++) gamma_[3].rep[i] = pow(-1.0,(double)i);
      for (i=g/2; i < g; i++) gamma_[3].rep[i] = pow(-1.0,(double)(i+1));

      for (i=4; i < nirrep_/2; i++) {
        if (nt==4) sprintf(label,"Eg");
        else sprintf(label,"E%dg",i-3);

        IrreducibleRepresentation ir(g,2,label);

        for (j=0; j < g/4 ; j++) ir.rep[j] = 2.0*cos((double)theta*j*(i-3));
        for (j=g/4; j < g/2; j++)
          ir.rep[j] = pow(-1.0,(double)(i+1))*2.0*cos((double)theta*j*(i-3));
        for (j=g/2; j < g; j++) ir.rep[j] = 0.0;

        gamma_[i] = ir;
        }

      { IrreducibleRepresentation ir(g,1,"A1u"); gamma_[i] = ir; i++; }
      { IrreducibleRepresentation ir(g,1,"A2u"); gamma_[i] = ir; i++; }
      { IrreducibleRepresentation ir(g,1,"B1u"); gamma_[i] = ir; i++; }
      { IrreducibleRepresentation ir(g,1,"B2u"); gamma_[i] = ir; i++; }

      for (; i < nirrep_; i++) {
        i0=i-(nirrep_/2);

        if (nt==4) sprintf(label,"Eu");
        else sprintf(label,"E%du",i0-3);

        IrreducibleRepresentation ir(g,2,label);
        gamma_[i] = ir;
        }

      for (i=nirrep_/2,i0=0; i < nirrep_; i0++,i++) {
        for (j=0; j < g/4; j++) gamma_[i].rep[j] =  gamma_[i0].rep[j];
        for (   ; j < g/2; j++) gamma_[i].rep[j] = -gamma_[i0].rep[j];
        for ( ; j < 3*g/4; j++) gamma_[i].rep[j] =  gamma_[i0].rep[j];
        for (   ; j < g  ; j++) gamma_[i].rep[j] = -gamma_[i0].rep[j];
        }
      }

    for (i=0; i < nt ; i++) {
      symop[i][0][0] = symop[i][1][1] =
      symop[i+nt][0][0] = symop[i+nt][1][1] = cos((double)theta*i);
      symop[i][0][1] = symop[i+nt][0][1] = sin((double)theta*i);
      symop[i][1][0] = symop[i+nt][1][0] = -sin((double)theta*i);
      symop[i][2][2] = 1.0;
      symop[i+nt][2][2] = -1.0;

      rot[i] = trans[i] = symop[i].trace();
      trans[i+nt] = symop[i+nt].trace();
      rot[i+nt] = -trans[i+nt];
      }

    for (i=0,j=2*nt; i < nt ; i++,j++) {
      symop[j][0][0] = cos((double)theta*i);
      symop[j][1][1] = -cos((double)theta*i);
      symop[j][1][0] = symop[j][0][1] = -sin((double)theta*i);
      symop[j][2][2] = -1.0;

      rot[j] = trans[j] = symop[j].trace();
      }

    for (i=0,j=3*nt; i < nt ; i++,j++) {
      symop[j][0][0] = cos((double)theta*i);
      symop[j][1][1] = -cos((double)theta*i);
      symop[j][1][0] = symop[j][0][1] = sin((double)theta*i);
      symop[j][2][2] = 1.0;

      trans[j] = symop[j].trace();
      rot[j] = -trans[j];
      }

    break;

  default:
    return -1;

    }
    
/* ok, we have the reducible representation of the rotations and translations,
 * now let's use projection operators to find out how many rotations and
 * translations there are in each irrep
 */

  if (pg != C1 && pg != CI && pg != CS ) {
    for (i=0; i < nirrep_; i++) {
      double nr=0; double nt=0;

      for (j=0; j < gamma_[i].g; j++) {
        nr += rot[j]*gamma_[i].rep[j];
        nt += trans[j]*gamma_[i].rep[j];
        }

      if (pg==CN || pg==CNH || pg==SN) {
        gamma_[i].nrot_ = (int) ((nr+0.5)/(gamma_[i].g*gamma_[i].degen));
        gamma_[i].ntrans_ = (int) ((nt+0.5)/(gamma_[i].g*gamma_[i].degen));
        }
      else {
        gamma_[i].nrot_ = (int) ((nr+0.5)/gamma_[i].g);
        gamma_[i].ntrans_ = (int) ((nt+0.5)/gamma_[i].g);
        }
      }
    }

  return 0;
}
