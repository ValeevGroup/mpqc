
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
    {
      IrreducibleRepresentation ir(1,1,"A");
      ir.nrot_ = 3;
      ir.ntrans_ = 3;
      ir.rep[0][0] = 1;

      gamma_[0] = ir;
    }

    symop[0](0,0) = symop[0](1,1) = symop[0](2,2) = 1;
    break;

  case CI:  // equivalent to S2 about the z axis
    {
      IrreducibleRepresentation ir(2,1,"Ag");
      ir.rep[0][0] = 1;
      ir.rep[0][1] = 1;
      ir.nrot_=3;

      gamma_[0] = ir;
    }

    symop[0](0,0) = symop[0](1,1) = symop[0](2,2) = 1;

    {
      IrreducibleRepresentation ir(2,1,"Au");
      ir.rep[0][0] = 1;
      ir.rep[0][1] = -1;
      ir.ntrans_=3;

      gamma_[1] = ir;
    }

    symop[1](0,0) = symop[1](1,1) = symop[1](2,2) = -1;

    break;

  case CS: // reflection through the xy plane
    {
      IrreducibleRepresentation ir(2,1,"A'");
      ir.rep[0][0] = 1;
      ir.rep[0][1] = 1;
      ir.nrot_=1;
      ir.ntrans_=2;

      gamma_[0] = ir;
    }

    symop[0](0,0) = symop[0](1,1) = symop[0](2,2) = 1;

    {
      IrreducibleRepresentation ir(2,1,"A\"");
      ir.rep[0][0] = 1;
      ir.rep[0][1] = -1;
      ir.nrot_=2;
      ir.ntrans_=1;

      gamma_[1] = ir;
    }

    symop[1](0,0) = symop[1](1,1) = 1;
    symop[1](2,2) = -1;

    break;

  case CN: // clockwise rotation about z axis by theta*i radians
    {
      IrreducibleRepresentation ir(g,1,"A");
      for (i=0; i < g; i++)
        ir.rep[0][i] = 1;

      gamma_[0] = ir;
    }

    if (nt%2) {
      for (i=1; i < nirrep_; i++) {
        if (nt==3)
          sprintf(label,"E");
        else
          sprintf(label,"E%d",i);

	IrreducibleRepresentation ir(g,2,label);

        ir.rep[0][0] = 1.0;
        ir.rep[1][0] = 0.0;
        for (j=1; j < g; j++) {
          ir.rep[0][j] = cos((double)theta*j*i);
          ir.rep[1][j] = sin((double)theta*j*i);
        }

        gamma_[i] = ir;
      }
    } else {
      {
        IrreducibleRepresentation ir(g,1,"B");
        for (i=0; i < g; i++)
          ir.rep[0][i] = pow(-1.0,(double)i);

        gamma_[1] = ir;
      }

      for (i=2; i < nirrep_; i++) {
        if (nt==4)
          sprintf(label,"E");
        else
          sprintf(label,"E%d",i-1);

	IrreducibleRepresentation ir(g,2,label);

        for (j=0; j < g; j++) {
          ir.rep[0][j] = cos((double)theta*j*(i-1));
          ir.rep[1][j] = sin((double)theta*j*(i-1));
        }

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
    {
      IrreducibleRepresentation ir(g,1,"A1");
      for (i=0; i < g; i++)
        ir.rep[0][i] = 1;

      gamma_[0] = ir;
    }

    {
      IrreducibleRepresentation ir(g,1,"A2");
      for (i=0; i < g/2; i++)
        ir.rep[0][i] = 1;
      for (; i < g; i++)
        ir.rep[0][i] = -1;

      gamma_[1] = ir;
    }

    if (nt%2) { // cnv with odd nt has only a1, a2, and e's
      for (i=2; i < nirrep_; i++) {
        char lab[4];
        if (nt==3)
          sprintf(lab,"E");
        else
          sprintf(lab,"E%d",i-1);

        IrreducibleRepresentation ir(g,2,lab);

        for (j=0; j < g/2; j++) {
          ir.rep[0][j] = 2.0*cos((double)theta*j*(i-1));
          ir.rep[1][j] = 2.0*sin((double)theta*j*(i-1));
        }

        for (; j < g; j++) {
          ir.rep[0][j] = 0;
          ir.rep[1][j] = 0;
        }

        gamma_[i] = ir;
        }
      }
    else {  // ok, we've got a1, a2, b1, b2, and e's
      {
        IrreducibleRepresentation ir(g,1,"B1");
        for (i=0; i < g; i++)
          ir.rep[0][i] = pow(-1.0,(double)i);

        gamma_[2] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"B2");
        for (i=0; i < g/2 ; i++)
          ir.rep[0][i] = pow(-1.0,(double)i);
        for (; i < g ; i++)
          ir.rep[0][i] = pow(-1.0,(double)(i+1));

        gamma_[3] = ir;
      }

      for (i=4; i < nirrep_; i++) {
        char lab[4];
        if (nt==4)
          sprintf(lab,"E");
        else
          sprintf(lab,"E%d",i-3);

        IrreducibleRepresentation ir(g,2,lab);

        for (j=0; j < g/2; j++) {
          ir.rep[0][j] = 2.0*cos((double)theta*j*(i-3));
          ir.rep[1][j] = 2.0*sin((double)theta*j*(i-3));
        }
        
        for (; j < g; j++) {
          ir.rep[0][j] = 0;
          ir.rep[1][j] = 0;
        }

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
      {
        IrreducibleRepresentation ir(g,1,"A'");
        for (i=0; i < g ; i++)
          ir.rep[0][i] = 1;

        gamma_[0] = ir;
      }

      for (i=1; i < nirrep_/2 ; i++) {
        if (nt==3)
          sprintf(label,"E'");
        else
          sprintf(label,"E%d'",i);

        IrreducibleRepresentation ir(g,2,label);

        ir.rep[0][0] = 1.0;
        ir.rep[1][0] = 0.0;

        for (j=1; j < g/2; j++) {
          ir.rep[0][j] = cos((double)theta*j*i);
          ir.rep[1][j] = sin((double)theta*j*i);
        }
        
        for (; j < g; j++) {
          ir.rep[0][j] = ir.rep[0][j-(g/2)];
          ir.rep[1][j] = ir.rep[1][j-(g/2)];
        }

        gamma_[i] = ir;
      }

      i=nirrep_/2;
      {
        IrreducibleRepresentation ir(g,1,"A\"");
        for (j=0; j < g/2; j++)
          ir.rep[0][j] = 1.0;
        for ( ; j < g; j++)
          ir.rep[0][j] = -1.0;

        gamma_[i] = ir;
      }

      for (i=i+1; i < nirrep_ ; i++) {
        i0 = i-(nirrep_/2);
        if (nt==3)
          sprintf(label,"E\"");
        else
          sprintf(label,"E%d\"",i0);

        IrreducibleRepresentation ir(g,2,label);

        for (j=0; j < g/2; j++) {
          ir.rep[0][j] = gamma_[i0].rep[0][j];
          ir.rep[1][j] = gamma_[i0].rep[1][j];
        }
        
        for (j=g/2; j < g ; j++) {
          ir.rep[0][j] = -gamma_[i0].rep[0][j];
          ir.rep[1][j] = -gamma_[i0].rep[1][j];
        }

        gamma_[i] = ir;
      }
    } else {
      double ineg = pow(-1.0,(double)nt/2);

      {
        IrreducibleRepresentation ir(g,1,"Ag");
        for (i=0; i < g; i++)
          ir.rep[0][i] = 1.0;

        gamma_[0] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"Bg");
        for (i=0; i < g/2; i++)
          ir.rep[0][i] = pow(-1.0,(double) i);
        for (; i < g; i++)
          ir.rep[0][i] = ineg*pow(-1.0,(double) i);

        gamma_[1] = ir;
      }

      for (i=2; i < nirrep_/2; i++) {
        ineg=pow(-1.0,(double) (i-1));

        if (nt==4)
          sprintf(label,"Eg");
        else
          sprintf(label,"E%dg",i-1);

        IrreducibleRepresentation ir(g,2,label);

        ir.rep[0][0] = 1.0;
        ir.rep[1][0] = 0.0;

        for (j=1; j < g/2; j++) {
          ir.rep[0][j] = cos((double)theta*j*(i-1));
          ir.rep[1][j] = sin((double)theta*j*(i-1));
        }
        
        for (; j < g; j++) {
          ir.rep[0][j] = ineg*ir.rep[0][j-(g/2)];
          ir.rep[1][j] = ineg*ir.rep[1][j-(g/2)];
        }

        gamma_[i] = ir;
      }

      i=nirrep_/2;
      {
        IrreducibleRepresentation ir(g,1,"Au");
        for (j=0; j < g/2; j++)
          ir.rep[0][j] = 1.0;
        for (j=g/2; j < g; j++)
          ir.rep[0][j] = -1.0;

        gamma_[i] = ir;
      }

      i++;
      {
        IrreducibleRepresentation ir(g,1,"Bu");
        for (j=0; j < g/2; j++)
          ir.rep[0][j] = gamma_[1].rep[0][j];
        for (j=g/2; j < g; j++)
          ir.rep[0][j] = -gamma_[1].rep[0][j];

        gamma_[i] = ir;
      }

      for (i=i+1; i < nirrep_ ; i++) {
        i0=i-(nirrep_/2);

        if (nt==4)
          sprintf(label,"Eu");
        else
          sprintf(label,"E%du",i0-1);

        IrreducibleRepresentation ir(g,2,label);

        for (j=0; j < g/2; j++) {
          ir.rep[0][j] =  gamma_[i0].rep[0][j];
          ir.rep[1][j] =  gamma_[i0].rep[1][j];
        }
        
        for (j=g/2; j < g; j++) {
          ir.rep[0][j] = -gamma_[i0].rep[0][j];
          ir.rep[1][j] = -gamma_[i0].rep[1][j];
        }

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
      {
        IrreducibleRepresentation ir(g,1,"Ag");
        for (i=0; i < g ; i++)
          ir.rep[0][i] = 1.0;

        gamma_[0] = ir;
      }

      for (i=1; i < nirrep_/2 ; i++) {
        if (nt==6)
          sprintf(label,"Eg");
        else
          sprintf(label,"E%dg",i);

        IrreducibleRepresentation ir(g,2,label);

        ir.rep[0][0] = 1.0;
        ir.rep[1][0] = 0.0;

        for (j=1; j < g ; j++) {
          ir.rep[0][j] = cos((double)2.0*theta*j*i);
          ir.rep[1][j] = sin((double)2.0*theta*j*i);
        }

        gamma_[i] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"Au");
        for (j=0; j < g ; j++)
          ir.rep[0][j] = pow(-1.0,(double)j);

        gamma_[i] = ir;
      }

      for (i=i+1; i < nirrep_ ; i++) {
        i0=i-nirrep_/2;

        if (nt==6)
          sprintf(label,"Eu");
        else
          sprintf(label,"E%du",i0);

        IrreducibleRepresentation ir(g,2,label);

        for (j=0; j < g; j++) {
          ir.rep[0][j] = pow(-1.0,(double) j)*gamma_[i0].rep[0][j];
          ir.rep[1][j] = pow(-1.0,(double) j)*gamma_[i0].rep[1][j];
        }

        gamma_[i] = ir;
      }
    } else {
      {
        IrreducibleRepresentation ir(g,1,"A");
        for (i=0; i < g; i++)
          ir.rep[0][i] = 1.0;

        gamma_[0] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"B");
        for (i=0; i < g; i++)
          ir.rep[0][i] = pow(-1.0,(double)i);

        gamma_[1] = ir;
      }

      for (i=2; i < nirrep_ ; i++) {
        if (nt==4)
          sprintf(label,"E");
        else
          sprintf(label,"E%d",i-1);

        IrreducibleRepresentation ir(g,2,label);

        for (j=0; j < g; j++) {
          ir.rep[0][j] = cos((double)theta*j*(i-1));
          ir.rep[1][j] = sin((double)theta*j*(i-1));
        }

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
      {
        IrreducibleRepresentation ir(g,1,"A");
        for (i=0; i < g; i++)
          ir.rep[0][i] = 1.0;

        gamma_[0] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"B1");
        ir.rep[0][0] = ir.rep[0][1] = 1.0;
        ir.rep[0][2] = ir.rep[0][3] = -1.0;
        gamma_[1] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"B2");
        ir.rep[0][0] = ir.rep[0][3] = 1.0;
        ir.rep[0][1] = ir.rep[0][2] = -1.0;
        gamma_[2] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"B3");
        ir.rep[0][0] = ir.rep[0][2] = 1.0;
        ir.rep[0][1] = ir.rep[0][3] = -1.0;
        gamma_[3] = ir;
      }
    } else if (nt%2) {
      {
        IrreducibleRepresentation ir(g,1,"A1");
        for (i=0; i < g; i++)
          ir.rep[0][i] = 1.0;

        gamma_[0] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"A2");
        for (i=0; i < g/2; i++)
          ir.rep[0][i] = 1.0;
        for (i=g/2; i < g; i++)
          ir.rep[0][i] = -1.0;
        gamma_[1] = ir;
      }

      for (i=2; i < nirrep_ ; i++) {
        if (nt==3)
          sprintf(label,"E");
        else
          sprintf(label,"E%d",i-1);

        IrreducibleRepresentation ir(g,2,label);

        ir.rep[0][0] = 2.0;
        ir.rep[1][0] = 0.0;

        for (j=1; j < g/2; j++) {
          ir.rep[0][j] = 2.0*cos((double)theta*j*(i-1));
          ir.rep[1][j] = 2.0*sin((double)theta*j*(i-1));
        }
        
        for (j=g/2; j < g; j++) {
          ir.rep[0][j] = 0.0;
          ir.rep[1][j] = 0.0;
        }

        gamma_[i] = ir;
      }
    } else {
      {
        IrreducibleRepresentation ir(g,1,"A1");
        for (i=0; i < g; i++)
          ir.rep[0][i] = 1.0;

        gamma_[0] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"A2");
        for (i=0; i < g/2; i++)
          ir.rep[0][i] = 1.0;
        for (i=g/2; i < g; i++)
          ir.rep[0][i] = -1.0;

        gamma_[1] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"B1");
        for (i=0; i < g; i++)
          ir.rep[0][i] = pow(-1.0,(double)i);

        gamma_[2] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"B2");
        for (i=0; i < g/2; i++)
          ir.rep[0][i] = pow(-1.0,(double)i);
        for (i=g/2; i < g; i++)
          ir.rep[0][i] = pow(-1.0,(double)(i+1));

        gamma_[3] = ir;
      }

      for (i=4; i < nirrep_; i++) {
        if (nt==4)
          sprintf(label,"E");
        else
          sprintf(label,"E%d",i-3);

        IrreducibleRepresentation ir(g,2,label);

        for (j=0; j < g/2; j++) {
          ir.rep[0][j] = 2.0*cos((double)theta*j*(i-3));
          ir.rep[1][j] = 2.0*sin((double)theta*j*(i-3));
        }
        
        for (j=g/2; j < g; j++) {
          ir.rep[0][j] = 0.0;
          ir.rep[1][j] = 0.0;
        }

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
      {
        IrreducibleRepresentation ir(g,1,"A1g");
        for (i=0; i < g; i++)
          ir.rep[0][i] = 1.0;

        gamma_[0] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"A2g");
        for (i=0; i < g/2; i++)
          ir.rep[0][i] = 1.0;
        for (i=g/2; i < g; i++)
          ir.rep[0][i] = -1.0;

        gamma_[1] = ir;
      }

      for (i=2; i < nirrep_/2 ; i++) {
        if (nt==3)
          sprintf(label,"Eg");
        else
          sprintf(label,"E%dg",i-1);

        IrreducibleRepresentation ir(g,2,label);

        ir.rep[0][0] = 2.0;
        ir.rep[1][0] = 0.0;

        for (j=1; j < g/2; j++) {
          ir.rep[0][j] = 2.0*cos((double)2.0*theta*j*(i-1));
          ir.rep[1][j] = 2.0*sin((double)2.0*theta*j*(i-1));
        }
        
        for (j=g/2; j < g; j++) {
          ir.rep[0][j] = 0.0;
          ir.rep[1][j] = 0.0;
        }

        gamma_[i] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"A1u");
        gamma_[i] = ir;
        i++;
      }
      {
        IrreducibleRepresentation ir(g,1,"A2u");
        gamma_[i] = ir;
        i++;
      }

      for (; i < nirrep_; i++) {
        i0=i-nirrep_/2;
        
        if (nt==3)
          sprintf(label,"Eu");
        else
          sprintf(label,"E%du",i0-1);

        IrreducibleRepresentation ir(g,2,label);
        gamma_[i] = ir;
      }

      for (i=nirrep_/2,i0=0; i < nirrep_; i++,i0++) {
        for (int d=0; d < gamma_[i].degen; d++) {
          for (j=0; j < g/2; j++)
            gamma_[i].rep[d][j] = pow(-1.0,(double) j)*gamma_[i0].rep[d][j];
          for (; j < 3*g/4; j++)
            gamma_[i].rep[d][j] =  gamma_[i0].rep[d][j];
          for (; j < g    ; j++)
            gamma_[i].rep[d][j] = -gamma_[i0].rep[d][j];
        }
      }
    } else {
      {
        IrreducibleRepresentation ir(g,1,"A1");
        for (i=0; i < g; i++)
          ir.rep[0][i] = 1.0;

        gamma_[0] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"A2");
        for (i=0; i < g/2; i++)
          ir.rep[0][i] = 1.0;
        for (i=g/2; i < g; i++)
          ir.rep[0][i] = -1.0;
        gamma_[1] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"B1");
        for (i=0; i < g/2; i++)
          ir.rep[0][i] = pow(-1.0,(double) i);
        for ( ; i < 3*g/4; i++)
          ir.rep[0][i] =  1.0;
        for ( ; i < g    ; i++)
          ir.rep[0][i] = -1.0;

        gamma_[2] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"B2");
        for (i=0; i < g/2; i++)
          ir.rep[0][i] = pow(-1.0,(double) i);
        for ( ; i < 3*g/4; i++)
          ir.rep[0][i] = -1.0;
        for ( ; i < g    ; i++)
          ir.rep[0][i] = 1.0;

        gamma_[3] = ir;
      }

      for (i=4; i < nirrep_; i++) {
        if (nt==2)
          sprintf(label,"E");
        else
          sprintf(label,"E%d",i-3);

        IrreducibleRepresentation ir(g,2,label);

        ir.rep[0][0] = 2.0;
        ir.rep[1][0] = 0.0;

        for (j=1; j < g/2; j++) {
          ir.rep[0][j] = 2.0*cos((double)theta*j*(i-3)*0.5);
          ir.rep[1][j] = 2.0*sin((double)theta*j*(i-3)*0.5);
        }
        
        for (j=g/2; j < g; j++) {
          ir.rep[0][j] = 0.0;
          ir.rep[1][j] = 0.0;
        }

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
      {
        IrreducibleRepresentation ir(g,1,"Ag");
        for (i=0; i < 8; i++)
          ir.rep[0][i] = 1.0;

        gamma_[0] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"B1g");
        for (i=0; i < 4; i++)
          ir.rep[0][i] = 1.0;
        for (i=4; i < 8; i++) 
          ir.rep[0][i] = -1.0;

        gamma_[1] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"B2g");
        ir.rep[0][0] = ir.rep[0][3] = ir.rep[0][5] = ir.rep[0][6] = 1.0;
        ir.rep[0][1] = ir.rep[0][2] = ir.rep[0][4] = ir.rep[0][7] = -1.0;

        gamma_[2] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"B3g");
        ir.rep[0][0] = ir.rep[0][3] = ir.rep[0][4] = ir.rep[0][7] = 1.0;
        ir.rep[0][1] = ir.rep[0][2] = ir.rep[0][5] = ir.rep[0][6] = -1.0;

        gamma_[3] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"Au");
        for (i=0; i < 2; i++)
          ir.rep[0][i] = 1.0;
        for (i=2; i < 4; i++)
          ir.rep[0][i] = -1.0;
        for (i=4; i < 6; i++)
          ir.rep[0][i] = 1.0;
        for (i=6; i < 8; i++)
          ir.rep[0][i] = -1.0;

        gamma_[4] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"B1u");
        for (i=0; i < 2; i++)
          ir.rep[0][i] = 1.0;
        for (i=2; i < 6; i++)
          ir.rep[0][i] = -1.0;
        for (i=6; i < 8; i++)
          ir.rep[0][i] = 1.0;

        gamma_[5] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"B2u");
        for (i=0; i < 4; i++)
          ir.rep[0][i] =  pow(-1.0,(double)i);
        for (i=4; i < 8; i++)
          ir.rep[0][i] = -pow(-1.0,(double)i);

        gamma_[6] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"B3u");
        for (i=0; i < 8; i++)
          ir.rep[0][i] = pow(-1.0,(double)i);

        gamma_[7] = ir;
      }
    }

    else if (nt%2) {
      {
        IrreducibleRepresentation ir(g,1,"A1'");
        for (i=0; i < g; i++)
          ir.rep[0][i] = 1.0;

        gamma_[0] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"A2'");
        for (i=0; i < g/2; i++)
          ir.rep[0][i] = 1.0;
        for (i=g/2; i < g; i++)
          ir.rep[0][i] = -1.0;

        gamma_[1] = ir;
      }

      for (i=2; i < nirrep_/2 ; i++) {
        if (nt==3)
          sprintf(label,"E'");
        else
          sprintf(label,"E%d'",i-1);

        IrreducibleRepresentation ir(g,2,label);

        ir.rep[0][0] = 2.0;
        ir.rep[1][0] = 0.0;

        for (j=1; j < g/4; j++) {
          ir.rep[0][j] = 2.0*cos((double)theta*j*(i-1));
          ir.rep[1][j] = 2.0*sin((double)theta*j*(i-1));
        }
        
        for (   ; j < g/2; j++) {
          ir.rep[0][j] = 2.0*cos((double)theta*j*(i-1));
          ir.rep[1][j] = 2.0*sin((double)theta*j*(i-1));
        }
        
        for (   ; j < g  ; j++) {
          ir.rep[0][j] = 0.0;
          ir.rep[1][j] = 0.0;
        }

        gamma_[i] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"A1\"");
        gamma_[i] = ir;
        i++;
      }
      {
        IrreducibleRepresentation ir(g,1,"A2\"");
        gamma_[i] = ir;
        i++;
      }

      for ( ; i < nirrep_; i++) {
        i0=i-(nirrep_/2);

        if (nt==3)
          sprintf(label,"E\"");
        else
          sprintf(label,"E%d\"",i0-1);

        IrreducibleRepresentation ir(g,2,label);

        gamma_[i] = ir;
      }

      for (i=nirrep_/2,i0=0; i < nirrep_; i0++,i++) {
        for (int d=0; d < gamma_[i].degen; d++) {
          for (j=0; j < g/4; j++) gamma_[i].rep[d][j] =  gamma_[i0].rep[d][j];
          for (   ; j < g/2; j++) gamma_[i].rep[d][j] = -gamma_[i0].rep[d][j];
          for ( ; j < 3*g/4; j++) gamma_[i].rep[d][j] =  gamma_[i0].rep[d][j];
          for ( ; j < g    ; j++) gamma_[i].rep[d][j] = -gamma_[i0].rep[d][j];
        }
      }

    } else {
      {
        IrreducibleRepresentation ir(g,1,"A1g");
        for (i=0; i < g; i++)
          ir.rep[0][i] = 1.0;

        gamma_[0] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"A2g");
        for (i=0; i < g/2; i++)
          ir.rep[0][i] = 1.0;
        for (i=g/2; i < g; i++)
          ir.rep[0][i] = -1.0;

        gamma_[1] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"B1g");
        for (i=0; i < g; i++)
          ir.rep[0][i] = pow(-1.0,(double)i);

        gamma_[2] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"B2g");
        for (i=0; i < g/2; i++)
          ir.rep[0][i] = pow(-1.0,(double)i);
        for (i=g/2; i < g; i++)
          ir.rep[0][i] = pow(-1.0,(double)(i+1));

        gamma_[3] = ir;
      }

      for (i=4; i < nirrep_/2; i++) {
        if (nt==4)
          sprintf(label,"Eg");
        else
          sprintf(label,"E%dg",i-3);

        IrreducibleRepresentation ir(g,2,label);

        for (j=0; j < g/4 ; j++) {
          ir.rep[0][j] = 2.0*cos((double)theta*j*(i-3));
          ir.rep[1][j] = 2.0*sin((double)theta*j*(i-3));
        }
        
        for (j=g/4; j < g/2; j++) {
          ir.rep[0][j] =
            pow(-1.0,(double)(i+1))*2.0*cos((double)theta*j*(i-3));
          ir.rep[1][j] =
            pow(-1.0,(double)(i+1))*2.0*sin((double)theta*j*(i-3));
        }

        for (j=g/2; j < g; j++) {
          ir.rep[0][j] = 0.0;
          ir.rep[1][j] = 0.0;
        }

        gamma_[i] = ir;
      }

      {
        IrreducibleRepresentation ir(g,1,"A1u");
        gamma_[i] = ir;
        i++;
      }
      {
        IrreducibleRepresentation ir(g,1,"A2u");
        gamma_[i] = ir;
        i++;
      }
      {
        IrreducibleRepresentation ir(g,1,"B1u");
        gamma_[i] = ir;
        i++;
      }
      {
        IrreducibleRepresentation ir(g,1,"B2u");
        gamma_[i] = ir;
        i++;
      }

      for (; i < nirrep_; i++) {
        i0=i-(nirrep_/2);

        if (nt==4)
          sprintf(label,"Eu");
        else
          sprintf(label,"E%du",i0-3);

        IrreducibleRepresentation ir(g,2,label);
        gamma_[i] = ir;
      }

      for (i=nirrep_/2,i0=0; i < nirrep_; i0++,i++) {
        for (int d=0; d < gamma_[i].degen; d++) {
        for (j=0; j < g/4; j++) gamma_[i].rep[d][j] =  gamma_[i0].rep[d][j];
        for (   ; j < g/2; j++) gamma_[i].rep[d][j] = -gamma_[i0].rep[d][j];
        for ( ; j < 3*g/4; j++) gamma_[i].rep[d][j] =  gamma_[i0].rep[d][j];
        for (   ; j < g  ; j++) gamma_[i].rep[d][j] = -gamma_[i0].rep[d][j];
        }
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
        nr += rot[j]*gamma_[i].rep[0][j];
        nt += trans[j]*gamma_[i].rep[0][j];
      }

      gamma_[i].nrot_ = (int) ((nr+0.5)/gamma_[i].g);
      gamma_[i].ntrans_ = (int) ((nt+0.5)/gamma_[i].g);
    }
  }

  delete[] rot;
  delete[] trans;
  
  return 0;
}
