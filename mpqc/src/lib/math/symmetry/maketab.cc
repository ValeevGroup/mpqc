
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
  int i,j,i0,ei,gi;
  char label[4];
  double itheta, jitheta, ctheta, stheta;

  if (!g) return 0;

  gamma_ = new IrreducibleRepresentation[nirrep_];

  for (i=0; i < nirrep_; i++)
    gamma_[i].init();

  symop = new SymmetryOperation[g];

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
      ir.rep[0] = 1;
      ir.proj[0][0] = 1;

      gamma_[0] = ir;
    }

    symop[0](0,0) = symop[0](1,1) = symop[0](2,2) = 1;
    break;

  case CI:  // equivalent to S2 about the z axis
    {
      IrreducibleRepresentation ir(2,1,"Ag");
      ir.rep[0] = 1;
      ir.rep[1] = 1;
      ir.proj[0][0] = 1;
      ir.proj[0][1] = 1;
      ir.nrot_=3;

      gamma_[0] = ir;
    }

    symop[0](0,0) = symop[0](1,1) = symop[0](2,2) = 1;

    {
      IrreducibleRepresentation ir(2,1,"Au");
      ir.rep[0] = 1;
      ir.rep[1] = -1;
      ir.proj[0][0] = 1;
      ir.proj[0][1] = -1;
      ir.ntrans_=3;

      gamma_[1] = ir;
    }

    symop[1](0,0) = symop[1](1,1) = symop[1](2,2) = -1;

    break;

  case CS: // reflection through the xy plane
    {
      IrreducibleRepresentation ir(2,1,"A'");
      ir.rep[0] = 1;
      ir.rep[1] = 1;
      ir.proj[0][0] = 1;
      ir.proj[0][1] = 1;
      ir.nrot_=1;
      ir.ntrans_=2;

      gamma_[0] = ir;
    }

    symop[0](0,0) = symop[0](1,1) = symop[0](2,2) = 1;

    {
      IrreducibleRepresentation ir(2,1,"A\"");
      ir.rep[0] = 1;
      ir.rep[1] = -1;
      ir.proj[0][0] = 1;
      ir.proj[0][1] = -1;
      ir.nrot_=2;
      ir.ntrans_=1;

      gamma_[1] = ir;
    }

    symop[1](0,0) = symop[1](1,1) = 1;
    symop[1](2,2) = -1;

    break;

  case CN:
    // clockwise rotation about z axis by theta*i radians
    //
    // for odd n, the irreps are A and E1...E(nir-1)
    // for even n, the irreps are A, B, and E1...E(nir-2)
    //
    {
      IrreducibleRepresentation ir(g,1,"A");
      for (gi=0; gi < g; gi++) {
        ir.rep[gi] = ir.proj[0][gi] = 1;
      }

      gamma_[0] = ir;
    }
    i=1;

    if (!(nt%2)) {
      IrreducibleRepresentation ir(g,1,"B");
      for (gi=0; gi < g; gi++) {
        ir.rep[gi] = ir.proj[0][gi] = (gi%2) ? -1.0 : 1.0;
      }

      gamma_[i] = ir;
      i++;
    }

    // for the E irreps, the projection operators are:
    //   Ei xx = cos(m*theta*i) m = 0-(nt-1)
    //      xy = -sin(m*theta*i)
    //      yx = sin(m*theta*i)
    //      yy = cos(m*theta*i)

    ei=1;
    itheta = theta;
    for (; i < nirrep_; i++, ei++, itheta += theta) {
      if (nt==3 || nt==4)
        sprintf(label,"E");
      else
        sprintf(label,"E%d",ei);

      IrreducibleRepresentation ir(g,2,label);

      jitheta = 0;
      for (j=0; j < g; j++, jitheta += itheta) {
        ctheta = cos(jitheta);
        stheta = sin(jitheta);
        
        ir.rep[j] = ctheta;

        ir.proj[0][j] = ctheta;
        ir.proj[1][j] = -stheta;
        ir.proj[2][j] = stheta;
        ir.proj[3][j] = ctheta;
      }

      gamma_[i] = ir;
    }

    itheta=0;
    for (i=0; i < nt ; i++, itheta += theta) {
      symop[i][0][0] = symop[i][1][1] = cos(itheta);
      symop[i][0][1] = sin(itheta);
      symop[i][1][0] = -sin(itheta);
      symop[i][2][2] = 1.0;

      rot[i] = trans[i] = symop[i].trace();
    }

    break;

  case CNV:
    // clockwise rotation about z axis by theta*i radians, then
    // reflect through the xz plane
    //
    // for odd n, the irreps are A1, A2, and E1...E(nir-2)
    // for even n, the irreps are A1, A2, B1, B2, and E1...E(nir-4)
    //
    {
      IrreducibleRepresentation ir1(g,1,"A1");
      IrreducibleRepresentation ir2(g,1,"A2");
      for (gi=0; gi < nt; gi++) {
        // Cn's
        ir1.rep[gi] = ir1.proj[0][gi] = 1;
        ir2.rep[gi] = ir2.proj[0][gi] = 1;

        // sigma's
        ir1.rep[gi+nt] = ir1.proj[0][gi+nt] = 1;
        ir2.rep[gi+nt] = ir2.proj[0][gi+nt] = -1;
      }

      gamma_[0] = ir1;
      gamma_[1] = ir2;
    }
    i=2;

    if (!(nt%2)) {
      IrreducibleRepresentation ir1(g,1,"B1");
      IrreducibleRepresentation ir2(g,1,"B2");
      for (gi=0; gi < nt ; gi++) {
        double ci = (gi%2) ? -1.0 : 1.0;
        
        // Cn's
        ir1.rep[gi] = ir1.proj[0][gi] = ci;
        ir2.rep[gi] = ir2.proj[0][gi] = ci;

        // sigma's
        ir1.rep[gi+nt] = ir1.proj[0][gi+nt] = ci;
        ir2.rep[gi+nt] = ir2.proj[0][gi+nt] = -ci;
      }

      gamma_[i] = ir1;
      i++;
      gamma_[i] = ir2;
      i++;
    }
    
    // for the E irreps, the projection operators are:
    // for the n Cn's:
    //   Ei xx = cos(m*theta*i) m = 0-(nt-1)
    //      xy = -sin(m*theta*i)
    //      yx = sin(m*theta*i)
    //      yy = cos(m*theta*i)
    //
    // for the n sigma's:
    //   Ei xx = cos(m*theta*i) m = 0-(nt-1)
    //      xy = sin(m*theta*i)
    //      yx = sin(m*theta*i)
    //      yy = -cos(m*theta*i)
    ei=1;
    itheta=theta;
    for (; i < nirrep_; i++, ei++, itheta += theta) {
      char lab[4];
      if (nt==3 || nt==4)
        sprintf(lab,"E");
      else
        sprintf(lab,"E%d",ei);

      IrreducibleRepresentation ir(g,2,lab);

      jitheta = 0;
      for (j=0; j < nt; j++, jitheta += itheta) {
        ctheta = cos(jitheta);
        stheta = sin(jitheta);
          
        // Cn's
        ir.rep[j] = 2.0*ctheta;

        ir.proj[0][j] = ctheta;
        ir.proj[1][j] = -stheta;
        ir.proj[2][j] = stheta;
        ir.proj[3][j] = ctheta;

        // sigma's
        ir.rep[j+nt] = 0;

        ir.proj[0][j+nt] = ctheta;
        ir.proj[1][j+nt] = stheta;
        ir.proj[2][j+nt] = stheta;
        ir.proj[3][j+nt] = -ctheta;
      }

      gamma_[i] = ir;
    }

    itheta = 0;
    for (i=0, j=nt; i < nt ; i++, j++, itheta += theta) {
      // Cn's
      symop[i][0][0] = symop[i][1][1] = cos(itheta);
      symop[i][0][1] = sin(itheta);
      symop[i][1][0] = -sin(itheta);
      symop[i][2][2] = 1.0;

      rot[i] = trans[i] = symop[i].trace();

      // sigma's
      symop[j][0][0] = cos(itheta);
      symop[j][1][1] = -cos(itheta);
      symop[j][1][0] = symop[i+nt][0][1] = sin(itheta);
      symop[j][2][2] = 1.0;

      trans[j] = symop[j].trace();
      rot[j] = -trans[j];
    }

    break;

  case CNH:
    // lockwise rotation about z axis by theta*i radians,
    // as well as rotation-reflection about same axis

    //
    // for odd n, the irreps are A', A", E1'...E(nir/2-1)', E1"...E(nir/2-1)''
    // for even n, the irreps are Ag, Bg, Au, Bu,
    //                            E1g...E(nir/2-1)g, E1u...E(nir/2-1)u
    //
    {
      IrreducibleRepresentation ir1(g,1, (nt%2) ? "A'" : "Ag");
      IrreducibleRepresentation ir2(g,1, (nt%2) ? "A\"" : "Au");
      for (gi=0; gi < nt; gi++) {
        // Cn's
        ir1.rep[gi] = ir1.proj[0][gi] = 1;
        ir2.rep[gi] = ir2.proj[0][gi] = 1;

        // Sn's
        ir1.rep[gi+nt] = ir1.proj[0][gi+nt] = 1;
        ir2.rep[gi+nt] = ir2.proj[0][gi+nt] = -1;
      }

      gamma_[0] = ir1;
      gamma_[nirrep_/2] = ir2;
    }
    i=1;

    if (!(nt%2)) {
      double ineg = ((nt/2)%2) ? -1.0 : 1.0;
      
      IrreducibleRepresentation ir1(g,1,"Bg");
      IrreducibleRepresentation ir2(g,1,"Bu");
      for (gi=0; gi < nt; gi++) {
        double ci = (gi%2) ? -1.0 : 1.0;

        // Cn's
        ir1.rep[gi] = ir1.proj[0][gi] = ci;
        ir2.rep[gi] = ir2.proj[0][gi] = ci;
      
        // Sn's
        ir1.rep[gi+nt] = ir1.proj[0][gi+nt] = ci*ineg;
        ir2.rep[gi+nt] = ir2.proj[0][gi+nt] = -ci*ineg;
      }

      gamma_[i] = ir1;
      gamma_[i+nirrep_/2] = ir2;
      i++;
    }

    // for the E irreps, the projection operators are:
    // for the n Cn's:
    //   Ei xx = cos(m*theta*i) m = 0-(nt-1)
    //      xy = -sin(m*theta*i)
    //      yx = sin(m*theta*i)
    //      yy = cos(m*theta*i)
    //
    // for the n Sn's:
    //   Ei xx = cos(m*theta*i) m = 0-(nt-1)
    //      xy = -sin(m*theta*i)
    //      yx = sin(m*theta*i)
    //      yy = cos(m*theta*i)
    ei=1;
    itheta=theta;
    for (; i < nirrep_/2 ; i++, ei++, itheta += theta) {
      if (nt==3 || nt==4)
        sprintf(label,(nt%2) ? "E'" : "Eg");
      else
        sprintf(label,"E%d%s", ei, (nt%2) ? "'" : "g");

      IrreducibleRepresentation ir1(g,2,label);

      if (nt==3 || nt==4)
        sprintf(label,(nt%2) ? "E\"" : "Eu");
      else
        sprintf(label,"E%d%s", ei, (nt%2) ? "\"" : "u");

      IrreducibleRepresentation ir2(g,2,label);

      jitheta = 0;
      double ineg = (nt%2) ? 1.0 : pow(-1.0,(double)ei);

      for (j=0; j < nt; j++, jitheta += itheta) {
        ctheta = cos(jitheta);
        stheta = sin(jitheta);
        
        ir1.rep[j] = ir2.rep[j] = ctheta;
        
        // Cn's
        ir1.proj[0][j] = ctheta;
        ir1.proj[1][j] = -stheta;
        ir1.proj[2][j] = stheta;
        ir1.proj[3][j] = ctheta;

        ir2.proj[0][j] = ctheta;
        ir2.proj[1][j] = -stheta;
        ir2.proj[2][j] = stheta;
        ir2.proj[3][j] = ctheta;

        // Sn's
        ir1.rep[j+nt] = ineg*ctheta;
        ir2.rep[j+nt] = -ineg*ctheta;
        
        ir1.proj[0][j+nt] = ineg*ctheta;
        ir1.proj[1][j+nt] = -ineg*stheta;
        ir1.proj[2][j+nt] = ineg*stheta;
        ir1.proj[3][j+nt] = ineg*ctheta;

        ir2.proj[0][j+nt] = -ineg*ctheta;
        ir2.proj[1][j+nt] = ineg*stheta;
        ir2.proj[2][j+nt] = -ineg*stheta;
        ir2.proj[3][j+nt] = -ineg*ctheta;
      }
        
      gamma_[i] = ir1;
      gamma_[i+nirrep_/2] = ir2;
    }

    itheta = 0;
    for (i=0; i < nt ; i++, itheta += theta) {
      symop[i][0][0] = symop[i][1][1] =
        symop[i+nt][0][0] = symop[i+nt][1][1] = cos(itheta);
      symop[i][0][1] = symop[i+nt][0][1] = sin(itheta);
      symop[i][1][0] = symop[i+nt][1][0] = -sin(itheta);
      symop[i][2][2] = 1.0;
      symop[i+nt][2][2] = -1.0;

      rot[i] = trans[i] = symop[i].trace();
      trans[i+nt] = symop[i+nt].trace();
      rot[i+nt] = -trans[i+nt];
    }

    break;

  case SN:
    // clockwise rotation-reflection by theta*i radians about z axis
    //
    // for odd n/2, the irreps are Ag, Au, E1g...E(nir/2-1)g,E1u...E(nir/2-1)u
    // for even n/2, the irreps are A, B, E1...E(nir-2)
    //
    i=0;
    {
      IrreducibleRepresentation ir1(g,1,((nt/2)%2) ? "Ag" : "A");
      IrreducibleRepresentation ir2(g,1,((nt/2)%2) ? "Au" : "B");
      for (gi=0; gi < nt; gi++) {
        ir1.rep[gi] = ir1.proj[0][gi] = 1.0;
        ir2.rep[gi] = ir2.proj[0][gi] = (gi%2) ? -1.0 : 1.0;
      }
      
      gamma_[i] = ir1;
      if ((nt/2)%2) {
        gamma_[i+nirrep_/2] = ir2;
      } else {
        i++;
        gamma_[i] = ir2;
      }
      i++;
    }

    // for the E irreps, the projection operators are:
    //   Ei xx = cos(m*theta*i) m = 0-(nt-1)
    //      xy = -sin(m*theta*i)
    //      yx = sin(m*theta*i)
    //      yy = cos(m*theta*i)
    ei=1;
    itheta=theta;
    if ((nt/2)%2) {
      for (; i < nirrep_/2 ; i++, ei++, itheta += theta) {
        if (nt==6)
          sprintf(label,"Eg");
        else
          sprintf(label,"E%dg",ei);

        IrreducibleRepresentation ir1(g,2,label);

        if (nt==6)
          sprintf(label,"Eu");
        else
          sprintf(label,"E%du", ei);

        IrreducibleRepresentation ir2(g,2,label);

        jitheta = 0;
        double ineg = -pow(-1.0,(double)ei);
        for (j=0; j < nt; j++, jitheta += itheta) {
          double ci1 = (j%2) ? -ineg : 1.0;
          double ci2 = (j%2) ? ineg : 1.0;
          
          ctheta = cos(jitheta);
          stheta = sin(jitheta);
        
          ir1.rep[j] = ci1*ctheta;
        
          ir1.proj[0][j] = ci1*ctheta;
          ir1.proj[1][j] = -ci1*stheta;
          ir1.proj[2][j] = ci1*stheta;
          ir1.proj[3][j] = ci1*ctheta;

          ir2.rep[j] = ci2*ctheta;
        
          ir2.proj[0][j] = ci2*ctheta;
          ir2.proj[1][j] = -ci2*stheta;
          ir2.proj[2][j] = ci2*stheta;
          ir2.proj[3][j] = ci2*ctheta;
        }

        gamma_[i] = ir1;
        gamma_[i+nirrep_/2] = ir2;
      }
    } else {
      for (; i < nirrep_; i++, ei++, itheta += theta) {
        if (nt==4)
          sprintf(label,"E");
        else
          sprintf(label,"E%d",ei);

        IrreducibleRepresentation ir(g,2,label);

        jitheta = 0;
        for (j=0; j < nt; j++, jitheta += itheta) {
          ctheta = cos(jitheta);
          stheta = sin(jitheta);
        
          ir.rep[j] = ctheta;
        
          ir.proj[0][j] = ctheta;
          ir.proj[1][j] = -stheta;
          ir.proj[2][j] = stheta;
          ir.proj[3][j] = ctheta;
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

  case DN:
    // clockwise rotation about z axis, followed by C2 about x axis

    // D2 is a special case
    if (nt==2) {
      IrreducibleRepresentation ir1(g,1,"A");
      IrreducibleRepresentation ir2(g,1,"B1");
      IrreducibleRepresentation ir3(g,1,"B2");
      IrreducibleRepresentation ir4(g,1,"B3");

      for (i=0; i < g; i++) {
        ir1.rep[i] = ir1.proj[0][i] = 1.0;
        ir2.rep[i] = ir2.proj[0][i] = (i < 2) ? 1.0 : -1.0;
        ir3.rep[i] = ir3.proj[0][i] = (i % 2) ? -1.0 : 1.0;
        ir4.rep[i] = ir4.proj[0][i] = (i < 2) ?
          ((i % 2) ? -1.0 : 1.0) : ((i%2) ? 1.0 : -1.0);
      }

      gamma_[0] = ir1;
      gamma_[1] = ir2;
      gamma_[2] = ir3;
      gamma_[3] = ir4;

    } else {
      // Dn is isomorphic with Cnv
      //
      // for odd n, the irreps are A1, A2, and E1...E(nir-2)
      // for even n, the irreps are A1, A2, B1, B2, and E1...E(nir-4)
      //
      {
        IrreducibleRepresentation ir1(g,1,"A1");
        IrreducibleRepresentation ir2(g,1,"A2");
        for (gi=0; gi < nt; gi++) {
          // Cn's
          ir1.rep[gi] = ir1.proj[0][gi] = 1.0;
          ir2.rep[gi] = ir2.proj[0][gi] = 1.0;

          // C2's
          ir1.rep[gi+nt] = ir1.proj[0][gi+nt] = 1.0;
          ir2.rep[gi+nt] = ir2.proj[0][gi+nt] = -1.0;
        }

        gamma_[0] = ir1;
        gamma_[1] = ir2;
      }
      i=2;

      if (!(nt%2)) {
        IrreducibleRepresentation ir1(g,1,"B1");
        IrreducibleRepresentation ir2(g,1,"B2");
        for (gi=0; gi < nt ; gi++) {
          double ci = (gi%2) ? -1.0 : 1.0;
        
          // Cn's
          ir1.rep[gi] = ir1.proj[0][gi] = ci;
          ir2.rep[gi] = ir2.proj[0][gi] = ci;

          // sigma's
          ir1.rep[gi+nt] = ir1.proj[0][gi+nt] = ci;
          ir2.rep[gi+nt] = ir2.proj[0][gi+nt] = -ci;
        }

        gamma_[i] = ir1;
        i++;
        gamma_[i] = ir2;
        i++;
      }

      // for the E irreps, the projection operators are:
      // for the n Cn's:
      //   Ei xx = cos(m*theta*i) m = 0-(nt-1)
      //      xy = -sin(m*theta*i)
      //      yx = sin(m*theta*i)
      //      yy = cos(m*theta*i)
      //
      // for the n C2's:
      //   Ei xx = cos(m*theta*i) m = 0-(nt-1)
      //      xy = -sin(m*theta*i)
      //      yx = -sin(m*theta*i)
      //      yy = -cos(m*theta*i)
      ei=1;
      itheta=theta;
      for (; i < nirrep_; i++, ei++, itheta += theta) {
        char lab[4];
        if (nt==3 || nt==4)
          sprintf(lab,"E");
        else
          sprintf(lab,"E%d",ei);

        IrreducibleRepresentation ir(g,2,lab);

        jitheta = 0;
        for (j=0; j < nt; j++, jitheta += itheta) {
          ctheta = cos(jitheta);
          stheta = sin(jitheta);
          
          // Cn's
          ir.rep[j] = 2.0*ctheta;

          ir.proj[0][j] = ctheta;
          ir.proj[1][j] = -stheta;
          ir.proj[2][j] = stheta;
          ir.proj[3][j] = ctheta;

          // C2's
          ir.rep[j+nt] = 0;

          ir.proj[0][j+nt] = ctheta;
          ir.proj[1][j+nt] = -stheta;
          ir.proj[2][j+nt] = -stheta;
          ir.proj[3][j+nt] = -ctheta;
        }

        gamma_[i] = ir;
      }
    }
    
    itheta=0;
    for (i=0, j=nt; i < nt ; i++, j++, itheta += theta) {
      // Cn's
      symop[i][0][0] = symop[i][1][1] = cos(itheta);
      symop[i][0][1] = sin(itheta);
      symop[i][1][0] = -sin(itheta);
      symop[i][2][2] = 1.0;

      rot[i] = trans[i] = symop[i].trace();

      // C2's
      symop[j][0][0] = cos(itheta);
      symop[j][1][1] = -cos(itheta);
      symop[j][1][0] = symop[j][0][1] = -sin(itheta);
      symop[j][2][2] = -1.0;

      rot[j] = trans[j] = symop[j].trace();
    }

    break;

  case DND:
    // rotation reflection about z axis by theta/2 radians, followed
    // by c2 about x axis, then reflection through yz plane
    //
    // for odd n, the irreps are A1g, A2g, A1u, A2u, E1g...E(nir/2-2)g,
    //                                               E1u...E(nir/2-2)u
    // for even n, the irreps are A1, A2, B1, B2, E1...E(nir-4)
    //
    {
      IrreducibleRepresentation ir1(g,1,(nt%2) ? "A1g" : "A1");
      IrreducibleRepresentation ir2(g,1,(nt%2) ? "A2g" : "A2");
      IrreducibleRepresentation ir3(g,1,(nt%2) ? "A1u" : "B1");
      IrreducibleRepresentation ir4(g,1,(nt%2) ? "A2u" : "B2");

      for (gi=0; gi < 2*nt; gi++) {
        // Sn
        ir1.rep[gi] = ir1.proj[0][gi] = 1.0;
        ir2.rep[gi] = ir2.proj[0][gi] = 1.0;
        ir3.rep[gi] = ir3.proj[0][gi] = (gi%2) ? -1.0 : 1.0;
        ir4.rep[gi] = ir4.proj[0][gi] = (gi%2) ? -1.0 : 1.0;

        // n C2's and n sigma's
        ir1.rep[gi+2*nt] = ir1.proj[0][gi+2*nt] = 1.0;
        ir2.rep[gi+2*nt] = ir2.proj[0][gi+2*nt] = -1.0;
        ir3.rep[gi+2*nt] = ir3.proj[0][gi+2*nt] = (gi < nt) ? 1.0 : -1.0;
        ir4.rep[gi+2*nt] = ir4.proj[0][gi+2*nt] = (gi < nt) ? -1.0 : 1.0;
      }

      gamma_[0] = ir1;
      gamma_[1] = ir2;

      if (nt%2) {
        gamma_[nirrep_/2] = ir3;
        gamma_[1+nirrep_/2] = ir4;
        i=2;
      } else {
        gamma_[2] = ir3;
        gamma_[3] = ir4;
        i=4;
      }
    }
    
    // for the E irreps, the projection operators are:
    // for the 2n Sn's
    //   Ei xx = cos(m*theta*i) m = 0-(nt-1)
    //      xy = -sin(m*theta*i)
    //      yx = sin(m*theta*i)
    //      yy = cos(m*theta*i)
    //
    // for the n C2's:
    //   Ei xx = cos(m*theta*i) m = 0-(nt-1)
    //      xy = -sin(m*theta*i)
    //      yx = -sin(m*theta*i)
    //      yy = -cos(m*theta*i)
    //
    // for the n sigma's:
    //   Ei xx = cos(m*theta*i) m = 0-(nt-1)
    //      xy = sin(m*theta*i)
    //      yx = sin(m*theta*i)
    //      yy = -cos(m*theta*i)
    ei=1;
    itheta=0.5*theta;

    if (nt%2) {
      for (; i < nirrep_/2 ; i++, ei++, itheta += 0.5*theta) {
        if (nt==3)
          sprintf(label,"Eg");
        else
          sprintf(label,"E%dg",ei);

        IrreducibleRepresentation ir1(g,2,label);

        if (nt==3)
          sprintf(label,"Eu");
        else
          sprintf(label,"E%du",ei);

        IrreducibleRepresentation ir2(g,2,label);

        jitheta=0;
        double ineg = -pow(-1.0,(double)ei);
        for (j=0; j < 2*nt; j++, jitheta += itheta) {
          double ci1 = (j%2) ? -ineg : 1.0;
          double ci2 = (j%2) ? ineg : 1.0;
          
          ctheta = cos(jitheta);
          stheta = sin(jitheta);
          
          // Sn's
          ir1.rep[j] = 2.0*ci1*ctheta;

          ir1.proj[0][j] = ci1*ctheta;
          ir1.proj[1][j] = -ci1*stheta;
          ir1.proj[2][j] = ci1*stheta;
          ir1.proj[3][j] = ci1*ctheta;

          ir2.rep[j] = 2.0*ci2*ctheta;

          ir2.proj[0][j] = ci2*ctheta;
          ir2.proj[1][j] = -ci2*stheta;
          ir2.proj[2][j] = ci2*stheta;
          ir2.proj[3][j] = ci2*ctheta;
        }
        
        jitheta = 0;
        for (j=0; j < nt; j++, jitheta += ei*theta) {
          ctheta = cos(jitheta);
          stheta = sin(jitheta);

          // C2's
          ir1.rep[j+2*nt] = 0.0;

          ir1.proj[0][j+2*nt] = ctheta;
          ir1.proj[1][j+2*nt] = -stheta;
          ir1.proj[2][j+2*nt] = -stheta;
          ir1.proj[3][j+2*nt] = -ctheta;

          ir2.rep[j+2*nt] = 0.0;

          ir2.proj[0][j+2*nt] = ctheta;
          ir2.proj[1][j+2*nt] = -stheta;
          ir2.proj[2][j+2*nt] = -stheta;
          ir2.proj[3][j+2*nt] = -ctheta;

          // sigma d's
          ctheta = cos(jitheta + ei*0.5*theta);
          stheta = sin(jitheta + ei*0.5*theta);

          ir1.rep[j+3*nt] = 0.0;

          ir1.proj[0][j+3*nt] = -ineg*ctheta;
          ir1.proj[1][j+3*nt] = -ineg*stheta;
          ir1.proj[2][j+3*nt] = -ineg*stheta;
          ir1.proj[3][j+3*nt] = ineg*ctheta;

          ir2.rep[j+3*nt] = 0.0;

          ir2.proj[0][j+3*nt] = ineg*ctheta;
          ir2.proj[1][j+3*nt] = ineg*stheta;
          ir2.proj[2][j+3*nt] = ineg*stheta;
          ir2.proj[3][j+3*nt] = -ineg*ctheta;
        }

        gamma_[i] = ir1;
        gamma_[i+nirrep_/2] = ir2;
      }
    } else {
      for (; i < nirrep_; i++, ei++, itheta += 0.5*theta) {
        if (nt==2)
          sprintf(label,"E");
        else
          sprintf(label,"E%d",ei);

        IrreducibleRepresentation ir(g,2,label);

        jitheta=0;
        for (j=0; j < 2*nt; j++, jitheta += itheta) {
          ctheta = cos(jitheta);
          stheta = sin(jitheta);
          
          ir.rep[j] = 2.0*ctheta;

          ir.proj[0][j] = ctheta;
          ir.proj[1][j] = -stheta;
          ir.proj[2][j] = stheta;
          ir.proj[3][j] = ctheta;
        }
        
        jitheta=0;
        for (j=0; j < nt; j++, jitheta += ei*theta) {
          ctheta = cos(jitheta);
          stheta = sin(jitheta);
          
          // C2's
          ir.rep[j+2*nt] = 0.0;

          ir.proj[0][j+2*nt] = ctheta;
          ir.proj[1][j+2*nt] = -stheta;
          ir.proj[2][j+2*nt] = -stheta;
          ir.proj[3][j+2*nt] = -ctheta;

          // sigma's
          ctheta = cos(jitheta + ei*0.5*theta);
          stheta = sin(jitheta + ei*0.5*theta);
          
          ir.rep[j+3*nt] = 0.0;

          ir.proj[0][j+3*nt] = ctheta;
          ir.proj[1][j+3*nt] = stheta;
          ir.proj[2][j+3*nt] = stheta;
          ir.proj[3][j+3*nt] = -ctheta;
        }

        gamma_[i] = ir;
      }
    }

    // Sn's
    itheta = 0;
    for (i=0; i < 2*nt ; i++, itheta += 0.5*theta) {
      symop[i][0][0] = symop[i][1][1] = cos(itheta);
      symop[i][0][1] = sin(itheta);
      symop[i][1][0] = -sin(itheta);
      symop[i][2][2] = (i%2) ? -1.0 : 1.0;

      trans[i] = symop[i].trace();
      rot[i] = (i%2) ? -trans[i] : trans[i];
    }

    // C2's
    itheta = 0;
    for (i=0,j=2*nt; i < nt ; i++, j++, itheta += theta) {
      symop[j][0][0] = cos(itheta);
      symop[j][1][1] = -cos(itheta);
      symop[j][1][0] = symop[j][0][1] = -sin(itheta);
      symop[j][2][2] = -1.0;

      rot[j] = trans[j] = symop[j].trace();
    }

    // sigma's
    itheta = 0.5*theta;
    for (i=0,j=3*nt; i < nt ; i++, j++, itheta += theta) {
      symop[j][0][0] = cos(itheta);
      symop[j][1][1] = -cos(itheta);
      symop[j][1][0] = symop[j][0][1] = sin(itheta);
      symop[j][2][2] = 1.0;

      trans[j] = symop[j].trace();
      rot[j] = -trans[j];
    }

    break;

  case DNH:
    // clockwise rotation and rotation-reflection about z axis,
    // followed by c2 about x axis and then reflection
    // through xz 

    // d2h is a special case
    if (nt==2) {
      IrreducibleRepresentation ir1(g,1,"Ag");
      IrreducibleRepresentation ir2(g,1,"B1g");
      IrreducibleRepresentation ir3(g,1,"B2g");
      IrreducibleRepresentation ir4(g,1,"B3g");

      IrreducibleRepresentation ir5(g,1,"Au");
      IrreducibleRepresentation ir6(g,1,"B1u");
      IrreducibleRepresentation ir7(g,1,"B2u");
      IrreducibleRepresentation ir8(g,1,"B3u");

      for (gi=0; gi < 8; gi++) {
        ir1.rep[gi] = ir1.proj[0][gi] = 1.0;
        ir2.rep[gi] = ir2.proj[0][gi] = (gi < 4) ? 1.0 : -1.0;

        ir3.rep[gi] = ir3.proj[0][gi] = (gi%4==0 || gi%4==3) ?
          ((gi < 4) ? 1.0 : -1.0) : ((gi < 4) ? -1.0 : 1.0);
        ir4.rep[gi] = ir4.proj[0][gi] = (gi%4==0 || gi%4==3) ? 1.0 : -1.0;
        
        ir5.rep[gi] = ir5.proj[0][gi] = (gi%4==0 || gi%4==1) ? 1.0 : -1.0;
        ir6.rep[gi] = ir6.proj[0][gi] = (gi%4==0 || gi%4==1) ?
          ((gi < 4) ? 1.0 : -1.0) : ((gi < 4) ? -1.0 : 1.0);

        ir7.rep[gi] = ir7.proj[0][gi] = (gi < 4) ?
          ((gi%2) ? -1.0 : 1.0) : ((gi%2) ? 1.0 : -1.0);
        ir8.rep[gi] = ir8.proj[0][gi] = (gi%2) ? -1.0 : 1.0;
      }

      gamma_[0] = ir1;
      gamma_[1] = ir2;
      gamma_[2] = ir3;
      gamma_[3] = ir4;
      gamma_[4] = ir5;
      gamma_[5] = ir6;
      gamma_[6] = ir7;
      gamma_[7] = ir8;

    } else {
      {
        IrreducibleRepresentation ir1(g,1, (nt%2) ? "A1'" : "A1g");
        IrreducibleRepresentation ir2(g,1, (nt%2) ? "A2'" : "A2g");
        IrreducibleRepresentation ir3(g,1, (nt%2) ? "A1\"" : "A1u");
        IrreducibleRepresentation ir4(g,1, (nt%2) ? "A2\"" : "A2u");

        for (gi=0; gi < nt; gi++) {
          // n Cn's
          ir1.rep[gi] = ir1.proj[0][gi] = 1.0;
          ir2.rep[gi] = ir2.proj[0][gi] = 1.0;
          ir3.rep[gi] = ir3.proj[0][gi] = 1.0;
          ir4.rep[gi] = ir4.proj[0][gi] = 1.0;

          // n Sn's
          ir1.rep[gi+nt] = ir1.proj[0][gi+nt] = 1.0;
          ir2.rep[gi+nt] = ir2.proj[0][gi+nt] = 1.0;
          ir3.rep[gi+nt] = ir3.proj[0][gi+nt] = -1.0;
          ir4.rep[gi+nt] = ir4.proj[0][gi+nt] = -1.0;

          // n C2's
          ir1.rep[gi+2*nt] = ir1.proj[0][gi+2*nt] = 1.0;
          ir2.rep[gi+2*nt] = ir2.proj[0][gi+2*nt] = -1.0;
          ir3.rep[gi+2*nt] = ir3.proj[0][gi+2*nt] = 1.0;
          ir4.rep[gi+2*nt] = ir4.proj[0][gi+2*nt] = -1.0;

          // n sigma's
          ir1.rep[gi+3*nt] = ir1.proj[0][gi+3*nt] = 1.0;
          ir2.rep[gi+3*nt] = ir2.proj[0][gi+3*nt] = -1.0;
          ir3.rep[gi+3*nt] = ir3.proj[0][gi+3*nt] = -1.0;
          ir4.rep[gi+3*nt] = ir4.proj[0][gi+3*nt] = 1.0;
        }

        gamma_[0] = ir1;
        gamma_[1] = ir2;
        gamma_[nirrep_/2] = ir3;
        gamma_[1+nirrep_/2] = ir4;
          
        i=2;
      }

      if (!(nt%2)) {
        IrreducibleRepresentation ir1(g,1,"B1g");
        IrreducibleRepresentation ir2(g,1,"B2g");
        IrreducibleRepresentation ir3(g,1,"B1u");
        IrreducibleRepresentation ir4(g,1,"B2u");

        for (gi=0; gi < nt; gi++) {
          // n Cn's
          ir1.rep[gi] = ir1.proj[0][gi] = (gi%2) ? -1.0 : 1.0;
          ir2.rep[gi] = ir2.proj[0][gi] = (gi%2) ? -1.0 : 1.0;
          ir3.rep[gi] = ir3.proj[0][gi] = (gi%2) ? -1.0 : 1.0;
          ir4.rep[gi] = ir4.proj[0][gi] = (gi%2) ? -1.0 : 1.0;

          // n Sn's
          ir1.rep[gi+nt] = ir1.proj[0][gi+nt] = (gi%2) ? -1.0 : 1.0;
          ir2.rep[gi+nt] = ir2.proj[0][gi+nt] = (gi%2) ? -1.0 : 1.0;
          ir3.rep[gi+nt] = ir3.proj[0][gi+nt] = (gi%2) ? 1.0 : -1.0;
          ir4.rep[gi+nt] = ir4.proj[0][gi+nt] = (gi%2) ? 1.0 : -1.0;

          // n C2's
          ir1.rep[gi+2*nt] = ir1.proj[0][gi+2*nt] = (gi%2) ? -1.0 : 1.0;
          ir2.rep[gi+2*nt] = ir2.proj[0][gi+2*nt] = (gi%2) ? 1.0 : -1.0;
          ir3.rep[gi+2*nt] = ir3.proj[0][gi+2*nt] = (gi%2) ? -1.0 : 1.0;
          ir4.rep[gi+2*nt] = ir4.proj[0][gi+2*nt] = (gi%2) ? 1.0 : -1.0;

          // n sigma's
          ir1.rep[gi+3*nt] = ir1.proj[0][gi+3*nt] = (gi%2) ? -1.0 : 1.0;
          ir2.rep[gi+3*nt] = ir2.proj[0][gi+3*nt] = (gi%2) ? 1.0 : -1.0;
          ir3.rep[gi+3*nt] = ir3.proj[0][gi+3*nt] = (gi%2) ? 1.0 : -1.0;
          ir4.rep[gi+3*nt] = ir4.proj[0][gi+3*nt] = (gi%2) ? -1.0 : 1.0;
        }
        
        gamma_[2] = ir1;
        gamma_[3] = ir2;
        gamma_[2+nirrep_/2] = ir3;
        gamma_[3+nirrep_/2] = ir4;

        i=4;
      }
      
      ei=1;
      itheta=theta;
      for (; i < nirrep_/2 ; i++, ei++, itheta += theta) {
        if (nt==3)
          sprintf(label,"E'");
        else if (nt==4)
          sprintf(label,"Eg");
        else
          sprintf(label,"E%d%s", ei, (nt%2) ? "'" : "g");

        IrreducibleRepresentation ir1(g,2,label);

        if (nt==3)
          sprintf(label,"E\"");
        else if (nt==4)
          sprintf(label,"Eu");
        else
          sprintf(label,"E%d%s", ei, (nt%2) ? "\"" : "u");

        IrreducibleRepresentation ir2(g,2,label);

        jitheta=0;
        double ineg = (nt%2) ? -1.0 : -pow(-1.0,(double)ei);

        for (j=0; j < nt; j++, jitheta += itheta) {
          ctheta = cos(jitheta);
          stheta = sin(jitheta);
          
          // Cn's
          ir1.rep[j] = 2.0*ctheta;

          ir1.proj[0][j] = ctheta;
          ir1.proj[1][j] = -stheta;
          ir1.proj[2][j] = stheta;
          ir1.proj[3][j] = ctheta;

          ir2.rep[j] = 2.0*ctheta;

          ir2.proj[0][j] = ctheta;
          ir2.proj[1][j] = -stheta;
          ir2.proj[2][j] = stheta;
          ir2.proj[3][j] = ctheta;

          // Sn's
          ir1.rep[j+nt] = -2.0*ineg*ctheta;

          ir1.proj[0][j+nt] = -ineg*ctheta;
          ir1.proj[1][j+nt] = ineg*stheta;
          ir1.proj[2][j+nt] = -ineg*stheta;
          ir1.proj[3][j+nt] = -ineg*ctheta;

          ir2.rep[j+nt] = 2.0*ineg*ctheta;

          ir2.proj[0][j+nt] = ineg*ctheta;
          ir2.proj[1][j+nt] = -ineg*stheta;
          ir2.proj[2][j+nt] = ineg*stheta;
          ir2.proj[3][j+nt] = ineg*ctheta;

          // C2's
          ir1.rep[j+2*nt] = 0.0;

          ir1.proj[0][j+2*nt] = ctheta;
          ir1.proj[1][j+2*nt] = -stheta;
          ir1.proj[2][j+2*nt] = -stheta;
          ir1.proj[3][j+2*nt] = -ctheta;

          ir2.rep[j+2*nt] = 0.0;

          ir2.proj[0][j+2*nt] = ctheta;
          ir2.proj[1][j+2*nt] = -stheta;
          ir2.proj[2][j+2*nt] = -stheta;
          ir2.proj[3][j+2*nt] = -ctheta;

          // sigma's
          ir1.rep[j+3*nt] = 0.0;

          ir1.proj[0][j+3*nt] = -ineg*ctheta;
          ir1.proj[1][j+3*nt] = -ineg*stheta;
          ir1.proj[2][j+3*nt] = -ineg*stheta;
          ir1.proj[3][j+3*nt] = ineg*ctheta;

          ir2.rep[j+3*nt] = 0.0;

          ir2.proj[0][j+3*nt] = ineg*ctheta;
          ir2.proj[1][j+3*nt] = ineg*stheta;
          ir2.proj[2][j+3*nt] = ineg*stheta;
          ir2.proj[3][j+3*nt] = -ineg*ctheta;
        }
        
        gamma_[i] = ir1;
        gamma_[i+nirrep_/2] = ir2;
      }
    }

    itheta=0;
    for (i=0; i < nt ; i++, itheta += theta) {
      ctheta = cos(itheta);
      stheta = sin(itheta);
      
      // Cn's
      symop[i][0][0] = symop[i][1][1] = ctheta;
      symop[i][0][1] = stheta;
      symop[i][1][0] = -stheta;
      symop[i][2][2] = 1.0;
      
      rot[i] = trans[i] = symop[i].trace();

      // Sn's
      symop[i+nt][0][0] = symop[i+nt][1][1] = ctheta;
      symop[i+nt][0][1] = stheta;
      symop[i+nt][1][0] = -stheta;
      symop[i+nt][2][2] = -1.0;

      trans[i+nt] = symop[i+nt].trace();
      rot[i+nt] = -trans[i+nt];

      // C2's
      symop[i+2*nt][0][0] = ctheta;
      symop[i+2*nt][1][1] = -ctheta;
      symop[i+2*nt][1][0] = symop[i+2*nt][0][1] = -stheta;
      symop[i+2*nt][2][2] = -1.0;

      rot[i+2*nt] = trans[i+2*nt] = symop[i+2*nt].trace();

      // sigma's
      symop[i+3*nt][0][0] = ctheta;
      symop[i+3*nt][1][1] = -ctheta;
      symop[i+3*nt][1][0] = symop[i+3*nt][0][1] = stheta;
      symop[i+3*nt][2][2] = 1.0;

      trans[i+3*nt] = symop[i+3*nt].trace();
      rot[i+3*nt] = -trans[i+3*nt];
    }

    break;

  case T:
    t();
    break;

  case TH:
    th();
    break;

  case TD:
    td();
    break;

  case O:
    o();
    break;

  case OH:
    oh();
    break;

  case I:
    this->i();
    break;

  case IH:
    ih();
    break;

  default:
    return -1;

  }
    
/* ok, we have the reducible representation of the rotations and translations,
 * now let's use projection operators to find out how many rotations and
 * translations there are in each irrep
 */

  if (pg != C1 && pg != CI && pg != CS && pg != T && pg != TD && pg != TH &&
      pg != O && pg != OH && pg != I && pg != IH) {
    for (i=0; i < nirrep_; i++) {
      double nr=0; double nt=0;

      for (j=0; j < gamma_[i].g; j++) {
        nr += rot[j]*gamma_[i].rep[j];
        nt += trans[j]*gamma_[i].rep[j];
      }

      gamma_[i].nrot_ = (int) ((nr+0.5)/gamma_[i].g);
      gamma_[i].ntrans_ = (int) ((nt+0.5)/gamma_[i].g);
    }
  }

  delete[] rot;
  delete[] trans;
  
  return 0;
}
