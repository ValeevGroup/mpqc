
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

  
  symop = new SymmetryOperation[g];
  SymmetryOperation so;

  // this array forms a reducible representation for rotations about x,y,z
  double *rot = new double[g];
  memset(rot,'\0',sizeof(double)*g);

  // this array forms a reducible representation for translations along x,y,z
  double *trans = new double[g];
  memset(trans,'\0',sizeof(double)*g);

  // the angle to rotate about the principal axis
  double theta = (nt) ? 2.0*M_PI/nt : 2.0*M_PI;

  switch (pg) {

  case C1:
    // no symmetry case
    gamma_[0].init(1,1,"A");
    gamma_[0].nrot_ = 3;
    gamma_[0].ntrans_ = 3;
    gamma_[0].rep[0][0][0] = 1.0;

    symop[0].unit();

    break;

  case CI:
    // equivalent to S2 about the z axis
    gamma_[0].init(2,1,"Ag");
    gamma_[0].rep[0][0][0] = 1.0;
    gamma_[0].rep[1][0][0] = 1.0;
    gamma_[0].nrot_=3;

    gamma_[1].init(2,1,"Au");
    gamma_[1].rep[0][0][0] =  1.0;
    gamma_[1].rep[1][0][0] = -1.0;
    gamma_[1].ntrans_=3;

    symop[0].unit();
    symop[1][0][0] = symop[1][1][1] = symop[1][2][2] = -1;

    break;

  case CS: // reflection through the xy plane
    gamma_[0].init(2,1,"A'");
    gamma_[0].rep[0][0][0] = 1.0;
    gamma_[0].rep[1][0][0] = 1.0;
    gamma_[0].nrot_=1;
    gamma_[0].ntrans_=2;

    gamma_[1].init(2,1,"A\"");
    gamma_[1].rep[0][0][0] =  1.0;
    gamma_[1].rep[1][0][0] = -1.0;
    gamma_[1].nrot_=2;
    gamma_[1].ntrans_=1;

    symop[0].unit();
    symop[1].unit();
    symop[1][2][2] = -1;

    break;

  case CN:
    // clockwise rotation about z axis by theta*i radians
    //
    // for odd n, the irreps are A and E1...E(nir-1)
    // for even n, the irreps are A, B, and E1...E(nir-2)
    //
    gamma_[0].init(g,1,"A");
    for (gi=0; gi < g; gi++)
      gamma_[0].rep[gi][0][0] = 1.0;

    i=1;

    if (!(nt%2)) {
      gamma_[1].init(g,1,"B");
      for (gi=0; gi < g; gi++)
        gamma_[1].rep[gi][0][0] = (gi%2) ? -1.0 : 1.0;

      i++;
    }

    // for the E irreps, the projection operators are:
    //   Ei xx = cos(m*theta*i) m = 0-(nt-1)
    //      xy = -sin(m*theta*i)
    //      yx = sin(m*theta*i)
    //      yy = cos(m*theta*i)

    ei=1;
    for (; i < nirrep_; i++, ei++) {
      IrreducibleRepresentation& ir = gamma_[i];

      if (nt==3 || nt==4)
        sprintf(label,"E");
      else
        sprintf(label,"E%d",ei);

      ir.init(g,2,label);
      ir.complex_=1;

      // identity
      ir.rep[0].unit();

      // Cn
      ir.rep[1].rotation(ei*theta);

      // the other n-1 Cn's
      for (j=2; j < g; j++)
        ir.rep[j] = ir.rep[j-1].operate(ir.rep[1]);
    }

    // identity
    symop[0].unit();

    // Cn
    symop[1].rotation(theta);
    
    // the other n-2 Cn's
    for (i=2; i < nt; i++)
      symop[i] = symop[i-1].operate(symop[1]);
    
    for (i=0; i < nt ; i++)
      rot[i] = trans[i] = symop[i].trace();

    break;

  case CNV:
    // clockwise rotation about z axis by theta*i radians, then
    // reflect through the xz plane
    //
    // for odd n, the irreps are A1, A2, and E1...E(nir-2)
    // for even n, the irreps are A1, A2, B1, B2, and E1...E(nir-4)
    //

    gamma_[0].init(g,1,"A1");
    gamma_[1].init(g,1,"A2");

    for (gi=0; gi < nt; gi++) {
      // Cn's
      gamma_[0].rep[gi][0][0] = 1.0;
      gamma_[1].rep[gi][0][0] = 1.0;

      // sigma's
      gamma_[0].rep[gi+nt][0][0] =  1.0;
      gamma_[1].rep[gi+nt][0][0] = -1.0;
    }

    i=2;

    if (!(nt%2)) {
      gamma_[2].init(g,1,"B1");
      gamma_[3].init(g,1,"B2");

      for (gi=0; gi < nt ; gi++) {
        double ci = (gi%2) ? -1.0 : 1.0;
        
        // Cn's
        gamma_[2].rep[gi][0][0] = ci;
        gamma_[3].rep[gi][0][0] = ci;

        // sigma's
        gamma_[2].rep[gi+nt][0][0] =  ci;
        gamma_[3].rep[gi+nt][0][0] = -ci;
      }

      i=4;
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
    for (; i < nirrep_; i++, ei++) {
      IrreducibleRepresentation& ir = gamma_[i];

      char lab[4];
      if (nt==3 || nt==4)
        sprintf(lab,"E");
      else
        sprintf(lab,"E%d",ei);

      ir.init(g,2,lab);

      // identity
      ir.rep[0].unit();

      // Cn
      ir.rep[1].rotation(ei*theta);

      // the other n-2 Cn's
      for (j=2; j < nt; j++)
        ir.rep[j] = ir.rep[j-1].operate(ir.rep[1]);
      
      // sigma xz
      ir.rep[nt][0][0] =  1.0;
      ir.rep[nt][1][1] = -1.0;

      SymRep sr(2);
      sr.rotation(ei*theta/2.0);
      
      // the other n-1 sigma's
      for (j=nt+1; j < g; j++)
        ir.rep[j] = ir.rep[j-1].sim_transform(sr);
    }

    // identity
    symop[0].unit();
    
    // Cn
    symop[1].rotation(theta);
    
    // the other n-2 Cn's
    for (i=2; i < nt; i++)
      symop[i] = symop[i-1].operate(symop[1]);
    
    // sigma xz
    symop[nt][0][0] =  1.0;
    symop[nt][1][1] = -1.0;
    symop[nt][2][2] =  1.0;

    so.rotation(theta/2.0);

    // the other n-1 sigma's
    for (j=nt+1; j < g; j++)
      symop[j] = symop[j-1].sim_transform(so);

    for (i=0; i < nt ; i++) {
      rot[i] = trans[i] = symop[i].trace();

      rot[i+nt] = -symop[i+nt].trace();
      trans[i+nt] = symop[i+nt].trace();
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
    gamma_[0].init(g,1, (nt%2) ? "A'" : "Ag");
    gamma_[nirrep_/2].init(g,1, (nt%2) ? "A\"" : "Au");

    for (gi=0; gi < nt; gi++) {
      // Cn's
      gamma_[0].rep[gi][0][0] = 1.0;
      gamma_[nirrep_/2].rep[gi][0][0] = 1.0;

      // Sn's
      gamma_[0].rep[gi+nt][0][0] = 1.0;
      gamma_[nirrep_/2].rep[gi+nt][0][0] = -1.0;
    }

    i=1;

    if (!(nt%2)) {
      double ineg = ((nt/2)%2) ? -1.0 : 1.0;
      
      gamma_[1].init(g,1,"Bg");
      gamma_[1+nirrep_/2].init(g,1,"Bu");

      for (gi=0; gi < nt; gi++) {
        double ci = (gi%2) ? -1.0 : 1.0;

        // Cn's
        gamma_[1].rep[gi][0][0] = ci;
        gamma_[1+nirrep_/2].rep[gi][0][0] = ci;
      
        // Sn's
        gamma_[1].rep[gi+nt][0][0] =  ci*ineg;
        gamma_[1+nirrep_/2].rep[gi+nt][0][0] = -ci*ineg;
      }

      i=2;
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
    for (; i < nirrep_/2 ; i++, ei++) {
      IrreducibleRepresentation& ir1 = gamma_[i];
      IrreducibleRepresentation& ir2 = gamma_[i+nirrep_/2];

      if (nt==3 || nt==4)
        sprintf(label,(nt%2) ? "E'" : "Eg");
      else
        sprintf(label,"E%d%s", ei, (nt%2) ? "'" : "g");

      ir1.init(g,2,label);

      if (nt==3 || nt==4)
        sprintf(label,(nt%2) ? "E\"" : "Eu");
      else
        sprintf(label,"E%d%s", ei, (nt%2) ? "\"" : "u");

      ir2.init(g,2,label);

      ir1.complex_=1;
      ir2.complex_=1;

      // identity
      ir1.rep[0].unit();
      ir2.rep[0].unit();

      // Cn
      ir1.rep[1].rotation(ei*theta);
      ir2.rep[1].rotation(ei*theta);

      for (j=2; j < nt; j++) {
        ir1.rep[j] = ir1.rep[j-1].operate(ir1.rep[1]);
        ir2.rep[j] = ir2.rep[j-1].operate(ir2.rep[1]);
      }

      double ineg = (nt%2) ? 1.0 : pow(-1.0,(double)ei);
      for (j=nt; j < g; j++) {
        for (int ri=0; ri < 2; ri++) {
          for (int rj=0; rj < 2; rj++) {
            ir1.rep[j][ri][rj] =  ineg*ir1.rep[j-nt][ri][rj];
            ir2.rep[j][ri][rj] = -ineg*ir2.rep[j-nt][ri][rj];
          }
        }
      }
    }

    // identity
    symop[0].unit();

    // Cn
    symop[1].rotation(theta);
    
    // the other n-2 Cn's
    for (i=2; i < nt; i++)
      symop[i] = symop[i-1].operate(symop[1]);

    for (i=0; i < nt ; i++) {
      symop[i+nt] = symop[i];
      symop[i+nt][2][2] = -1.0;
      
      rot[i] = trans[i] = symop[i].trace();
      trans[i+nt] = symop[i+nt].trace();
      rot[i+nt] = -trans[i+nt];
    }

    break;

#if 0
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

#endif
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
        nr += rot[j]*gamma_[i].character(j);
        nt += trans[j]*gamma_[i].character(j);
      }

      gamma_[i].nrot_ = (int) ((nr+0.5)/gamma_[i].g);
      gamma_[i].ntrans_ = (int) ((nt+0.5)/gamma_[i].g);
    }
  }

  delete[] rot;
  delete[] trans;
  
  return 0;
}
