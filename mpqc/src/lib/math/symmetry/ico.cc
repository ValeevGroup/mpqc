
#include <math.h>
#include <string.h>

#include <math/symmetry/pointgrp.h>
#include <math/symmetry/tform.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// returns r1xr2
static SymmetryOperation
operate(SymmetryOperation& r1, SymmetryOperation& r2)
{
  SymmetryOperation ret;

  for (int i=0; i < 3; i++) {
    for (int j=0; j < 3; j++) {
      double t=0;
      for (int k=0; k < 3; k++)
        t += r1[i][k] * r2[k][j];
      ret[i][j] = t;
    }
  }
  
  return ret;
}

// given matrix R in so and E in frame, perform similarity transform
// R' = ~E * R * E
static SymmetryOperation
sim_transform(SymmetryOperation& so, const SymmetryOperation& frame)
{
  int i,j,k;
  SymmetryOperation foo,ret;

  // foo = ~E*R 
  for (i=0; i < 3; i++) {
    for (j=0; j < 3; j++) {
      double t=0;
      for (k=0; k < 3; k++) t += frame(i,k)*so(k,j);
      foo(i,j) = t;
    }
  }

  // R' = ~E*R*E = foo*E
  for (i=0; i < 3; i++) {
    for (j=0; j < 3; j++) {
      double t=0;
      for (k=0; k < 3; k++) t += foo(i,k)*frame(j,k);
      ret(i,j) = t;
    }
  }

  return ret;
}

// these are the operations which make up T
static void
i_ops(SymmetryOperation *symop, SymmetryOperation *rep)
{
  int i;
  
  // identity
  symop[0][0][0] = symop[0][1][1] = symop[0][2][2] = 1.0;
  rep[0] = symop[0];

  //
  // 12 C5's
  //
  // first the 2 C5's about the z axis
  symop[1][0][0] =  cos(2.0*M_PI/5.0);
  symop[1][0][1] =  sin(2.0*M_PI/5.0);
  symop[1][1][0] = -sin(2.0*M_PI/5.0);
  symop[1][1][1] =  cos(2.0*M_PI/5.0);
  symop[1][2][2] =  1.0;
  
  symop[2][0][0] =  cos(2.0*M_PI/5.0);
  symop[2][0][1] = -sin(2.0*M_PI/5.0);
  symop[2][1][0] =  sin(2.0*M_PI/5.0);
  symop[2][1][1] =  cos(2.0*M_PI/5.0);
  symop[2][2][2] =  1.0;
  
  rep[1] = operate(symop[1],symop[1]);
  rep[2] = operate(symop[2],symop[2]);
  rep[13] = symop[2];
  rep[14] = symop[1];

  // now rotate the first C5's by 2pi/3 degrees about the zx axis (sort of)
  SymmetryOperation so;

  double cosd = sin(2.0*M_PI/5.0)/((1-cos(2.0*M_PI/5.0))*sqrt(3.0));
  double cosd2 = cosd*cosd;
  double sind2 = 1 - cosd2;
  double sind = sqrt(sind2);

  so.zero();
  so[0][0] =  1 - 1.5*cosd2;
  so[1][0] =  0.5*sqrt(3.0)*cosd;
  so[2][0] =  1.5*cosd*sind;
  so[0][1] = -0.5*sqrt(3.0)*cosd;
  so[1][1] = -0.5;
  so[2][1] =  0.5*sqrt(3.0)*sind;
  so[0][2] =  1.5*cosd*sind;
  so[1][2] = -0.5*sqrt(3.0)*sind;
  so[2][2] =  1 - 1.5*sind2;
  
  symop[3] = sim_transform(symop[1],so);
  symop[4] = sim_transform(symop[2],so);

  // rotate twice to get the first one aligned along the x axis
  symop[3] = sim_transform(symop[3],symop[1]);
  symop[3] = sim_transform(symop[3],symop[1]);
  symop[4] = sim_transform(symop[4],symop[1]);
  symop[4] = sim_transform(symop[4],symop[1]);

  rep[3] = operate(symop[4],symop[4]);
  rep[4] = operate(symop[3],symop[3]);
  rep[15] = symop[3];
  rep[16] = symop[4];
  
  // and then rotate those by 2pi/5 about the z axis 4 times
  for (i=5; i < 13; i++) {
    symop[i] = sim_transform(symop[i-2], symop[1]);
    rep[i] = sim_transform(rep[i-2], rep[1]);
    rep[i+12] = sim_transform(rep[i+10], rep[1]);
  }

  //
  // 12 C5^2's
  //
  // get these from operating on each of the C5's with itself
  for (i=13; i < 25; i++)
    symop[i] = operate(symop[i-12],symop[i-12]);

  //
  // 20 C3's
  //
  // first we have 2 C3's about the zx axis
  symop[25] = so;
  symop[26] = operate(so,so);
  
  // and then rotate those by 2pi/5 about the z axis 4 times
  for (i=27; i < 35; i++)
    symop[i] = sim_transform(symop[i-2], symop[1]);

  // now rotate one of the above C3's by 2pi/3 about the zx axis
  symop[35] = sim_transform(symop[27],so);
  symop[36] = sim_transform(symop[28],so);

  // and then rotate those by 2pi/5 about the z axis 4 times
  for (i=37; i < 45; i++)
    symop[i] = sim_transform(symop[i-2], symop[1]);

  rep[25] = symop[35];
  rep[26] = symop[36];
  
  for (i=27; i < 35; i++)
    rep[i] = sim_transform(rep[i-2],rep[1]);
  
  rep[35] = symop[26];
  rep[36] = symop[25];
  
  for (i=37; i < 45; i++)
    rep[i] = sim_transform(rep[i-2],rep[1]);

  //
  // 15 C2's
  //
  // first we have a C2 about the y axis
  symop[45][0][0] = -1.0;
  symop[45][1][1] =  1.0;
  symop[45][2][2] = -1.0;

  rep[45] = symop[45];
  
  // and rotate that by 2pi/5 about the z axis 4 times
  for (i=46; i < 50; i++) {
    symop[i] = sim_transform(symop[i-1], symop[1]);
    rep[i] = sim_transform(rep[i-1],rep[1]);
  }

  // now take the C2 about the y axis and rotate it by 2pi/3 about the zx axis
  symop[50] = symop[45];
  symop[50] = sim_transform(symop[50],so);

  // align this c2 along the x axis
  symop[50] = sim_transform(symop[50],symop[2]);
  symop[50] = sim_transform(symop[50],symop[2]);

  // and rotate that by 2pi/5 about the z axis 4 times
  for (i=51; i < 55; i++)
    symop[i] = sim_transform(symop[i-1], symop[1]);

  // finally, take a C2 about the y axis, and rotate it by 2pi/3 about the
  // xz axis
  symop[55] = symop[45];
  symop[55] = sim_transform(symop[55],symop[35]);
  symop[55] = sim_transform(symop[55],symop[1]);

  // and then rotate that by 2pi/5 about the z axis 4 times
  for (i=56; i < 60; i++)
    symop[i] = sim_transform(symop[i-1], symop[1]);

  rep[50] = symop[55];
  rep[55] = symop[50];
  
  for (i=51; i < 55; i++) {
    rep[i] = sim_transform(rep[i-1], rep[1]);
    rep[i+5] = sim_transform(rep[i+4], rep[1]);
  }
}

// this gives us the operations in Ih which come from ixI (ie, the inverse
// operating on all the symmetry operations from I).
static void
ih_ops(SymmetryOperation *symop)
{
  //i_ops(symop);

  for (int i=0; i < 60; i++)
    for (int j=0; j < 3; j++)
      for (int k=0; k < 3; k++)
        symop[i][j][k] *= -1.0;
}

void CharacterTable::i()
{
  SymmetryOperation *t2rep = new SymmetryOperation[60];
  
  // t_ops gives us all the symmetry operations we need
  i_ops(symop,t2rep);

  int i,j,k;

  {
    IrreducibleRepresentation ir(g,1,"A");
    for (i=0; i < g; i++) {
      ir.rep[i] = 1;
      ir.proj[0][i] = 1;
    }

    gamma_[0] = ir;
  }

  // the symmetry operation matrices give us a basis for irrep T1.
  {
    IrreducibleRepresentation ir1(g,3,"T1");
    IrreducibleRepresentation ir2(g,3,"T2");
    IrreducibleRepresentation irg(g,4,"G");
    
    ir1.nrot_ = 1;
    ir1.ntrans_ = 1;

    for (i=0; i < g; i++) {
      //SymRotation r(3,symop[i],1);
      SymRotation r(2,symop[i],1);
      
      ir1.rep[i]=symop[i].trace();
      ir2.rep[i]=t2rep[i].trace();
      //irg.rep[i]=r.trace()-t2rep[i].trace();
      irg.rep[i]=r.trace()-r[0][0];

      for (j=0; j < 3; j++) {
        for (k=0; k < 3; k++) {
          ir1.proj[3*j+k][i] = symop[i][k][j];
          ir2.proj[3*j+k][i] = t2rep[i][k][j];
        }
      }
    }

    gamma_[1] = ir1;
    gamma_[2] = ir2;
    gamma_[3] = irg;
  }

  // the pure d rotation matrices give us a basis for H
  {
    IrreducibleRepresentation ir(g,5,"H");

    for (i=0; i < g; i++) {
      SymRotation r(2,symop[i],1);
      
      ir.rep[i]=r.trace();
      
      for (j=0; j < 5; j++) {
        for (k=0; k < 5; k++) {
          ir.proj[5*j+k][i] = r[k][j];
        }
      }
    }
    
    gamma_[4] = ir;
  }
    
}


void CharacterTable::ih()
{
  // first get the ExT operations, then the ixT operations
  //i_ops(symop);
  ih_ops(&symop[60]);
  
  int i,j,k;
  {
    IrreducibleRepresentation ir1(g,1,"Ag");
    IrreducibleRepresentation ir2(g,1,"Au");
    for (i=0; i < 60; i++) {
      ir1.rep[i] = ir1.proj[0][i] = 1;
      ir2.rep[i] = ir2.proj[0][i] = 1;

      ir1.rep[i+60] = ir1.proj[0][i+60] = 1;
      ir2.rep[i+60] = ir2.proj[0][i+60] = -1;
    }

    gamma_[0] = ir1;
    gamma_[5] = ir2;
  }

  // the symmetry operation matrices form a basis for T1u.
  {
    IrreducibleRepresentation ir1(g,3,"T1g");
    IrreducibleRepresentation ir2(g,3,"T2g");
    IrreducibleRepresentation ir3(g,3,"T1u");
    IrreducibleRepresentation ir4(g,3,"T2u");

    ir1.nrot_ = 1;
    ir3.ntrans_ = 1;

    for (i=0; i < g; i++) {
      ir3.rep[i]=0;
      
      for (j=0; j < 3; j++) {
        for (k=0; k < 3; k++)
          ir3.proj[3*j+k][i] = symop[i][k][j];
        ir3.rep[i] += ir3.proj[3*j+j][i];
      }
    }

    for (i=0; i < g/2; i++) {
      ir1.rep[i] = ir3.rep[i];
      ir1.rep[i+60] = -ir3.rep[i+60];
      
      for (j=0; j < 9; j++) {
        ir1.proj[j][i] = ir3.proj[j][i];
        ir1.proj[j][i+60] = -ir3.proj[j][i+60];
      }
    }

    gamma_[1] = ir1;
    gamma_[2] = ir2;
    gamma_[6] = ir3;
    gamma_[7] = ir4;
  }

  {
    IrreducibleRepresentation ir1(g,4,"Gg");
    IrreducibleRepresentation ir2(g,4,"Gu");

    gamma_[3] = ir1;
    gamma_[8] = ir2;
  }

  {
    IrreducibleRepresentation ir1(g,5,"Hg");
    IrreducibleRepresentation ir2(g,5,"Hu");

    gamma_[4] = ir1;
    gamma_[9] = ir2;
  }
}
