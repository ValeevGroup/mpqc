
#include <math.h>
#include <string.h>

#include <math/symmetry/pointgrp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


// these are the operations which make up T
static void
t_ops(SymmetryOperation *symop)
{
  // identity
  symop[0][0][0] = symop[0][1][1] = symop[0][2][2] = 1.0;

  // a = ( 1, 1, 1)
  // b = (-1,-1, 1)
  // c = ( 1,-1,-1)
  // d = (-1, 1,-1)
  // C3 (a)
  symop[1][0][2] =  1.0;
  symop[1][1][0] =  1.0;
  symop[1][2][1] =  1.0;

  // C3 (b)
  symop[2][0][2] = -1.0;
  symop[2][1][0] =  1.0;
  symop[2][2][1] = -1.0;

  // C3 (c)
  symop[3][0][2] = -1.0;
  symop[3][1][0] = -1.0;
  symop[3][2][1] =  1.0;

  // C3 (d)
  symop[4][0][2] =  1.0;
  symop[4][1][0] = -1.0;
  symop[4][2][1] = -1.0;

  // C3^2 (a)
  symop[5][0][1] =  1.0;
  symop[5][1][2] =  1.0;
  symop[5][2][0] =  1.0;

  // C3^2 (b)
  symop[6][0][1] =  1.0;
  symop[6][1][2] = -1.0;
  symop[6][2][0] = -1.0;

  // C3^2 (c)
  symop[7][0][1] = -1.0;
  symop[7][1][2] =  1.0;
  symop[7][2][0] = -1.0;

  // C3^2 (d)
  symop[8][0][1] = -1.0;
  symop[8][1][2] = -1.0;
  symop[8][2][0] =  1.0;

  // C2(x)
  symop[9][0][0] =  1.0;
  symop[9][1][1] = -1.0;
  symop[9][2][2] = -1.0;

  // C2(y)
  symop[10][0][0] = -1.0;
  symop[10][1][1] =  1.0;
  symop[10][2][2] = -1.0;

  // C2(z)
  symop[11][0][0] = -1.0;
  symop[11][1][1] = -1.0;
  symop[11][2][2] =  1.0;
}

// this gives us the operations in Th which come from ixT (ie, the inverse
// operating on all the symmetry operations from T).
static void
th_ops(SymmetryOperation *symop)
{
  t_ops(symop);

  for (int i=0; i < 12; i++)
    for (int j=0; j < 3; j++)
      for (int k=0; k < 3; k++)
        symop[i][j][k] *= -1.0;
}

// this gives us the operations in Td which aren't in T.
static void
td_ops(SymmetryOperation *symop)
{
  // S4 (x)
  symop[0][0][0] = -1.0;
  symop[0][1][2] = -1.0;
  symop[0][2][1] =  1.0;

  // S4^3 (x)
  symop[1][0][0] = -1.0;
  symop[1][1][2] =  1.0;
  symop[1][2][1] = -1.0;

  // S4 (y)
  symop[2][0][2] =  1.0;
  symop[2][1][1] = -1.0;
  symop[2][2][0] = -1.0;

  // S4^3 (y)
  symop[3][0][2] = -1.0;
  symop[3][1][1] = -1.0;
  symop[3][2][0] =  1.0;

  // S4 (z)
  symop[4][0][1] = -1.0;
  symop[4][1][0] =  1.0;
  symop[4][2][2] = -1.0;

  // S4^3 (z)
  symop[5][0][1] =  1.0;
  symop[5][1][0] = -1.0;
  symop[5][2][2] = -1.0;

  // a = ( 1, 1, 1)
  // b = (-1,-1, 1)
  // c = ( 1,-1,-1)
  // d = (-1, 1,-1)
  // sigma (ac)
  symop[6][0][0] =  1.0;
  symop[6][1][2] =  1.0;
  symop[6][2][1] =  1.0;

  // sigma (bd)
  symop[7][0][0] =  1.0;
  symop[7][1][2] = -1.0;
  symop[7][2][1] = -1.0;

  // sigma (ad)
  symop[8][0][2] =  1.0;
  symop[8][1][1] =  1.0;
  symop[8][2][0] =  1.0;

  // sigma (bc)
  symop[9][0][2] = -1.0;
  symop[9][1][1] =  1.0;
  symop[9][2][0] = -1.0;

  // sigma (ab)
  symop[10][0][1] =  1.0;
  symop[10][1][0] =  1.0;
  symop[10][2][2] =  1.0;

  // sigma (dc)
  symop[11][0][1] = -1.0;
  symop[11][1][0] = -1.0;
  symop[11][2][2] =  1.0;
}

// this gives us the operations in O which aren't in T.  They are much
// like the Td operations, except they are inverted
static void
o_ops(SymmetryOperation *symop)
{
  td_ops(symop);

  for (int i=0; i < 12; i++)
    for (int j=0; j < 3; j++)
      for (int k=0; k < 3; k++)
        symop[i][j][k] *= -1.0;
}

void CharacterTable::t()
{
#if 0
  // t_ops gives us all the symmetry operations we need
  t_ops(symop);

  int i,j,k;

  {
    IrreducibleRepresentation ir(g,1,"A");
    for (i=0; i < g; i++) {
      ir.rep[i] = 1;
      ir.proj[0][i] = 1;
    }

    gamma_[0] = ir;
  }

  {
    IrreducibleRepresentation ir(g,2,"E");

    // identity
    ir.rep[0] = 1.0;
    ir.proj[0][0] = 1.0;
    ir.proj[1][0] = 0.0;
    ir.proj[2][0] = 0.0;
    ir.proj[3][0] = 1.0;

    // 4 C3's
    for (j=1; j < 5; j++) {
      ir.rep[j] = -0.5;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] = -sqrt(3.0)*0.5;
      ir.proj[2][j] =  sqrt(3.0)*0.5;
      ir.proj[3][j] = -0.5;
    }

    // 4 C3^2's
    for (j=5; j < 9; j++) {
      ir.rep[j] = -0.5;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] =  sqrt(3.0)*0.5;
      ir.proj[2][j] = -sqrt(3.0)*0.5;
      ir.proj[3][j] = -0.5;
    }

    // 3 C2's
    for (j=9; j < 12; j++) {
      ir.rep[j] = 1.0;
      ir.proj[0][j] = 1.0;
      ir.proj[1][j] = 0.0;
      ir.proj[2][j] = 0.0;
      ir.proj[3][j] = 1.0;
    }

    gamma_[1] = ir;
  }

  // the symmetry operation matrices give us a basis for irrep T
  {
    IrreducibleRepresentation ir(g,3,"T");
    
    ir.nrot_ = 1;
    ir.ntrans_ = 1;

    for (i=0; i < g; i++) {
      ir.rep[i]=0;
      
      for (j=0; j < 3; j++) {
        for (k=0; k < 3; k++)
          ir.proj[3*j+k][i] = symop[i][k][j];
        ir.rep[i] += ir.proj[3*j+j][i];
      }
    }

    gamma_[2] = ir;
  }

#endif
}


void CharacterTable::th()
{
#if 0
  // first get the ExT operations, then the ixT operations
  t_ops(symop);
  th_ops(&symop[12]);
  
  int i,j,k;
  {
    IrreducibleRepresentation ir1(g,1,"Ag");
    IrreducibleRepresentation ir2(g,1,"Au");
    for (i=0; i < 12; i++) {
      ir1.rep[i] = ir1.proj[0][i] = 1;
      ir2.rep[i] = ir2.proj[0][i] = 1;

      ir1.rep[i+12] = ir1.proj[0][i+12] = 1;
      ir2.rep[i+12] = ir2.proj[0][i+12] = -1;
    }

    gamma_[0] = ir1;
    gamma_[1] = ir2;
  }

  {
    IrreducibleRepresentation ir1(g,2,"Eg");
    IrreducibleRepresentation ir2(g,2,"Eu");

    // identity
    ir1.rep[0] = 1.0;
    ir1.proj[0][0] = 1.0;
    ir1.proj[1][0] = 0.0;
    ir1.proj[2][0] = 0.0;
    ir1.proj[3][0] = 1.0;

    // 4 C3's
    for (j=1; j < 5; j++) {
      ir1.rep[j] = -0.5;

      ir1.proj[0][j] = -0.5;
      ir1.proj[1][j] = -sqrt(3.0)*0.5;
      ir1.proj[2][j] =  sqrt(3.0)*0.5;
      ir1.proj[3][j] = -0.5;
    }

    // 4 C3^2's
    for (j=5; j < 9; j++) {
      ir1.rep[j] = -0.5;

      ir1.proj[0][j] = -0.5;
      ir1.proj[1][j] =  sqrt(3.0)*0.5;
      ir1.proj[2][j] = -sqrt(3.0)*0.5;
      ir1.proj[3][j] = -0.5;
    }

    // 3 C2's
    for (j=9; j < 12; j++) {
      ir1.rep[j] = 1.0;
      ir1.proj[0][j] = 1.0;
      ir1.proj[1][j] = 0.0;
      ir1.proj[2][j] = 0.0;
      ir1.proj[3][j] = 1.0;
    }

    for (j=0; j < 12; j++) {
      ir2.rep[j] = ir1.rep[j];
      ir1.rep[j+12] = ir1.rep[j];
      ir2.rep[j+12] = -ir1.rep[j];

      for (k=0; k < 4; k++) {
        ir2.proj[k][j] = ir1.proj[k][j];
        ir1.proj[k][j+12] = ir1.proj[k][j];
        ir2.proj[k][j+12] = -ir1.proj[k][j];
      }
    }

    gamma_[2] = ir1;
    gamma_[3] = ir2;
  }

  // the symmetry operation matrices form a basis for Tu.  Tg(g)=Tu(g) for
  // the proper rotations, and = -Tu(g) for the improper ones
  {
    IrreducibleRepresentation ir1(g,3,"Tg");
    IrreducibleRepresentation ir2(g,3,"Tu");

    ir1.nrot_ = 1;
    ir2.ntrans_ = 1;

    for (i=0; i < g; i++) {
      ir2.rep[i]=0;
      
      for (j=0; j < 3; j++) {
        for (k=0; k < 3; k++)
          ir2.proj[3*j+k][i] = symop[i][k][j];
        ir2.rep[i] += ir2.proj[3*j+j][i];
      }
    }

    for (i=0; i < g/2; i++) {
      ir1.rep[i] = ir2.rep[i];
      ir1.rep[i+12] = -ir2.rep[i+12];
      
      for (j=0; j < 9; j++) {
        ir1.proj[j][i] = ir2.proj[j][i];
        ir1.proj[j][i+12] = -ir2.proj[j][i+12];
      }
    }

    gamma_[4] = ir1;
    gamma_[5] = ir2;
  }
#endif
}

void CharacterTable::td()
{
#if 0
  // first get the T operations, then the Td operations
  t_ops(symop);
  td_ops(&symop[12]);
  
  int i,j,k;
  
  {
    IrreducibleRepresentation ir1(g,1,"A1");
    IrreducibleRepresentation ir2(g,1,"A2");
    for (i=0; i < 12; i++) {
      ir1.rep[i] = ir1.proj[0][i] = 1;
      ir2.rep[i] = ir2.proj[0][i] = 1;

      ir1.rep[i+12] = ir1.proj[0][i+12] = 1;
      ir2.rep[i+12] = ir2.proj[0][i+12] = -1;
    }

    gamma_[0] = ir1;
    gamma_[1] = ir2;
  }

  {
    IrreducibleRepresentation ir(g,2,"E");

    // identity
    ir.rep[0] = 2.0;

    ir.proj[0][0] = 1.0;
    ir.proj[1][0] = 0.0;
    ir.proj[2][0] = 0.0;
    ir.proj[3][0] = 1.0;

    // 4 C3's
    for (j=1; j < 5; j++) {
      ir.rep[j] = -1.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] = -0.5*sqrt(3.0);
      ir.proj[2][j] =  0.5*sqrt(3.0);
      ir.proj[3][j] = -0.5;
    }

    // 4 C3^2's
    for (j=5; j < 9; j++) {
      ir.rep[j] = -1.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] =  0.5*sqrt(3.0);
      ir.proj[2][j] = -0.5*sqrt(3.0);
      ir.proj[3][j] = -0.5;
    }

    // 3 C2's
    for (j=9; j < 12; j++) {
      ir.rep[j] = 2.0;

      ir.proj[0][j] = 1.0;
      ir.proj[1][j] = 0.0;
      ir.proj[2][j] = 0.0;
      ir.proj[3][j] = 1.0;
    }

    // 6 S4's
    for (j=12; j < 14; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] =  1.0;
      ir.proj[1][j] =  0.0;
      ir.proj[2][j] =  0.0;
      ir.proj[3][j] = -1.0;
    }

    for (j=14; j < 16; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] = -0.5*sqrt(3.0);
      ir.proj[2][j] = -0.5*sqrt(3.0);
      ir.proj[3][j] =  0.5;
    }

    for (j=16; j < 18; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] =  0.5*sqrt(3.0);
      ir.proj[2][j] =  0.5*sqrt(3.0);
      ir.proj[3][j] =  0.5;
    }

    // 6 sigmas
    for (j=18; j < 20; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] =  1.0;
      ir.proj[1][j] =  0.0;
      ir.proj[2][j] =  0.0;
      ir.proj[3][j] = -1.0;
    }

    for (j=20; j < 22; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] = -0.5*sqrt(3.0);
      ir.proj[2][j] = -0.5*sqrt(3.0);
      ir.proj[3][j] =  0.5;
    }

    for (j=22; j < 24; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] =  0.5*sqrt(3.0);
      ir.proj[2][j] =  0.5*sqrt(3.0);
      ir.proj[3][j] =  0.5;
    }


    gamma_[2] = ir;
  }

  // the symmetry operation matrices form a basis for T2.  T1(g)=T2(g) for
  // the proper rotations, and = -T2(g) for the improper ones
  {
    IrreducibleRepresentation ir1(g,3,"T1");
    IrreducibleRepresentation ir2(g,3,"T2");

    ir1.nrot_ = 1;
    ir2.ntrans_ = 1;

    for (i=0; i < g; i++) {
      ir2.rep[i]=0;
      
      for (j=0; j < 3; j++) {
        for (k=0; k < 3; k++)
          ir2.proj[3*j+k][i] = symop[i][k][j];
        ir2.rep[i] += ir2.proj[3*j+j][i];
      }
    }

    for (i=0; i < g/2; i++) {
      ir1.rep[i] = ir2.rep[i];
      ir1.rep[i+12] = -ir2.rep[i+12];
      
      for (j=0; j < 9; j++) {
        ir1.proj[j][i] = ir2.proj[j][i];
        ir1.proj[j][i+12] = -ir2.proj[j][i+12];
      }
    }

    gamma_[3] = ir1;
    gamma_[4] = ir2;
  }
#endif
}

void CharacterTable::o()
{
#if 0
  // first get the T operations, then the O operations
  t_ops(symop);
  o_ops(&symop[12]);
  
  int i,j,k;
  
  {
    IrreducibleRepresentation ir1(g,1,"A1");
    IrreducibleRepresentation ir2(g,1,"A2");
    for (i=0; i < 12; i++) {
      ir1.rep[i] = ir1.proj[0][i] = 1;
      ir2.rep[i] = ir2.proj[0][i] = 1;

      ir1.rep[i+12] = ir1.proj[0][i+12] = 1;
      ir2.rep[i+12] = ir2.proj[0][i+12] = -1;
    }

    gamma_[0] = ir1;
    gamma_[1] = ir2;
  }

  {
    IrreducibleRepresentation ir(g,2,"E");

    // identity
    ir.rep[0] = 2.0;

    ir.proj[0][0] = 1.0;
    ir.proj[1][0] = 0.0;
    ir.proj[2][0] = 0.0;
    ir.proj[3][0] = 1.0;

    // 4 C3's
    for (j=1; j < 5; j++) {
      ir.rep[j] = -1.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] = -0.5*sqrt(3.0);
      ir.proj[2][j] =  0.5*sqrt(3.0);
      ir.proj[3][j] = -0.5;
    }

    // 4 C3^2's
    for (j=5; j < 9; j++) {
      ir.rep[j] = -1.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] =  0.5*sqrt(3.0);
      ir.proj[2][j] = -0.5*sqrt(3.0);
      ir.proj[3][j] = -0.5;
    }

    // 3 C2's
    for (j=9; j < 12; j++) {
      ir.rep[j] = 2.0;

      ir.proj[0][j] = 1.0;
      ir.proj[1][j] = 0.0;
      ir.proj[2][j] = 0.0;
      ir.proj[3][j] = 1.0;
    }

    // 6 C4's
    for (j=12; j < 14; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] =  1.0;
      ir.proj[1][j] =  0.0;
      ir.proj[2][j] =  0.0;
      ir.proj[3][j] = -1.0;
    }

    for (j=14; j < 16; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] = -0.5*sqrt(3.0);
      ir.proj[2][j] = -0.5*sqrt(3.0);
      ir.proj[3][j] =  0.5;
    }

    for (j=16; j < 18; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] =  0.5*sqrt(3.0);
      ir.proj[2][j] =  0.5*sqrt(3.0);
      ir.proj[3][j] =  0.5;
    }

    // 6 C2's
    for (j=18; j < 20; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] =  1.0;
      ir.proj[1][j] =  0.0;
      ir.proj[2][j] =  0.0;
      ir.proj[3][j] = -1.0;
    }

    for (j=20; j < 22; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] = -0.5*sqrt(3.0);
      ir.proj[2][j] = -0.5*sqrt(3.0);
      ir.proj[3][j] =  0.5;
    }

    for (j=22; j < 24; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] =  0.5*sqrt(3.0);
      ir.proj[2][j] =  0.5*sqrt(3.0);
      ir.proj[3][j] =  0.5;
    }


    gamma_[2] = ir;
  }

  // the symmetry operation matrices form a basis for T1.  T2(g)=T1(g) for
  // the proper rotations, and = -T1(g) for the improper ones
  {
    IrreducibleRepresentation ir2(g,3,"T1");
    IrreducibleRepresentation ir1(g,3,"T2");

    ir2.nrot_ = 1;
    ir2.ntrans_ = 1;

    for (i=0; i < g; i++) {
      ir2.rep[i]=0;
      
      for (j=0; j < 3; j++) {
        for (k=0; k < 3; k++)
          ir2.proj[3*j+k][i] = symop[i][k][j];
        ir2.rep[i] += ir2.proj[3*j+j][i];
      }
    }

    for (i=0; i < g/2; i++) {
      ir1.rep[i] = ir2.rep[i];
      ir1.rep[i+12] = -ir2.rep[i+12];
      
      for (j=0; j < 9; j++) {
        ir1.proj[j][i] = ir2.proj[j][i];
        ir1.proj[j][i+12] = -ir2.proj[j][i+12];
      }
    }

    gamma_[3] = ir2;
    gamma_[4] = ir1;
  }
#endif
}

void CharacterTable::oh()
{
#if 0
  // first get the T operations, then the O operations, then the Th
  // operations, then the Td operations
  t_ops(symop);
  o_ops(&symop[12]);
  th_ops(&symop[24]);
  td_ops(&symop[36]);

  int i,j,k;
  
  {
    IrreducibleRepresentation ir1(g,1,"A1g");
    IrreducibleRepresentation ir2(g,1,"A2g");
    IrreducibleRepresentation ir3(g,1,"A1u");
    IrreducibleRepresentation ir4(g,1,"A2u");

    for (i=0; i < 12; i++) {
      ir1.rep[i] = ir1.proj[0][i] = 1;
      ir2.rep[i] = ir2.proj[0][i] = 1;
      ir3.rep[i] = ir3.proj[0][i] = 1;
      ir4.rep[i] = ir4.proj[0][i] = 1;

      ir1.rep[i+12] = ir1.proj[0][i+12] = 1;
      ir2.rep[i+12] = ir2.proj[0][i+12] = -1;
      ir3.rep[i+12] = ir3.proj[0][i+12] = 1;
      ir4.rep[i+12] = ir4.proj[0][i+12] = -1;

      ir1.rep[i+24] = ir1.proj[0][i+24] = 1;
      ir2.rep[i+24] = ir2.proj[0][i+24] = 1;
      ir3.rep[i+24] = ir3.proj[0][i+24] = -1;
      ir4.rep[i+24] = ir4.proj[0][i+24] = -1;

      ir1.rep[i+36] = ir1.proj[0][i+36] = 1;
      ir2.rep[i+36] = ir2.proj[0][i+36] = -1;
      ir3.rep[i+36] = ir3.proj[0][i+36] = -1;
      ir4.rep[i+36] = ir4.proj[0][i+36] = 1;
    }

    gamma_[0] = ir1;
    gamma_[1] = ir2;
    gamma_[5] = ir3;
    gamma_[6] = ir4;
  }

  {
    IrreducibleRepresentation ir(g,2,"Eg");
    IrreducibleRepresentation ir2(g,2,"Eu");

    // identity
    ir.rep[0] = 2.0;

    ir.proj[0][0] = 1.0;
    ir.proj[1][0] = 0.0;
    ir.proj[2][0] = 0.0;
    ir.proj[3][0] = 1.0;

    // 4 C3's
    for (j=1; j < 5; j++) {
      ir.rep[j] = -1.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] = -0.5*sqrt(3.0);
      ir.proj[2][j] =  0.5*sqrt(3.0);
      ir.proj[3][j] = -0.5;
    }

    // 4 C3^2's
    for (j=5; j < 9; j++) {
      ir.rep[j] = -1.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] =  0.5*sqrt(3.0);
      ir.proj[2][j] = -0.5*sqrt(3.0);
      ir.proj[3][j] = -0.5;
    }

    // 3 C2's
    for (j=9; j < 12; j++) {
      ir.rep[j] = 2.0;

      ir.proj[0][j] = 1.0;
      ir.proj[1][j] = 0.0;
      ir.proj[2][j] = 0.0;
      ir.proj[3][j] = 1.0;
    }

    // 6 C4's
    for (j=12; j < 14; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] =  1.0;
      ir.proj[1][j] =  0.0;
      ir.proj[2][j] =  0.0;
      ir.proj[3][j] = -1.0;
    }

    for (j=14; j < 16; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] = -0.5*sqrt(3.0);
      ir.proj[2][j] = -0.5*sqrt(3.0);
      ir.proj[3][j] =  0.5;
    }

    for (j=16; j < 18; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] =  0.5*sqrt(3.0);
      ir.proj[2][j] =  0.5*sqrt(3.0);
      ir.proj[3][j] =  0.5;
    }

    // 6 C2's
    for (j=18; j < 20; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] =  1.0;
      ir.proj[1][j] =  0.0;
      ir.proj[2][j] =  0.0;
      ir.proj[3][j] = -1.0;
    }

    for (j=20; j < 22; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] = -0.5*sqrt(3.0);
      ir.proj[2][j] = -0.5*sqrt(3.0);
      ir.proj[3][j] =  0.5;
    }

    for (j=22; j < 24; j++) {
      ir.rep[j] = 0.0;

      ir.proj[0][j] = -0.5;
      ir.proj[1][j] =  0.5*sqrt(3.0);
      ir.proj[2][j] =  0.5*sqrt(3.0);
      ir.proj[3][j] =  0.5;
    }


    for (j=0; j < 24; j++) {
      ir.rep[j+24] = ir.rep[j];
      ir2.rep[j] = ir.rep[j];
      ir2.rep[j+24] = -ir.rep[j];

      for (k=0; k < 4; k++) {
        ir2.proj[k][j] = ir.proj[k][j];
        ir.proj[k][j+24] = ir.proj[k][j];
        ir2.proj[k][j+24] = -ir.proj[k][j];
      }
    }
    
    gamma_[2] = ir;
    gamma_[7] = ir2;
  }

  // the symmetry operation matrices form a basis for T1u.  T2u(g)=T1u(g) for
  // the proper rotations, and = -T1(g) for the improper ones.
  // T1g(g)=T1u(g) for the O part, and = -T1u(g) for the ixO part.
  // T2g(g)=T1g(g) for proper rotations and =-T1g(g) for improper
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

    for (i=0; i < g/4; i++) {
      ir4.rep[i] = ir3.rep[i];
      ir4.rep[i+12] = -ir3.rep[i+12];
      ir4.rep[i+24] =  ir3.rep[i+24];
      ir4.rep[i+36] = -ir3.rep[i+36];

      ir1.rep[i] = ir3.rep[i];
      ir1.rep[i+12] =  ir3.rep[i+12];
      ir1.rep[i+24] = -ir3.rep[i+24];
      ir1.rep[i+36] = -ir3.rep[i+36];
      
      ir2.rep[i] = ir3.rep[i];
      ir2.rep[i+12] = -ir3.rep[i+12];
      ir2.rep[i+24] = -ir3.rep[i+24];
      ir2.rep[i+36] =  ir3.rep[i+36];
      
      for (j=0; j < 9; j++) {
        ir4.proj[j][i] = ir3.proj[j][i];
        ir4.proj[j][i+12] = -ir3.proj[j][i+12];
        ir4.proj[j][i+24] =  ir3.proj[j][i+24];
        ir4.proj[j][i+36] = -ir3.proj[j][i+36];

        ir1.proj[j][i] = ir3.proj[j][i];
        ir1.proj[j][i+12] =  ir3.proj[j][i+12];
        ir1.proj[j][i+24] = -ir3.proj[j][i+24];
        ir1.proj[j][i+36] = -ir3.proj[j][i+36];

        ir2.proj[j][i] = ir3.proj[j][i];
        ir2.proj[j][i+12] = -ir3.proj[j][i+12];
        ir2.proj[j][i+24] = -ir3.proj[j][i+24];
        ir2.proj[j][i+36] =  ir3.proj[j][i+36];
      }
    }

    gamma_[3] = ir1;
    gamma_[4] = ir2;
    gamma_[8] = ir3;
    gamma_[9] = ir4;
  }
#endif
}
