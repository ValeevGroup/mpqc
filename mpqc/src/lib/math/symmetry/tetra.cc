
#include <math.h>
#include <string.h>

#include <math/symmetry/pointgrp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


void CharacterTable::t()
{
  int i,j;

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

  {
    IrreducibleRepresentation ir(g,3,"T");
    
    ir.nrot_ = 1;
    ir.ntrans_ = 1;

    // identity
    ir.rep[0] = 3.0;

    // 8 C3's
    for (j=1; j < 9; j++)
      ir.rep[j] = 0.0;

    // 3 C2's
    for (j=9; j < 12; j++)
      ir.rep[j] = -1.0;

    // identity
    ir.proj[0][0] = 1.0;
    ir.proj[1][0] = 0.0;
    ir.proj[2][0] = 0.0;

    ir.proj[3][0] = 0.0;
    ir.proj[4][0] = 1.0;
    ir.proj[5][0] = 0.0;

    ir.proj[6][0] = 0.0;
    ir.proj[7][0] = 0.0;
    ir.proj[8][0] = 1.0;

    // 8 C3's
    ir.proj[0][1] =  0.0;
    ir.proj[1][1] =  1.0;
    ir.proj[2][1] =  0.0;
    ir.proj[3][1] =  0.0;
    ir.proj[4][1] =  0.0;
    ir.proj[5][1] =  1.0;
    ir.proj[6][1] =  1.0;
    ir.proj[7][1] =  0.0;
    ir.proj[8][1] =  0.0;

    ir.proj[0][2] =  0.0;
    ir.proj[1][2] =  1.0;
    ir.proj[2][2] =  0.0;
    ir.proj[3][2] =  0.0;
    ir.proj[4][2] =  0.0;
    ir.proj[5][2] = -1.0;
    ir.proj[6][2] = -1.0;
    ir.proj[7][2] =  0.0;
    ir.proj[8][2] =  0.0;

    ir.proj[0][3] =  0.0;
    ir.proj[1][3] = -1.0;
    ir.proj[2][3] =  0.0;
    ir.proj[3][3] =  0.0;
    ir.proj[4][3] =  0.0;
    ir.proj[5][3] =  1.0;
    ir.proj[6][3] = -1.0;
    ir.proj[7][3] =  0.0;
    ir.proj[8][3] =  0.0;

    ir.proj[0][4] =  0.0;
    ir.proj[1][4] = -1.0;
    ir.proj[2][4] =  0.0;
    ir.proj[3][4] =  0.0;
    ir.proj[4][4] =  0.0;
    ir.proj[5][4] = -1.0;
    ir.proj[6][4] =  1.0;
    ir.proj[7][4] =  0.0;
    ir.proj[8][4] =  0.0;

    ir.proj[0][5] =  0.0;
    ir.proj[1][5] =  0.0;
    ir.proj[2][5] =  1.0;
    ir.proj[3][5] =  1.0;
    ir.proj[4][5] =  0.0;
    ir.proj[5][5] =  0.0;
    ir.proj[6][5] =  0.0;
    ir.proj[7][5] =  1.0;
    ir.proj[8][5] =  0.0;

    ir.proj[0][6] =  0.0;
    ir.proj[1][6] =  0.0;
    ir.proj[2][6] = -1.0;
    ir.proj[3][6] =  1.0;
    ir.proj[4][6] =  0.0;
    ir.proj[5][6] =  0.0;
    ir.proj[6][6] =  0.0;
    ir.proj[7][6] = -1.0;
    ir.proj[8][6] =  0.0;

    ir.proj[0][7] =  0.0;
    ir.proj[1][7] =  0.0;
    ir.proj[2][7] = -1.0;
    ir.proj[3][7] = -1.0;
    ir.proj[4][7] =  0.0;
    ir.proj[5][7] =  0.0;
    ir.proj[6][7] =  0.0;
    ir.proj[7][7] =  1.0;
    ir.proj[8][7] =  0.0;

    ir.proj[0][8] =  0.0;
    ir.proj[1][8] =  0.0;
    ir.proj[2][8] =  1.0;
    ir.proj[3][8] = -1.0;
    ir.proj[4][8] =  0.0;
    ir.proj[5][8] =  0.0;
    ir.proj[6][8] =  0.0;
    ir.proj[7][8] = -1.0;
    ir.proj[8][8] =  0.0;

    // 3 C2's
    ir.proj[0][9] =  1.0;
    ir.proj[1][9] =  0.0;
    ir.proj[2][9] =  0.0;
    ir.proj[3][9] =  0.0;
    ir.proj[4][9] = -1.0;
    ir.proj[5][9] =  0.0;
    ir.proj[6][9] =  0.0;
    ir.proj[7][9] =  0.0;
    ir.proj[8][9] = -1.0;

    ir.proj[0][10] = -1.0;
    ir.proj[1][10] =  0.0;
    ir.proj[2][10] =  0.0;
    ir.proj[3][10] =  0.0;
    ir.proj[4][10] =  1.0;
    ir.proj[5][10] =  0.0;
    ir.proj[6][10] =  0.0;
    ir.proj[7][10] =  0.0;
    ir.proj[8][10] = -1.0;

    ir.proj[0][11] = -1.0;
    ir.proj[1][11] =  0.0;
    ir.proj[2][11] =  0.0;
    ir.proj[3][11] =  0.0;
    ir.proj[4][11] = -1.0;
    ir.proj[5][11] =  0.0;
    ir.proj[6][11] =  0.0;
    ir.proj[7][11] =  0.0;
    ir.proj[8][11] =  1.0;

    gamma_[2] = ir;
  }

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


void CharacterTable::th()
{
  fprintf(stderr,"no Th yet\n");
  abort();
}

void CharacterTable::td()
{
  int i,j;
  
  {
    IrreducibleRepresentation ir(g,1,"A1");
    for (i=0; i < g; i++) {
      ir.rep[i] = 1;
      ir.proj[0][i] = 1;
    }

    gamma_[0] = ir;
  }

  {
    IrreducibleRepresentation ir(g,1,"A2");
    for (i=0; i < 12; i++) {
      ir.rep[i] = 1;
      ir.proj[0][i] = 1;
    }
    for (i=12; i < 24; i++) {
      ir.rep[i] = -1;
      ir.proj[0][i] = -1;
    }

    gamma_[1] = ir;
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

  {
    IrreducibleRepresentation ir1(g,3,"T1");
    IrreducibleRepresentation ir2(g,3,"T2");

    ir1.nrot_ = 1;
    ir2.ntrans_ = 1;

    // identity
    ir1.rep[0] = 3.0;
    ir2.rep[0] = 3.0;

    // 8 C3's
    for (j=1; j < 9; j++) {
      ir1.rep[j] = 0.0;
      ir2.rep[j] = 0.0;
    }

    // 3 C2's
    for (j=9; j < 12; j++) {
      ir1.rep[j] = -1.0;
      ir2.rep[j] = -1.0;
    }

    // 6 S4's
    for (j=12; j < 18; j++) {
      ir1.rep[j] =  1.0;
      ir2.rep[j] = -1.0;
    }

    // 6 sigma's
    for (j=18; j < 24; j++) {
      ir1.rep[j] = -1.0;
      ir2.rep[j] =  1.0;
    }

    // identity
    ir1.proj[0][0] = 1.0;
    ir1.proj[1][0] = 0.0;
    ir1.proj[2][0] = 0.0;

    ir1.proj[3][0] = 0.0;
    ir1.proj[4][0] = 1.0;
    ir1.proj[5][0] = 0.0;

    ir1.proj[6][0] = 0.0;
    ir1.proj[7][0] = 0.0;
    ir1.proj[8][0] = 1.0;

    // 8 C3's
    ir1.proj[0][1] =  0.0;
    ir1.proj[1][1] =  1.0;
    ir1.proj[2][1] =  0.0;
    ir1.proj[3][1] =  0.0;
    ir1.proj[4][1] =  0.0;
    ir1.proj[5][1] =  1.0;
    ir1.proj[6][1] =  1.0;
    ir1.proj[7][1] =  0.0;
    ir1.proj[8][1] =  0.0;

    ir1.proj[0][2] =  0.0;
    ir1.proj[1][2] =  1.0;
    ir1.proj[2][2] =  0.0;
    ir1.proj[3][2] =  0.0;
    ir1.proj[4][2] =  0.0;
    ir1.proj[5][2] = -1.0;
    ir1.proj[6][2] = -1.0;
    ir1.proj[7][2] =  0.0;
    ir1.proj[8][2] =  0.0;

    ir1.proj[0][3] =  0.0;
    ir1.proj[1][3] = -1.0;
    ir1.proj[2][3] =  0.0;
    ir1.proj[3][3] =  0.0;
    ir1.proj[4][3] =  0.0;
    ir1.proj[5][3] =  1.0;
    ir1.proj[6][3] = -1.0;
    ir1.proj[7][3] =  0.0;
    ir1.proj[8][3] =  0.0;

    ir1.proj[0][4] =  0.0;
    ir1.proj[1][4] = -1.0;
    ir1.proj[2][4] =  0.0;
    ir1.proj[3][4] =  0.0;
    ir1.proj[4][4] =  0.0;
    ir1.proj[5][4] = -1.0;
    ir1.proj[6][4] =  1.0;
    ir1.proj[7][4] =  0.0;
    ir1.proj[8][4] =  0.0;

    ir1.proj[0][5] =  0.0;
    ir1.proj[1][5] =  0.0;
    ir1.proj[2][5] =  1.0;
    ir1.proj[3][5] =  1.0;
    ir1.proj[4][5] =  0.0;
    ir1.proj[5][5] =  0.0;
    ir1.proj[6][5] =  0.0;
    ir1.proj[7][5] =  1.0;
    ir1.proj[8][5] =  0.0;

    ir1.proj[0][6] =  0.0;
    ir1.proj[1][6] =  0.0;
    ir1.proj[2][6] = -1.0;
    ir1.proj[3][6] =  1.0;
    ir1.proj[4][6] =  0.0;
    ir1.proj[5][6] =  0.0;
    ir1.proj[6][6] =  0.0;
    ir1.proj[7][6] = -1.0;
    ir1.proj[8][6] =  0.0;

    ir1.proj[0][7] =  0.0;
    ir1.proj[1][7] =  0.0;
    ir1.proj[2][7] = -1.0;
    ir1.proj[3][7] = -1.0;
    ir1.proj[4][7] =  0.0;
    ir1.proj[5][7] =  0.0;
    ir1.proj[6][7] =  0.0;
    ir1.proj[7][7] =  1.0;
    ir1.proj[8][7] =  0.0;

    ir1.proj[0][8] =  0.0;
    ir1.proj[1][8] =  0.0;
    ir1.proj[2][8] =  1.0;
    ir1.proj[3][8] = -1.0;
    ir1.proj[4][8] =  0.0;
    ir1.proj[5][8] =  0.0;
    ir1.proj[6][8] =  0.0;
    ir1.proj[7][8] = -1.0;
    ir1.proj[8][8] =  0.0;

    // 3 C2's
    ir1.proj[0][9] =  1.0;
    ir1.proj[1][9] =  0.0;
    ir1.proj[2][9] =  0.0;
    ir1.proj[3][9] =  0.0;
    ir1.proj[4][9] = -1.0;
    ir1.proj[5][9] =  0.0;
    ir1.proj[6][9] =  0.0;
    ir1.proj[7][9] =  0.0;
    ir1.proj[8][9] = -1.0;

    ir1.proj[0][10] = -1.0;
    ir1.proj[1][10] =  0.0;
    ir1.proj[2][10] =  0.0;
    ir1.proj[3][10] =  0.0;
    ir1.proj[4][10] =  1.0;
    ir1.proj[5][10] =  0.0;
    ir1.proj[6][10] =  0.0;
    ir1.proj[7][10] =  0.0;
    ir1.proj[8][10] = -1.0;

    ir1.proj[0][11] = -1.0;
    ir1.proj[1][11] =  0.0;
    ir1.proj[2][11] =  0.0;
    ir1.proj[3][11] =  0.0;
    ir1.proj[4][11] = -1.0;
    ir1.proj[5][11] =  0.0;
    ir1.proj[6][11] =  0.0;
    ir1.proj[7][11] =  0.0;
    ir1.proj[8][11] =  1.0;

    // 6 S4's
    ir1.proj[0][12] =  1.0;
    ir1.proj[1][12] =  0.0;
    ir1.proj[2][12] =  0.0;
    ir1.proj[3][12] =  0.0;
    ir1.proj[4][12] =  0.0;
    ir1.proj[5][12] = -1.0;
    ir1.proj[6][12] =  0.0;
    ir1.proj[7][12] =  1.0;
    ir1.proj[8][12] =  0.0;

    ir1.proj[0][13] =  1.0;
    ir1.proj[1][13] =  0.0;
    ir1.proj[2][13] =  0.0;
    ir1.proj[3][13] =  0.0;
    ir1.proj[4][13] =  0.0;
    ir1.proj[5][13] =  1.0;
    ir1.proj[6][13] =  0.0;
    ir1.proj[7][13] = -1.0;
    ir1.proj[8][13] =  0.0;

    ir1.proj[0][14] =  0.0;
    ir1.proj[1][14] =  0.0;
    ir1.proj[2][14] =  1.0;
    ir1.proj[3][14] =  0.0;
    ir1.proj[4][14] =  1.0;
    ir1.proj[5][14] =  0.0;
    ir1.proj[6][14] = -1.0;
    ir1.proj[7][14] =  0.0;
    ir1.proj[8][14] =  0.0;

    ir1.proj[0][15] =  0.0;
    ir1.proj[1][15] =  0.0;
    ir1.proj[2][15] = -1.0;
    ir1.proj[3][15] =  0.0;
    ir1.proj[4][15] =  1.0;
    ir1.proj[5][15] =  0.0;
    ir1.proj[6][15] =  1.0;
    ir1.proj[7][15] =  0.0;
    ir1.proj[8][15] =  0.0;

    ir1.proj[0][16] =  0.0;
    ir1.proj[1][16] = -1.0;
    ir1.proj[2][16] =  0.0;
    ir1.proj[3][16] =  1.0;
    ir1.proj[4][16] =  0.0;
    ir1.proj[5][16] =  0.0;
    ir1.proj[6][16] =  0.0;
    ir1.proj[7][16] =  0.0;
    ir1.proj[8][16] =  1.0;

    ir1.proj[0][17] =  0.0;
    ir1.proj[1][17] =  1.0;
    ir1.proj[2][17] =  0.0;
    ir1.proj[3][17] = -1.0;
    ir1.proj[4][17] =  0.0;
    ir1.proj[5][17] =  0.0;
    ir1.proj[6][17] =  0.0;
    ir1.proj[7][17] =  0.0;
    ir1.proj[8][17] =  1.0;

    // 6 sigma's
    ir1.proj[0][18] = -1.0;
    ir1.proj[1][18] =  0.0;
    ir1.proj[2][18] =  0.0;
    ir1.proj[3][18] =  0.0;
    ir1.proj[4][18] =  0.0;
    ir1.proj[5][18] = -1.0;
    ir1.proj[6][18] =  0.0;
    ir1.proj[7][18] = -1.0;
    ir1.proj[8][18] =  0.0;

    ir1.proj[0][19] = -1.0;
    ir1.proj[1][19] =  0.0;
    ir1.proj[2][19] =  0.0;
    ir1.proj[3][19] =  0.0;
    ir1.proj[4][19] =  0.0;
    ir1.proj[5][19] =  1.0;
    ir1.proj[6][19] =  0.0;
    ir1.proj[7][19] =  1.0;
    ir1.proj[8][19] =  0.0;

    ir1.proj[0][20] =  0.0;
    ir1.proj[1][20] =  0.0;
    ir1.proj[2][20] = -1.0;
    ir1.proj[3][20] =  0.0;
    ir1.proj[4][20] = -1.0;
    ir1.proj[5][20] =  0.0;
    ir1.proj[6][20] = -1.0;
    ir1.proj[7][20] =  0.0;
    ir1.proj[8][20] =  0.0;

    ir1.proj[0][21] =  0.0;
    ir1.proj[1][21] =  0.0;
    ir1.proj[2][21] =  1.0;
    ir1.proj[3][21] =  0.0;
    ir1.proj[4][21] = -1.0;
    ir1.proj[5][21] =  0.0;
    ir1.proj[6][21] =  1.0;
    ir1.proj[7][21] =  0.0;
    ir1.proj[8][21] =  0.0;

    ir1.proj[0][22] =  0.0;
    ir1.proj[1][22] = -1.0;
    ir1.proj[2][22] =  0.0;
    ir1.proj[3][22] = -1.0;
    ir1.proj[4][22] =  0.0;
    ir1.proj[5][22] =  0.0;
    ir1.proj[6][22] =  0.0;
    ir1.proj[7][22] =  0.0;
    ir1.proj[8][22] = -1.0;

    ir1.proj[0][23] =  0.0;
    ir1.proj[1][23] =  1.0;
    ir1.proj[2][23] =  0.0;
    ir1.proj[3][23] =  1.0;
    ir1.proj[4][23] =  0.0;
    ir1.proj[5][23] =  0.0;
    ir1.proj[6][23] =  0.0;
    ir1.proj[7][23] =  0.0;
    ir1.proj[8][23] = -1.0;

    for (int k=0; k < 9; k++) {
      for (j=0; j < 12; j++)
        ir2.proj[k][j] = ir1.proj[k][j];

      for (j=12; j < 24; j++)
        ir2.proj[k][j] = -ir1.proj[k][j];
    }

    gamma_[3] = ir1;
    gamma_[4] = ir2;
  }

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

  // S4 (x)
  symop[12][0][0] = -1.0;
  symop[12][1][2] = -1.0;
  symop[12][2][1] =  1.0;

  // S4^3 (x)
  symop[13][0][0] = -1.0;
  symop[13][1][2] =  1.0;
  symop[13][2][1] = -1.0;

  // S4 (y)
  symop[14][0][2] =  1.0;
  symop[14][1][1] = -1.0;
  symop[14][2][0] = -1.0;

  // S4^3 (y)
  symop[15][0][2] = -1.0;
  symop[15][1][1] = -1.0;
  symop[15][2][0] =  1.0;

  // S4 (z)
  symop[16][0][1] = -1.0;
  symop[16][1][0] =  1.0;
  symop[16][2][2] = -1.0;

  // S4^3 (z)
  symop[17][0][1] =  1.0;
  symop[17][1][0] = -1.0;
  symop[17][2][2] = -1.0;

  // a = ( 1, 1, 1)
  // b = (-1,-1, 1)
  // c = ( 1,-1,-1)
  // d = (-1, 1,-1)
  // sigma (ac)
  symop[18][0][0] =  1.0;
  symop[18][1][2] =  1.0;
  symop[18][2][1] =  1.0;

  // sigma (bd)
  symop[19][0][0] =  1.0;
  symop[19][1][2] = -1.0;
  symop[19][2][1] = -1.0;

  // sigma (ad)
  symop[20][0][2] =  1.0;
  symop[20][1][1] =  1.0;
  symop[20][2][0] =  1.0;

  // sigma (bc)
  symop[21][0][2] = -1.0;
  symop[21][1][1] =  1.0;
  symop[21][2][0] = -1.0;

  // sigma (ab)
  symop[22][0][1] =  1.0;
  symop[22][1][0] =  1.0;
  symop[22][2][2] =  1.0;

  // sigma (dc)
  symop[23][0][1] = -1.0;
  symop[23][1][0] = -1.0;
  symop[23][2][2] =  1.0;
}
