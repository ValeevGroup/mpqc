
#include <math.h>
#include <string.h>

#include <math/symmetry/pointgrp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// these are the operations which make up T
static void
i_ops(SymRep *t1rep, SymRep *t2rep, SymRep *grep, SymRep *hrep)
{
#if 0
  int i;
  
  // identity
  for (i=0; i < 3; i++) {
    t1rep[0][i][i] = 1.0;
    t2rep[0][i][i] = 1.0;
    grep[0][i][i] = 1.0;
    hrep[0][i][i] = 1.0;
  }
  grep[0][3][3] = 1.0;
  hrep[0][3][3] = 1.0;
  hrep[0][4][4] = 1.0;
    
  //
  // 12 C5's
  //
  // first the 2 C5's about the z axis
  double c2p5 = cos(2.0*M_PI/5.0);
  double s2p5 = sin(2.0*M_PI/5.0);
  double c4p5 = cos(4.0*M_PI/5.0);
  double s4p5 = sin(4.0*M_PI/5.0);

  t1rep[1][0][0] =  c2p5;
  t1rep[1][0][1] =  s2p5;
  t1rep[1][1][0] = -s2p5;
  t1rep[1][1][1] =  c2p5;
  t1rep[1][2][2] =  1.0;
  
  t1rep[2][0][0] =  c2p5;
  t1rep[2][0][1] = -s2p5;
  t1rep[2][1][0] =  s2p5;
  t1rep[2][1][1] =  c2p5;
  t1rep[2][2][2] =  1.0;
  
  t2rep[1] = t1rep[1].operate(t1rep[1]);
  t2rep[2] = t1rep[2].operate(t1rep[2]);

  grep[1][0][0] =  c2p5;
  grep[1][0][1] =  s2p5;
  grep[1][1][0] = -s2p5;
  grep[1][1][1] =  c2p5;
  grep[1][2][2] =  c4p5;
  grep[1][2][3] = -s4p5;
  grep[1][3][2] =  s4p5;
  grep[1][3][3] =  c4p5;
  
  grep[2][0][0] =  c2p5;
  grep[2][0][1] = -s2p5;
  grep[2][1][0] =  s2p5;
  grep[2][1][1] =  c2p5;
  grep[2][2][2] =  c4p5;
  grep[2][2][3] =  s4p5;
  grep[2][3][2] = -s4p5;
  grep[2][3][3] =  c4p5;
  
  hrep[1][0][0] = 1.0;
  hrep[1][1][1] =  c4p5;
  hrep[1][1][2] =  s4p5;
  hrep[1][2][1] = -s4p5;
  hrep[1][2][2] =  c4p5;
  hrep[1][3][3] =  c2p5;
  hrep[1][3][4] = -s2p5;
  hrep[1][4][3] =  s2p5;
  hrep[1][4][4] =  c2p5;
  
  hrep[2][0][0] = 1.0;
  hrep[2][1][1] =  c4p5;
  hrep[2][1][2] = -s4p5;
  hrep[2][2][1] =  s4p5;
  hrep[2][2][2] =  c4p5;
  hrep[2][3][3] =  c2p5;
  hrep[2][3][4] =  s2p5;
  hrep[2][4][3] = -s2p5;
  hrep[2][4][4] =  c2p5;
  
  // form rotation matrices for the C3 axis about the zx axis (these were
  // taken from turbomole version 2, which claims they were sort of inherited
  // from hondo
  Rep t1so(3);
  Rep gso(4);
  Rep hso(5);

  double cosd = s2p5/((1.0-c2p5)*sqrt(3.0));
  double cosd2 = cosd*cosd;
  double sind2 = 1.0 - cosd2;
  double sind = sqrt(sind2);

  t1so[0][0] =  1.0 - 1.5*cosd2;
  t1so[1][0] =  0.5*sqrt(3.0)*cosd;
  t1so[2][0] =  1.5*cosd*sind;
  t1so[0][1] = -0.5*sqrt(3.0)*cosd;
  t1so[1][1] = -0.5;
  t1so[2][1] =  0.5*sqrt(3.0)*sind;
  t1so[0][2] =  1.5*cosd*sind;
  t1so[1][2] = -0.5*sqrt(3.0)*sind;
  t1so[2][2] =  1.0 - 1.5*sind2;

  gso[0][0] = (3.0*sqrt(5.0)+5.0)/20.0;
  gso[0][1] = cosd*sqrt(3.0)*(sqrt(5.0)-1.0)/4.0;
  gso[0][2] = 3.0*sqrt(5.0)/10.0;
  gso[0][3] = -sqrt(5.0-2.0*sqrt(5.0))*sqrt(5.0)/10.0;
  gso[1][0] = -gso[0][1];
  gso[1][1] = (1-sqrt(5))/4.0;
  gso[1][2] = cosd*sqrt(3)/2.0;
  gso[1][3] = cosd*sqrt(5-2*sqrt(5))*sqrt(3)/2.0;
  gso[2][0] = gso[0][2];
  gso[2][1] = -gso[1][2];
  gso[2][2] = (5-3*sqrt(5))/20.0;
  gso[2][3] = sqrt(5-2*sqrt(5))*(sqrt(5)+5)/20;
  gso[3][0] = -gso[0][3];
  gso[3][1] = gso[1][3];
  gso[3][2] = -gso[2][3];
  gso[3][3] = (sqrt(5)+1)/4.0;

  hso[0][0] = -1.0/5.0;
  hso[0][4] = sqrt(3)*(sqrt(5)+1)/10.0;
  hso[0][3] = 3.0*cosd*(3.0*sqrt(5.0)-5.0)/10.0;
  hso[0][2] = 3.0*cosd*(5.0-sqrt(5.0))/10.0;
  hso[0][1] = sqrt(3.0)*(sqrt(5.0)-1.0)/10.0;
  hso[4][0] = hso[0][4];
  hso[4][4] = (2.0*sqrt(5.0)+1.0)/10.0;
  hso[4][3] = sqrt(3.0)*cosd*(5.0-2.0*sqrt(5.0))/10.0;
  hso[4][2] = sqrt(3.0)*cosd*(5.0-3.0*sqrt(5.0))/5.0;
  hso[4][1] = 2.0/5.0;
  hso[3][0] = -hso[0][3];
  hso[3][4] = -hso[4][3];
  hso[3][3] = -1.0/2.0;
  hso[3][2] = 0.0;
  hso[3][1] = sqrt(3.0)*cosd*(5.0-sqrt(5.0))/5.0;
  hso[2][0] = -hso[0][2];
  hso[2][4] = -hso[4][2];
  hso[2][3] = 0.0;
  hso[2][2] = -1.0/2.0;
  hso[2][1] = -sqrt(3.0)*sqrt(5.0)*cosd/10.0;
  hso[1][0] = hso[0][1];
  hso[1][4] = hso[4][1];
  hso[1][3] = -hso[3][1];
  hso[1][2] = -hso[2][1];
  hso[1][1] = (1.0-2.0*sqrt(5.0))/10.0;
  
  // now rotate the first C5's by 2pi/3 degrees about the zx axis (sort of)
  t1rep[3] = t1rep[1].sim_transform(t1so);
  t1rep[4] = t1rep[2].sim_transform(t1so);

  grep[3] = grep[1].sim_transform(gso);
  grep[4] = grep[2].sim_transform(gso);

  hrep[3] = hrep[1].sim_transform(hso);
  hrep[4] = hrep[2].sim_transform(hso);

  // rotate twice to get the first one aligned along the x axis
  t1rep[3] = t1rep[3].sim_transform(t1rep[1]).sim_transform(t1rep[1]);
  t1rep[4] = t1rep[4].sim_transform(t1rep[1]).sim_transform(t1rep[1]);

  grep[3] = grep[3].sim_transform(grep[1]).sim_transform(grep[1]);
  grep[4] = grep[4].sim_transform(grep[1]).sim_transform(grep[1]);

  hrep[3] = hrep[3].sim_transform(hrep[1]).sim_transform(hrep[1]);
  hrep[4] = hrep[4].sim_transform(hrep[1]).sim_transform(hrep[1]);

  t2rep[3] = t1rep[4].operate(t1rep[4]);
  t2rep[4] = t1rep[3].operate(t1rep[3]);

  t2rep[13] = t1rep[2];
  t2rep[14] = t1rep[1];

  t2rep[15] = t1rep[3];
  t2rep[16] = t1rep[4];
  
  // and then rotate those by 2pi/5 about the z axis 4 times
  for (i=5; i < 13; i++) {
    t1rep[i] = t1rep[i-2].sim_transform(t1rep[1]);
    grep[i] = grep[i-2].sim_transform(grep[1]);
    hrep[i] = hrep[i-2].sim_transform(hrep[1]);

    t2rep[i] = t2rep[i-2].sim_transform(t2rep[1]);
    t2rep[i+12] = t2rep[i+10].sim_transform(t2rep[1]);
  }

  //
  // 12 C5^2's
  //
  // get these from operating on each of the C5's with itself
  for (i=13; i < 25; i++) {
    t1rep[i] = t1rep[i-12].operate(t1rep[i-12]);
    grep[i] = grep[i-12].operate(grep[i-12]);
    hrep[i] = hrep[i-12].operate(hrep[i-12]);
  }

  //
  // 20 C3's
  //
  // first we have 2 C3's about the zx axis
  t1rep[25] = t1so;
  t1rep[26] = t1so.operate(t1so);
  
  grep[25] = gso;
  grep[26] = gso.operate(gso);
  
  hrep[25] = hso;
  hrep[26] = hso.operate(hso);
  
  // and then rotate those by 2pi/5 about the z axis 4 times
  for (i=27; i < 35; i++) {
    t1rep[i] = t1rep[i-2].sim_transform(t1rep[1]);
    grep[i] = grep[i-2].sim_transform(grep[1]);
    hrep[i] = hrep[i-2].sim_transform(hrep[1]);
  }

  // now rotate one of the above C3's by 2pi/3 about the zx axis
  t1rep[35] = t1rep[27].sim_transform(t1so);
  t1rep[36] = t1rep[28].sim_transform(t1so);

  grep[35] = grep[27].sim_transform(gso);
  grep[36] = grep[28].sim_transform(gso);

  hrep[35] = hrep[27].sim_transform(hso);
  hrep[36] = hrep[28].sim_transform(hso);

  // and then rotate those by 2pi/5 about the z axis 4 times
  for (i=37; i < 45; i++) {
    t1rep[i] = t1rep[i-2].sim_transform(t1rep[1]);
    grep[i] = grep[i-2].sim_transform(grep[1]);
    hrep[i] = hrep[i-2].sim_transform(hrep[1]);
  }

  t2rep[25] = t1rep[35];
  t2rep[26] = t1rep[36];
  
  for (i=27; i < 35; i++)
    t2rep[i] = t2rep[i-2].sim_transform(t2rep[1]);
  
  t2rep[35] = t1rep[26];
  t2rep[36] = t1rep[25];
  
  for (i=37; i < 45; i++)
    t2rep[i] = t2rep[i-2].sim_transform(t2rep[1]);

  //
  // 15 C2's
  //
  // first we have a C2 about the y axis
  t1rep[45][0][0] = -1.0;
  t1rep[45][1][1] =  1.0;
  t1rep[45][2][2] = -1.0;

  t2rep[45] = t1rep[45];
  
  grep[45][0][0] = -1.0;
  grep[45][1][1] =  1.0;
  grep[45][2][2] = -1.0;
  grep[45][3][3] =  1.0;
  
  hrep[45][0][0] =  1.0;
  hrep[45][1][1] =  1.0;
  hrep[45][2][2] = -1.0;
  hrep[45][3][3] = -1.0;
  hrep[45][4][4] =  1.0;
  
  // and rotate that by 2pi/5 about the z axis 4 times
  for (i=46; i < 50; i++) {
    t1rep[i] = t1rep[i-1].sim_transform(t1rep[1]);
    t2rep[i] = t2rep[i-1].sim_transform(t2rep[1]);
    grep[i] = grep[i-1].sim_transform(grep[1]);
    hrep[i] = hrep[i-1].sim_transform(hrep[1]);
  }

  // now take the C2 about the y axis and rotate it by 2pi/3 about the zx axis
  t1rep[50] = t1rep[45].sim_transform(t1so);
  grep[50] = grep[45].sim_transform(gso);
  hrep[50] = hrep[45].sim_transform(hso);

  // align this c2 along the x axis
  t1rep[50] = t1rep[50].sim_transform(t1rep[2]).sim_transform(t1rep[2]);
  grep[50] = grep[50].sim_transform(grep[2]).sim_transform(grep[2]);
  hrep[50] = hrep[50].sim_transform(hrep[2]).sim_transform(hrep[2]);

  // and rotate that by 2pi/5 about the z axis 4 times
  for (i=51; i < 55; i++) {
    t1rep[i] = t1rep[i-1].sim_transform(t1rep[1]);
    grep[i] = grep[i-1].sim_transform(grep[1]);
    hrep[i] = hrep[i-1].sim_transform(hrep[1]);
  }

  // finally, take a C2 about the y axis, and rotate it by 2pi/3 about the
  // xz axis, and align it along the x axis
  t1rep[55] = t1rep[45].sim_transform(t1rep[35]).sim_transform(t1rep[1]);
  grep[55] = grep[45].sim_transform(grep[35]).sim_transform(grep[1]);
  hrep[55] = hrep[45].sim_transform(hrep[35]).sim_transform(hrep[1]);

  // and then rotate that by 2pi/5 about the z axis 4 times
  for (i=56; i < 60; i++) {
    t1rep[i] = t1rep[i-1].sim_transform(t1rep[1]);
    grep[i] = grep[i-1].sim_transform(grep[1]);
    hrep[i] = hrep[i-1].sim_transform(hrep[1]);
  }

  t2rep[50] = t1rep[55];
  t2rep[55] = t1rep[50];
  
  for (i=51; i < 55; i++) {
    t2rep[i] = t2rep[i-1].sim_transform(t2rep[1]);
    t2rep[i+5] = t2rep[i+4].sim_transform(t2rep[1]);
  }
#endif
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
#if 0
  int i,j,k;

  Rep *t1rep = new Rep[60];
  Rep *t2rep = new Rep[60];
  Rep *grep = new Rep[60];
  Rep *hrep = new Rep[60];
  
  for (i=0; i < 60; i++) {
    t1rep[i].set_dim(3);
    t2rep[i].set_dim(3);
    grep[i].set_dim(4);
    hrep[i].set_dim(5);
  }
  
  // i_ops gives us all the symmetry operations we need
  i_ops(t1rep,t2rep,grep,hrep);

  IrreducibleRepresentation ira(g,1,"A");
  IrreducibleRepresentation ir1(g,3,"T1");
  IrreducibleRepresentation ir2(g,3,"T2");
  IrreducibleRepresentation irg(g,4,"G");
  IrreducibleRepresentation irh(g,5,"H");
    
  ir1.nrot_ = 1;
  ir1.ntrans_ = 1;

  for (i=0; i < g; i++) {
    ira.rep[i] = ira.proj[0][i] = 1;

    ir1.rep[i] = t1rep[i].trace();
    ir2.rep[i] = t2rep[i].trace();
    irg.rep[i] = grep[i].trace();
    irh.rep[i] = hrep[i].trace();

    for (j=0; j < 3; j++) {
      for (k=0; k < 3; k++) {
        ir1.proj[3*j+k][i] = t1rep[i][k][j];
        ir2.proj[3*j+k][i] = t2rep[i][k][j];
      }
    }

    for (j=0; j < 4; j++)
      for (k=0; k < 4; k++)
        irg.proj[4*j+k][i] = grep[i][k][j];

    for (j=0; j < 5; j++)
      for (k=0; k < 5; k++)
        irh.proj[5*j+k][i] = hrep[i][k][j];

    symop[i] = t1rep[i];
  }

  gamma_[0] = ira;
  gamma_[1] = ir1;
  gamma_[2] = ir2;
  gamma_[3] = irg;
  gamma_[4] = irh;

  delete[] t1rep;
  delete[] t2rep;
  delete[] grep;
  delete[] hrep;
#endif
}


void CharacterTable::ih()
{
#if 0
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
#endif
}
