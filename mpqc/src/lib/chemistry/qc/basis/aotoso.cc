
#include <util/container/array.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/rot.h>

#include <chemistry/qc/integral/integralv2.h>
#include <util/group/picl.h>
extern "C" {
#include <util/misc/timer.h>
}

////////////////////////////////////////////////////////////////////////////

struct lin_comb {
    int ns;
    int f0;
    int mapf0;
    double **c;

    lin_comb(int ins, int if0, int imf0) : ns(ins), f0(if0), mapf0(imf0) {
      int i;
      
      c = new double*[ns];
      for (i=0; i < ns; i++) {
        c[i] = new double[ns];
        memset(c[i],0,sizeof(double)*ns);
      }
    }

    ~lin_comb() {
      if (c) {
        for (int i=0; i < ns; i++)
          if (c[i])
            delete[] c[i];
        delete[] c;
        c=0;
      }
    }

    void print() const {
      int i;
      for (i=0; i < ns; i++)
        printf(" %10d",mapf0+i);
      printf("\n");
      
      for (i=0; i < ns; i++) {
        printf("%2d",f0+i);
        for (int j=0; j < ns; j++)
          printf(" %10.7f",c[i][j]);
        printf("\n");
      }
    }
};

struct contribution {
    int bfn;
    double coef;

    contribution() {}
    contribution(int b, double c) : bfn(b), coef(c) {}
};
ARRAY_dec(contribution);
ARRAY_def(contribution);

struct SO {
    Arraycontribution cont;
};
ARRAY_dec(SO);
ARRAY_def(SO);

struct SO_block {
    ArraySO so;

    void add(SO& s, int i) {
      if (i >= so.length())
        so.reset_length(i+1);
      so[i] = s;
    }
};


RefSCMatrix
PetiteList::aotoso()
{
  tim_enter("aotoso");
  
  int i, d, ii, jj, g, fn, s, c, ir, f;

  GaussianBasisSet& gbs = *_gbs.pointer();
  Molecule& mol = *gbs.molecule().pointer();

  // create the character table for the point group
  CharacterTable ct = mol.point_group().char_table();
  SymmetryOperation so;

  // ncomp is the number of symmetry blocks we have
  int ncomp=0;
  for (i=0; i < _nirrep; i++)
    ncomp += ct.gamma(i).degeneracy();
  
  // saoelem is the current SO in a block
  int *saoelem = new int[ncomp];
  memset(saoelem,0,sizeof(int)*ncomp);

  int *whichir = new int[ncomp];
  int *whichcmp = new int[ncomp];
  for (i=ii=0; i < _nirrep; i++) {
    for (int j=0; j < ct.gamma(i).degeneracy(); j++,ii++) {
      whichir[ii] = i;
      whichcmp[ii] = j;
    }
  }
  
  // SOs is an array of SO_blocks which holds the redundant SO's
  SO_block *SOs = new SO_block[ncomp];
  for (i=0; i < ncomp; i++) {
    ir = whichir[i];
    SOs[i].so.set_length(_nbf_in_ir[ir]);
  }
  
  double *blc = new double[gbs.nbasis()];
  lin_comb **lc = new lin_comb*[_ng];

  // loop over all unique shells
  tim_enter("doit");
  for (i=0; i < _natom; i++) {
    for (s=0; s < gbs.nshell_on_center(i); s++) {
      int shell_i = gbs.shell_on_center(i,s);
      
      // here I'm looking for the first shell in a group of equivalent
      // shells
      for (g=0; g < _ng; g++)
        if (_shell_map[shell_i][g] < shell_i)
          break;
      
      // means I broke out of the above loop
      if (g != _ng)
        continue;
      
      // we now have a unique shell.  now operate on this shell with all
      // _ng symmetry operations.  the basis function that each basis
      // function in this shell is mapped to by symmetry operation g
      // is stored in lc[g]
      
      tim_enter("lin comb");
      for (g=0; g < _ng; g++) {
        so = ct.symm_operation(g);
        int j = _atom_map[i][g];

        int func_i = gbs.shell_to_function(gbs.shell_on_center(i,s));
        int func_j = gbs.shell_to_function(gbs.shell_on_center(j,s));

        lc[g] = new lin_comb(gbs(i,s).nfunction(),func_i,func_j);
        lin_comb& lcg = *lc[g];
        
        int fi=0;
        for (c=0; c < gbs(i,s).ncontraction(); c++) {
          int am=gbs(i,s).am(c);

          if (am==0) {
            lcg.c[fi][fi] = 1.0;
          } else {
            Rotation rr(am,so,gbs(i,s).is_pure(c));
            for (ii=0; ii < rr.dim(); ii++) {
              for (jj=0; jj < rr.dim(); jj++)
                lcg.c[fi+ii][fi+jj] = rr(ii,jj);
            }
          }

          fi += gbs(i,s).nfunction(c);
          func_i += gbs(i,s).nfunction(c);
          func_j += gbs(i,s).nfunction(c);
        }
#if 0
        lcg.print();
#endif
      }
      tim_exit("lin comb");

      
      // now operate on the linear combinations with the projection operators
      // for each irrep to form the SO's
      tim_enter("form");
      int irnum=0;
      for (ir=0; ir < ct.nirrep(); ir++) {
        for (fn=0; fn < gbs(i,s).nfunction(); fn++) {
          for (d=0; d < ct.gamma(ir).degeneracy(); d++) {
            // if this is a point group with a complex E representation,
            // then only do the first set of projections for E
            if (d && ct.complex() && ct.gamma(ir).degeneracy()==2)
              break;
            
            for (int comp=0; comp < ct.gamma(ir).degeneracy(); comp++) {

              // form the projection for this irrep
              tim_enter("project");
              memset(blc,0,sizeof(double)*gbs.nbasis());

              for (g=0; g < _ng; g++) {
                double ccdg = ct.gamma(ir).p(comp,d,g);

                if (fabs(ccdg) < 1.0e-5)
                  continue;
                
                double *coef_fn = lc[g]->c[fn];
                double *blcp = &blc[lc[g]->mapf0];

                for (f=0; f < lc[g]->ns; f++)
                  *blcp++ += ccdg * *coef_fn++;
              }
              tim_exit("project");

              // find out how many nonzero elements there are
              tim_enter("zero");
              int nonzero=0;
              double *blcp = blc;
              for (ii=gbs.nbasis(); ii; ii--,blcp++)
                if ((*blcp* *blcp) > 0.0009)
                  nonzero++;
              tim_exit("zero");

              if (!nonzero)
                continue;
              
              tim_enter("tso");
              SO tso;
              tso.cont.set_length(nonzero);
              jj=0;
              blcp = blc;
              for (int iii=gbs.nbasis(), ii=0; iii; iii--, ii++, blcp++) {
                if ((*blcp* *blcp) > 0.0009) {
                  tso.cont[jj].bfn = ii;
                  tso.cont[jj].coef = *blcp;
                  jj++;
                }
              }
              tim_exit("tso");

              // normalize the linear combination if c1 is non-zero
              tim_enter("norm");
              double c1=0;
              for (ii=0; ii < tso.cont.length(); ii++)
                c1 += tso.cont[ii].coef*tso.cont[ii].coef;

              c1 = 1.0/sqrt(c1);
              
              for (ii=0; ii < tso.cont.length(); ii++)
                tso.cont[ii].coef *= c1;
              tim_exit("norm");

              // add this SO to the appropriate block
              tim_enter("insert");
              SOs[irnum+comp].add(tso,saoelem[irnum+comp]);
              tim_exit("insert");
            
              saoelem[irnum+comp]++;
            }
          }
        }

        irnum += ct.gamma(ir).degeneracy();
      }
      tim_exit("form");

      for (g=0; g < _ng; g++)
        delete lc[g];
    }
  }
  tim_exit("doit");

  delete[] lc;

  // ok, now let's actually stick all the neato SO information into a
  // matrix
  tim_enter("finish");

  RefSCMatrix ret = AO_basisdim()->create_matrix(SO_basisdim());
  ret.assign(0.0);

  BlockedSCMatrix *retp = BlockedSCMatrix::castdown(ret);
  if (!retp) {
    fprintf(stderr,"PetiteList::aotoso: bad things happening <shudder>\n");
    fprintf(stderr,"ret is not a blocked matrix\n");
    abort();
  }
  
  for (i=0; i < ncomp; i++) {
    ir = whichir[i];
    int cmp = whichcmp[i];

    if (saoelem[i] == _nbf_in_ir[ir]) {
      // if we found the right number, do nothing else

      for (ii=0; ii < saoelem[i]; ii++) {
        Arraycontribution& cont = SOs[i].so[ii].cont;
        
        for (jj=0; jj < cont.length(); jj++)
          retp->mats_[i].set_element(cont[jj].bfn, ii, cont[jj].coef);
      }

    } else if (saoelem[i] < _nbf_in_ir[ir]) {
      // if we found too few, there are big problems
      
      fprintf(stderr,"trouble making SO's for irrep %s\n",
              ct.gamma(ir).symbol());
      fprintf(stderr,"  only found %d out of %d SO's\n",
              saoelem[i], _nbf_in_ir[ir]);
      abort();

    } else {
      for (ii=0; ii < _nbf_in_ir[ir]; ii++) {
        Arraycontribution& cont = SOs[i].so[ii].cont;
        
        for (jj=0; jj < cont.length(); jj++)
          retp->mats_[i].set_element(cont[jj].bfn, ii, cont[jj].coef);
      }
    }
  }
  tim_exit("finish");

#if 0
  tim_enter("orthog");
  RefSymmSCMatrix S(gbs.basisdim());
  S.assign(0.0);

  tim_enter("overlap");
  RefSCElementOp op = new GaussianOverlapIntv2(_gbs);
  S.element_op(op);
  op=0;
  tim_exit("overlap");

  tim_enter("schmidt");
  RefSCMatrix bs = ret.clone();
  BlockedSCMatrix::castdown(bs)->block(S);
  S=0;
  ret->schmidt_orthog(bs.pointer(), gbs.nbasis());
  bs=0;
  tim_exit("schmidt");
  tim_exit("orthog");
#endif

  delete[] blc;
  delete[] saoelem;
  delete[] whichir;

  tim_exit("aotoso");

  {
    int nproc,me,host;
    int top,ord,dir;
    RefMessageGrp grp;

    grp = MessageGrp::get_default_messagegrp();

    open0_messagegrp(&nproc,&me,&host,grp);
    setarc0(&nproc,&top,&ord,&dir);

    tim_print(0);
  }

  return ret;
}
