
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/rot.h>

#include <chemistry/qc/integral/integralv2.h>

////////////////////////////////////////////////////////////////////////////

contribution::contribution(int b, double c) : bfn(b), coef(c)
{
}
ARRAY_def(contribution);

int
SO::equiv(const SO& so)
{
  int i;
      
  if (so.cont.length() != cont.length())
    return 0;

  double c=0;
  for (i=0; i < cont.length(); i++) {
    if (cont[i].bfn != so.cont[i].bfn)
      return 0;
    c += cont[i].coef*so.cont[i].coef;
  }
      
  // if the overlap == 1.0, they're equal (SO's should have been
  // normalized by now)
  if (fabs(fabs(c)-1.0) < 1.0e-3)
    return 1;

  return 0;
}
ARRAY_def(SO);

int
SO_block::add(SO& s, int i)
{
  // first check to see if s is already here
  for (int j=0; j < ((i < so.length()) ? i : so.length()); j++)
    if (so[j].equiv(s))
      return 0;
      
  if (i >= so.length())
    so.reset_length(i+1);
  so[i] = s;

  return 1;
}

void
SO_block::print(const char *title)
{
  printf("SO block %s\n",title);
  for (int i=0; i < so.length(); i++) {
    printf("SO %d\n",i+1);
    for (int j=0; j < so[i].cont.length(); j++)
      printf(" %10d",so[i].cont[j].bfn);
    printf("\n");
    for (int j=0; j < so[i].cont.length(); j++)
      printf(" %10.7f",so[i].cont[j].coef);
    printf("\n");
  }
}

////////////////////////////////////////////////////////////////////////////

SymmetryOrbitals::SymmetryOrbitals() :
  sos_(0), som_(0)
{
}

SymmetryOrbitals::SymmetryOrbitals(const RefGaussianBasisSet& gbs) :
  gbs_(gbs)
{
  PetiteList pl(gbs);
  sos_ = pl.aotoso();

  CharacterTable ct = gbs->molecule()->point_group().char_table();
  
  som_ = new SO_block*[ct.nirrep()];

  SO_block *sp = sos_;
  for (int i=0; i < ct.nirrep(); i++) {
    som_[i] = sp;
    sp += ct.gamma(i).degeneracy();
  }
}

SymmetryOrbitals::~SymmetryOrbitals()
{
  gbs_=0;
  if (sos_) {
    delete[] sos_;
    sos_=0;
  }
  if (som_) {
    delete[] som_;
    som_=0;
  }
}

void
SymmetryOrbitals::print(FILE *out)
{
  char label[80];
  
  fprintf(out,"SymmetryOrbitals:\n");
  CharacterTable ct = gbs_->molecule()->point_group().char_table();

  for (int i=0; i < ct.nirrep(); i++) {
    fprintf(out,"  irrep %s:\n",ct.gamma(i).symbol());
    for (int j=0; j < ct.gamma(i).degeneracy(); j++) {
      sprintf(label,"%d",j+1);
      som_[i][j].print(label);
    }
    fprintf(out,"\n");
  }
}
    
RefBlockedSCDimension
SymmetryOrbitals::AO_basisdim()
{
  RefBlockedSCDimension ret =
    new BlockedSCDimension(gbs_->matrixkit(),gbs_->nbasis());

  return ret;
}

RefBlockedSCDimension
SymmetryOrbitals::SO_basisdim()
{
  PetiteList pl(gbs_);
  return pl.SO_basisdim();
}
  
////////////////////////////////////////////////////////////////////////////

AOSO_Transformation::AOSO_Transformation(
  const RefGaussianBasisSet& gbs)
  : sos(gbs)
{
}

void
AOSO_Transformation::process(SCMatrixBlockIter&)
{
  fprintf(stderr,"AOSO_Transformation::process(SCMatrixBlockIter&):"
          " can only handle RectBlocks\n");
  abort();
}

void
AOSO_Transformation::process(SCMatrixRectBlock* blk)
{
  SO_block& sob = sos.sos_[current_block()];

  if (blk->jend > sob.so.length()) {
    fprintf(stderr,"AOSO_Transformation::process(SCMatrixRectBlock*):"
            " blk is wrong size\n");
    abort();
  }
  
  int isize = blk->iend - blk->istart;
  int jsize = blk->jend - blk->jstart;
  
  memset(blk->data,0,isize*jsize*sizeof(double));

  for (int j=blk->jstart; j < blk->jend; j++) {
    Arraycontribution& tc = sob.so[j].cont;

    for (int i=0; i < tc.length(); i++) {
      int bfn = tc[i].bfn;
      if (bfn >= blk->istart && bfn < blk->iend)
        blk->data[(bfn-blk->istart)*jsize+j] = tc[i].coef;
    }
  }
}

////////////////////////////////////////////////////////////////////////////

AOSO_Unit::AOSO_Unit(const RefBlockedSCDimension& rd1,
                     const RefBlockedSCDimension& rd2)
  : d1(rd1), d2(rd2)
{
}

void
AOSO_Unit::process(SCMatrixBlockIter&)
{
  fprintf(stderr,"AOSO_Unit::process(SCMatrixBlockIter&):"
          " can only handle RectBlocks\n");
  abort();
}

void
AOSO_Unit::process(SCMatrixRectBlock* blk)
{
  int fi = (d1->nblocks()==1) ? d1->first(0) : d1->first(current_block());
  int fj = (d2->nblocks()==1) ? d2->first(0) : d2->first(current_block());

  int isize = blk->iend - blk->istart;
  int jsize = blk->jend - blk->jstart;
  
  memset(blk->data,0,isize*jsize*sizeof(double));

  for (int i=blk->istart; i < blk->iend; i++)
    for (int j=blk->jstart; j < blk->jend; j++)
      if (i+fi==j+fj)
        blk->data[(i-blk->istart)*jsize+j-blk->jstart]=1;
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

struct shell_overlap {
    int nf;
    double **c;

    void init(GaussianShell& gs) {
      int i,j;

      if (c) {
        for (i=0; i < nf; i++)
          if (c[i]) delete[] c[i];
        delete[] c;
      }

      nf = gs.nfunction();
      c = new double*[nf];
      for (i=0; i < nf; i++) {
        c[i] = new double[nf];
        memset(c[i],0,sizeof(double)*nf);
      }

      i=j=0;
      for (int n=0; n < gs.ncontraction(); n++) {
        for (int fi=0; fi < gs.nfunction(n); fi++)
          for (int fj=0; fj < gs.nfunction(n); fj++)
            c[i+fi][j+fj] = gs.relative_overlap(n,fi,fj);
        i += gs.nfunction(n);
        j += gs.nfunction(n);
      }
    }

    shell_overlap() : nf(0), c(0) {}

    shell_overlap(GaussianShell& gs) : nf(0), c(0) {
      init(gs);
    }

    ~shell_overlap() {
      if (c) {
        for (int i=0; i < nf; i++)
          if (c[i]) delete[] c[i];
        delete[] c;
        c=0;
      }
    }

    void print() const {
      for (int fi=0; fi < nf; fi++) {
        printf(" %2d ",fi+1);
        for (int fj=0; fj < nf; fj++)
          printf(" %10.7f",c[fi][fj]);
        printf("\n");
      }
    }
};


////////////////////////////////////////////////////////////////////////////

SO_block *
PetiteList::aotoso()
{
  int i, d, ii, jj, g, fn, s, c, ir, f;

  GaussianBasisSet& gbs = *gbs_.pointer();
  Molecule& mol = *gbs.molecule().pointer();

  // create the character table for the point group
  CharacterTable ct = mol.point_group().char_table();
  SymmetryOperation so;

  // ncomp is the number of symmetry blocks we have. for point groups with
  // complex E representations, this will be cut in two eventually
  int ncomp=0;
  for (i=0; i < nirrep_; i++)
    ncomp += ct.gamma(i).degeneracy();
  
  // saoelem is the current SO in a block
  int *saoelem = new int[ncomp];
  memset(saoelem,0,sizeof(int)*ncomp);

  int *whichir = new int[ncomp];
  int *whichcmp = new int[ncomp];
  for (i=ii=0; i < nirrep_; i++) {
    for (int j=0; j < ct.gamma(i).degeneracy(); j++,ii++) {
      whichir[ii] = i;
      whichcmp[ii] = j;
    }
  }
  
  // SOs is an array of SO_blocks which holds the redundant SO's
  SO_block *SOs = new SO_block[ncomp];
  for (i=0; i < ncomp; i++) {
    ir = whichir[i];
    int len = (ct.gamma(ir).complex()) ? nbf_in_ir_[ir]/2 : nbf_in_ir_[ir];
    SOs[i].so.set_length(len);
  }
  
  double *blc = new double[gbs.nbasis()];
  lin_comb **lc = new lin_comb*[ng_];

  shell_overlap sov;

  // loop over all unique shells
  for (i=0; i < natom_; i++) {
    for (s=0; s < gbs.nshell_on_center(i); s++) {
      int shell_i = gbs.shell_on_center(i,s);
      
      // here I'm looking for the first shell in a group of equivalent
      // shells
      for (g=0; g < ng_; g++)
        if (shell_map_[shell_i][g] < shell_i)
          break;
      
      // means I broke out of the above loop
      if (g != ng_)
        continue;
      
      // we now have a unique shell.  now operate on this shell with all
      // ng_ symmetry operations.  the basis function that each basis
      // function in this shell is mapped to by symmetry operation g
      // is stored in lc[g]
      
      // test to see if there are any high am cartesian functions in this
      // shell
      int cartfunc=0;
      for (c=0; c < gbs(i,s).ncontraction(); c++) {
        if (gbs(i,s).am(c) > 1 && gbs(i,s).is_cartesian(c)) {
          cartfunc=1;
          sov.init(gbs(i,s));
          break;
        }
      }

      // for now don't allow symmetry with cartesian functions...I just can't
      // seem to get them working.
      if (cartfunc) {
        fprintf(stderr,"PetiteList::aotoso:  cannot yet handle symmetry for"
                " angular momentum >= 2\n");
        abort();
      }

      for (g=0; g < ng_; g++) {
        so = ct.symm_operation(g);
        int j = atom_map_[i][g];

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
      }

      // now operate on the linear combinations with the projection operators
      // for each irrep to form the SO's
      int irnum=0;
      for (ir=0; ir < ct.nirrep(); ir++) {
        int cmplx = (ct.complex() && ct.gamma(ir).complex());
        
        for (fn=0; fn < gbs(i,s).nfunction(); fn++) {
          for (d=0; d < ct.gamma(ir).degeneracy(); d++) {
            // if this is a point group with a complex E representation,
            // then only do the first set of projections for E
            if (d && cmplx)
              break;
            
            for (int comp=0; comp < ct.gamma(ir).degeneracy(); comp++) {

              // form the projection for this irrep
              memset(blc,0,sizeof(double)*gbs.nbasis());

              for (g=0; g < ng_; g++) {
                double ccdg = ct.gamma(ir).p(comp,d,g);

                if (fabs(ccdg) < 1.0e-5)
                  continue;
                
                double *coef_fn = lc[g]->c[fn];
                double *blcp = &blc[lc[g]->mapf0];

                for (f=0; f < lc[g]->ns; f++)
                  *blcp++ += ccdg * *coef_fn++;
              }

              // find out how many nonzero elements there are
              int nonzero=0;
              double *blcp = blc;
              for (ii=gbs.nbasis(); ii; ii--,blcp++)
                if ((*blcp* *blcp) > 0.0009)
                  nonzero++;

              if (!nonzero)
                continue;
              
              // copy the nonzero bits to tso
              SO tso;
              tso.cont.set_length(nonzero);
              ii=jj=0; blcp = blc;
              for (int iii=gbs.nbasis(); iii; iii--, ii++, blcp++) {
                if ((*blcp* *blcp) > 0.0009) {
                  tso.cont[jj].bfn = ii;
                  tso.cont[jj].coef = *blcp;
                  jj++;
                }
              }

              // normalize the linear combination
              double c1=0;
              for (ii=0; ii < tso.cont.length(); ii++)
                c1 += tso.cont[ii].coef*tso.cont[ii].coef;

              c1 = 1.0/sqrt(c1);
              
              for (ii=0; ii < tso.cont.length(); ii++)
                tso.cont[ii].coef *= c1;

              // add this SO to the appropriate block
              if (SOs[irnum+comp].add(tso,saoelem[irnum+comp]))
                saoelem[irnum+comp]++;
            }
          }
        }

        irnum += ct.gamma(ir).degeneracy();
      }

      for (g=0; g < ng_; g++)
        delete lc[g];
    }
  }

  delete[] lc;
  delete[] blc;

  for (i=0; i < ncomp; i++) {
    ir = whichir[i];
    int cmp = whichcmp[i];
    int scal = ct.gamma(ir).complex() ? 2 : 1;

    if (saoelem[i] < nbf_in_ir_[ir]/scal) {
      // if we found too few, there are big problems
      
      fprintf(stderr,"trouble making SO's for irrep %s\n",
              ct.gamma(ir).symbol());
      fprintf(stderr,"  only found %d out of %d SO's\n",
              saoelem[i], nbf_in_ir_[ir]/scal);
      SOs[i].print("");

      abort();

    } else if (saoelem[i] > nbf_in_ir_[ir]/scal) {
      // there are some redundant so's left...need to do something to get
      // the elements we want
      
      fprintf(stderr,"trouble making SO's for irrep %s\n",
              ct.gamma(ir).symbol());
      fprintf(stderr,"  found %d SO's, but there should only be %d\n",
              saoelem[i], nbf_in_ir_[ir]/scal);
      SOs[i].print("");
      abort();
    }
  }

  if (ct.complex()) {
    SO_block *nSOs = new SO_block[nblocks_];

    int in=0;

    for (i=ir=0; ir < nirrep_; ir++) {
      if (ct.gamma(ir).complex()) {
        nSOs[in].so.set_length(nbf_in_ir_[ir]);
        int k;
        for (k=0; k < SOs[i].so.length(); k++)
          nSOs[in].add(SOs[i].so[k],k);
        i++;

        for (int j=0; j < SOs[i].so.length(); j++,k++)
          nSOs[in].add(SOs[i].so[j],k);

        i++;
        in++;
      } else {
        for (int j=0; j < ct.gamma(ir).degeneracy(); j++,i++,in++) {
          nSOs[in].so.set_length(nbf_in_ir_[ir]);
          for (int k=0; k < SOs[i].so.length(); k++)
            nSOs[in].add(SOs[i].so[k],k);
        }
      }
    }

    SO_block *tmp= SOs;
    SOs = nSOs;
    delete[] tmp;
  }

  delete[] saoelem;
  delete[] whichir;
  delete[] whichcmp;

  return SOs;
}
