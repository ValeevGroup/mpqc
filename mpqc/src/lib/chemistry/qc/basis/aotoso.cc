
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/rot.h>

////////////////////////////////////////////////////////////////////////////

contribution::contribution()
{
}

contribution::~contribution()
{
}

contribution::contribution(int b, double c) : bfn(b), coef(c)
{
}

////////////////////////////////////////////////////////////////////////////

SO::SO() : cont(0), len(0)
{
}

SO::SO(int l) : cont(0), len(0)
{
  set_length(l);
}

SO::~SO()
{
  set_length(0);
}

SO&
SO::operator=(const SO& so)
{
  set_length(so.len);
  for (int i=0; i < len; i++)
    cont[i] = so.cont[i];
  return *this;
}

void
SO::set_length(int l)
{
  len=l;
  if (cont) {
    delete[] cont;
    cont=0;
  }

  if (l)
    cont = new contribution[l];
}

void
SO::reset_length(int l)
{
  if (l <= len)
    return;
  
  l=l+10;
  
  contribution *newcont = new contribution[l];
  
  if (cont) {
    for (int i=0; i < len; i++)
      newcont[i] = cont[i];
    
    delete[] cont;
  }

  cont=newcont;
  len=l;
}

int
SO::equiv(const SO& so)
{
  int i;
      
  if (so.len != len)
    return 0;

  double c=0;
  for (i=0; i < len; i++) {
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

////////////////////////////////////////////////////////////////////////////

SO_block::SO_block() : so(0), len(0)
{
}

SO_block::SO_block(int l) : so(0), len(0)
{
  set_length(l);
}

SO_block::~SO_block()
{
  set_length(0);
}

void
SO_block::set_length(int l)
{
  len=l;
  if (so) {
    delete[] so;
    so=0;
  }

  if (l)
    so = new SO[l];
}

void
SO_block::reset_length(int l)
{
  SO *newso = new SO[l];
  
  if (so) {
    for (int i=0; i < len; i++)
      newso[i] = so[i];
    
    delete[] so;
  }

  so=newso;
  len=l;
}

int
SO_block::add(SO& s, int i)
{
  // first check to see if s is already here
  for (int j=0; j < ((i < len) ? i : len); j++)
    if (so[j].equiv(s))
      return 0;
      
  if (i >= len)
    reset_length(i+1);
  so[i] = s;

  return 1;
}

void
SO_block::print(const char *title)
{
  int i,j;
  printf("SO block %s\n",title);
  for (i=0; i < len; i++) {
    printf("SO %d\n",i+1);
    for (j=0; j < so[i].len; j++)
      printf(" %10d",so[i].cont[j].bfn);
    printf("\n");
    for (j=0; j < so[i].len; j++)
      printf(" %10.7f",so[i].cont[j].coef);
    printf("\n");
  }
}

static int
soblock_length(SO_block *sob, int nb)
{
  int lt=0;
  for (int b=0; b < nb; b++) {
    int lb = 0;
    for (int i=0; i < sob[b].len; i++) {
      lb += sizeof(contribution)*sob[b].so[i].len;
    }
    printf("sizeof block %3d = %12d\n",b+1,lb);
    lt += lb;
  }
  printf("total size = %12d\n",lt);
  return lt;
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

  Molecule& mol = *gbs_.molecule().pointer();

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
      int wi = ct.which_irrep(ii);
      int wc = ct.which_comp(ii);
      whichir[ii] = i;
      whichcmp[ii] = j;
    }
  }
  
  // SOs is an array of SO_blocks which holds the redundant SO's
  SO_block *SOs = new SO_block[ncomp];
  for (i=0; i < ncomp; i++) {
    ir = whichir[i];
    int len = (ct.gamma(ir).complex()) ? nbf_in_ir_[ir]/2 : nbf_in_ir_[ir];
    SOs[i].set_length(len);
  }
  
  double *blc = new double[gbs_.nbasis()];
  lin_comb **lc = new lin_comb*[ng_];

  shell_overlap sov;

  // loop over all unique shells
  for (i=0; i < natom_; i++) {
    for (s=0; s < gbs_.nshell_on_center(i); s++) {
      int shell_i = gbs_.shell_on_center(i,s);
      
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
      for (c=0; c < gbs_(i,s).ncontraction(); c++) {
        if (gbs_(i,s).am(c) > 1 && gbs_(i,s).is_cartesian(c)) {
          cartfunc=1;
          sov.init(gbs_(i,s));
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

        int func_i = gbs_.shell_to_function(gbs_.shell_on_center(i,s));
        int func_j = gbs_.shell_to_function(gbs_.shell_on_center(j,s));

        lc[g] = new lin_comb(gbs_(i,s).nfunction(),func_i,func_j);
        lin_comb& lcg = *lc[g];
        
        int fi=0;
        for (c=0; c < gbs_(i,s).ncontraction(); c++) {
          int am=gbs_(i,s).am(c);

          if (am==0) {
            lcg.c[fi][fi] = 1.0;
          } else {
            Rotation rr(am,so,gbs_(i,s).is_pure(c));
            for (ii=0; ii < rr.dim(); ii++) {
              for (jj=0; jj < rr.dim(); jj++)
                lcg.c[fi+ii][fi+jj] = rr(ii,jj);
            }
          }

          fi += gbs_(i,s).nfunction(c);
          func_i += gbs_(i,s).nfunction(c);
          func_j += gbs_(i,s).nfunction(c);
        }
      }

      // now operate on the linear combinations with the projection operators
      // for each irrep to form the SO's
      int irnum=0;
      for (ir=0; ir < ct.nirrep(); ir++) {
        int cmplx = (ct.complex() && ct.gamma(ir).complex());
        
        for (fn=0; fn < gbs_(i,s).nfunction(); fn++) {
          for (d=0; d < ct.gamma(ir).degeneracy(); d++) {
            // if this is a point group with a complex E representation,
            // then only do the first set of projections for E
            if (d && cmplx)
              break;
            
            for (int comp=0; comp < ct.gamma(ir).degeneracy(); comp++) {

              // form the projection for this irrep
              memset(blc,0,sizeof(double)*gbs_.nbasis());

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
              for (ii=gbs_.nbasis(); ii; ii--,blcp++)
                if ((*blcp* *blcp) > 0.0009)
                  nonzero++;

              if (!nonzero)
                continue;
              
              // copy the nonzero bits to tso
              SO tso;
              tso.set_length(nonzero);
              ii=jj=0; blcp = blc;
              for (int iii=gbs_.nbasis(); iii; iii--, ii++, blcp++) {
                if ((*blcp* *blcp) > 0.0009) {
                  tso.cont[jj].bfn = ii;
                  tso.cont[jj].coef = *blcp;
                  jj++;
                }
              }

              // normalize the linear combination
              double c1=0;
              for (ii=0; ii < tso.len; ii++)
                c1 += tso.cont[ii].coef*tso.cont[ii].coef;

              c1 = 1.0/sqrt(c1);
              
              for (ii=0; ii < tso.len; ii++)
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
        nSOs[in].set_length(nbf_in_ir_[ir]);
        int k;
        for (k=0; k < SOs[i].len; k++)
          nSOs[in].add(SOs[i].so[k],k);
        i++;

        for (int j=0; j < SOs[i].len; j++,k++)
          nSOs[in].add(SOs[i].so[j],k);

        i++;
        in++;
      } else {
        for (int j=0; j < ct.gamma(ir).degeneracy(); j++,i++,in++) {
          nSOs[in].set_length(nbf_in_ir_[ir]);
          for (int k=0; k < SOs[i].len; k++)
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

static void
do_transform(const RefSymmSCMatrix& skel, const RefSymmSCMatrix& sym,
             const GaussianBasisSet& gbs_, PetiteList& pl,
             CharacterTable& ct)
{
  int b,c,i,j,k,l;
  
  BlockedSymmSCMatrix *lsym = BlockedSymmSCMatrix::castdown(sym.pointer());
  BlockedSymmSCMatrix *lskl = BlockedSymmSCMatrix::castdown(skel.pointer());

  SO_block *sos = pl.aotoso();

  SCMatrixLTriBlock *sklblk, *symblk;
  SCMatrixLTriSubBlock *symsblk;

  RefSCMatrixSubblockIter skliter =
    skel->all_blocks(SCMatrixSubblockIter::Read);
  RefSCMatrixSubblockIter symiter;

#if 1
  //printf("\n sizeof sos:\n");
  //soblock_length(sos, lsym->nblocks());
  
  SO_block *SU = new SO_block[lsym->nblocks()];
  
  for (skliter->begin(); skliter->ready(); skliter->next()) {
    sklblk = SCMatrixLTriBlock::castdown(skliter->block());

    double *skldata = sklblk->data;

    int kstart = sklblk->start;
    int kend = sklblk->end;
    
    for (b=0; b < lsym->nblocks(); b++) {
      if (lsym->block(b).null())
        continue;

      // form first transform SU
      SO_block& sob = sos[b];
      SO_block& sub = SU[b];
      
      sub.set_length(sob.len);

      for (j=0; j < sub.len; j++) {
        SO& soj = sob.so[j];
        SO& suj = sub.so[j];
        
        int sojl = soj.len;
        contribution *coj = soj.cont;

        int kk=0, dk=0, idk=0;
        for (k=kstart; k < kend; k++, dk++) {
          idk += dk;

          double sukj=0;
          double *skdata = skldata+idk;
          contribution *sojt = coj;

          int bfl;
          for (l=sojl; l && (bfl=(*sojt).bfn) < kend; l--,sojt++)
            if (bfl >= kstart)
              break;
          
          int dl;
          for (; l && (bfl=(*sojt).bfn) < kend && (dl=bfl-kstart) <= dk;
                 l--,sojt++) {
            sukj += (*sojt).coef * skdata[dl];
          }

          skdata = skldata+dk;
          
          for (; l && (bfl=(*sojt).bfn) < kend; l--,sojt++) {
            dl = bfl-kstart;
            sukj += (*sojt).coef * skdata[(dl*(dl+1))>>1];
          }

          if (fabs(sukj) > 1.0e-12) {
            suj.reset_length(kk+1);
            suj.cont[kk].bfn = k;
            suj.cont[kk].coef = sukj;
            kk++;
          }
        }
      }
    }
  }

  //printf("\n sizeof SU:\n");
  //soblock_length(SU, lsym->nblocks());

  // now form Sym = U~ * SU
  for (b=0; b < lsym->nblocks(); b++) {
    if (lsym->block(b).null())
      continue;
    int ir = ct.which_irrep(b);
    double skal = (double)ct.order()/(double)ct.gamma(ir).degeneracy();

    SO_block& sob = sos[b];
    SO_block& sub = SU[b];

    symiter = lsym->block(b)->local_blocks(SCMatrixSubblockIter::Write);

    for (symiter->begin(); symiter->ready(); symiter->next()) {
      // symblk can either be an LTri block or an LTriSub block
      SCMatrixBlock* blk = symiter->block();
        
      if (symblk = SCMatrixLTriBlock::castdown(blk)) {
        int ij=0;
        for (i=symblk->start; i < symblk->end; i++) {
          contribution *ci = sob.so[i].cont;
          int cilen = sob.so[i].len;

          for (j=symblk->start; j <= i; j++,ij++) {
            contribution *cj = sub.so[j].cont;
            int cjlen = sub.so[j].len;

            if (!cjlen)
              continue;
        
            int ii=0,jj=0;
            double tij = 0;

            for (k=0; k < gbs_.nbasis(); k++) {
              int subf = cj[jj].bfn;
              int ubf = ci[ii].bfn;

              if (k < ubf && k < subf)
                continue;
              else if (k < subf) {
                ii++;
              } else if (k < ubf) {
                jj++;
              } else {
                tij += cj[jj].coef*ci[ii].coef*skal;
                ii++;
                jj++;
              }

              if (ii >= cilen || jj >= cjlen)
                break;
            }

            symblk->data[ij] = tij;
          }
        }
      }
    }
  }

#else

  sym.assign(0.0);
  for (skliter->begin(); skliter->ready(); skliter->next()) {
    sklblk = SCMatrixLTriBlock::castdown(skliter->block());

    int kstart = sklblk->start;
    int kend = sklblk->end;
    
    for (b=0; b < lsym->nblocks(); b++) {
      if (lsym->block(b).null())
        continue;
      int ir = ct.which_irrep(b);
      double skal = (double)ct.order()/(double)ct.gamma(ir).degeneracy();

      symiter = lsym->block(b)->local_blocks();

      SO_block& sob = sos.sos_[b];

      for (symiter->begin(); symiter->ready(); symiter->next()) {
        // symblk can either be an LTri block or an LTriSub block
        SCMatrixBlock* blk = symiter->block();
        
        if (symblk = SCMatrixLTriBlock::castdown(blk)) {
          int isize = symblk->end - symblk->start;
          int jsize = isize;
          int ij=0;

          for (i=symblk->start; i < symblk->end; i++) {
            SO& soi = sob.so[i];
            contribution *ci = soi.cont;

            for (j=symblk->start; j <= i; j++,ij++) {
              SO& soj = sob.so[j];
              contribution *cj = soj.cont;
              double tij = 0;

              contribution *cit = ci;
              
              int kk;
              for (kk=0; kk < soi.len; kk++,cit++)
                if ((*cit).bfn >= kstart)
                  break;
              
              for (; kk < soi.len && (k=(*cit).bfn) < kend; kk++,cit++) {
                
                int dk = k - kstart;
                int idk = (dk*(dk+1))>>1;
                
                double *skdata = sklblk->data+dk;
                double *skidata = sklblk->data+idk;
                
                double cik = skal*(*cit).coef;
                double ttij=0;
                
                contribution *cjt = cj;
                
                int ll;
                for (ll=soj.len; ll ; ll--,cjt++)
                  if ((*cjt).bfn >= kstart)
                    break;
                
                int dl;
                for (; ll && (l=(*cjt).bfn) < kend && (dl=l-kstart) <= dk;
                       ll--,cjt++) {
                  ttij += (*cjt).coef*skidata[dl];
                }

                for (; ll && (l=(*cjt).bfn) < kend; ll--,cjt++) {
                  dl = l-kstart;
                  ttij += (*cjt).coef*skdata[(dl*(dl+1))>>1];
                }

                tij += cik*ttij;
              }

              symblk->data[ij] += tij;
            }
          }
        }
      }
    }
  }
#endif
}

void
PetiteList::symmetrize(const RefSymmSCMatrix& skel,
                       const RefSymmSCMatrix& sym)
{
  int b,c;

  CharacterTable ct = gbs_.molecule()->point_group().char_table();

#if 0

  RefSCMatrix aotoso(AO_basisdim(), SO_basisdim());
  RefSCElementOp sotrans = new AOSO_Transformation(gbs_);
  aotoso.element_op(sotrans);
  sotrans=0;

  BlockedSCMatrix *lu = BlockedSCMatrix::castdown(aotoso);

  for (b=0; b < lu->nblocks(); b++) {
    if (lu->block(b).null())
      continue;
    
    int ir = ct.which_irrep(b);
  
    double skal = (double)ct.order()/(double)ct.gamma(ir).degeneracy();
    skal = sqrt(skal);
    lu->block(b).scale(skal);
  }

  sym.assign(0.0);
  sym.accumulate_transform(aotoso.t(),skel);

#else
  do_transform(skel,sym,gbs_,*this,ct);
#endif

  BlockedSymmSCMatrix *la = BlockedSymmSCMatrix::castdown(sym.pointer());
  
  // loop through blocks and finish symmetrizing degenerate blocks
  for (b=0; b < la->nblocks(); b++) {
    if (la->block(b).null())
      continue;

    int ir=ct.which_irrep(b);

    if (ct.gamma(ir).degeneracy()==1)
      continue;

    if (ct.gamma(ir).complex()) {
      int nbf = nbf_in_ir_[ir]/2;
      
      RefSymmSCMatrix irrep = la->block(b).get_subblock(0, nbf-1);
      irrep.accumulate(la->block(b).get_subblock(nbf, 2*nbf-1));

      RefSCMatrix sub = la->block(b).get_subblock(nbf, 2*nbf-1, 0, nbf-1);
      RefSCMatrix subt = sub.t();
      subt.scale(-1.0);
      sub.accumulate(subt);
      subt=0;

      la->block(b).assign_subblock(irrep,  0, nbf-1);
      la->block(b).assign_subblock(irrep,nbf, 2*nbf-1);
      la->block(b).assign_subblock(sub, nbf, 2*nbf-1, 0, nbf-1);

    } else {
      RefSymmSCMatrix irrep = la->block(b).copy();
      for (c=1; c < ct.gamma(ir).degeneracy(); c++)
        irrep.accumulate(la->block(b+c));
      
      for (c=0; c < ct.gamma(ir).degeneracy(); c++)
        la->block(b+c).assign(irrep);

      b += ct.gamma(ir).degeneracy()-1;
    }
  }
}
