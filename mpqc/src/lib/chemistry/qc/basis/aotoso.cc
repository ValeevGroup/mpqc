
#include <util/misc/formio.h>

#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/shellrot.h>
#include <chemistry/qc/basis/petite.h>

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

SO::SO() : cont(0), len(0), length(0)
{
}

SO::SO(int l) : cont(0), len(0), length(0)
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
  length = so.length;
  for (int i=0; i < length; i++)
    cont[i] = so.cont[i];
  return *this;
}

void
SO::set_length(int l)
{
  len=l;
  length=l;
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
  length=l;

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
      
  if (so.length != length)
    return 0;

  double c=0;
  for (i=0; i < length; i++) {
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
  cout << node0 << indent << "SO block " << title << endl;
  for (i=0; i < len; i++) {
    cout << node0 << indent << "SO " << i+1 << endl << indent;
    for (j=0; j < so[i].length; j++)
      cout << node0 << scprintf(" %10d",so[i].cont[j].bfn);
    cout << node0 << endl << indent;
    for (j=0; j < so[i].length; j++)
      cout << node0 << scprintf(" %10.7f",so[i].cont[j].coef);
    cout << node0 << endl;
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
    cout << node0 << indent << scprintf("sizeof block %3d = %12d\n",b+1,lb);
    lt += lb;
  }
  cout << node0 << indent << scprintf("total size = %12d\n",lt);
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
      cout << node0 << indent;
      for (i=0; i < ns; i++)
        cout << node0 << scprintf(" %10d",mapf0+i);
      cout << node0 << endl;
      
      for (i=0; i < ns; i++) {
        cout << node0 << indent << scprintf("%2d",f0+i);
        for (int j=0; j < ns; j++)
          cout << node0 << scprintf(" %10.7f",c[i][j]);
        cout << node0 << endl;
      }
    }
};

////////////////////////////////////////////////////////////////////////////

SO_block *
PetiteList::aotoso_info()
{
  int i, d, ii, jj, g, fn, s, c, ir, f;

  GaussianBasisSet& gbs = *gbs_.pointer();
  Molecule& mol = *gbs.molecule().pointer();

  // create the character table for the point group
  CharacterTable ct = mol.point_group().char_table();
  SymmetryOperation so;

  if (c1_) {
    SO_block *SOs = new SO_block[1];
    SOs[0].set_length(gbs.nbasis());
    for (i=0; i < gbs.nbasis(); i++) {
      SOs[0].so[i].set_length(1);
      SOs[0].so[i].cont[0].bfn=i;
      SOs[0].so[i].cont[0].coef=1.0;
    }
    return SOs;
  }

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
  
  double *blc = new double[gbs.nbasis()];
  lin_comb **lc = new lin_comb*[ng_];

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
          break;
        }
      }

      // for now don't allow symmetry with cartesian functions...I just can't
      // seem to get them working.
      if (cartfunc && ng_ != nirrep_) {
        cerr << node0 << indent
             << "PetiteList::aotoso: cannot yet handle symmetry for "
             << " angular momentum >= 2\n";
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
            ShellRotation rr =
              ints_->shell_rotation(am,so,gbs(i,s).is_pure(c));
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
              tso.set_length(nonzero);
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
              for (ii=0; ii < tso.length; ii++)
                c1 += tso.cont[ii].coef*tso.cont[ii].coef;

              c1 = 1.0/sqrt(c1);
              
              for (ii=0; ii < tso.length; ii++)
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
      
      cerr << node0 << indent
           << scprintf("trouble making SO's for irrep %s\n",
                       ct.gamma(ir).symbol());
      cerr << node0 << indent
           << scprintf("  only found %d out of %d SO's\n",
                       saoelem[i], nbf_in_ir_[ir]/scal);
      SOs[i].print("");
      abort();

    } else if (saoelem[i] > nbf_in_ir_[ir]/scal) {
      // there are some redundant so's left...need to do something to get
      // the elements we want
      
      cerr << node0 << indent
           << scprintf("trouble making SO's for irrep %s\n",
                       ct.gamma(ir).symbol());
      cerr << node0 << indent
           << scprintf("  found %d SO's, but there should only be %d\n",
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

RefSCMatrix
PetiteList::aotoso()
{
  RefSCMatrix aoso(AO_basisdim(), SO_basisdim(), gbs_->so_matrixkit());
  aoso.assign(0.0);

  if (c1_) {
    aoso->unit();
    return aoso;
  }
  
  SO_block *sos = aotoso_info();
  
  BlockedSCMatrix *aosop = BlockedSCMatrix::castdown(aoso.pointer());

  for (int b=0; b < aosop->nblocks(); b++) {
    RefSCMatrix aosb = aosop->block(b);

    if (aosb.null())
      continue;
    
    SO_block& sob = sos[b];
    
    RefSCMatrixSubblockIter iter =
      aosb->local_blocks(SCMatrixSubblockIter::Write);

    for (iter->begin(); iter->ready(); iter->next()) {
      if (SCMatrixRectBlock::castdown(iter->block())) {
        SCMatrixRectBlock *blk = SCMatrixRectBlock::castdown(iter->block());

        int jlen = blk->jend-blk->jstart;
    
        for (int j=0; j < sob.len; j++) {
          if (j < blk->jstart || j >= blk->jend)
            continue;
      
          SO& soj = sob.so[j];
      
          for (int i=0; i < soj.len; i++) {
            int ii=soj.cont[i].bfn;
            
            if (ii < blk->istart || ii >= blk->iend)
              continue;

            blk->data[(ii-blk->istart)*jlen+(j-blk->jstart)] =
              soj.cont[i].coef;
          }
        }
      } else {
        SCMatrixRectSubBlock *blk =
          SCMatrixRectSubBlock::castdown(iter->block());

        for (int j=0; j < sob.len; j++) {
          if (j < blk->jstart || j >= blk->jend)
            continue;
      
          SO& soj = sob.so[j];
      
          for (int i=0; i < soj.len; i++) {
            int ii=soj.cont[i].bfn;
        
            if (ii < blk->istart || ii >= blk->iend)
              continue;

            blk->data[ii*blk->istride+j] = soj.cont[i].coef;
          }
        }
      }
    }
  }
  
  return aoso;
}

RefSCMatrix
PetiteList::sotoao()
{
  if (c1_)
    return aotoso();
  else
    return aotoso().i();
}

RefSymmSCMatrix
PetiteList::to_SO_basis(const RefSymmSCMatrix& a)
{
  // if C1, then do nothing
  if (c1_)
    return a;
  
  RefSymmSCMatrix aomatrix = BlockedSymmSCMatrix::castdown(a.pointer());
  if (aomatrix.null()) {
    aomatrix = gbs_->so_matrixkit()->symmmatrix(AO_basisdim());
    aomatrix->convert(a);
  }

  RefSymmSCMatrix somatrix(SO_basisdim(), gbs_->so_matrixkit());
  somatrix.assign(0.0);
  somatrix->accumulate_transform(aotoso().t(), aomatrix);

  return somatrix;
}

RefSymmSCMatrix
PetiteList::to_AO_basis(const RefSymmSCMatrix& somatrix)
{
  // if C1, then do nothing
  if (c1_)
    return somatrix;
  
  RefSymmSCMatrix aomatrix(AO_basisdim(), gbs_->so_matrixkit());
  aomatrix.assign(0.0);
  aomatrix->accumulate_transform(sotoao().t(), somatrix);

  RefSymmSCMatrix aom(gbs_->basisdim(), gbs_->matrixkit());
  aom->convert(aomatrix);

  return aom;
}

RefSCMatrix
PetiteList::evecs_to_SO_basis(const RefSCMatrix& aoev)
{
  cerr << node0 << indent
       << "PetiteList::evecs_to_SO_basis: don't work yet\n";
  abort();
  
  RefSCMatrix aoevecs = BlockedSCMatrix::castdown(aoev.pointer());
  if (aoevecs.null()) {
    aoevecs = gbs_->so_matrixkit()->matrix(AO_basisdim(), AO_basisdim());
    aoevecs->convert(aoev);
  }

  RefSCMatrix soev =  aotoso().t() * aoevecs;
  soev.print("soev");

  RefSCMatrix soevecs(SO_basisdim(), SO_basisdim(), gbs_->so_matrixkit());
  soevecs->convert(soev);

  return soevecs;
}

RefSCMatrix
PetiteList::evecs_to_AO_basis(const RefSCMatrix& soevecs)
{
  // if C1, then do nothing
  if (c1_)
    return soevecs;
  
  RefSCMatrix aoev = aotoso() * soevecs;

  RefSCMatrix aoevecs(gbs_->basisdim(), gbs_->basisdim(), gbs_->matrixkit());
  aoevecs->convert(aoev);

  return aoevecs;
}

/////////////////////////////////////////////////////////////////////////////

static void
do_transform(const RefSymmSCMatrix& skel, const RefSymmSCMatrix& sym,
             const GaussianBasisSet& gbs_, PetiteList& pl,
             CharacterTable& ct)
{
  int b,c,i,j,k,l;
  
  BlockedSymmSCMatrix *lsym = BlockedSymmSCMatrix::castdown(sym.pointer());
  BlockedSymmSCMatrix *lskl = BlockedSymmSCMatrix::castdown(skel.pointer());

  SO_block *sos = pl.aotoso_info();

  SCMatrixLTriBlock *sklblk, *symblk;
  SCMatrixLTriSubBlock *symsblk;

  RefSCMatrixSubblockIter skliter =
    skel->all_blocks(SCMatrixSubblockIter::Read);
  RefSCMatrixSubblockIter symiter;

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
        
        int sojl = soj.length;
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
          int cilen = sob.so[i].length;

          for (j=symblk->start; j <= i; j++,ij++) {
            contribution *cj = sub.so[j].cont;
            int cjlen = sub.so[j].length;

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
}

void
PetiteList::symmetrize(const RefSymmSCMatrix& skel,
                       const RefSymmSCMatrix& sym)
{
  // if C1, then do nothing
  if (c1_) {
    sym->convert(skel);
    return;
  }
  
  int b,c;

  GaussianBasisSet& gbs = *gbs_.pointer();
  CharacterTable ct = gbs.molecule()->point_group().char_table();

#if 1
  RefSymmSCMatrix bskel = BlockedSymmSCMatrix::castdown(skel.pointer());
  if (bskel.null()) {
    bskel = gbs.so_matrixkit()->symmmatrix(AO_basisdim());
    bskel->convert(skel);
  }
  
  RefSCMatrix aoso = aotoso().t();
  BlockedSCMatrix *lu = BlockedSCMatrix::castdown(aoso.pointer());

  for (b=0; b < lu->nblocks(); b++) {
    if (lu->block(b).null())
      continue;
    
    int ir = ct.which_irrep(b);
  
    double skal = (double)ct.order()/(double)ct.gamma(ir).degeneracy();
    skal = sqrt(skal);
    lu->block(b).scale(skal);
  }

  sym.assign(0.0);
  sym.accumulate_transform(aoso,bskel);
  aoso=0;
#else
  do_transform(skel,sym,gbs,*this,ct);
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

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
