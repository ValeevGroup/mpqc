//
// aotoso.cc --- more symmetry stuff
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <float.h>

#include <math/scmat/lapack.h>
#include <util/misc/formio.h>

#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/shellrot.h>
#include <chemistry/qc/basis/petite.h>

using namespace std;
using namespace sc;

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

SO::SO() : len(0), length(0), cont(0)
{
}

SO::SO(int l) : len(0), length(0), cont(0)
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
  set_length(so.length);
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

SO_block::SO_block() : len(0), so(0)
{
}

SO_block::SO_block(int l) : len(0), so(0)
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
  if (len == l) return;

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
  ExEnv::out0() << indent << "SO block " << title << endl;
  for (i=0; i < len; i++) {
    ExEnv::out0() << indent << "SO " << i+1 << endl << indent;
    for (j=0; j < so[i].length; j++)
      ExEnv::out0() << scprintf(" %10d",so[i].cont[j].bfn);
    ExEnv::out0() << endl << indent;
    for (j=0; j < so[i].length; j++)
      ExEnv::out0() << scprintf(" %10.7f",so[i].cont[j].coef);
    ExEnv::out0() << endl;
  }
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
      ExEnv::out0() << indent;
      for (i=0; i < ns; i++)
        ExEnv::out0() << scprintf(" %10d",mapf0+i);
      ExEnv::out0() << endl;
      
      for (i=0; i < ns; i++) {
        ExEnv::out0() << indent << scprintf("%2d",f0+i);
        for (int j=0; j < ns; j++)
          ExEnv::out0() << scprintf(" %10.7f",c[i][j]);
        ExEnv::out0() << endl;
      }
    }
};

////////////////////////////////////////////////////////////////////////////

SO_block *
PetiteList::aotoso_info()
{
  int iuniq, i, j, d, ii, jj, g, s, c, ir;

  GaussianBasisSet& gbs = *gbs_.pointer();
  Molecule& mol = *gbs.molecule().pointer();

  // create the character table for the point group
  CharacterTable ct = mol.point_group()->char_table();
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
  
  // loop over all unique shells
  for (iuniq=0; iuniq < gbs.molecule()->nunique(); iuniq++) {
    int nequiv = gbs.molecule()->nequivalent(iuniq);
    i = gbs.molecule()->unique(iuniq);
    for (s=0; s < gbs.nshell_on_center(i); s++) {
      int shell_i = gbs.shell_on_center(i,s);
      
      // test to see if there are any high am cartesian functions in this
      // shell.  for now don't allow symmetry with cartesian functions...I
      // just can't seem to get them working.
      for (c=0; c < gbs(i,s).ncontraction(); c++) {
        if (gbs(i,s).am(c) > 1 && gbs(i,s).is_cartesian(c)) {
          if (ng_ != nirrep_) {
            ExEnv::err0() << indent
                         << "PetiteList::aotoso: cannot yet handle"
                         << " symmetry for angular momentum >= 2\n";
            abort();
          }
        }
      }

      // the functions do not mix between contractions
      // so the contraction loop can be done outside the symmetry
      // operation loop
      int bfn_offset_in_shell = 0;
      for (c=0; c < gbs(i,s).ncontraction(); c++) {
        const blasint nfuncuniq = gbs(i,s).nfunction(c);
        const blasint nfuncall = nfuncuniq * nequiv;

        // allocate an array to store linear combinations of orbitals
        // this is destroyed by the SVD routine
        double **linorb = new double*[nfuncuniq];
        linorb[0] = new double[nfuncuniq*nfuncall];
        for (j=1; j<nfuncuniq; j++) {
          linorb[j] = &linorb[j-1][nfuncall];
        }

        // a copy of linorb
        double **linorbcop = new double*[nfuncuniq];
        linorbcop[0] = new double[nfuncuniq*nfuncall];
        for (j=1; j<nfuncuniq; j++) {
          linorbcop[j] = &linorbcop[j-1][nfuncall];
        }

        // allocate an array for the SVD routine
        double **u = new double*[nfuncuniq];
        u[0] = new double[nfuncuniq*nfuncuniq];
        for (j=1; j<nfuncuniq; j++) {
          u[j] = &u[j-1][nfuncuniq];
        }

        // loop through each irrep to form the linear combination
        // of orbitals of that symmetry
        int irnum = 0;
        for (ir=0; ir<ct.nirrep(); ir++) {
          int cmplx = (ct.complex() && ct.gamma(ir).complex());
          for (int comp=0; comp < ct.gamma(ir).degeneracy(); comp++) {
            memset(linorb[0], 0, nfuncuniq*nfuncall*sizeof(double));
            for (d=0; d < ct.gamma(ir).degeneracy(); d++) {
              // if this is a point group with a complex E representation,
              // then only do the first set of projections for E
              if (d && cmplx) break;

              // operate on each function in this contraction with each
              // symmetry operation to form symmetrized linear combinations
              // of orbitals

              for (g=0; g<ng_; g++) {
                // the character
                double ccdg = ct.gamma(ir).p(comp,d,g);

                so = ct.symm_operation(g);
                int equivatom = atom_map_[i][g];
                int atomoffset
                  = gbs.molecule()->atom_to_unique_offset(equivatom);
        
                ShellRotation rr
                  = ints_->shell_rotation(gbs(i,s).am(c),
                                          so,gbs(i,s).is_pure(c));

                for (ii=0; ii < rr.dim(); ii++) {
                  for (jj=0; jj < rr.dim(); jj++) {
                    linorb[ii][atomoffset*nfuncuniq+jj] += ccdg * rr(ii,jj);
                  }
                }
              }
            }
            // find the linearly independent SO's for this irrep/component
            memcpy(linorbcop[0], linorb[0], nfuncuniq*nfuncall*sizeof(double));
            double *singval = new double[nfuncuniq];
            double djunk=0; blasint ijunk=1;
            blasint lwork = 5*nfuncall;
            double *work = new double[lwork];
            blasint info;
            // solves At = V SIGMA Ut (since FORTRAN array layout is used)
            F77_DGESVD("N","A",&nfuncall,&nfuncuniq,linorb[0],&nfuncall,
                       singval, &djunk, &ijunk, u[0], &nfuncuniq,
                       work, &lwork, &info);
            // put the lin indep symm orbs into the so array
            for (j=0; j<nfuncuniq; j++) {
              if (singval[j] > 1.0e-6) {
                SO tso;
                tso.set_length(nfuncall);
                int ll = 0, llnonzero = 0;
                for (int k=0; k<nequiv; k++) {
                  for (int l=0; l<nfuncuniq; l++,ll++) {
                    double tmp = 0.0;
                    for (int m=0; m<nfuncuniq; m++) {
                      tmp += u[m][j] * linorbcop[m][ll] / singval[j];
                    }
                    if (fabs(tmp) > DBL_EPSILON) {
                      int equivatom = gbs.molecule()->equivalent(iuniq,k);
                      tso.cont[llnonzero].bfn
                        = l
                        + bfn_offset_in_shell
                        + gbs.shell_to_function(gbs.shell_on_center(equivatom,
                                                                    s));
                      tso.cont[llnonzero].coef = tmp;
                      llnonzero++;
                    }
                  }
                }
                tso.reset_length(llnonzero);
                if (llnonzero == 0) {
                  ExEnv::errn() << "aotoso: internal error: no bfns in SO"
                               << endl;
                  abort();
                }
                if (SOs[irnum+comp].add(tso,saoelem[irnum+comp])) {
                  saoelem[irnum+comp]++;
                }
                else {
                  ExEnv::errn() << "aotoso: internal error: "
                               << "impossible duplicate SO"
                               << endl;
                  abort();
                }
              }
            }
            delete[] singval;
            delete[] work;
          }
          irnum += ct.gamma(ir).degeneracy();
        }
        bfn_offset_in_shell += gbs(i,s).nfunction(c);

        delete[] linorb[0];
        delete[] linorb;
        delete[] linorbcop[0];
        delete[] linorbcop;
        delete[] u[0];
        delete[] u;
      }
    }
  }

  // Make sure all the nodes agree on what the symmetry orbitals are.
  // (All the above work for me > 0 is ignored.)
  Ref<MessageGrp> grp = MessageGrp::get_default_messagegrp();
  for (i=0; i<ncomp; i++) {
    int len = SOs[i].len;
    grp->bcast(len);
    SOs[i].reset_length(len);
    for (j=0; j<len; j++) {
      int solen = SOs[i].so[j].length;
      grp->bcast(solen);
      SOs[i].so[j].reset_length(solen);
      for (int k=0; k<solen; k++) {
        grp->bcast(SOs[i].so[j].cont[k].bfn);
        grp->bcast(SOs[i].so[j].cont[k].coef);
      }
    }
  }

  for (i=0; i < ncomp; i++) {
    ir = whichir[i];
    int scal = ct.gamma(ir).complex() ? 2 : 1;

    if (saoelem[i] < nbf_in_ir_[ir]/scal) {
      // if we found too few, there are big problems
      
      ExEnv::err0() << indent
           << scprintf("trouble making SO's for irrep %s\n",
                       ct.gamma(ir).symbol());
      ExEnv::err0() << indent
           << scprintf("  only found %d out of %d SO's\n",
                       saoelem[i], nbf_in_ir_[ir]/scal);
      SOs[i].print("");
      abort();

    } else if (saoelem[i] > nbf_in_ir_[ir]/scal) {
      // there are some redundant so's left...need to do something to get
      // the elements we want
      
      ExEnv::err0() << indent
           << scprintf("trouble making SO's for irrep %s\n",
                       ct.gamma(ir).symbol());
      ExEnv::err0() << indent
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

        for (j=0; j < SOs[i].len; j++,k++)
          nSOs[in].add(SOs[i].so[j],k);

        i++;
        in++;
      } else {
        for (j=0; j < ct.gamma(ir).degeneracy(); j++,i++,in++) {
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
  
  BlockedSCMatrix *aosop = dynamic_cast<BlockedSCMatrix*>(aoso.pointer());

  for (int b=0; b < aosop->nblocks(); b++) {
    RefSCMatrix aosb = aosop->block(b);

    if (aosb == 0)
      continue;
    
    SO_block& sob = sos[b];
    
    Ref<SCMatrixSubblockIter> iter =
      aosb->local_blocks(SCMatrixSubblockIter::Write);

    for (iter->begin(); iter->ready(); iter->next()) {
      if (dynamic_cast<SCMatrixRectBlock*>(iter->block())) {
        SCMatrixRectBlock *blk = dynamic_cast<SCMatrixRectBlock*>(iter->block());

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
          dynamic_cast<SCMatrixRectSubBlock*>(iter->block());

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
  
  delete[] sos;

  return aoso;
}

RefSCMatrix
PetiteList::sotoao()
{
  if (c1_)
    return aotoso();
  else if (nirrep_ == ng_) // subgroup of d2h
    return aotoso().t();
  else
    return aotoso().i();
}

RefSymmSCMatrix
PetiteList::to_SO_basis(const RefSymmSCMatrix& a)
{
  // SO basis is always blocked, so first make sure a is blocked
  RefSymmSCMatrix aomatrix = dynamic_cast<BlockedSymmSCMatrix*>(a.pointer());
  if (aomatrix == 0) {
    aomatrix = gbs_->so_matrixkit()->symmmatrix(AO_basisdim());
    aomatrix->convert(a);
  }

  // if C1, then do nothing
  if (c1_)
    return aomatrix.copy();
  
  RefSymmSCMatrix somatrix(SO_basisdim(), gbs_->so_matrixkit());
  somatrix.assign(0.0);
  somatrix->accumulate_transform(aotoso(), aomatrix,
                                 SCMatrix::TransposeTransform);

  return somatrix;
}

RefSymmSCMatrix
PetiteList::to_AO_basis(const RefSymmSCMatrix& somatrix)
{
  // if C1, then do nothing
  if (c1_)
    return somatrix.copy();
  
  RefSymmSCMatrix aomatrix(AO_basisdim(), gbs_->so_matrixkit());
  aomatrix.assign(0.0);

  if (nirrep_ == ng_) // subgroup of d2h
    aomatrix->accumulate_transform(aotoso(), somatrix);
  else
    aomatrix->accumulate_transform(aotoso().i(), somatrix,
                                   SCMatrix::TransposeTransform);

  return aomatrix;
}

RefSCMatrix
PetiteList::evecs_to_SO_basis(const RefSCMatrix& aoev)
{
  // if C1, then do nothing
  if (c1_)
    return aoev.copy();

  // given aoev may be non-blocked or its dimensions may have different sub-blocking
  // in that case copy into a blocked matrix of the desired dimensions
  RefSCMatrix aoevecs = dynamic_cast<BlockedSCMatrix*>(aoev.pointer());
  const bool need_to_copy = aoevecs == 0 || (aoevecs && !AO_basisdim()->equiv(aoev.rowdim()));
  if (need_to_copy) {
    aoevecs = gbs_->so_matrixkit()->matrix(AO_basisdim(), aoev.coldim());
    aoevecs->convert(aoev);
  }

  RefSCMatrix soev =  sotoao() * aoevecs;

  // redimension
  RefSCMatrix soevecs(SO_basisdim(), SO_basisdim(), gbs_->so_matrixkit());
  // this convert messes it all up!
  // because Blocked ops are not implemented!
  // copy manually
  // TODO replace with convert when it's fixed to work with Blocked matrices properly
#define SCMATRIX_CONVERT_BROKEN 1
#if SCMATRIX_CONVERT_BROKEN
  const int nrow = soevecs.nrow();
  const int ncol = soevecs.ncol();
  for(int r=0; r<nrow; ++r) {
    for(int c=0; c<ncol; ++c) {
      soevecs.set_element(r, c, soev.get_element(r, c));
    }
  }
#else
  soevecs->convert(soev);
#endif

#if 0
  aoev.print("PetiteList::evecs_to_SO_basis: input C_ao");
  aoevecs.print("PetiteList::evecs_to_SO_basis: input C_ao (blocked)");
  soev.print("PetiteList::evecs_to_SO_basis: output C_so");
  soevecs.print("PetiteList::evecs_to_SO_basis: output C_so (blocked)");
#endif

  return soevecs;
}

RefSCMatrix
PetiteList::evecs_to_AO_basis(const RefSCMatrix& soevecs)
{
  // if C1, then do nothing
  if (c1_)
    return soevecs.copy();
  
  RefSCMatrix aoev = aotoso() * soevecs;

  return aoev;
}

/////////////////////////////////////////////////////////////////////////////

void
PetiteList::symmetrize(const RefSymmSCMatrix& skel,
                       const RefSymmSCMatrix& sym)
{
  GaussianBasisSet& gbs = *gbs_.pointer();

  // SO basis is always blocked, so first make sure skel is blocked
  RefSymmSCMatrix bskel = dynamic_cast<BlockedSymmSCMatrix*>(skel.pointer());
  if (bskel == 0) {
    bskel = gbs.so_matrixkit()->symmmatrix(AO_basisdim());
    bskel->convert(skel);
  }
  
  // if C1, then do nothing
  if (c1_) {
    sym.assign(bskel);
    return;
  }
  
  int b,c;

  CharacterTable ct = gbs.molecule()->point_group()->char_table();

  RefSCMatrix aoso = aotoso();
  BlockedSCMatrix *lu = dynamic_cast<BlockedSCMatrix*>(aoso.pointer());

  for (b=0; b < lu->nblocks(); b++) {
    if (lu->block(b) == 0)
      continue;
    
    int ir = ct.which_irrep(b);
  
    double skal = (double)ct.order()/(double)ct.gamma(ir).degeneracy();
    skal = sqrt(skal);
    lu->block(b).scale(skal);
  }

  sym.assign(0.0);
  sym.accumulate_transform(aoso,bskel,SCMatrix::TransposeTransform);
  aoso=0;

  BlockedSymmSCMatrix *la = dynamic_cast<BlockedSymmSCMatrix*>(sym.pointer());
  
  // loop through blocks and finish symmetrizing degenerate blocks
  for (b=0; b < la->nblocks(); b++) {
    if (la->block(b) == 0)
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
// c-file-style: "ETS"
// End:
