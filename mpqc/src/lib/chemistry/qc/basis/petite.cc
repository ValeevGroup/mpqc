
#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/molecule/localdef.h>

#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/shellrot.h>

////////////////////////////////////////////////////////////////////////////

PetiteList::PetiteList(const RefGaussianBasisSet &gbs,
                       const RefIntegral& ints) :
  gbs_(gbs),
  ints_(ints)
{
  init();
}

PetiteList::~PetiteList()
{
  if (p1_)
    delete[] p1_;

  if (lamij_)
    delete[] lamij_;

  if (nbf_in_ir_)
    delete[] nbf_in_ir_;
  
  if (atom_map_) {
    for (int i=0; i < natom_; i++)
      delete[] atom_map_[i];
    delete[] atom_map_;
  }

  if (shell_map_) {
    for (int i=0; i < nshell_; i++)
      delete[] shell_map_[i];
    delete[] shell_map_;
  }

  natom_=0;
  nshell_=0;
  ng_=0;
  nblocks_=0;
  nirrep_=0;
  p1_=0;
  atom_map_=0;
  shell_map_=0;
  lamij_=0;
  nbf_in_ir_=0;
}

void
PetiteList::init()
{
  int i;

  // grab references to the Molecule and BasisSet for convenience
  GaussianBasisSet& gbs = *gbs_.pointer();
  Molecule& mol = *gbs.molecule().pointer();

  // create the character table for the point group
  CharacterTable ct = mol.point_group().char_table();
  
  // initialize private members
  ng_ = ct.order();
  natom_ = mol.natom();
  nshell_ = gbs.nshell();
  nirrep_ = ct.nirrep();

  // allocate storage for arrays
  p1_ = new char[nshell_];
  lamij_ = new char[i_offset(nshell_)];

  atom_map_ = new int*[natom_];
  for (i=0; i < natom_; i++)
    atom_map_[i] = new int[ng_];
  
  shell_map_ = new int*[nshell_];
  for (i=0; i < nshell_; i++)
    shell_map_[i] = new int[ng_];
  
  // set up atom and shell mappings
  double np[3];
  SymmetryOperation so;
  
  // loop over all centers
  for (i=0; i < natom_; i++) {
    AtomicCenter ac = mol.atom(i);

    // then for each symop in the pointgroup, transform the coordinates of
    // center "i" and see which atom it maps into
    for (int g=0; g < ng_; g++) {
      so = ct.symm_operation(g);

      for (int ii=0; ii < 3; ii++) {
        np[ii] = 0;
        for (int jj=0; jj < 3; jj++)
          np[ii] += so(ii,jj) * ac[jj];
      }

      atom_map_[i][g] = mol.atom_at_position(np, 0.05);
    }

    // hopefully, shells on equivalent centers will be numbered in the same
    // order
    for (int s=0; s < gbs.nshell_on_center(i); s++) {
      int shellnum = gbs.shell_on_center(i,s);
      for (int g=0; g < ng_; g++) {
        shell_map_[shellnum][g] = gbs.shell_on_center(atom_map_[i][g],s);
      }
    }
  }

  memset(p1_,0,nshell_);
  memset(lamij_,0,i_offset(nshell_));
  
  // now we do p1_ and lamij_
  for (i=0; i < nshell_; i++) {
    int g;

    // we want the highest numbered shell in a group of equivalent shells
    for (g=0; g < ng_; g++)
      if (shell_map_[i][g] > i)
        break;
    
    if (g < ng_)
      continue;
    
    // i is in the group P1
    p1_[i] = 1;

    for (int j=0; j <= i; j++) {
      int ij = i_offset(i)+j;
      int nij = 0;

      // test to see if IJ is in the group P2, if it is, then set lambda(ij)
      // equal to the number of equivalent shell pairs.  This number is
      // just the order of the group divided by the number of times ij is
      // mapped into itself
      int gg;
      for (gg=0; gg < ng_; gg++) {
        int gi = shell_map_[i][gg];
        int gj = shell_map_[j][gg];
        int gij = ij_offset(gi,gj);
        if (gij > ij)
          break;
        else if (gij == ij)
          nij++;
      }

      if (gg < ng_)
        continue;

      lamij_[ij] = (char) (ng_/nij);
    }
  }

  // form reducible representation of the basis functions
  double *red_rep = new double[ng_];
  memset(red_rep,0,sizeof(double)*ng_);
  
  for (i=0; i < natom_; i++) {
    for (int g=0; g < ng_; g++) {
      so = ct.symm_operation(g);
      int j= atom_map_[i][g];

      if (i!=j)
        continue;
      
      for (int s=0; s < gbs.nshell_on_center(i); s++) {
        for (int c=0; c < gbs(i,s).ncontraction(); c++) {
          int am=gbs(i,s).am(c);

          if (am==0)
            red_rep[g] += 1.0;
          else {
            ShellRotation r(am,so,ints_,gbs(i,s).is_pure(c));
            red_rep[g] += r.trace();
          }
        }
      }
    }
  }

  // and then use projection operators to figure out how many SO's of each
  // symmetry type there will be
  nblocks_ = 0;
  nbf_in_ir_ = new int[nirrep_];
  for (i=0; i < nirrep_; i++) {
    double t=0;
    for (int g=0; g < ng_; g++)
      t += ct.gamma(i).character(g)*red_rep[g];

    nbf_in_ir_[i] = ((int) (t+0.5))/ng_;
    if (ct.gamma(i).complex()) {
      nblocks_++;
      nbf_in_ir_[i] *= 2;
    } else {
      nblocks_ += ct.gamma(i).degeneracy();
    }
  }

  delete[] red_rep;
}

RefSCDimension
PetiteList::AO_basisdim()
{
  RefSCDimension dim = new SCDimension(gbs_->nbasis(),1);
  dim->blocks()->set_subdim(0, gbs_->basisdim());
  return dim;
}

RefSCDimension
PetiteList::SO_basisdim()
{
  int i, j, ii;
  
  // grab a reference to the basis set
  GaussianBasisSet& gbs = *gbs_.pointer();
  
  // create the character table for the point group
  CharacterTable ct = gbs.molecule()->point_group().char_table();

  // ncomp is the number of symmetry blocks we have
  int ncomp=nblocks();
  
  // saoelem is the current SO in a block
  int *nao = new int [ncomp];
  memset(nao,0,sizeof(int)*ncomp);

  for (i=ii=0; i < nirrep_; i++) {
    int je = ct.gamma(i).complex() ? 1 : ct.gamma(i).degeneracy();
    for (j=0; j < je; j++,ii++)
      nao[ii] = nbf_in_ir_[i];
  }

  RefSCDimension ret = new SCDimension(gbs.nbasis(),ncomp,nao);
  delete[] nao;

  for (i=ii=0; i < nirrep_; i++) {
    int je = ct.gamma(i).complex() ? 1 : ct.gamma(i).degeneracy();
    for (j=0; j < je; j++,ii++) {
      char lab[24];
      sprintf(lab,"irrep %s comp %d", ct.gamma(i).symbol(), j);
      ret->blocks()->set_subdim(ii, new SCDimension(nbf_in_ir_[i]));
    }
  }

  return ret;
}

void
PetiteList::print(FILE *o, int verbose)
{
  int i;

  fprintf(o,"PetiteList:\n");

  if (verbose) {
    fprintf(o,"  natom_ = %d\n",natom_);
    fprintf(o,"  nshell_ = %d\n",nshell_);
    fprintf(o,"  ng_ = %d\n",ng_);
    fprintf(o,"  nirrep_ = %d\n",nirrep_);

    fprintf(o,"\n");
    fprintf(o,"  atom_map_ = \n");

    for (i=0; i < natom_; i++) {
      fprintf(o,"    ");
      for (int g=0; g < ng_; g++)
        fprintf(o,"%5d ",atom_map_[i][g]);
      fprintf(o,"\n");
    }

    fprintf(o,"\n");
    fprintf(o,"  shell_map_ = \n");
    for (i=0; i < nshell_; i++) {
      fprintf(o,"    ");
      for (int g=0; g < ng_; g++)
      fprintf(o,"%5d ",shell_map_[i][g]);
      fprintf(o,"\n");
    }

    fprintf(o,"\n");
    fprintf(o,"  p1_ = \n");
    for (i=0; i < nshell_; i++)
      fprintf(o,"    %5d\n",p1_[i]);

    fprintf(o,"  lamij_ = \n");
    for (i=0; i < nshell_; i++) {
      fprintf(o,"    ");
      for (int j=0; j <= i; j++)
        fprintf(o,"%5d ",lamij_[i_offset(i)+j]);
      fprintf(o,"\n");
    }
    fprintf(o,"\n");
  }

  CharacterTable ct = gbs_->molecule()->point_group().char_table();
  for (i=0; i < nirrep_; i++)
    fprintf(o,"  %5d functions of %s symmetry\n",
            nbf_in_ir_[i], ct.gamma(i).symbol());
}

// forms the basis function rotation matrix for the g'th symmetry operation
// in the point group
RefSCMatrix
PetiteList::r(int g)
{
  // grab a reference to the basis set
  GaussianBasisSet& gbs = *gbs_.pointer();
  
  SymmetryOperation so =
    gbs.molecule()->point_group().char_table().symm_operation(g);

  RefSCMatrix ret = gbs.matrixkit()->matrix(gbs.basisdim(), gbs.basisdim());
  ret.assign(0.0);
  
  // this should be replaced with an element op at some point
  for (int i=0; i < natom_; i++) {
    int j = atom_map_[i][g];

    for (int s=0; s < gbs.nshell_on_center(i); s++) {
      int func_i = gbs.shell_to_function(gbs.shell_on_center(i,s));
      int func_j = gbs.shell_to_function(gbs.shell_on_center(j,s));
      
      for (int c=0; c < gbs(i,s).ncontraction(); c++) {
        int am=gbs(i,s).am(c);

        if (am==0) {
          ret.set_element(func_j,func_i,1.0);
        } else {
          ShellRotation rr(am,so,ints_,gbs(i,s).is_pure(c));
          for (int ii=0; ii < rr.dim(); ii++)
            for (int jj=0; jj < rr.dim(); jj++)
              ret.set_element(func_j+jj,func_i+ii,rr(ii,jj));
        }

        func_i += gbs(i,s).nfunction(c);
        func_j += gbs(i,s).nfunction(c);
      }
    }
  }
  return ret;
}
