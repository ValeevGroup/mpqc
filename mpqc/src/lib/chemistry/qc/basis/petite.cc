
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/scmat/local.h>
#include <chemistry/molecule/localdef.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/rot.h>

////////////////////////////////////////////////////////////////////////////

PetiteList::PetiteList()
{
  _gbs=0;
  _natom=0;
  _nshell=0;
  _nirrep=0;
  _ng=0;
  _p1=0;
  _atom_map=0;
  _shell_map=0;
  _lamij=0;
  _nbf_in_ir=0;
}

PetiteList::PetiteList(const RefGaussianBasisSet &gbs)
{
  init(gbs);
}

PetiteList::~PetiteList()
{
  if (_p1)
    delete[] _p1;

  if (_lamij)
    delete[] _lamij;

  if (_nbf_in_ir)
    delete[] _nbf_in_ir;
  
  if (_atom_map) {
    for (int i=0; i < _natom; i++)
      delete[] _atom_map[i];
    delete[] _atom_map;
  }

  if (_shell_map) {
    for (int i=0; i < _nshell; i++)
      delete[] _shell_map[i];
    delete[] _shell_map;
  }

  _gbs=0;
  _natom=0;
  _nshell=0;
  _ng=0;
  _nirrep=0;
  _p1=0;
  _atom_map=0;
  _shell_map=0;
  _lamij=0;
  _nbf_in_ir=0;
}

static int
atom_num(Point& p, Molecule& mol)
{
  for (int i=0; i < mol.natom(); i++) {
    if (dist(p,mol.atom(i).point()) < 0.05)
      return i;
  }
  return -1;
}

void
PetiteList::init(const RefGaussianBasisSet &gb)
{
  int i;

  _gbs = gb;
    
  // grab references to the Molecule and BasisSet for convenience
  GaussianBasisSet& gbs = *_gbs.pointer();
  Molecule& mol = *gbs.molecule().pointer();

  // create the character table for the point group
  CharacterTable ct = mol.point_group().char_table();
  
  // initialize private members
  _ng = ct.order();
  _natom = mol.natom();
  _nshell = gbs.nshell();
  _nirrep = ct.nirrep();

  // allocate storage for arrays
  _p1 = new char[_nshell];
  _lamij = new char[ioff(_nshell)];

  _atom_map = new int*[_natom];
  for (i=0; i < _natom; i++)
    _atom_map[i] = new int[_ng];
  
  _shell_map = new int*[_nshell];
  for (i=0; i < _nshell; i++)
    _shell_map[i] = new int[_ng];
  
  // set up atom and shell mappings
  Point np;
  SymmetryOperation so;
  
  // loop over all centers
  for (i=0; i < _natom; i++) {
    AtomicCenter ac = mol.atom(i);

    // then for each symop in the pointgroup, transform the coordinates of
    // center "i" and see which atom it maps into
    for (int g=0; g < _ng; g++) {
      so = ct.symm_operation(g);

      for (int ii=0; ii < 3; ii++) {
        np[ii] = 0;
        for (int jj=0; jj < 3; jj++)
          np[ii] += so(ii,jj) * ac[jj];
      }

      _atom_map[i][g] = atom_num(np,mol);
    }

    // hopefully, shells on equivalent centers will be numbered in the same
    // order
    for (int s=0; s < gbs.nshell_on_center(i); s++) {
      int shellnum = gbs.shell_on_center(i,s);
      for (int g=0; g < _ng; g++) {
        _shell_map[shellnum][g] = gbs.shell_on_center(_atom_map[i][g],s);
      }
    }
  }

  memset(_p1,0,_nshell);
  memset(_lamij,0,ioff(_nshell));
  
  // now we do _p1 and _lamij
  for (i=0; i < _nshell; i++) {
    int g;

    // we want the highest numbered shell in a group of equivalent shells
    for (g=0; g < _ng; g++)
      if (_shell_map[i][g] > i)
        break;
    
    if (g < _ng)
      continue;
    
    // i is in the group P1
    _p1[i] = 1;

    for (int j=0; j <= i; j++) {
      int ij = ioff(i)+j;
      int nij = 0;

      // test to see if IJ is in the group P2, if it is, then set lambda(ij)
      // equal to the number of equivalent shell pairs.  This number is
      // just the order of the group divided by the number of times ij is
      // mapped into itself
      int gg;
      for (gg=0; gg < _ng; gg++) {
        int gi = _shell_map[i][gg];
        int gj = _shell_map[j][gg];
        int gij = ioff(gi,gj);
        if (gij > ij)
          break;
        else if (gij == ij)
          nij++;
      }

      if (gg < _ng)
        continue;

      _lamij[ij] = (char) (_ng/nij);
    }
  }

  // form reducible representation of the basis functions
  double *red_rep = new double[_ng];
  memset(red_rep,0,sizeof(double)*_ng);
  
  for (i=0; i < _natom; i++) {
    for (int g=0; g < _ng; g++) {
      so = ct.symm_operation(g);
      int j= _atom_map[i][g];

      if (i!=j)
        continue;
      
      for (int s=0; s < gbs.nshell_on_center(i); s++) {
        for (int c=0; c < gbs(i,s).ncontraction(); c++) {
          int am=gbs(i,s).am(c);

          if (am==0)
            red_rep[g] += 1.0;
          else {
            Rotation r(am,so,gbs(i,s).is_pure(c));
            red_rep[g] += r.trace();
          }
        }
      }
    }
  }

  // and then use projection operators to figure out how many SO's of each
  // symmetry type there will be
  _nbf_in_ir = new int[_nirrep];
  for (i=0; i < _nirrep; i++) {
    double t=0;
    for (int g=0; g < _ng; g++)
      t += ct[i][g]*red_rep[g];

    _nbf_in_ir[i] = ((int) (t+0.5))/_ng;
  }

  delete[] red_rep;
}

// forms the basis function rotation matrix for the g'th symmetry operation
// in the point group
RefSCMatrix
PetiteList::r(int g)
{
  SymmetryOperation so =
    _gbs->molecule()->point_group().char_table().symm_operation(g);
  GaussianBasisSet& gbs = *_gbs.pointer();

  RefSCMatrix ret = gbs.basisdim()->create_matrix(gbs.basisdim());
  ret.assign(0.0);
  
  // this should be replaced with an element op at some point
  for (int i=0; i < _natom; i++) {
    int j = _atom_map[i][g];

    for (int s=0; s < gbs.nshell_on_center(i); s++) {
      int func_i = gbs.shell_to_function(gbs.shell_on_center(i,s));
      int func_j = gbs.shell_to_function(gbs.shell_on_center(j,s));
      
      for (int c=0; c < gbs(i,s).ncontraction(); c++) {
        int am=gbs(i,s).am(c);

        if (am==0) {
          ret.set_element(func_j,func_i,1.0);
        } else {
          Rotation rr(am,so,gbs(i,s).is_pure(c));
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

void
PetiteList::print(FILE *o)
{
  int i;

  fprintf(o,"PetiteList:\n");
  fprintf(o,"  _natom = %d\n",_natom);
  fprintf(o,"  _nshell = %d\n",_nshell);
  fprintf(o,"  _ng = %d\n",_ng);
  fprintf(o,"  _nirrep = %d\n",_nirrep);

  fprintf(o,"\n");
  fprintf(o,"  _atom_map = \n");

  for (i=0; i < _natom; i++) {
    fprintf(o,"    ");
    for (int g=0; g < _ng; g++)
      fprintf(o,"%5d ",_atom_map[i][g]);
    fprintf(o,"\n");
  }

  fprintf(o,"\n");
  fprintf(o,"  _shell_map = \n");
  for (i=0; i < _nshell; i++) {
    fprintf(o,"    ");
    for (int g=0; g < _ng; g++)
      fprintf(o,"%5d ",_shell_map[i][g]);
    fprintf(o,"\n");
  }

  fprintf(o,"\n");
  fprintf(o,"  _p1 = \n");
  for (i=0; i < _nshell; i++)
    fprintf(o,"    %5d\n",_p1[i]);
    
  fprintf(o,"  _lamij = \n");
  for (i=0; i < _nshell; i++) {
    fprintf(o,"    ");
    for (int j=0; j <= i; j++)
      fprintf(o,"%5d ",_lamij[ioff(i)+j]);
    fprintf(o,"\n");
  }
  
  fprintf(o,"\n");
  CharacterTable ct = _gbs->molecule()->point_group().char_table();
  for (i=0; i < _nirrep; i++)
    fprintf(o,"  %5d functions of %s symmetry\n",_nbf_in_ir[i],ct[i].symbol());
}

////////////////////////////////////////////////////////////////////////////

class lin_comb {
  private:
    int _nsh;
    int _nbf;
    int *fn;
    double **c;
  public:
    lin_comb(int nsh, int nbf, int f0) : _nsh(nsh), _nbf(nbf) {
      int i;
      fn = new int[_nsh];
      for (i=0; i < _nsh; i++)
        fn[i] = f0+i;
      
      c = new double*[_nsh];
      for (i=0; i < _nsh; i++) {
        c[i] = new double[_nbf];
        memset(c[i],0,sizeof(double)*_nbf);
      }
    }
    ~lin_comb() {
      if (fn) delete[] fn; fn=0;
      if (c) {
        for (int i=0; i < _nsh; i++)
          if (c[i])
            delete[] c[i];
        delete[] c;
      }
      c=0; _nsh=_nbf=0;
    }

    int numbf() const { return _nbf; }
    int numsh() const { return _nsh; }
    int bfnum(int i) const { return fn[i]; }
    double& coef(int i, int j) { return c[i][j]; }
    void print() const {
      for (int i=0; i < _nsh; i++) {
        printf("%2d",fn[i]);
        for (int j=0; j < _nbf; j++)
          printf(" %10.7f",c[i][j]);
        printf("\n");
      }
    }
};
    
RefSCMatrix
PetiteList::aotoso()
{
  int i, d, ii, jj, g, fn, s, c, ir, f;

  GaussianBasisSet& gbs = *_gbs.pointer();
  Molecule& mol = *gbs.molecule().pointer();

  RefSCMatrix ret = gbs.basisdim()->create_matrix(gbs.basisdim());
  ret.assign(0.0);
  
  // create the character table for the point group
  CharacterTable ct = mol.point_group().char_table();
  SymmetryOperation so;

  // ncomp is the number of symmetry blocks we have
  int ncomp=0;
  for (i=0; i < _nirrep; i++)
    ncomp += ct[i].degeneracy();
  
  // fao is the first SO in a block, saoelem is the current SO in a block,
  // and maxsao is the last SO in a block
  int *fao = new int[ncomp];
  int *saoelem = new int[ncomp];
  int *maxsao = new int[ncomp];

  saoelem[0]=0;
  maxsao[0]=_nbf_in_ir[0];
  fao[0]=0;
  
  ii=1;
  for (i=1; i < _nirrep; i++) {
    for (d=0; d < ct[i].degeneracy(); d++,ii++) {
      maxsao[ii] = maxsao[ii-1]+_nbf_in_ir[i];
      saoelem[ii] = maxsao[ii-1];
      fao[ii] = saoelem[ii];
    }
  }
    
  double *blc = new double[gbs.nbasis()];

  // loop over all unique shells
  for (i=0; i < _natom; i++) {
    for (s=0; s < gbs.nshell_on_center(i); s++) {
      int shell_i = gbs.shell_on_center(i,s);
      
      if (!_p1[shell_i])
        continue;
      
      // form linear combinations of the basis functions in this shell
      lin_comb **lc = new lin_comb*[_ng];
      
      for (g=0; g < _ng; g++) {
        so = ct.symm_operation(g);
        int j = _atom_map[i][g];

        int func_i = gbs.shell_to_function(gbs.shell_on_center(i,s));
        int func_j = gbs.shell_to_function(gbs.shell_on_center(j,s));

        lc[g] = new lin_comb(gbs(i,s).nfunction(),gbs.nbasis(),func_i);
        lin_comb& lcg = *lc[g];
        
        int fi=0;
        for (c=0; c < gbs(i,s).ncontraction(); c++) {
          int am=gbs(i,s).am(c);

          if (am==0) {
            lcg.coef(fi,func_j) = 1.0;
          } else {
            Rotation rr(am,so,gbs(i,s).is_pure(c));
            for (ii=0; ii < rr.dim(); ii++)
              for (jj=0; jj < rr.dim(); jj++)
                lcg.coef(fi+ii,func_j+jj) = rr(ii,jj);
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
        for (fn=0; fn < gbs(i,s).nfunction(); fn++) {
          for (d=0; d < ct[ir].nproj(); d++) {

            int comp = d % ct[ir].degeneracy();

            // form the projection for this irrep
            memset(blc,0,sizeof(double)*gbs.nbasis());
            for (g=0; g < _ng; g++) {
              lin_comb& lcg = *lc[g];

              for (f=0; f < gbs.nbasis(); f++)
                blc[f] += ct[ir](d,g)*lcg.coef(fn,f);
            }

            // normalize the linear combination if c1 is non-zero
            double c1=0;
            for (ii=0; ii < gbs.nbasis(); ii++)
              c1 += blc[ii]*blc[ii];

            if (c1 > 1.0e-3) {
              c1 = 1.0/sqrt(c1);
              for (ii=0; ii < gbs.nbasis(); ii++)
                blc[ii] *= c1;
            } else {
              continue;
            }

            // check to see if this SO has already been generated
            int saw_so=0;
            for (jj=fao[irnum];
                 jj < saoelem[irnum+ct[ir].degeneracy()-1]; jj++) {
              double t=0;
              for (ii=0; ii < gbs.nbasis(); ii++)
                t += blc[ii]*ret.get_element(ii,jj);

              if (fabs(t) > 0.1) {
                saw_so=1;
                break;
              }
            }

            if (saw_so)
              continue;

            int tir = saoelem[irnum+comp];
            if (tir >= maxsao[irnum+comp])
              continue;
            
            for (ii=0; ii < gbs.nbasis(); ii++)
              ret.set_element(ii,tir,blc[ii]);

            saoelem[irnum+comp]++;
          }
        }

        irnum += ct[ir].degeneracy();
      }

      for (g=0; g < _ng; g++)
        delete lc[g];
      delete lc;
    }
  }

  delete[] blc;
  delete[] saoelem;
  delete[] maxsao;

  return ret;
}
