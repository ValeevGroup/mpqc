
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
  int i;
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
    int leave=0;

    // we want the lowest numbered shell in a group of equivalent shells
    for (int g=0; g < _ng; g++)
      if (_shell_map[i][g] < i) {
        leave=1;
        break;
      }
    
    if (leave)
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
      leave=0;
      for (int g=0; g < _ng; g++) {
        int gi = _shell_map[i][g];
        int gj = _shell_map[j][g];
        int gij = ioff(gi,gj);
        if (gij < ij) {
          leave=1;
          break;
        }
        else if (gij == ij)
          nij++;
      }

      if (leave)
        continue;

      _lamij[ij] = (char) (_ng/nij);
    }
  }

  // form reducible representation of the basis functions
  double *red_rep = new double[_ng];
  memset(red_rep,0,sizeof(double)*_ng);
  
  for (int i=0; i < _natom; i++) {
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

  _nbf_in_ir = new int[_nirrep];
  for (int i=0; i < _nirrep; i++) {
    double t=0;
    for (int g=0; g < _ng; g++)
      t += ct[i][g]*red_rep[g];

    _nbf_in_ir[i] = ((int) (t+0.5))/_ng;
    if (!ct.complex())
      _nbf_in_ir[i] *= (int)(ct[i].degeneracy()+0.5);
  }

  delete[] red_rep;
}

RefSCMatrix
PetiteList::r(int g)
{
  SymmetryOperation so =
    _gbs->molecule()->point_group().char_table().symm_operation(g);
  GaussianBasisSet& gbs = *_gbs.pointer();

  RefSCMatrix ret = gbs.basisdim()->create_matrix(gbs.basisdim());
  ret.assign(0.0);
  
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
              ret.set_element(func_j+jj,func_i+ii,rr(jj,ii));
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
  fprintf(o,"PetiteList:\n");
  fprintf(o,"  _natom = %d\n",_natom);
  fprintf(o,"  _nshell = %d\n",_nshell);
  fprintf(o,"  _ng = %d\n",_ng);
  fprintf(o,"  _nirrep = %d\n",_nirrep);

  fprintf(o,"\n");
  fprintf(o,"  _atom_map = \n");
  int i;
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
  for (int i=0; i < _nirrep; i++)
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
      fn = new int[_nsh];
      for (int i=0; i < _nsh; i++)
        fn[i] = f0+i;
      
      c = new double*[_nsh];
      for (int i=0; i < _nsh; i++) {
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
};
    
RefSCMatrix
PetiteList::aotoso()
{
  GaussianBasisSet& gbs = *_gbs.pointer();
  Molecule& mol = *gbs.molecule().pointer();

  RefSCMatrix ret = gbs.basisdim()->create_matrix(gbs.basisdim());
  ret.assign(0.0);
  
  // create the character table for the point group
  CharacterTable ct = mol.point_group().char_table();
  SymmetryOperation so;
  ct.print();
  
  double *red_rep = new double[_ng];
  int *ninir = new int[_nirrep];

  int *saoelem = new int[_nirrep];
  saoelem[0]=0;
  for (int i=1; i < _nirrep; i++)
    saoelem[i] = saoelem[i-1]+_nbf_in_ir[i-1];

  // loop over unique shells
  for (int i=0; i < _natom; i++) {
    for (int s=0; s < gbs.nshell_on_center(i); s++) {
      int shell_i = gbs.shell_on_center(i,s);
      
      if (!_p1[shell_i])
        continue;
        
      // find out how many shells are equivalent to this one
      int neqs=0;
      for (int g=0; g < _ng; g++)
        if (shell_i==_shell_map[shell_i][g])
          neqs++;

      neqs = _ng/neqs;
      
      // loop over contractions now to get a reducible representation for
      // the shell
      memset(red_rep,0,sizeof(double)*_ng);

      for (int g=0; g < _ng; g++) {
        int j = _atom_map[i][g];
        if (i!=j && _atom_map[j][g]!=j)
          continue;
        
        so = ct.symm_operation(g);
        
        for (int c=0; c < gbs(i,s).ncontraction(); c++) {
          int am=gbs(i,s).am(c);
          
          if (am==0)
            red_rep[g] += neqs*1.0;
          else {
            Rotation r(am,so,gbs(i,s).is_pure(c));
            red_rep[g] += neqs*r.trace();
          }
        }
      }

      // now extract number of functions of each symmetry that we can expect
      memset(ninir,0,sizeof(int)*_nirrep);
      for (int ir=0; ir < _nirrep; ir++) {
        double t=0;
        for (int g=0; g < _ng; g++)
          t += ct[ir][g]*red_rep[g];
        
        ninir[ir] = ((int) (t+0.5))/_ng;
      }

      // ok, we know how many functions we're looking for now, so start
      // forming linear combinations of basis functions

      lin_comb **lc = new lin_comb*[_ng];
      
      for (int g=0; g < _ng; g++) {
        so = ct.symm_operation(g);
        int j = _atom_map[i][g];

        int func_i = gbs.shell_to_function(gbs.shell_on_center(i,s));
        int func_j = gbs.shell_to_function(gbs.shell_on_center(j,s));

        lc[g] = new lin_comb(gbs(i,s).nfunction(),gbs.nbasis(),func_i);
        lin_comb& lcg = *lc[g];
        
        int fi=0;
        for (int c=0; c < gbs(i,s).ncontraction(); c++) {
          int am=gbs(i,s).am(c);

          if (am==0) {
            lcg.coef(fi,func_j) = 1.0;
          } else {
            Rotation rr(am,so,gbs(i,s).is_pure(c));
            for (int ii=0; ii < rr.dim(); ii++)
              for (int jj=0; jj < rr.dim(); jj++)
                lcg.coef(fi+ii,func_j+jj) = rr(ii,jj);
          }

          fi += gbs(i,s).nfunction(c);
          func_i += gbs(i,s).nfunction(c);
          func_j += gbs(i,s).nfunction(c);
        }
      }

      // form the combinations
      double *blc = new double[gbs.nbasis()];
      
      for (int ir=0; ir < ct.nirrep(); ir++) {
        //printf("irrep %s\n",ct[ir].symbol());
        
        for (int fn=0; fn < gbs(i,s).nfunction(); fn++) {
          memset(blc,0,sizeof(double)*gbs.nbasis());

          for (int g=0; g < _ng; g++) {
            lin_comb& lcg = *lc[g];

            for (int f=0; f < gbs.nbasis(); f++)
              blc[f] += ct[ir][g]*lcg.coef(fn,f);
          }

          double c1=0;
          for (int ii=0; ii < gbs.nbasis(); ii++)
            c1 += blc[ii]*blc[ii];

          if (c1 < 1.0e-3)
            continue;
          
          c1 = 1.0/sqrt(c1);
          
          for (int ii=0; ii < gbs.nbasis(); ii++)
            blc[ii] *= c1;

          // check to see if we already have this SO (it happens sometimes,
          // so sue me).
          int break_this=0;
          for (int jj=0;  jj < saoelem[ir]; jj++) {
            double t=0;
            for (int ii=0; ii < gbs.nbasis(); ii++)
              t += blc[ii]*ret.get_element(ii,jj);
            if (fabs(t) > .9) {
              break_this=1;
              break;
            }
          }

          if (break_this)
            break;

          for (int ii=0; ii < gbs.nbasis(); ii++) {
            ret.set_element(ii,saoelem[ir],blc[ii]);
          }
          saoelem[ir]++;
          
          // if this is a degenerate irrep and this shell is not centered at
          // the origin, then look for the second component (T's are right
          // out!).
          if (neqs==1 || ct[ir].degeneracy() < 1.5)
            continue;
          
          // find a symop which doesn't map i into itself
          int g=0;
          while (_atom_map[i][g]==i) g++;
          int j=_atom_map[i][g];
          
          saoelem[ir]++;
#if 0
          printf("  %d",lc[0]->bfnum(fn));
          for (int ii=0; ii < gbs.nbasis(); ii++)
            printf(" %10.7f",blc[ii]);
          printf("\n");
#endif
        }
      }

      delete[] blc;
      for (int g=0; g < _ng; g++)
        delete lc[g];
      delete lc;
    }
  }

  delete[] red_rep;
  delete[] ninir;
  delete[] saoelem;

  return ret;
}
    
////////////////////////////////////////////////////////////////////////////
