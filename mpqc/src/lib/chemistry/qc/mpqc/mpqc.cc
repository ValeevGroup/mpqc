
extern "C" {
#include <stdio.h>
#include <stdlib.h>
}

#include <math.h>
#include <math/newmat7/newmat.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include "mpqc.h"

static const char *first_line = "%MPQC generated input";

int MPQC::do_exchange_energy(int f)
{
  int old = _do_exchange_energy;
  _do_exchange_energy = f;
  return old;
}

int MPQC::do_eigenvectors(int f)
{
  int old = _do_eigenvectors;
  _do_eigenvectors = f;
  return old;
}

int MPQC::do_eigenvectors_on_disk(int f)
{
  int old = _do_eigenvectors_on_disk;
  _do_eigenvectors_on_disk = f;
  return old;
}

int MPQC::do_reuse_old_eigenvectors(int f)
{
  int old = _do_reuse_old_eigenvectors;
  _do_reuse_old_eigenvectors = f;
  return old;
}

void MPQC::init(KeyVal&keyval)
{

  // if mpqc.in is around make sure that is is from MPQC
  FILE*input = fopen("mpqc.in","r");
  if (input) {
      char line[200];
      fgets(line,200,input);
      if (strncmp(first_line,line,strlen(first_line))) {
          fprintf(stderr,"MPQC::init: \"mpqc.in\" exists--remove it\n");
          abort();
        }
    }

  _have_eigenvectors_on_disk = 0;
  x_changed();
  bs_values = 0;
  bsg_values = 0;
  _do_exchange_energy = 0;
  _do_eigenvectors = 0;
  _do_eigenvectors_on_disk = 1;
  _do_reuse_old_eigenvectors = 1;

  _maxiter = keyval.intvalue("maxiter");
  if (keyval.error() != KeyVal::OK) _maxiter = 30;

  _nocc = keyval.intvalue("nocc");
  if (keyval.error() != KeyVal::OK) {
      _nocc = 0;
      int i;
      for (i=0; i<_mol.natom(); i++) {
          _nocc += (int)(_mol[i].element().charge()+0.5);
        }
      if ((_nocc/2)*2 != _nocc) {
          fprintf(stderr,"MPQC::can only do closed shells for now\n");
          abort();
        }
      else {
          _nocc = _nocc/2;
        }
    }

  _do_exchange_energy = keyval.intvalue("exchange");
  if (keyval.error() != KeyVal::OK) _do_exchange_energy = 0;
}

MPQC::MPQC(KeyVal&keyval,Molecule&mol,GaussianBasisSet&gbs):
OneBodyWavefunction(keyval,mol,gbs)
{
  init(keyval);
}

MPQC::MPQC(KeyVal&keyval,Molecule&mol,GaussianBasisSet&gbs,MolecularCoor&mc):
OneBodyWavefunction(keyval,mol,gbs,mc)
{
  init(keyval);
}

MPQC::~MPQC()
{
  if (bs_values) delete[] bs_values;
  if (bsg_values) delete[] bsg_values;
}

static void seek_to_line(FILE*fp,const char* match)
{
  char line[200];
  char* s;
  line[0] = '\0';
  while ((s = fgets(line,200,fp)) && strncmp(line,match,strlen(match))) {
      //printf("match = \"%s\"\n",match);
      //printf("line  = \"%s\"\n",line);
    }
  if (!s) {
      fprintf(stderr,"MPQC: got a end of file while parsing the output\n");
      abort();
    }
}

void MPQC::compute()
{
  GaussianBasisSet&_gbs = basis();

  // transfer the current postion into the molecule
  X_to_molecule();

  // create the input file for mpqc
  FILE* input;
  input = fopen("mpqc.in","w");
  fprintf(input,"%s\n",first_line);
  fprintf(input,"default:(\n");
  fprintf(input,"  restart = %d\n",_have_old_eigenvectors_on_disk);
  fprintf(input,"  warmrestart = no\n");
  fprintf(input,"  nocc = [ %d ]\n",_nocc);
  fprintf(input,"  local_P = yes\n");
  fprintf(input,"  maxiter = %d\n",_maxiter);
  fprintf(input,"  exchange = %d\n",_do_exchange_energy);
  fprintf(input,"  wfn = scf\n");
  fprintf(input,"  dertype = %s\n",do_gradient()?"first":"none");
  //fprintf(input,"  basis = basis\n");
  fprintf(input,"  optimize_geometry = no\n");
  fprintf(input,"  print_geometry = yes\n");
  fprintf(input,"  properties = no\n");
  fprintf(input,"  { basis atoms geometry } = {\n");
  Pix center;
  int i=0;
  for (center=_mol.first(); center!=0; _mol.next(center)) {
      fprintf(input,"    basis_%d %s [ %f %f %f ]\n",
              i,
	      _mol(center).element().symbol(),
	      _mol(center)[0],
	      _mol(center)[1],
	      _mol(center)[2]);
      i++;
    }
  fprintf(input,"    }\n");
  fprintf(input,"  )\n");

  // the basis set info
  // the centers in _gbs must correspond to those in _mol
  fprintf(input,"basis:(\n");
  for (i=0,center=_mol.first(); center!=0; _mol.next(center)) {
    fprintf(input,"  %s:basis_%d: [\n",_mol(center).element().name(),i);
    for (Pix j=_gbs.first_shell_on_center(center);
         j!=0;
         _gbs.next_shell_on_center(center,j)) {
        int k;

        // print out the angular momentum
        fprintf(input,"    (type: [");
        for (k=0; k < _gbs[j].ncontraction(); k++) {
            fprintf(input," am = %c",_gbs[j].amchar(k));
          }
        fprintf(input,"]\n");

        // print out the exp and coef array heading
        fprintf(input,"     {exp");
        for (k=0; k < _gbs[j].ncontraction(); k++) {
            fprintf(input," coef:%d",k);
          }
        fprintf(input,"}={\n");

        const char *fmt = " %16.12f";

        // loop thru each primitive and print out the exponent
        for (k=0; k<_gbs[j].nprimitive(); k++) {
            fprintf(input,fmt,_gbs[j].exponent(k));
            // print out the contraction coefficients
            for (int l=0; l<_gbs[j].ncontraction(); l++) {
                fprintf(input,fmt,_gbs[j].coefficient_norm(l,k));
              }
            fprintf(input,"\n");
          }
        fprintf(input,"      })\n");
      }
    fprintf(input,"    ]\n");
    i++;
    }
  fprintf(input,"  )\n");
  fclose(input);

  // run mpqc
  system("mpqcproc > mpqc.out");

  // parse the output file from mpqc
  FILE*output = fopen("mpqc.out","r");

  // find the energy
  double energy;
  seek_to_line(output,"  iter       total energy");
  char line[200];
  while (fgets(line,200,output)) {
      char first[30];
      double lastenergy = energy;
      if (sscanf(line,"%s %lf",first,&energy) == 1) {
	  energy= lastenergy;
	  if (!strcmp(first,"converged")) break;
	  else if (!(!strcmp(first,"checkpointing")
		     ||!strcmp(first,"resetting"))) {
	      fprintf(stderr,"MPQC: unexpected data in output: aborting\n");
	      abort();
	    }
	}
    }
  set_energy(energy);
  _have_eigenvectors_on_disk = 1;
  //printf("MPQC got energy %16.10f\n",energy);

  // find the exchange energy
  if (_do_exchange_energy) {
      fscanf(output,"%*s %*s %*s %*s %lf",&_exchange_energy);
      _have_exchange_energy = 1;
    }

  if (do_gradient()) {
      seek_to_line(output,"The total gradient:");
      int atom;
      double x,y,z;
      ColumnVector g(3*_mol.natom());
      while(fscanf(output,"%d %lf %lf %lf",&atom,&x,&y,&z) == 4) {
	  g(atom*3 + 1) = x;
	  g(atom*3 + 2) = y;
	  g(atom*3 + 3) = z;
	}
      set_gradient(g);
      //printf("MPQC got gradient:\n");
      //Print(g);
    }

}

double MPQC::exchange_energy()
{
  if (!_have_exchange_energy) {
      int old = do_exchange_energy(1);
      compute();
      do_exchange_energy(old);
    }
  return _exchange_energy;
}

const Matrix& MPQC::eigenvectors()
{
  const GaussianBasisSet&_gbs=basis();
  if (!_have_eigenvectors) {
      if (!_have_eigenvectors_on_disk) {
          int old = do_eigenvectors_on_disk(1);
          compute();
          do_eigenvectors_on_disk(old);
        }
      int nbasis = _gbs.nbasis();
      _eigenvectors.ReDimension(nbasis,nbasis);
      read_vector("mpqc.scfvec",nbasis,_eigenvectors);
    }
  return _eigenvectors;
}

void MPQC::x_changed()
{
  OneBodyWavefunction::x_changed();
  _have_exchange_energy = 0;
  _have_eigenvectors = 0;
  if (_do_reuse_old_eigenvectors) {
      _have_old_eigenvectors_on_disk = _have_eigenvectors_on_disk;
    }
  _have_eigenvectors_on_disk = 0;
}

void MPQC::maxiter(int m)
{
  _maxiter = m;
}

double MPQC::occupation(int i)
{
  if (i<_nocc) return 2.0;
  else return 0.0;
}

void MPQC::print(FILE*fp)
{
  fprintf(fp,"MPQC:\n");
  fprintf(fp,"  maxiter = %d\n",_maxiter);
  fprintf(fp,"  exchange = %d\n",_do_exchange_energy);
  _mol.print(fp);
  //_gbs.print(fp);
}
