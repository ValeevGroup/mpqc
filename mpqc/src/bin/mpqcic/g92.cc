
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

extern "C" {
#include <tmpl.h>
#include <util/sgen/sgen.h>
#include <math/array/math_lib.h>
}

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <util/misc/libmisc.h>
#include <math/scmat/local.h>
#include <chemistry/molecule/molecule.h>

#include "mpqc_int.h"
#include "g92_int.h"

//////////////////////////////////////////////////////////////////////////////

static void
read_geometry(RefMolecule& mol, const RefKeyVal& keyval, FILE *outfp)
{
  StateInBinXDR si("geom.dat","r+");

  int iter;
  si.get(iter);

  RefMolecule molecule;
  molecule.restore_state(si);
  mol = molecule;
}

//////////////////////////////////////////////////////////////////////////////

int
main(int argc, char *argv[])
{
  int i,j,ij;
  int errcod, geom_code=-1;
  
  int read_geom, opt_geom, nopt;
  int print_geometry, make_pdb=0;
  int geometry_converged;
  int do_freq;

  RefMolecule mol;
  
  struct stat stbuf;
  RefSCVector gradient, geometry;

 // initialize timing for mpqc

  tim_enter("g92");
  tim_enter("input");

  FILE *outfile = stdout;
  
  char *filename = (argv[1]) ? argv[1] : "mpqc.in";

  RefKeyVal keyval;

  fprintf(outfile,"\n       SC swallows G92 whole\n\n\n");

   // initialize keyval
  RefKeyVal pkv(new ParsedKeyVal(filename));
  RefKeyVal ppkv(new PrefixKeyVal(":mpqc :default",pkv));
  pkv = new ParsedKeyVal("input",ppkv);
  keyval = new AggregateKeyVal(ppkv,pkv);

  pkv = ppkv = 0;

  mol = keyval->describedclassvalue("molecule");
    
  print_geometry = keyval->booleanvalue("print_geometry");
  make_pdb = keyval->booleanvalue("write_pdb");
  read_geom = keyval->booleanvalue("read_geometry");
  opt_geom = keyval->booleanvalue("optimize_geometry");
  do_freq = keyval->booleanvalue("frequencies");
  if (keyval->error() != KeyVal::OK)
    do_freq = 1;
  
  char *molname = keyval->pcharvalue("filename");
  if (keyval->error() != KeyVal::OK)
    molname = "g92foo";

  nopt = keyval->intvalue("nopt");
  if (keyval->error() != KeyVal::OK)
    nopt=1;

  if (keyval->exists("filename"))
    filename = keyval->pcharvalue("filename");

  fprintf(outfile,"\n  mpqc options:\n");
  fprintf(outfile,"    optimize_geometry  = %s\n",(opt_geom)?"YES":"NO");
  fprintf(outfile,"    nopt               = %d\n",nopt);
  fprintf(outfile,"    print_geometry     = %s\n",(print_geometry)?"YES":"NO");
  fprintf(outfile,"    frequencies        = %s\n",(do_freq)?"YES":"NO");

  if (print_geometry)
    mol->print();
    
  // initialize the geometry optimization stuff
  if (opt_geom) {
    fprintf(outfile,"\n");
    geom_code = Geom_init_mpqc(mol,keyval);
  } else if (read_geom && stat("geom.dat",&stbuf)==0 && stbuf.st_size!=0) {
    read_geometry(mol,keyval,outfile);
  }

  // write pdb file if requested
  if (make_pdb)
    Geom_write_pdb(keyval,mol,"initial geometry");

  gradient = Geom_dim_natom3()->create_vector();
  
  RefKeyVal g92kv = new PrefixKeyVal(":g92 :default",keyval);

  tim_exit("input");

  int iter=0;
  while (geom_code != GEOM_DONE && geom_code != GEOM_ABORT && iter < nopt) {

    gradient.assign(0.0);
    double energy=0;

    // calculate new scf_vector and gradient
    tim_enter("scf_vect");
    // call g92 interface from mike
    if (run_g92(molname, g92kv, mol, energy, gradient) < 0) {
      fprintf(stderr,"g92 puked\n");
      exit(1);
    }
    tim_exit("scf_vect");

    // g92 gives us forces, not the gradient
    gradient.scale(-1.0);
    
    printf("\n  energy from g92 = %20.10f\n\n",energy);

    printf("  gradient from g92:\n");
    for (i=ij=0; i < mol->natom(); i++) {
      printf("  %5d",i+1);
      for (j=0; j < 3; j++,ij++) {
        printf(" %20.13f",gradient.get_element(ij));
      }
      printf("\n");
    }

    fprintf(outfile,"\n");
    geom_code = Geom_update_mpqc(gradient, keyval);
    iter++;
  }

  if (opt_geom && geom_code != GEOM_DONE && geom_code != GEOM_ABORT &&
      iter == nopt) {
    geometry_converged = 0;
    fprintf(outfile,"  Too many geometry iterations: quitting\n");
  } else {
    geometry_converged = 1;
  }

  if (opt_geom) {
    Geom_done_mpqc(keyval, geometry_converged);

   // write pdb file if requested
    if (make_pdb && geometry_converged)
      Geom_write_pdb(keyval,mol,"final geometry");
    else if (make_pdb)
      Geom_write_pdb(keyval,mol,"converged geometry");
  }

  if (do_freq) {
    double energy;
    int nmodes, nimag;

    RefSCDimension dim1 = Geom_dim_natom3();
    RefSCDimension dim2 = new LocalSCDimension(dim1->n()-6);
    RefSCDimension dim3 = new LocalSCDimension(dim1->n() * dim2->n());
    RefSCVector normals = dim3->create_vector();
    RefSCVector freqs = dim2->create_vector();
    RefSCVector fc = dim3->create_vector();
    
    if (g92_freq_driver(molname, mol, g92kv, 0, energy, gradient, freqs,
                        normals, fc, nmodes, nimag) < 0) {
      fprintf(stderr,"could not do g92 freqs\n");
    }

    freqs->print("frequencies");
    
  }

  tim_print(0);

  fflush(outfile);

  return 0;
}

