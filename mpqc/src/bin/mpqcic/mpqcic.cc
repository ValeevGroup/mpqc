
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

extern "C" {
#include <tmpl.h>
#include <util/sgen/sgen.h>
#include <math/dmt/libdmt.h>
#include <math/array/math_lib.h>
#include <comm/picl/picl.h>
#include <comm/picl/ext/piclext.h>
}

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <util/misc/libmisc.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include <chemistry/qc/dmtscf/scf_dmt.h>
#include <chemistry/qc/force/libforce.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/molecule/molecule.h>

extern "C" {
#include "scf_ffo.gbl"
#include "opt2_fock.gbl"
}

#include "mpqc_int.h"
#include "opt2.h"

//////////////////////////////////////////////////////////////////////////////

static void read_geometry(centers_t&, const RefKeyVal&, FILE*);
static void mkcostvec(centers_t*, sym_struct_t*, dmt_cost_t*);

///////////////////////////////////////////////////////////////////////////

static void
clean_and_exit()
{
  picl_prober();

  close0(0);
  exit(0);
}

static void
init_mp()
{
  int nproc,me,host;
  int top,ord,dir;

  open0(&nproc,&me,&host);
  setarc0(&nproc,&top,&ord,&dir);
}

static void
init_dmt(centers_t *centers, sym_struct_t *sym_info)
{
  int i;

  int_initialize_1e(0,0,centers,centers);
  int_initialize_offsets1(centers,centers);

  int nshell = centers->nshell;
  int nshtr = nshell * (nshell+1) / 2;

  int *shellmap = new int[nshell];
  memset(shellmap,0,sizeof(int)*nshell);

  for (i=0; i < nshell; i++) shellmap[i] = INT_SH_NFUNC(centers,i);

  dmt_cost_t *costvec = new dmt_cost_t[nshtr];

  if (costvec) {
    mkcostvec(centers,sym_info,costvec);
    dmt_def_map2(centers->nfunc, centers->nshell, shellmap, costvec, 0);
    delete[] costvec;
  } else {
    dmt_def_map2(centers->nfunc, centers->nshell, shellmap, 0, 1);
  }

  delete[] shellmap;  

  if (mynode0() == 0) printf("\n");
  dmt_map_examine();

  int_done_offsets1(centers,centers);
  int_done_1e();
}

static void
reset_centers(centers_t& centers, RefMolecule& mol)
{
  for (int i=0; i < mol->natom(); i++) {
    centers.center[i].r[0] = mol->atom(i)[0];
    centers.center[i].r[1] = mol->atom(i)[1];
    centers.center[i].r[2] = mol->atom(i)[2];
  }
}

////////////////////////////////////////////////////////////////////////////

int
main(int argc, char *argv[])
{
  int errcod, geom_code=-1;
  
  int do_scf, do_grad, do_mp2, do_opt2_v1, do_opt2_v2;
  int read_geom, opt_geom, nopt, proper;
  int save_fock, save_vector, print_geometry, make_pdb=0;
  int localp, throttle, sync_loop, node_timings;
  int geometry_converged;

  char *dertype="none";
  double dens;

  RefMolecule mol;
  RefGaussianBasisSet gbs;
  
  centers_t centers, oldcenters, *tcenters;
  scf_struct_t scf_info;
  sym_struct_t sym_info;

  dmt_matrix Scf_Vec,Fock,FockO;
  double_matrix_t gradient;
  
  struct stat stbuf;

 // initialize the picl routines
  init_mp();

 // initialize timing for mpqc

  tim_enter("mpqcnode");
  tim_enter("input");

  FILE *outfile = stdout;
  
  char *filename = (argv[1]) ? argv[1] : "mpqc.in";

  RefKeyVal keyval;

  int nfzc, nfzv, mem_alloc;

  if (mynode0() == 0) {
    fprintf(outfile,
        "\n       MPQC: Massively Parallel Quantum Chemistry\n\n\n");

    fprintf(outfile,"  Running on a %s with %d nodes.\n",
            machine_type(), numnodes0());
    fflush(outfile);

   // initialize keyval
    RefKeyVal pkv(new ParsedKeyVal(filename));
    RefKeyVal ppkv(new PrefixKeyVal(":mpqc :default",pkv));
    pkv = new ParsedKeyVal("input",ppkv);
    keyval = new AggregateKeyVal(ppkv,pkv);

    pkv = ppkv = 0;

    mol = keyval->describedclassvalue("molecule");
    
    gbs = keyval->describedclassvalue("basis");
    tcenters = gbs->convert_to_centers_t(mol.pointer());

    init_centers(&centers);
    init_centers(&oldcenters);
    assign_centers(&centers,tcenters);
    free_centers(tcenters);

    if (sym_struct_from_pg(mol->point_group(), centers, sym_info) < 0) {
      fprintf(stderr,"mpqcic:  could not form sym_info\n");
      exit(1);
    }

    RefKeyVal scfkv = new PrefixKeyVal(":scf :default",keyval);
    if (scf_init_scf_struct(scfkv, centers, scf_info) < 0) {
      fprintf(stderr,"mpqcic:  could not form scf_info\n");
      exit(1);
    }
    scfkv=0;

    // Read in opt2 stuff

    RefKeyVal opt2input(new PrefixKeyVal(":opt2 :mpqc :default", keyval));
    nfzc = opt2input->intvalue("frozen_docc");
    nfzv = opt2input->intvalue("frozen_uocc");
    mem_alloc = opt2input->intvalue("mem");
    if (!mem_alloc) mem_alloc = 8000000;

    opt2input = 0;

   // read input, and initialize various structs

    do_grad=0;
    if (keyval->exists("dertype")) dertype = keyval->pcharvalue("dertype");
    if (strcmp(dertype,"none")) do_grad=1;

    do_scf = 1;
    if (keyval->exists("do_scf")) do_scf = keyval->booleanvalue("do_scf");

    proper = 1;
    if (keyval->exists("properties"))
      proper = keyval->booleanvalue("properties");

    localp = 0;
    if (keyval->exists("local_P"))
      localp = keyval->booleanvalue("local_P");

    print_geometry = 0;
    if (keyval->exists("print_geometry"))
      print_geometry = keyval->booleanvalue("print_geometry");

    make_pdb = keyval->booleanvalue("write_pdb");

    do_mp2 = 0;
    if (keyval->exists("mp2"))
      do_mp2 = keyval->booleanvalue("mp2");

    do_opt2_v1 = 0;
    if (keyval->exists("opt2_v1"))
      do_opt2_v1 = keyval->booleanvalue("opt2_v1");

    do_opt2_v2 = 0;
    if (keyval->exists("opt2_v2"))
      do_opt2_v2 = keyval->booleanvalue("opt2_v2");

    read_geom = 0;
    if (keyval->exists("read_geometry"))
      read_geom = keyval->booleanvalue("read_geometry");

    opt_geom = 0;
    if (keyval->exists("optimize_geometry"))
      opt_geom = keyval->booleanvalue("optimize_geometry");
    if (opt_geom) do_grad=1;

    nopt=1;
    if(opt_geom) {
      if (keyval->exists("nopt")) nopt = keyval->intvalue("nopt");
      }

    throttle = 0;
    if (keyval->exists("throttle")) throttle = keyval->intvalue("throttle");

    sync_loop = 1;
    if (keyval->exists("sync_loop")) sync_loop = keyval->intvalue("sync_loop");

    save_fock=0;
    if (keyval->exists("save_fock"))
      save_fock = keyval->booleanvalue("save_fock");

    save_vector=1;
    if (keyval->exists("save_vector"))
      save_vector = keyval->booleanvalue("save_vector");

    node_timings=0;
    if (keyval->exists("node_timings"))
      node_timings = keyval->booleanvalue("node_timings");

    dens=2.5;
    if (keyval->exists("points_per_ang"))
      dens = keyval->doublevalue("points_per_ang");

    if (keyval->exists("filename"))
      filename = keyval->pcharvalue("filename");

    fprintf(outfile,"\n  mpqc options:\n");
    fprintf(outfile,"    do_scf             = %s\n",(do_scf)?"YES":"NO");
    fprintf(outfile,"    dertype            = %s\n",dertype);
    fprintf(outfile,"    optimize_geometry  = %s\n",(opt_geom)?"YES":"NO");
    fprintf(outfile,"    nopt               = %d\n",nopt);
    fprintf(outfile,"    properties         = %s\n",(proper)?"YES":"NO");
    fprintf(outfile,"    points_per_ang     = %f\n",dens);
    fprintf(outfile,"    save_vector        = %s\n",(save_vector)?"YES":"NO");
    fprintf(outfile,"    save_fock          = %s\n",(save_fock)?"YES":"NO");
    fprintf(outfile,"    node_timings       = %s\n",(node_timings)?"YES":"NO");
    fprintf(outfile,"    throttle           = %d\n",throttle);
    fprintf(outfile,"    sync_loop          = %d\n",sync_loop);
    fprintf(outfile,"    print_geometry     = %s\n",
                                                  (print_geometry)?"YES":"NO");
    fprintf(outfile,"    mp2                = %s\n\n", (do_mp2)?"YES":"NO");
    fprintf(outfile,"    opt2_v1            = %s\n\n", (do_opt2_v1)?"YES":"NO");
    fprintf(outfile,"    opt2_v2            = %s\n\n", (do_opt2_v2)?"YES":"NO");

    if (save_vector) {
      fprintf(outfile,"  scf vector will be written to file %s.scfvec\n",
              scf_info.fname);
    }
    
    if (save_fock) {
      fprintf(outfile,
              "  fock matrices will be written to file(s) %s.fock (%s.fock)\n",
              scf_info.fname, scf_info.fname);
    }
    
   // pretty print the scf options
    scf_print_options(stdout, scf_info);

    if (print_geometry) {
      mol->print();
    }
    
   // initialize the geometry optimization stuff
    if (opt_geom) {
      fprintf(outfile,"\n");
      geom_code = Geom_init_mpqc(mol,keyval);
    } else if (read_geom && stat("geom.dat",&stbuf)==0 && stbuf.st_size!=0) {
      read_geometry(centers,keyval,outfile);
    }

   // we may have changed the geometry in mol, so reform centers
    reset_centers(centers,mol);

   // write pdb file if requested
    if (make_pdb)
      Geom_write_pdb(keyval,mol,"initial geometry");
  }

  sgen_reset_bcast0();

  bcast0_scf_struct(&scf_info,0,0);
  bcast0_sym_struct(&sym_info,0,0);
  bcast0_centers(&centers,0,0);

  bcast0(&do_scf,sizeof(int),mtype_get(),0);
  bcast0(&do_grad,sizeof(int),mtype_get(),0);
  bcast0(&nopt,sizeof(int),mtype_get(),0);
  bcast0(&opt_geom,sizeof(int),mtype_get(),0);
  bcast0(&throttle,sizeof(int),mtype_get(),0);
  bcast0(&sync_loop,sizeof(int),mtype_get(),0);
  bcast0(&save_fock,sizeof(int),mtype_get(),0);
  bcast0(&save_vector,sizeof(int),mtype_get(),0);
  bcast0(&node_timings,sizeof(int),mtype_get(),0);
  bcast0(&localp,sizeof(int),mtype_get(),0);
  bcast0(&proper,sizeof(int),mtype_get(),0);
  bcast0(&do_mp2,sizeof(int),mtype_get(),0);
  bcast0(&do_opt2_v1,sizeof(int),mtype_get(),0);
  bcast0(&do_opt2_v2,sizeof(int),mtype_get(),0);
  bcast0(&dens,sizeof(double),mtype_get(),0);
  bcast0(&nfzc,sizeof(int),mtype_get(),0);
  bcast0(&nfzv,sizeof(int),mtype_get(),0);
  bcast0(&mem_alloc,sizeof(int),mtype_get(),0);

 // if we're using a projected guess vector, then initialize oldcenters
  if (scf_info.proj_vector) {
    if (mynode0()==0) {
      RefGaussianBasisSet gbs = keyval->describedclassvalue("pbasis");
      tcenters = gbs->convert_to_centers_t(mol.pointer());

      assign_centers(&oldcenters,tcenters);
      free_centers(tcenters);
      
      int_normalize_centers(&oldcenters);
    }

    bcast0_centers(&oldcenters,0,0);
  }

 // initialize the dmt library
  init_dmt(&centers,&sym_info);
  
 // initialize force and geometry routines

  if (do_grad) {
    RefKeyVal fkv;
    
    if (mynode0()==0) {
      fprintf(outfile,"\n");
    
      fkv = new PrefixKeyVal(":force :default",keyval);
    }
    
    if (scf_info.iopen)
      dmt_force_osscf_keyval_init(fkv.pointer(),outfile);
    else
      dmt_force_csscf_keyval_init(fkv.pointer(),outfile);

    allocbn_double_matrix(&gradient,"n1 n2",3,centers.n);
  }

  if (opt_geom) {
    bcast0(&geom_code,sizeof(int),mtype_get(),0);

    if (geom_code==GEOM_ABORT || geom_code==GEOM_DONE) {
      fprintf(outfile,"mpqcnode: geom_code says you are done or in trouble\n");
      clean_and_exit();
    }
  }

 // set the throttle for libdmt loops
  dmt_set_throttle(throttle);

 // set the sync_loop for libdmt loops
  dmt_set_sync_loop(sync_loop);

 // allocate memory for vector and fock matrices

  Scf_Vec = dmt_create("scf vector",scf_info.nbfao,COLUMNS);
  Fock = dmt_create("fock matrix",scf_info.nbfao,SCATTERED);
  if (scf_info.iopen)
    FockO = dmt_create("open fock matrix",scf_info.nbfao,SCATTERED);
  else
    FockO = dmt_nil();

 // read the scf vector if it's there
  char vecfile[512];
  sprintf(vecfile,"./%s.scfvec",scf_info.fname);
  
  if (scf_info.restart && stat(vecfile,&stbuf)==0 && stbuf.st_size!=0) {
    dmt_read(vecfile,Scf_Vec);
    
    if (mynode0()==0)
      fprintf(outfile,"\n  read vector from file %s\n\n",vecfile);
  } else {
    scf_info.restart = 0;
  }
  
  tim_exit("input");

 // if we need vector, get one
  if (do_scf) {

    int iter=0;
    while(geom_code != GEOM_DONE && geom_code != GEOM_ABORT && iter < nopt) {

      // broadcast new geometry information */
      for (int i=0; i < centers.n; i++)
        bcast0(centers.center[i].r,sizeof(double)*3,mtype_get(),0);

      // calculate new scf_vector

      tim_enter("scf_vect");
      errcod = scf_vector(&scf_info, &sym_info, &centers, Fock, FockO, Scf_Vec,
                          &oldcenters, outfile);
      tim_exit("scf_vect");

      if (errcod != 0) {
        fprintf(outfile,"trouble forming scf vector\n");
        clean_and_exit();
      }
    
      scf_info.restart=1;

      if (save_vector)
        dmt_write(vecfile,Scf_Vec);

      // get new geometry

      if (mynode0() == 0) fprintf(outfile,"\n");

      if (geom_code == GEOM_COMPUTE_GRADIENT) {
        if (!scf_info.iopen) {
          dmt_force_csscf(outfile, Fock, Scf_Vec,
                          &centers, &sym_info, scf_info.nclosed, &gradient);
        } else {
          dmt_force_osscf(outfile, Fock, FockO, Scf_Vec,
                          &centers, &sym_info, scf_info.nclosed,
                          scf_info.nopen, &gradient);
        }

        if (mynode0()==0) {
          fprintf(outfile,"\n");

          RefSCVector gradv(mol->dim_natom3());
          
          int i,j,ij;
          for (j=ij=0; j < gradient.n2; j++) {
            for (i=0; i < gradient.n1; i++,ij++) {
              gradv->set_element(ij,gradient.d[i][j]);
            }
          }
          
          geom_code = Geom_update_mpqc(gradv, keyval);
          reset_centers(centers,mol);
        }
      }
    
      else if (do_grad) {
        if (!scf_info.iopen) {
          dmt_force_csscf(outfile, Fock, Scf_Vec,
                          &centers, &sym_info, scf_info.nclosed, &gradient);
        } else {
          dmt_force_osscf(outfile, Fock, FockO, Scf_Vec,
                          &centers, &sym_info, scf_info.nclosed,
                          scf_info.nopen, &gradient);
        }
      }

      bcast0(&geom_code,sizeof(geom_code),mtype_get(),0);
      iter++;
    }

    if (opt_geom && geom_code != GEOM_DONE && geom_code != GEOM_ABORT &&
        mynode0() == 0 && iter == nopt) {
      geometry_converged = 0;
      fprintf(outfile,"  Too many geometry iterations: quitting\n");
    } else {
      geometry_converged = 1;
    }

    if (save_fock) {
      char fockfile[512];
      sprintf(fockfile,"%s.fock",scf_info.fname);
      dmt_write(fockfile,Fock);
      
      if (scf_info.iopen) {
        sprintf(fockfile,"%s.focko",scf_info.fname);
        dmt_write(fockfile,FockO);
      }
    }
  }

  if (opt_geom && mynode0()==0) {
    Geom_done_mpqc(keyval, geometry_converged);

   // write pdb file if requested
    if (make_pdb && geometry_converged)
      Geom_write_pdb(keyval,mol,"final geometry");
    else if (make_pdb)
      Geom_write_pdb(keyval,mol,"converged geometry");
  }

  if (do_grad) {
    if (scf_info.iopen)
      dmt_force_osscf_done();
    else
      dmt_force_csscf_done();
  }

 /* print out some useful stuff */
#if 0
  if (proper) {
    tim_enter("properties");
    int nat=centers.n;
    int natr=nat*(nat+1)/2;
    double_vector_t dipole,charge,mcharge,lcharge;
    double_matrix_t bond_pops,bond_indx;
    expts_t mulpts;

    allocbn_expts(&mulpts,"n",nat+natr);
    for(i=0; i < mulpts.n; i++) allocbn_exmul(&mulpts.p[i],"charge",0.0);

    allocbn_double_matrix(&bond_pops,"n1 n2",nat,nat);
    allocbn_double_matrix(&bond_indx,"n1 n2",nat,nat);
    allocbn_double_vector(&charge,"n",nat);
    allocbn_double_vector(&mcharge,"n",nat);
    allocbn_double_vector(&lcharge,"n",nat);
    allocbn_double_vector(&dipole,"n",3);

    scf_mulliken(&centers,&scf_info,&irreps,SCF_VEC,&bond_pops,
       &bond_indx,&mcharge,outfile);

    scf_lowdin(&centers,&scf_info,&irreps,SCF_VEC,
       &bond_indx,&lcharge,outfile);

    scf_dipole_and_ex_mulliken(&centers,&scf_info,&irreps,SCF_VEC,
        &bond_pops,&mulpts,&dipole,outfile);

    if(me==0) {
      FILE *dp = fopen("points.dat","w");
      FILE *mp = fopen("mopac.dat","w");
      if (dp)  {
	fprintf(dp,"%d\n",mulpts.n);
	for(i=0; i < mulpts.n; i++) {
	  fprintf(dp,"%20.13e %20.13e %20.13e %20.13e\n",
	               mulpts.p[i].charge,
                       mulpts.p[i].r[0],mulpts.p[i].r[1],mulpts.p[i].r[2]);
	  }
        fclose(dp);
	}
      if (mp) {
        expts_t at;
        at.n = centers.n;
        at.p = &mulpts.p[mulpts.n-centers.n];

        fprintf(mp," 1SCF AVSWRT ESP CONNOLLY GRADIENT DEPVAR=1.0\n");
        fprintf(mp," no comment\n\n");
        for (i=0; i < centers.n; i++)
          fprintf(mp,"%4s%13.4f 1 %13.4f 1 %13.4f 1\n",
            centers.center[i].atom,
            at.p[i].r[0]*0.529177,
            at.p[i].r[1]*0.529177,
            at.p[i].r[2]*0.529177
            );
        fclose(mp);
	}
      }

    Scf_charges_from_esp(&centers,SCF_VEC,&charge,&dipole,
                                              &mulpts,dens,-1,outfile,keyval);

    if(me==0) {
      char *atomlb;

      fprintf(outfile,"\n  Net Charges\n");
      fprintf(outfile,
   "      Atom            Mulliken              Lowdin                 ESP\n");
      for(i=0; i < centers.n; i++) {
        if(atom_labels) atomlb=atom_labels[i];
        else atomlb=centers.center[i].atom;

        fprintf(outfile,"%5d %5s %20.10f %20.10f %20.10f\n",
            i+1,atomlb,mcharge.d[i],lcharge.d[i],charge.d[i]);
        }
      printf("\n");
      }

    free_double_vector(&charge);
    free_double_vector(&mcharge);
    free_double_vector(&lcharge);
    free_double_vector(&dipole);
    free_double_matrix(&bond_indx);
    free_double_matrix(&bond_pops);
    free_expts(&mulpts);
    tim_exit("properties");
  }

  if (do_mp2) {
    if (!do_scf) dmt_read(fockfile,FOCK);
    tim_enter("mp2");
    mp2_hah(&centers,&scf_info,SCF_VEC,FOCK,outfile,keyval);
    tim_exit("mp2");
  }
#endif

  if (do_opt2_v1 || do_opt2_v2) {
    dmt_matrix S = dmt_create("libscfv3 overlap matrix",scf_info.nbfao,SCATTERED);
    dmt_matrix SAHALF;
    dmt_matrix SC;
    dmt_matrix EVECS;
    dmt_matrix SCR;
    double_vector_t occ_num;
    double_vector_t evals;
    dmt_matrix SCR1,SCR2,SCR3;
    // this got free'ed somewhere
    if(scf_info.iopen)
      FockO = dmt_create("opt2 open fock matrix",scf_info.nbfao,SCATTERED);
    SCR1 = dmt_create("opt2:scr1",scf_info.nbfao,COLUMNS);
    SCR2 = dmt_create("opt2:scr2",scf_info.nbfao,COLUMNS);
    SCR3 = dmt_create("opt2:scr3",scf_info.nbfao,COLUMNS);

    if (!do_scf) {

      if(mynode0()==0) fprintf(stderr,"Must do scf before opt2. Program exits\n");
      clean_and_exit();
      }

    tim_enter("opt2");

    scf_ffo(S, &scf_info, &sym_info, &centers, Scf_Vec, Fock, FockO);

    if (scf_info.iopen) {
      scf_make_opt2_fock(&scf_info,Fock,FockO,Scf_Vec,SCR1,SCR2,SCR3);
      dmt_free(SCR3);
      dmt_copy(Fock,SCR1); /* dmt_diag needs a columns dist. matrix */
      allocbn_double_vector(&evals,"n",scf_info.nbfao);
      dmt_diag(SCR1,SCR2,evals.d); /*SCR2 transforms from old to new mo basis */
      dmt_copy(Scf_Vec,SCR1);
      dmt_transpose(SCR1);
      dmt_mult(SCR1,SCR2,Scf_Vec);
      dmt_free(SCR1);
      dmt_free(SCR2);
      }
    else {
      /* form 'sahalf' matrix sahalf = u*ei^-0.5*u~ */
      allocbn_double_vector(&evals,"n",scf_info.nbfao);
      SAHALF= dmt_create("libscfv3 scf_core_guess scr4",scf_info.nbfao,COLUMNS);      EVECS = dmt_create("libscfv3 scf_core_guess scr3",scf_info.nbfao,COLUMNS);      SCR = dmt_create("libscfv3 scf_core_guess scr5",scf_info.nbfao,COLUMNS);
      SC = dmt_columns("libscfv3 scf_core_guess scr1",S);
      /* diagonalize overlap matrix */
      dmt_diag(SC,EVECS,evals.d);
      /* form SAHALF matrix (s^(-1/2), Sz&Ostl p. 143) */
      for(int i=0; i < scf_info.nbfao; i++) evals.d[i] = 1.0/sqrt(evals.d[i]);
      dmt_fill(SAHALF,0.0);
      dmt_set_diagonal(SAHALF,evals.d);
      /* form the orthogonalization matrix S^(-1/2) (Szabo&Ostlund p. 143)
       * (called SAHALF here) */
      dmt_transpose(EVECS);
      dmt_mult(SAHALF,EVECS,SCR);
      dmt_mult(EVECS,SCR,SAHALF);
      dmt_free(EVECS);
      dmt_free(SC);
      dmt_free(SCR);
      dmt_free(S);
      dmt_copy(Fock,SCR1); /* need a columns distr. matrix */
      dmt_mult(SCR1,SAHALF,SCR2);
      dmt_mult(SAHALF,SCR2,SCR3);
      /* SCR3 is now the Fock matrix in the orthogonalized ao basis */
      dmt_diag(SCR3,Scf_Vec,evals.d);
      dmt_copy(Scf_Vec,SCR1);
      dmt_mult(SAHALF,SCR1,Scf_Vec); /* Sz&Ostl p.146 point 9 */
      dmt_free(SCR1);
      dmt_free(SCR2);
      dmt_free(SCR3);
      }

    if (do_opt2_v1) {
        sync0();
        tim_enter("opt2_v1");
        opt2_v1(&centers,&scf_info,Scf_Vec,&evals,nfzc,nfzv,mem_alloc,
                    outfile);
        tim_exit("opt2_v1");
      }
    if (do_opt2_v2) {
        sync0();
        tim_enter("opt2_v2");
        opt2_v2(&centers,&scf_info,Scf_Vec,&evals,nfzc,nfzv,mem_alloc,
                    outfile);
        tim_exit("opt2_v2");
      }
    free_double_vector(&evals);

    tim_exit("opt2");
    }

  tim_print(node_timings);

  fflush(outfile);

  clean_and_exit();
}


static void
read_geometry(centers_t& centers, const RefKeyVal& keyval, FILE *outfp)
{
#if 0  
  StateInBinXDR si("geom.dat","r+");

  int iter;
  si.get(iter);

  fprintf(outfp,"\n using geometry from iteration %d\n",iter);

  RefSymmCoList symm_coords;
  symm_coords = SymmCoList::restore_state(si);
  symm_coords = 0;

  si.get(iter);
  if (iter) {
    symm_coords = SymmCoList::restore_state(si);
    symm_coords = 0;
  }

  Molecule mol(si);

  for(int i=0; i < centers.n; i++) {
    centers.center[i].r[0] = mol[i][0];
    centers.center[i].r[1] = mol[i][1];
    centers.center[i].r[2] = mol[i][2];
    }

  mol.print(outfp);

  RefSimpleCoList list = Geom_form_simples(mol);

  int nadd;
  if(nadd=keyval->count("add_simp")) {
    for(int i=0; i < nadd; i++) {
      char *val = keyval->pcharvalue("add_simp",i,0);

      if (!strcmp("stre",val))
        list->add(new Stre(keyval.pointer(),"add_simp",i));
      else if (!strcmp("bend",val))
        list->add(new Bend(keyval.pointer(),"add_simp",i));
      else if (!strcmp("tors",val))
        list->add(new Tors(keyval.pointer(),"add_simp",i));
      else if (!strcmp("out",val))
        list->add(new Out(keyval.pointer(),"add_simp",i));
      else if (!strcmp("linip",val))
        list->add(new LinIP(keyval.pointer(),"add_simp",i));
      else if (!strcmp("linop",val))
        list->add(new LinOP(keyval.pointer(),"add_simp",i));
      delete[] val;
      }
    }

  Geom_calc_simples(list,mol);

  fprintf(outfp,"\n  internal coordinates\n");
  Geom_print_pretty(list);

  list=0;
#endif
}

static void
mkcostvec(centers_t *centers,sym_struct_t *sym_info,dmt_cost_t *costvec)
{
  int flags;
  int i,j,k,l;
  int ij,kl,ijkl;
  int ioffi,ioffij;
  int g,gi,gj,gk,gl,gij,gkl,gijkl;
  int Qvecij,bound,cost;
  int nb;
  int leavel;
  int use_symmetry=(sym_info->g>1);
  int nproc=numnodes0();
  int me=mynode0();
  double *intbuf;
  extern signed char *scf_bnd_Qvec;

 /* free these up for now */
  int_done_offsets1(centers,centers);
  int_done_1e();

  int_initialize_offsets2(centers,centers,centers,centers);

  flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2;

  intbuf =
    int_initialize_erep(flags,0,centers,centers,centers,centers);

  scf_init_bounds(centers,intbuf);

  for (i=ij=0; i<centers->nshell; i++) {
    int nconi = centers->center[centers->center_num[i]].
                       basis.shell[centers->shell_num[i]].ncon;
    int nprimi = centers->center[centers->center_num[i]].
                       basis.shell[centers->shell_num[i]].nprim;
    int ami = centers->center[centers->center_num[i]].
                       basis.shell[centers->shell_num[i]].type[0].am;
    for (j=0; j<=i; j++,ij++) {
      int nconj = centers->center[centers->center_num[j]].
                         basis.shell[centers->shell_num[j]].ncon;
      int nprimj = centers->center[centers->center_num[j]].
                         basis.shell[centers->shell_num[j]].nprim;
      int amj = centers->center[centers->center_num[j]].
                         basis.shell[centers->shell_num[j]].type[0].am;

      costvec[ij].i = i;
      costvec[ij].j = j;
      costvec[ij].magnitude = (int) scf_bnd_Qvec[ij];
      costvec[ij].ami = ami;
      costvec[ij].amj = amj;
      costvec[ij].nconi = nconi;
      costvec[ij].nconj = nconj;
      costvec[ij].nprimi = nprimi;
      costvec[ij].nprimj = nprimj;
      costvec[ij].dimi = INT_SH_NFUNC(centers,i);;
      costvec[ij].dimj = INT_SH_NFUNC(centers,j);;
    }
  }

  int_done_erep();
  int_done_offsets2(centers,centers,centers,centers);
  int_initialize_1e(0,0,centers,centers);
  int_initialize_offsets1(centers,centers);
  scf_done_bounds();
}
