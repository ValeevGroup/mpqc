
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

extern "C" {
#include <tmpl.h>
#include <util/ipv2/ip_libv2.h>
#include <util/sgen/sgen.h>
#include <math/dmt/libdmt.h>
#include <util/misc/libmisc.h>
#include <util/bio/libbio.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include <chemistry/qc/dmtscf/scf_dmt.h>
#include <chemistry/qc/force/libforce.h>
#include <chemistry/qc/geom/libgeom.h>
#include <comm/picl/picl.h>
#include <comm/picl/ext/piclext.h>


#if defined(I860) || defined(SGI)
int led(int);
void bzero(void*,int);
#endif
}

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/simple.h>
#include <chemistry/molecule/symm.h>
#include <chemistry/molecule/simpleQCList.h>
#include <chemistry/molecule/symmQCList.h>
#include <math/nihmatrix/nihmatrix.h>
#include <math/nihmatrix/lmath.h>

#include "mpqc_int.h"

#define ioff(i) ((i)*((i)+1)/2)
#define IOFF(i,j) ((i)>(j)) ? ioff((i))+(j) : ioff((j))+(i)

static void clean_and_exit(int);
static void read_geometry(centers_t&,RefKeyVal,FILE*);

int host;
char* argv0;

dmt_cost_t *costvec;
static void mkcostvec(centers_t*, sym_struct_t*, dmt_cost_t*);

main(int argc, char *argv[])
{
  int i;
  int errcod;
  int nproc,me,top,ord,dir;
  int do_scf,do_grad,read_geom,opt_geom,nopt,iter,proper;
  int localp;
  int nshell,nshtr;
  int throttle,geom_code=-1,sync_loop;
  int save_fock,save_vector,restart,print_geometry;
  int node_timings;
  int size,count;
  int *shellmap;
  double energy,dens;
  int mp2;

  char *dertype="none";
  char *filename="mpqc";
  char **atom_labels;
  char vecfile[512],fockfile[512],fockofile[512],oldvecfile[512];
  FILE *infile,*outfile;

  scf_struct_t scf_info;
  sym_struct_t sym_info;
  scf_irreps_t irreps;
  centers_t centers;
  struct stat stbuf;

  dmt_matrix SCF_VEC,FOCK,FOCKO;
  double_matrix_t gradient;

/* turn out the lights */

#if defined(I860) && !defined(PARAGON)
  led(0);
#endif

/* some PICL emulators need argv0 */
  argv0 = argv[0];

/* get some info from the host program */

  open0(&nproc,&me,&host);
  setarc0(&nproc,&top,&ord,&dir);

/* prepare for potential debugging */
#if (defined(SGI) || defined(SUN4)) && defined(USE_DEBUG)
  debug_init(argv[0]);
#endif

#if ((defined(SGI) || defined(SUN4)) && (!defined(SABER))) && defined(USE_DEBUG)
  malloc_debug_on_error();
#endif

/* initialize timing for mpqc */

  tim_enter("mpqcnode");
  tim_enter("input");

/* set up i/o. node 0 will do all the reading and pass what's needed on to
 * the other nodes 
 */

  outfile = stdout;

  RefKeyVal keyval;

  if(me==0) {
    if(host0()!=0) {
      fprintf(outfile,"\n  loaded from host program\n");
      }
    else {
      fprintf(outfile,
        "\n       MPQC: Massively Parallel Quantum Chemistry\n\n\n");
      fprintf(outfile,"  loaded from workstation or SRM\n");
      }
    fprintf(outfile,"  Running on a %s with %d nodes.\n",machine_type(),nproc);
    fflush(outfile);

// various things still use the ipv2 routines to get input
// so this hack is needed for now
#define NEED_IPV2
#ifdef NEED_IPV2
    infile = fopen("mpqc.in","r+");
    if (!infile) {
      perror("infile");
      exit(1);
      }
    ip_initialize(infile,outfile);

    ip_cwk_add(":default");
    ip_cwk_add(":mpqc");

    ip_append_from_input("input",stdout);
    ip_cwk_clear();
    ip_cwk_add(":default");
    ip_cwk_add(":mpqc");
#endif
// If ipv2 also must read the input, then keyval specific things
// will break it.  In this case the keyval specific input will
// be read from keyval.in.
#ifdef NEED_IPV2
    ParsedKeyVal *pparsed = new ParsedKeyVal("mpqc.in");
    if (stat("keyval.in",&stbuf)==0 && stbuf.st_size!=0)
      pparsed->read("keyval.in");
    RefKeyVal parsed(pparsed);
#else
    RefKeyVal parsed(new ParsedKeyVal("mpqc.in"));
#endif
    keyval = new PrefixKeyVal(":intco :mpqc :default",*parsed.pointer());
    RefKeyVal extension(new ParsedKeyVal("input",*keyval.pointer()));
    RefKeyVal agg(new AggregateKeyVal(*parsed.pointer(),*extension.pointer()));
    keyval = new PrefixKeyVal(":intco :mpqc :default",*agg.pointer());

    // drop references to the temporary keyvals
    parsed = 0;
    extension = 0;
    agg = 0;


/* read input, and initialize various structs */

    do_grad=0;
    if (keyval->exists("dertype")) dertype = keyval->pcharvalue("dertype");
    if(strcmp(dertype,"none")) do_grad=1;

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

    mp2 = 0;
    if (keyval->exists("mp2"))
      mp2 = keyval->booleanvalue("mp2");

    read_geom = 0;
    if (keyval->exists("read_geometry"))
      read_geom = keyval->booleanvalue("read_geometry");

    opt_geom = 0;
    if (keyval->exists("optimize_geometry"))
      opt_geom = keyval->booleanvalue("optimize_geometry");
    if(opt_geom) do_grad=1;

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

    restart=1;
    if (keyval->exists("restart"))
      restart = keyval->booleanvalue("restart");

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
    fprintf(outfile,"    restart            = %s\n",(restart)?"YES":"NO");
    fprintf(outfile,"    save_vector        = %s\n",(save_vector)?"YES":"NO");
    fprintf(outfile,"    save_fock          = %s\n",(save_fock)?"YES":"NO");
    fprintf(outfile,"    node_timings       = %s\n",(node_timings)?"YES":"NO");
    fprintf(outfile,"    throttle           = %d\n",throttle);
    fprintf(outfile,"    sync_loop          = %d\n",sync_loop);
    fprintf(outfile,"    print_geometry     = %s\n",
                                                  (print_geometry)?"YES":"NO");
    fprintf(outfile,"    mp2                = %s\n\n", (mp2)?"YES":"NO");

    sprintf(vecfile,"%s.scfvec",filename);
    sprintf(oldvecfile,"%s.oldvec",filename);
    sprintf(fockfile,"%s.fock",filename);
    sprintf(fockofile,"%s.focko",filename);

    if(save_vector)
      fprintf(outfile,"  scf vector will be written to file %s\n",vecfile);
    if(save_fock)
      fprintf(outfile,"  fock matrices will be written to file(s) %s (%s)\n",
         fockfile,fockofile);

    errcod = scf_get_input(&scf_info,&sym_info,&irreps,&centers,outfile);
    if(errcod != 0) {
      fprintf(outfile,"trouble reading input\n");
      clean_and_exit(host);
      }

  /* read in user-defined labels for each atom */
  
    count = keyval->count("atom_labels");
    if(count==centers.n) {
      atom_labels = (char **) malloc(sizeof(char *)*centers.n);
      check_alloc(atom_labels,"mpqcnode: atom labels");

      char * tmp;
      for(i=0; i < count; i++) {
        tmp = keyval->pcharvalue("atom_labels",i);
        atom_labels[i]=tmp;
        }
      }
    else atom_labels=NULL;

    if(opt_geom) {
      fprintf(outfile,"\n");
      geom_code = Geom_init_mpqc(outfile,outfile,&centers,keyval);
      }
    else if (read_geom && stat("geom.dat",&stbuf)==0 && stbuf.st_size!=0) {
      read_geometry(centers,keyval,outfile);
      }
    }

  sgen_reset_bcast0();

  bcast0_scf_struct(&scf_info,0,0);
  bcast0_sym_struct(&sym_info,0,0);
  bcast0_centers(&centers,0,0);
  bcast0_scf_irreps(&irreps,0,0);

  bcast0(&do_scf,sizeof(int),mtype_get(),0);
  bcast0(&do_grad,sizeof(int),mtype_get(),0);
  bcast0(&nopt,sizeof(int),mtype_get(),0);
  bcast0(&opt_geom,sizeof(int),mtype_get(),0);
  bcast0(&throttle,sizeof(int),mtype_get(),0);
  bcast0(&sync_loop,sizeof(int),mtype_get(),0);
  bcast0(&restart,sizeof(int),mtype_get(),0);
  bcast0(&save_fock,sizeof(int),mtype_get(),0);
  bcast0(&save_vector,sizeof(int),mtype_get(),0);
  bcast0(&node_timings,sizeof(int),mtype_get(),0);
  bcast0(&localp,sizeof(int),mtype_get(),0);
  bcast0(&proper,sizeof(int),mtype_get(),0);
  bcast0(&mp2,sizeof(int),mtype_get(),0);
  bcast0(&dens,sizeof(double),mtype_get(),0);

  if(me==0) size=strlen(vecfile)+1;
  bcast0(&size,sizeof(int),mtype_get(),0);
  bcast0(vecfile,size,mtype_get(),0);

  if(me==0) size=strlen(oldvecfile)+1;
  bcast0(&size,sizeof(int),mtype_get(),0);
  bcast0(oldvecfile,size,mtype_get(),0);

  if(me==0) size=strlen(fockfile)+1;
  bcast0(&size,sizeof(int),mtype_get(),0);
  bcast0(fockfile,size,mtype_get(),0);

  if(me==0) size=strlen(fockofile)+1;
  bcast0(&size,sizeof(int),mtype_get(),0);
  bcast0(fockofile,size,mtype_get(),0);

/* initialize centers */

  int_initialize_1e(0,0,&centers,&centers);
  int_initialize_offsets1(&centers,&centers);

/* initialize force and geometry routines */

  if(do_grad) {
    if(me==0) fprintf(outfile,"\n");
    if(scf_info.iopen) dmt_force_osscf_init(outfile);
    else dmt_force_csscf_init(outfile);

    allocbn_double_matrix(&gradient,"n1 n2",3,centers.n);
    }

  if(opt_geom) {
    bcast0(&geom_code,sizeof(int),mtype_get(),0);

    if(geom_code==GEOM_ABORT || geom_code==GEOM_DONE) {
      fprintf(outfile,"mpqcnode: geom_code says you are done or in trouble\n");
      clean_and_exit(host);
      }
    }

/* set up shell map for libdmt */

  nshell=centers.nshell;
  nshtr=nshell*(nshell+1)/2;

  shellmap = (int *) malloc(sizeof(int)*nshell);
  bzero(shellmap,sizeof(int)*nshell);

  for(i=0; i < centers.nshell ; i++) shellmap[i] = INT_SH_NFUNC(&centers,i);

  costvec = (dmt_cost_t*) malloc(sizeof(dmt_cost_t)*nshtr);

  mkcostvec(&centers,&sym_info,costvec);

  dmt_def_map2(scf_info.nbfao,centers.nshell,shellmap,costvec,0);

  free(costvec);

  free(shellmap);
  if(me==0) printf("\n");
  dmt_map_examine();

  int_done_offsets1(&centers,&centers);
  int_done_1e();

/* set the throttle for libdmt loops */

  dmt_set_throttle(throttle);

/* set the sync_loop for libdmt loops */

  dmt_set_sync_loop(sync_loop);

/* allocate memory for vector and fock matrices */

  SCF_VEC = dmt_create("scf vector",scf_info.nbfao,COLUMNS);
  FOCK = dmt_create("fock matrix",scf_info.nbfao,SCATTERED);
  if(scf_info.iopen)
    FOCKO = dmt_create("open fock matrix",scf_info.nbfao,SCATTERED);
  else
    FOCKO = dmt_nil();

/* if restart, the read in old scf vector */
  
  if(restart) {
    dmt_read(vecfile,SCF_VEC);
    if(me==0) fprintf(outfile,"\n  read vector from file %s\n\n",vecfile);
    scf_info.restart=1;
    }

  tim_exit("input");

/* release memory used by the input parser */
#ifdef NEED_IPV2
  if(!scf_info.proj_vector) ip_done();
#endif
  // keyval is needed by the Geom routines
  //if (!scf_info.proj_vector) keyval = 0;

/* skip scf if we already have a vector and just want properties */
  if (!do_scf) goto THERE;

/* begin iterations here */

  iter=0;
  while( geom_code != GEOM_DONE && geom_code != GEOM_ABORT &&
         iter < nopt) {

  /* broadcast new geometry information */
    for (i=0; i<centers.n; i++) {
      bcast0(centers.center[i].r,sizeof(double)*3,mtype_get(),0);
      }

  /* print geometry */

    if(me==0 && print_geometry) {
      char *atomlb;

      int_initialize_1e(0,0,&centers,&centers);
      int_initialize_offsets1(&centers,&centers);
      fprintf(outfile,"\n\n         molecular geometry in atomic units");
      fprintf(outfile,"        basis set\n\n");
      for(i=0; i < centers.n ; i++) {
        if(atom_labels) atomlb=atom_labels[i];
        else atomlb=centers.center[i].atom;

        fprintf(outfile,"%5d %5s %12.7f %12.7f %12.7f      %s\n",i+1,
                atomlb,
                centers.center[i].r[0],
                centers.center[i].r[1],
                centers.center[i].r[2],
                centers.center[i].basis.name);
        }
      fprintf(outfile,"\n");
      int_done_offsets1(&centers,&centers);
      int_done_1e();
      }

  /* calculate new scf_vector */

    tim_enter("scf_vect");
    errcod = scf_vector(&scf_info,&sym_info,&irreps,&centers,
                          FOCK,FOCKO,SCF_VEC,outfile);
    tim_exit("scf_vect");

    if(errcod != 0) {
      fprintf(outfile,"trouble forming scf vector\n");
      clean_and_exit(host);
      }
    scf_info.restart=1;

    if(save_vector) dmt_write(vecfile,SCF_VEC);

  /* release memory used by the input parser */
#ifdef NEED_IPV2
    if(iter==0 && scf_info.proj_vector) ip_done();
#endif
    // keyval is needed by the Geom routines
    //if(iter==0 && scf_info.proj_vector) keyval = 0;

    energy = scf_info.nuc_rep+scf_info.e_elec;

 /* get new geometry */

    if(geom_code == GEOM_COMPUTE_ENERGY) {
      int_initialize_1e(0,0,&centers,&centers);
      int_initialize_offsets1(&centers,&centers);
      if (mynode0()==0) geom_code = Geom_update_mpqc(energy,NULL,NULL,keyval);
      int_done_offsets1(&centers,&centers);
      int_done_1e();
      }
    else if(geom_code == GEOM_COMPUTE_GRADIENT) {
      /* tim_enter("gradient"); */
      if(!scf_info.iopen) {
        dmt_force_csscf(outfile,FOCK,SCF_VEC,
            &centers,&sym_info,irreps.ir[0].nclosed,&gradient);
        }
      else {
        int ndoc=irreps.ir[0].nclosed;
        int nsoc=irreps.ir[0].nopen;
        dmt_force_osscf(outfile,FOCK,FOCKO,SCF_VEC,
            &centers,&sym_info,ndoc,nsoc,&gradient);
        }
      /* tim_exit("gradient"); */

      int_initialize_1e(0,0,&centers,&centers);
      int_initialize_offsets1(&centers,&centers);
      if (mynode0()==0)
        geom_code = Geom_update_mpqc(energy,&gradient,NULL,keyval);
      int_done_offsets1(&centers,&centers);
      int_done_1e();
      }
    else if (geom_code == GEOM_COMPUTE_HESSIAN) {
      if (mynode0()==0) {
        printf("A hessian was requested, but cannot be computed\n");
        geom_code = GEOM_ABORT;
        }
      }
    else if(do_grad) {
      /* tim_enter("gradient"); */
      if(!scf_info.iopen) {
        dmt_force_csscf(outfile,FOCK,SCF_VEC,
            &centers,&sym_info,irreps.ir[0].nclosed,&gradient);
        }
      else {
        int ndoc=irreps.ir[0].nclosed;
        int nsoc=irreps.ir[0].nopen;
        dmt_force_osscf(outfile,FOCK,FOCKO,SCF_VEC,
            &centers,&sym_info,ndoc,nsoc,&gradient);
        }
      /* tim_exit("gradient"); */
      }

    bcast0(&geom_code,sizeof(geom_code),mtype_get(),0);
    iter++;
    }

  int geometry_converged;
  if(opt_geom && geom_code != GEOM_DONE && geom_code != GEOM_ABORT && me==0 &&
     iter==nopt) {
      geometry_converged = 0;
      fprintf(outfile,"  Too many geometry iterations: quitting\n");
    }
  else {
      geometry_converged = 1;
    }

  if(save_fock) {
    dmt_write(fockfile,FOCK);
    if(scf_info.iopen) dmt_write(fockofile,FOCKO);
    }

THERE:

  if(opt_geom && me==0) Geom_done_mpqc(keyval,geometry_converged);

  if(do_grad) {
    if(scf_info.iopen) dmt_force_osscf_done();
    else dmt_force_csscf_done();
    }

 /* print out some usefull stuff */

  if(proper) {
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

  if (mp2) {
    if (!do_scf) dmt_read(fockfile,FOCK);
    //tim_enter("mp2 loop");
    //mp2_loop(&centers,&scf_info,SCF_VEC,FOCK,outfile,keyval);
    //tim_exit("mp2 loop");
    tim_enter("mp2");
    mp2_hah(&centers,&scf_info,SCF_VEC,FOCK,outfile,keyval);
    tim_exit("mp2");
    //tim_print(node_timings);
    //tim_enter("mp2 big");
    //mp2_hah_big(&centers,&scf_info,SCF_VEC,FOCK,outfile,keyval);
    //tim_exit("mp2 big");
    }

  tim_print(node_timings);

  clean_and_exit(host);
  }

static void
clean_and_exit(int host)
{
  int junk;

#if !defined(I860)
  picl_prober();

  /* tell host that we're done */
  if(mynode0()==0 && host0()) send0(&junk,sizeof(int),1,host);
#endif

 /* this one too */
  close0(0);
  exit(0);
  }

#if 0
double ftn_i_dsign(double a, double b)
{
  return (((b)>=0.0)? (fabs(a)): (-fabs(a)));
}
#endif

static void
read_geometry(centers_t& centers, RefKeyVal keyval, FILE *outfp)
{
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
}

#if defined(SUNMOS)
extern "C" {
void * sbrk(int) { int *foo = (int*) 0x7fffffff;  return (void*) foo; }

int openlog() { return 0; }
int syslog() { return 0; }

int atexit(void(*func)(void)) { return 0; }
}
#endif

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
  extern signed char *Qvec;

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
      costvec[ij].magnitude = (int) Qvec[ij];
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
