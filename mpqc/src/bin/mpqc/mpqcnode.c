
static char rcsid[]="$Id$";


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

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

#define ioff(i) ((i)*((i)+1)/2)
#define IOFF(i,j) ((i)>(j)) ? ioff((i))+(j) : ioff((j))+(i)

static void clean_and_exit();
static void mkcostvec(centers_t *centers,sym_struct_t *sym_info,int *costvec);

int host;
char* argv0;

main(argc,argv)
int argc;
char *argv[];
{
  int i;
  int errcod;
  int nproc,me,top,ord,dir;
  int do_grad,opt_geom,nopt,iter,proper;
  int localp;
  int nshell,nshtr;
  int throttle,geom_code,sync_loop;
  int save_fock,save_vector,restart,print_geometry;
  int node_timings;
  int size,count,use_dip;
  int *shellmap;
  double energy,dens;
  int_vector_t costvec;
  int mp2;
  int ij=0;

  char *dertype="none";
  char *filename="mpqc";
  char **atom_labels;
  char vecfile[512],fockfile[512],fockofile[512],oldvecfile[512];
  FILE *infile,*outfile;

  scf_struct_t scf_info;
  sym_struct_t sym_info;
  scf_irreps_t irreps;
  centers_t centers;

  dmt_matrix SCF_VEC,FOCK,FOCKO;
  double_matrix_t gradient;

/* turn out the lights */

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

  if(me==0) {
    if(host0()!=0) {
      extern char* compile_date();
      fprintf(outfile,"\n  loaded from host program\n");
      fprintf(outfile,"\n  compiled on %s\n\n",compile_date());
      }
    else {
      extern char* compile_date();
      fprintf(outfile,
        "\n       MPQC: Massively Parallel Quantum Chemistry\n\n\n");
      fprintf(outfile,"  loaded from workstation or SRM\n");
      fprintf(outfile,"\n  compiled on %s\n\n",compile_date());
      }
    fprintf(outfile,"  Running on a %s with %d nodes.\n",machine_type(),nproc);
    fflush(outfile);

    infile = fopen("mpqc.in","r+");
    ip_initialize(infile,outfile);

    ip_cwk_add(":default");
    ip_cwk_add(":mpqc");

    ip_append_from_input("input",stdout);
    ip_cwk_clear();
    ip_cwk_add(":default");
    ip_cwk_add(":mpqc");

/* read input, and initialize various structs */

    do_grad=0;
    ip_string("dertype",&dertype,0);
    if(strcmp(dertype,"none")) do_grad=1;

    proper = 1;
    errcod = ip_boolean("properties",&proper,0);

    use_dip = 1;
    errcod = ip_boolean("use_dipole",&use_dip,0);

    localp = 0;
    errcod = ip_boolean("local_P",&localp,0);

    print_geometry = 0;
    errcod = ip_boolean("print_geometry",&print_geometry,0);

    mp2 = 0;
    errcod = ip_boolean("mp2",&mp2,0);

    opt_geom = 0;
    errcod = ip_boolean("optimize_geometry",&opt_geom,0);
    if(opt_geom) do_grad=1;

    nopt=1;
    if(opt_geom) ip_data("nopt","%d",&nopt,0);

    throttle = 0;
    ip_data("throttle","%d",&throttle,0);

    sync_loop = 1;
    ip_data("sync_loop","%d",&sync_loop,0);

    save_fock=0;
    errcod = ip_boolean("save_fock",&save_fock,0);

    save_vector=1;
    errcod = ip_boolean("save_vector",&save_vector,0);

    node_timings=0;
    errcod = ip_boolean("node_timings",&node_timings,0);

    restart=1;
    errcod = ip_boolean("restart",&restart,0);

    dens=5.0;
    errcod = ip_data("points_per_ang","%lf",&dens,0);

    errcod = ip_string("filename",&filename,0);

    fprintf(outfile,"\n  mpqc options:\n");
    fprintf(outfile,"    dertype            = %s\n",dertype);
    fprintf(outfile,"    optimize_geometry  = %s\n",(opt_geom)?"YES":"NO");
    fprintf(outfile,"    nopt               = %d\n",nopt);
    fprintf(outfile,"    properties         = %s\n",(proper)?"YES":"NO");
    fprintf(outfile,"    use_dipole         = %s\n",(use_dip)?"YES":"NO");
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
  
    errcod = ip_count("labels",&count,0);
    if(errcod == 0 && count==centers.n) {
      atom_labels = (char **) malloc(sizeof(char *)*centers.n);
      check_alloc(atom_labels,"mpqcnode: atom labels");

      for(i=0; i < count; i++) ip_string("labels",&atom_labels[i],1,i);
      }
    else atom_labels=NULL;
    }

  sgen_reset_bcast0();

  bcast0_scf_struct(&scf_info,0,0);
  bcast0_sym_struct(&sym_info,0,0);
  bcast0_centers(&centers,0,0);
  bcast0_scf_irreps(&irreps,0,0);

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
  bcast0(&use_dip,sizeof(int),mtype_get(),0);
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

  if(opt_geom && me==0) {
    fprintf(outfile,"\n");
    geom_code = geom_initialize(outfile,outfile,&centers);
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

  if(!localp) {
    allocbn_int_vector(&costvec,"n",nshtr);
    zero_int_vector(&costvec);

#if 1
    mkcostvec(&centers,&sym_info,costvec.i);

#if 0
    gsum0(costvec.i,costvec.n,2,mtype_get(),0);
    bcast0_int_vector(&costvec,0,0);
#endif

    for(i=0; i < nshtr ; i++) costvec.i[i]+=1;
#else
    for (i=0; i < nshell; i++) {
      int j;
      int nconi = centers.center[centers.center_num[i]].
                         basis.shell[centers.shell_num[i]].ncon;
      int nprimi = centers.center[centers.center_num[i]].
                         basis.shell[centers.shell_num[i]].nprim;
      int ami = centers.center[centers.center_num[i]].
                         basis.shell[centers.shell_num[i]].type[0].am;

      for (j=0; j <= i; j++,ij++) {
        int nconj = centers.center[centers.center_num[j]].
                           basis.shell[centers.shell_num[j]].ncon;
        int nprimj = centers.center[centers.center_num[j]].
                           basis.shell[centers.shell_num[j]].nprim;
        int amj = centers.center[centers.center_num[j]].
                           basis.shell[centers.shell_num[j]].type[0].am;

        /* best so far
         * costvec[ij] = shellmap[i]*shellmap[j] + ami*amj + nprimi + nprimj;
         */

        /* close second */
        costvec.i[ij] = shellmap[i]*shellmap[j] + (1+ami*amj)*(nprimi + nprimj);
        }
      }
#endif

    dmt_def_map2(scf_info.nbfao,centers.nshell,shellmap,costvec.i,0);
    free_int_vector(&costvec);
    }
  else {
    costvec.n=0;
    costvec.i=NULL;
    dmt_def_map2(scf_info.nbfao,centers.nshell,shellmap,costvec.i,1);
    }

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
    if(!scf_info.proj_vector) ip_done();

/* begin iterations here */

  iter=0;
  while(  (    (!opt_geom)
            || (geom_code != GEOM_DONE && geom_code != GEOM_ABORT))
        && iter < nopt) {

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

  /* release memory used by the input parser */
    if(iter==0 && scf_info.proj_vector) ip_done();

    energy = scf_info.nuc_rep+scf_info.e_elec;

 /* get new geometry */

    if(opt_geom && geom_code == GEOM_COMPUTE_ENERGY) {
      int_initialize_1e(0,0,&centers,&centers);
      int_initialize_offsets1(&centers,&centers);
      if (mynode0()==0) geom_code = geom_update(energy,NULL,NULL);
      int_done_offsets1(&centers,&centers);
      int_done_1e();
      }
    else if(opt_geom && geom_code == GEOM_COMPUTE_GRADIENT) {
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
      if (mynode0()==0) geom_code = geom_update(energy,&gradient,NULL);
      int_done_offsets1(&centers,&centers);
      int_done_1e();
      }
    else if (opt_geom && geom_code == GEOM_COMPUTE_HESSIAN) {
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

  if(opt_geom && geom_code != GEOM_DONE && geom_code != GEOM_ABORT && me==0 &&
     iter==nopt) {
    fprintf(outfile,"  Too many geometry iterations: quitting\n");
    }

  if(save_vector) dmt_write(vecfile,SCF_VEC);
  if(save_fock) {
    dmt_write(fockfile,FOCK);
    if(scf_info.iopen) dmt_write(fockofile,FOCKO);
    }

  if(opt_geom && me==0) geom_done();

  if(do_grad) {
    if(scf_info.iopen) dmt_force_osscf_done();
    else dmt_force_csscf_done();
    }

 /* print out some usefull stuff */

  if(proper) {
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

#if 0
    if(me==0) {
      fprintf(outfile,"\n\n  Results from Extended Mulliken Analysis\n");
      fprintf(outfile,"\n       Charge        x          y          z\n");
      for(i=0; i < mulpts.n; i++) {
        fprintf(outfile,"%3d %10.5f %10.5f %10.5f %10.5f\n",i+1,
        mulpts.p[i].charge,mulpts.p[i].r[0],mulpts.p[i].r[1],mulpts.p[i].r[2]);
        }
      }
#endif

    scf_charges_from_esp(&centers,&scf_info,&charge,&dipole,
                                              &mulpts,dens,use_dip,outfile);

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
    }

  tim_print(node_timings);

  clean_and_exit(host);
  }

static void
clean_and_exit(host)
int host;
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


static void
mkcostvec(centers_t *centers,sym_struct_t *sym_info,int *costvec)
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

#if 1
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

#if 0
      costvec[ij] = Qvec[ij];
      costvec[ij] *= 2;
      costvec[ij] += 26;
#endif

      if (ami+amj==0) costvec[ij]=1;
      else if (ami+amj==1) costvec[ij]=15;
      else if (ami+amj==2) costvec[ij]=180;
      else costvec[ij]=680;

      if (costvec[ij]<0) costvec[ij]=1;
      }
    }
#else
 /* loop over shell blocks and figure out where they fit in */
  for (i=me; i<centers->nshell; i+=nproc) {
    if(use_symmetry && !sym_info->p1[i]) continue;
    ioffi=ioff(i);

    for (j=0; j<=i; j++) {
      ij=ioffi+j;
      if(use_symmetry) {
        if (!sym_info->lamij[ij]) continue;
        else ioffij=ioff(ij);
        leavel=0;
        }
      Qvecij=(int)Qvec[ij];

      for (k=0; k<=i; k++) {
        kl=ioff(k);
        for (l=0; l<=(k==i?j:k); l++,kl++) {
          if(use_symmetry) {
            ijkl=ioffij+kl;
            leavel=0;
            for(g=0; g < sym_info->g ; g++) {
              gi = sym_info->shell_map[i][g];
              gj = sym_info->shell_map[j][g];
              gk = sym_info->shell_map[k][g];
              gl = sym_info->shell_map[l][g];
              gij = IOFF(gi,gj);
              gkl = IOFF(gk,gl);
              gijkl = IOFF(gij,gkl);
              if(gijkl > ijkl) leavel=1;
              }
            if(leavel) continue;
            }

        /* this tell us how likely it is we'll calculate the integral */
          bound = Qvecij + (int) Qvec[kl];
          bound += bound;
          bound += 26;

        /* this tells us how many times we'll have to do the bugger */
          if(j==k && (i==j || k==l)) nb=1;
          else if (i==j||i==k||j==k||j==l||k==l) nb=2;
          else nb=3;

          costvec[ij] += nb;
          }
        }
      }
    }
#endif

  int_done_erep();
  int_done_offsets2(centers,centers,centers,centers);
  int_initialize_1e(0,0,centers,centers);
  int_initialize_offsets1(centers,centers);
  scf_done_bounds();
  }
