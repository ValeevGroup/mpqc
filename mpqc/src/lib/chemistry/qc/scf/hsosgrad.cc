
#define BOUNDS 1

#include <util/misc/newstring.h>
#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/scf/hsosscf.h>

static void
set_scale(double& coulombscale, double& exchangescale,
          int i, int j, int k, int l)
{
  double scale = 1.0;

  if ((i!=k)||(j!=l))
    scale *= 2.0;

  if (i!=j)
    scale *= 2.0;

  coulombscale = 0.5*scale;
  exchangescale = 0.25*scale;

  if (k!=l)
    coulombscale *= 2.0;

  if ((k!=l)&&(i==j))
    exchangescale *= 2.0;
}

static void
gr_density(const RefSCMatrix& vec, const RefSymmSCMatrix& dens,
           const RefSymmSCMatrix& opdens,
           int ndocc, int nsocc, double& pmax)
{
  int k;

  pmax=0.0;

  for (int i=0; i < vec->nrow(); i++) {
    for (int j=0; j <= i; j++) {
      double pt=0;
      for (k=0; k < ndocc; k++)
        pt += vec->get_element(i,k)*vec->get_element(j,k);
      
      double po=0;
      for (k=ndocc; k < ndocc+nsocc; k++)
        po += vec->get_element(i,k)*vec->get_element(j,k);
      
      dens->set_element(i,j,pt);
      opdens->set_element(i,j,po);

      if (fabs(pt)>pmax)
        pmax=fabs(pt);
    }
  }
  dens->scale(2.0);
  opdens->scale(2.0);
}

void
HSOSSCF::do_gradient(const RefSCVector& gradient)
{
  int i,x;

  double alpha[3][3], beta[3][3];

  memset(alpha,0,sizeof(double)*9);

  alpha[0][0] = 1.0;
  alpha[1][0] = alpha[0][1] = 0.5;
  alpha[1][1] = 0.25;

  memset(beta,0,sizeof(double)*9);

  beta[0][0] = -1.0;
  beta[1][0] = beta[0][1] = beta[1][1] = -0.5;
  
  // grab a reference to the scf_vector, presumably it is current
  _gr_vector = _eigenvectors.result_noupdate();
  
  // allocate storage for the temp arrays
  _gr_dens = _fock.clone();
  _gr_op_dens = _fock.clone();
  
  // form energy weighted density
  // first form MO fock matrices
  RefSymmSCMatrix mofock = _fock.clone();
  mofock.assign(0.0);
  mofock.accumulate_transform(_gr_vector.t(),_fock);
  
  RefSymmSCMatrix moofock = _op_fock.clone();
  moofock.assign(0.0);
  moofock.accumulate_transform(_gr_vector.t(),_op_fock);

  // now form the MO lagrangian
  //       c    o   v
  //  c  |2*FC|2*FC|0|
  //     -------------
  //  o  |2*FC| FO |0|
  //     -------------
  //  v  | 0  |  0 |0|
  //

  mofock.scale(2.0);

  for (i=_ndocc; i < _ndocc+_nsocc; i++)
    for (int j=_ndocc; j <= i; j++)
      mofock.set_element(i,j,moofock.get_element(i,j));

  for (i=_ndocc+_nsocc; i < basis()->nbasis(); i++)
    for (int j=0; j <= i; j++)
      mofock.set_element(i,j,0.0);

  moofock.assign(0.0);
  moofock.accumulate_transform(_gr_vector,mofock);
  moofock.scale(-1.0);

  // zero out gradient
  gradient.assign(0.0);

  // grab the centers struct
  centers_t *centers = basis()->convert_to_centers_t();
    
  // calculate nuclear repulsion contribution to the gradient
  double_vector_t dv;
  allocbn_double_vector(&dv,"n",3);

  int natom=centers->n;
  int nbasis=basis()->nbasis();

  for (int iatom=0; iatom < natom; iatom++) {
    int_nuclear_repulsion_1der(centers,centers,&dv,centers,iatom);
    for (int xyz=0; xyz < 3; xyz++)
      gradient.accumulate_element(iatom*3+xyz,dv.d[xyz]);
  }
  
  // gradient->print("nuclear repulsion terms");
  
  // now do the overlap contribution
  int_initialize_offsets1(centers,centers);
  double *oneebuff = int_initialize_1e(0,1,centers,centers);
  
  RefSCVector ovlp = gradient.clone();
  ovlp.assign(0.0);

  for (x=0; x < centers->n; x++) {
    for (int ish=0; ish < centers->nshell; ish++) {
      int istart = centers->func_num[ish];
      int iend = istart + INT_SH_NFUNC((centers),ish);
      
      for (int jsh=0; jsh <= ish; jsh++) {
        int jstart = centers->func_num[jsh];
        int jend = jstart + INT_SH_NFUNC((centers),jsh);

        int_shell_overlap_1der(centers,centers,oneebuff,ish,jsh,centers,x);

        zero_double_vector(&dv);
        
        int index=0;
        for (i=istart; i < iend; i++) {
          for (int j=jstart; j < jend; j++) {
            for (int k=0; k < 3; k++) {
              dv.d[k] += oneebuff[index] * moofock.get_element(i,j);
              index++;
            }
          }
        }
        if (ish!=jsh) {
          for (int k=0; k < 3; k++)
            dv.d[k] *= 2.0;
        }

        for (int k=0; k < 3; k++)
          ovlp.accumulate_element(x*3+k, dv.d[k]);
      }
    }
  }
    
  moofock=0;
  // ovlp.print("overlap contribution");
  gradient.accumulate(ovlp);
  
  // and now the one-electron contributions
  RefSCVector oneelec = ovlp;
  ovlp=0;
  oneelec.assign(0.0);

  // form density
  double pmax;
  gr_density(_gr_vector,_gr_dens,_gr_op_dens,_ndocc,_nsocc,pmax);

  mofock.assign(_gr_op_dens);
  mofock.scale(0.5);
  mofock.accumulate(_gr_dens);
  
  for (x=0; x < centers->n; x++) {
    for (int ish=0; ish < centers->nshell; ish++) {
      int istart = centers->func_num[ish];
      int iend = istart + INT_SH_NFUNC((centers),ish);
      
      for (int jsh=0; jsh <= ish; jsh++) {
        int jstart = centers->func_num[jsh];
        int jend = jstart + INT_SH_NFUNC((centers),jsh);

        int_shell_kinetic_1der(centers,centers,oneebuff,ish,jsh,centers,x);
        int_accum_shell_nuclear_hf_1der(centers,centers,oneebuff,
                                        ish,jsh,centers,x);
        int_accum_shell_nuclear_nonhf_1der(centers,centers,oneebuff,
                                           ish,jsh,centers,x);

        zero_double_vector(&dv);
        
        int index=0;
        for (i=istart; i < iend; i++) {
          for (int j=jstart; j < jend; j++) {
            for (int k=0; k < 3; k++) {
              dv.d[k] += oneebuff[index] * mofock.get_element(i,j);
              index++;
            }
          }
        }
        if (ish!=jsh) {
          for (int k=0; k < 3; k++)
            dv.d[k] *= 2.0;
        }

        for (int k=0; k < 3; k++)
          oneelec.accumulate_element(x*3+k, dv.d[k]);
      }
    }
  }
  mofock=0;

  // oneelec.print("one electron contribution");
  gradient.accumulate(oneelec);
  // gradient.print("gradient sans two electron contribution");
  
  // done with the one-electron stuff
  int_done_offsets1(centers,centers);
  int_done_1e();

  // now for the tricky part
  int flags = INT_EREP|INT_NOSTR1|INT_NOSTR2|INT_NOSTRB;
#if !BOUNDS
  flags |= INT_NODERB;
#endif

  int_initialize_offsets2(centers,centers,centers,centers);
  double *ints = int_initialize_erep(flags,1,centers,centers,centers,centers);
#if BOUNDS
  int_init_bounds_1der();
  int Pmax = int_bound_log(pmax);
  int PPmax = int_bound_log(6.0*pmax*pmax);
  int threshold = int_bound_log(_gradient.desired_accuracy());
#endif
  
  RefSCVector twoelec = oneelec;
  oneelec=0;
  twoelec.assign(0.0);
  
  double tnint=0;

  for (i=0; i < centers->nshell; i++) {
    for (int j=0; j <= i; j++) {

#if BOUNDS
      if (int_erep_2bound_1der(i,j)+PPmax < threshold)
        continue;
#endif

      for (int k=0; k <= i; k++) {
        for (int l=0; l <= ((i==k)?j:k); l++) {

#if BOUNDS
          if (int_erep_4bound_1der(i,j,k,l)+PPmax < threshold)
            continue;
#endif
          int sh[4], sz[4];
          double coulombscale,exchangescale;
          der_centers_t dercenters;
          
          sh[0]=i; sh[1]=j; sh[2]=k; sh[3]=l;

          int_erep_all1der_v(INT_EREP|INT_REDUND|INT_NOPERM,
                             sh,sz,&dercenters);

          tnint += (double) sz[0]*sz[1]*sz[2]*sz[3];

          set_scale(coulombscale,exchangescale,i,j,k,l);
          
          int indexijkl=0;
          for (int derset=0; derset < dercenters.n; derset++) {
            for (int xyz=0; xyz < 3; xyz++) {
              for (int ip=0; ip < sz[0]; ip++) {
                int io = ip + centers->func_num[i];
                for (int jp=0; jp < sz[1]; jp++) {
                  int jo = jp + centers->func_num[j];
                  for (int kp=0; kp < sz[2]; kp++) {
                    int ko = kp + centers->func_num[k];
                    for (int lp=0; lp < sz[3]; lp++) {
                      int lo = lp + centers->func_num[l];

                      if (fabs(ints[indexijkl]) < 1.0e-14) {
                        indexijkl++;
                        continue;
                      }

                      double contrib,contmp;

                      contrib=0;
                      contmp = coulombscale*ints[indexijkl];
                      contrib = alpha[0][0] * _gr_dens.get_element(io,jo)*
                                              _gr_dens.get_element(ko,lo)
                              + alpha[1][0] * _gr_op_dens.get_element(io,jo)*
                                              _gr_dens.get_element(ko,lo)
                              + alpha[0][1] * _gr_dens.get_element(io,jo)*
                                              _gr_op_dens.get_element(ko,lo)
                              + alpha[1][1] * _gr_op_dens.get_element(io,jo)*
                                              _gr_op_dens.get_element(ko,lo);
                      contrib *= contmp;

                      twoelec.accumulate_element(xyz+dercenters.num[derset]*3,
                                                 contrib);
                      twoelec.accumulate_element(xyz+dercenters.onum*3,
                                                 -contrib);
                      
                      contrib=0;
                      contmp = exchangescale*ints[indexijkl];
                      contrib = beta[0][0] * _gr_dens.get_element(io,ko)*
                                             _gr_dens.get_element(jo,lo)
                              + beta[1][0] * _gr_op_dens.get_element(io,ko)*
                                             _gr_dens.get_element(jo,lo)
                              + beta[0][1] * _gr_dens.get_element(io,ko)*
                                             _gr_op_dens.get_element(jo,lo)
                              + beta[1][1] * _gr_op_dens.get_element(io,ko)*
                                             _gr_op_dens.get_element(jo,lo);
                      contrib *= contmp;

                      twoelec.accumulate_element(xyz+dercenters.num[derset]*3,
                                                 contrib);
                      twoelec.accumulate_element(xyz+dercenters.onum*3,
                                                 -contrib);

                      if (i!=j && k!=l) {
                        contrib=0;
                        contrib = beta[0][0] * _gr_dens.get_element(io,lo)*
                                               _gr_dens.get_element(jo,ko)
                                + beta[1][0] * _gr_op_dens.get_element(io,lo)*
                                               _gr_dens.get_element(jo,ko)
                                + beta[0][1] * _gr_dens.get_element(io,lo)*
                                               _gr_op_dens.get_element(jo,ko)
                                + beta[1][1] * _gr_op_dens.get_element(io,lo)*
                                               _gr_op_dens.get_element(jo,ko);
                        contrib *= contmp;

                        twoelec.accumulate_element(
                                                 xyz+dercenters.num[derset]*3,
                                                 contrib);
                        twoelec.accumulate_element(xyz+dercenters.onum*3,
                                                 -contrib);
                      }

                      indexijkl++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // twoelec.print("two electron contribution");
  gradient.accumulate(twoelec);
  //gradient.print("cartesian gradient");

  printf("%20.0f derivative integrals\n",tnint);

#if BOUNDS
  int_done_bounds_1der();
#endif
  int_done_offsets2(centers,centers,centers,centers);
  int_done_erep();
  
  // clean up some things
  _gr_dens = 0;
  _gr_op_dens = 0;
  _gr_vector = 0;

  free_double_vector(&dv);
  free_centers(centers);
  free(centers);
}
