
#define BOUNDS 1

#include <util/misc/newstring.h>
#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/scf/clscf.h>

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
  exchangescale = -0.25*scale;

  if (k!=l)
    coulombscale *= 2.0;

  if ((k!=l)&&(i==j))
    exchangescale *= 2.0;
}

static void
gr_density(const RefSCMatrix& vec, const RefSymmSCMatrix& dens, int ndocc,
           double& pmax)
{
  pmax=0.0;

  for (int i=0; i < vec->nrow(); i++) {
    for (int j=0; j <= i; j++) {
      double pt=0;
      for (int k=0; k < ndocc; k++)
        pt += vec->get_element(i,k)*vec->get_element(j,k);
      
      dens->set_element(i,j,pt);
      if (fabs(pt)>pmax)
        pmax=fabs(pt);
    }
  }
  dens->scale(2.0);
}

static void
ew_density(const RefSCMatrix& vec, const RefDiagSCMatrix& evals,
           const RefSymmSCMatrix& ewdens, int ndocc)
{
  for (int i=0; i < vec->nrow(); i++) {
    for (int j=0; j <= i; j++) {
      double pt=0;
      for (int k=0; k < ndocc; k++)
        pt += vec->get_element(i,k)*vec->get_element(j,k)*
          evals->get_element(k);
      
      ewdens->set_element(i,j,pt);
    }
  }
  ewdens->scale(-2.0);
}

void
CLSCF::do_gradient(const RefSCVector& gradient)
{
  // grab a reference to the scf_vector, presumably it is current
  _gr_vector = _eigenvectors.result_noupdate();
  
  // allocate storage for the temp arrays
  _gr_dens = _fock.clone();
  
  // form energy weighted density
  ew_density(_gr_vector,_fock_evals,_gr_dens,_ndocc);
  
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
  
  //gradient->print("nuclear repulsion terms");
  
  // now do the overlap contribution
  int_initialize_offsets1(centers,centers);
  double *oneebuff = int_initialize_1e(0,1,centers,centers);
  
  RefSCVector ovlp = gradient.clone();
  ovlp.assign(0.0);

  for (int x=0; x < centers->n; x++) {
    for (int ish=0; ish < centers->nshell; ish++) {
      int istart = centers->func_num[ish];
      int iend = istart + INT_SH_NFUNC((centers),ish);
      
      for (int jsh=0; jsh <= ish; jsh++) {
        int jstart = centers->func_num[jsh];
        int jend = jstart + INT_SH_NFUNC((centers),jsh);

        int_shell_overlap_1der(centers,centers,oneebuff,ish,jsh,centers,x);

        zero_double_vector(&dv);
        
        int index=0;
        for (int i=istart; i < iend; i++) {
          for (int j=jstart; j < jend; j++) {
            for (int k=0; k < 3; k++) {
              dv.d[k] += oneebuff[index] * _gr_dens.get_element(i,j);
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
    
  //ovlp.print("overlap contribution");
  gradient.accumulate(ovlp);
  
  // and now the one-electron contributions
  RefSCVector oneelec = ovlp;
  oneelec.assign(0.0);

  // form density
  double pmax;
  gr_density(_gr_vector,_gr_dens,_ndocc,pmax);

  for (int x=0; x < centers->n; x++) {
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
        for (int i=istart; i < iend; i++) {
          for (int j=jstart; j < jend; j++) {
            for (int k=0; k < 3; k++) {
              dv.d[k] += oneebuff[index] * _gr_dens.get_element(i,j);
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

  //oneelec.print("one electron contribution");
  gradient.accumulate(oneelec);
  //gradient.print("gradient sans two electron contribution");
  
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
  twoelec.assign(0.0);
  
  double tnint=0;

  PetiteList pl(basis());

  for (int i=0; i < centers->nshell; i++) {
    if (!pl.in_p1(i))
      continue;
    
    for (int j=0; j <= i; j++) {
      int ij = (i*(i+1) >> 1) + j;
      if (!pl.in_p2(ij))
        continue;

#if BOUNDS
      if (int_erep_2bound_1der(i,j)+PPmax < threshold)
        continue;
#endif

      for (int k=0; k <= i; k++) {
        for (int l=0; l <= ((i==k)?j:k); l++) {
          int kl = (k*(k+1)>>1) + l;
          int qijkl;
          if (!(qijkl=pl.in_p4(ij,kl,i,j,k,l)))
            continue;
          
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

                      double contrib;

                      contrib = coulombscale*ints[indexijkl]*qijkl*
                                             _gr_dens.get_element(io,jo)*
                                             _gr_dens.get_element(ko,lo);

                      twoelec.accumulate_element(xyz+dercenters.num[derset]*3,
                                                 contrib);
                      twoelec.accumulate_element(xyz+dercenters.onum*3,
                                                 -contrib);
                      
                      contrib = exchangescale*ints[indexijkl]*qijkl*
                                              _gr_dens.get_element(io,ko)*
                                              _gr_dens.get_element(jo,lo);
                      twoelec.accumulate_element(xyz+dercenters.num[derset]*3,
                                                 contrib);
                      twoelec.accumulate_element(xyz+dercenters.onum*3,
                                                 -contrib);

                      if (i!=j && k!=l) {
                        contrib = exchangescale*ints[indexijkl]*qijkl*
                                              _gr_dens.get_element(io,lo)*
                                              _gr_dens.get_element(jo,ko);
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

  twoelec.print("two electron contribution");

  RefSCVector sym2ei = twoelec.copy();
  CharacterTable ct = _mol->point_group().char_table();
  SymmetryOperation so;
  
  for (int alpha=0; alpha < centers->n; alpha++) {
    for (int g=1; g < ct.order(); g++) {
      so = ct.symm_operation(g);
      int ap = pl.atom_map(alpha,g);

      sym2ei.accumulate_element(alpha*3+0,
                               twoelec.get_element(ap*3+0)*so(0,0) +
                               twoelec.get_element(ap*3+1)*so(1,0) +
                               twoelec.get_element(ap*3+2)*so(2,0));
      sym2ei.accumulate_element(alpha*3+1,
                               twoelec.get_element(ap*3+0)*so(0,1) +
                               twoelec.get_element(ap*3+1)*so(1,1) +
                               twoelec.get_element(ap*3+2)*so(2,1));
      sym2ei.accumulate_element(alpha*3+2,
                               twoelec.get_element(ap*3+0)*so(0,2) +
                               twoelec.get_element(ap*3+1)*so(1,2) +
                               twoelec.get_element(ap*3+2)*so(2,2));
    }
  }
    
  sym2ei.scale(1.0/ct.order());
  sym2ei.print("symmetrized two electron contribution");

  gradient.accumulate(sym2ei);
  //gradient.print("cartesian gradient");

  printf("%20.0f derivative integrals\n",tnint);

#if BOUNDS
  int_done_bounds_1der();
#endif
  int_done_offsets2(centers,centers,centers,centers);
  int_done_erep();
  
  // clean up some things
  _gr_dens = 0;

  free_double_vector(&dv);
  free_centers(centers);
  free(centers);
}
