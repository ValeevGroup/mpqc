
#ifndef _chemistry_qc_scf_ltbgrad_h
#define _chemistry_qc_scf_ltbgrad_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math.h>

#include <util/misc/timer.h>
#include <math/scmat/offset.h>

#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/scf/tbgrad.h>
  
template<class T>
class LocalTBGrad : public TBGrad<T> {
  protected:
    RefMessageGrp grp_;
    RefTwoBodyDerivInt tbi_;
    RefIntegral integral_;
    RefGaussianBasisSet gbs_;

  public:
    LocalTBGrad(T& t, const RefIntegral& ints, const RefGaussianBasisSet& bs,
                const RefMessageGrp& g) :
      TBGrad<T>(t), grp_(g), integral_(ints), gbs_(bs)
    {
      tbi_ = integral_->electron_repulsion_deriv();
    }

    ~LocalTBGrad() {}
    
    void build_tbgrad(double * tbgrad, double pmax, double accuracy) {
      tim_enter("two electron gradient");

      int me = grp_->me();
      int nproc = grp_->n();
      
      // grab ref for convenience
      GaussianBasisSet& gbs = *gbs_.pointer();
      Molecule& mol = *gbs.molecule().pointer();
      RefPetiteList rpl = integral_->petite_list();
      PetiteList& pl = *rpl.pointer();
      TwoBodyDerivInt& tbi = *tbi_.pointer();
      
      // create vector to hold skeleton gradient
      double *tbint = new double[mol.natom()*3];
      memset(tbint, 0, sizeof(double)*mol.natom()*3);

      // for bounds checking
      int PPmax = (int) (log(6.0*pmax*pmax)/log(2.0));
      int threshold = (int) (log(accuracy)/log(2.0));
  
      int kindex=0;
      for (int i=0; i < gbs.nshell(); i++) {
        if (!pl.in_p1(i))
          continue;
    
        int ni=gbs(i).nfunction();
        int fi=gbs.shell_to_function(i);
    
        for (int j=0; j <= i; j++) {
          int ij=i_offset(i)+j;
          if (!pl.in_p2(ij))
            continue;
      
          if (tbi.log2_shell_bound(i,j,-1,-1)+PPmax < threshold)
            continue;
      
          int nj=gbs(j).nfunction();
          int fj=gbs.shell_to_function(j);
    
          for (int k=0; k <= i; k++,kindex++) {
            if (kindex%nproc != me)
              continue;
            
            int nk=gbs(k).nfunction();
            int fk=gbs.shell_to_function(k);
    
            for (int l=0; l <= ((i==k)?j:k); l++) {
              if (tbi.log2_shell_bound(i,j,k,l)+PPmax < threshold)
                continue;
          
              int kl=i_offset(k)+l;
              int qijkl;
              if (!(qijkl=pl.in_p4(ij,kl,i,j,k,l)))
                continue;
          
              int nl=gbs(l).nfunction();
              int fl=gbs.shell_to_function(l);

              DerivCenters cent;
              tim_enter("quartet");
              tbi.compute_shell(i,j,k,l,cent);
              tim_exit("quartet");

              const double * buf = tbi.buffer();
          
              double cscl, escl;

              set_scale(cscl, escl, i, j, k, l);

              int indijkl=0;
              int nx=cent.n();
              if (cent.has_omitted_center()) nx--;
              for (int x=0; x < nx; x++) {
                int ix=cent.atom(x);
                int io=cent.omitted_atom();
                for (int ixyz=0; ixyz < 3; ixyz++) {
                  for (int ip=0, ii=fi; ip < ni; ip++, ii++) {
                    for (int jp=0, jj=fj; jp < nj; jp++, jj++) {
                      for (int kp=0, kk=fk; kp < nk; kp++, kk++) {
                        for (int lp=0, ll=fl; lp < nl; lp++, ll++, indijkl++) {
                          double contrib;
                          double qint = buf[indijkl]*qijkl;

                          contrib = cscl*qint*
                            contribution.cont1(ij_offset(ii,jj),
                                               ij_offset(kk,ll));

                          tbint[ixyz+ix*3] += contrib;
                          tbint[ixyz+io*3] -= contrib;

                          contrib = escl*qint*
                            contribution.cont2(ij_offset(ii,kk),
                                               ij_offset(jj,ll));

                          tbint[ixyz+ix*3] += contrib;
                          tbint[ixyz+io*3] -= contrib;

                          if (i!=j && k!=l) {
                            contrib = escl*qint*
                              contribution.cont2(ij_offset(ii,ll),
                                                 ij_offset(jj,kk));

                            tbint[ixyz+ix*3] += contrib;
                            tbint[ixyz+io*3] -= contrib;
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
      }
      
      grp_->sum(tbint, mol.natom()*3);
      
      CharacterTable ct = mol.point_group().char_table();
      SymmetryOperation so;

      for (int alpha=0; alpha < mol.natom(); alpha++) {
        double tbx = tbint[alpha*3+0];
        double tby = tbint[alpha*3+1];
        double tbz = tbint[alpha*3+2];
        
        for (int g=1; g < ct.order(); g++) {
          so = ct.symm_operation(g);
          int ap = pl.atom_map(alpha,g);

          tbx += tbint[ap*3+0]*so(0,0) + tbint[ap*3+1]*so(1,0) +
                 tbint[ap*3+2]*so(2,0);
          tby += tbint[ap*3+0]*so(0,1) + tbint[ap*3+1]*so(1,1) +
                 tbint[ap*3+2]*so(2,1);
          tbz += tbint[ap*3+0]*so(0,2) + tbint[ap*3+1]*so(1,2) +
                 tbint[ap*3+2]*so(2,2);
        }
        double scl = 1.0/(double)ct.order();
        tbgrad[alpha*3+0] += tbx*scl;
        tbgrad[alpha*3+1] += tby*scl;
        tbgrad[alpha*3+2] += tbz*scl;
      }
    
      delete[] tbint;
      
      tim_exit("two electron gradient");
    }
};

#endif
