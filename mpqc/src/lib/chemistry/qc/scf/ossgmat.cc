
#include <math/array/math_lib.h>
#include <math/scmat/local.h>
#include <chemistry/qc/integral/integralv2.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/scf/ossscf.h>

#define ioff(i) (((i)*((i)+1))>>1)
#define IOFF(a,b) (((a)>(b))?(ioff(a)+(b)):(ioff(b)+(a)))

static signed char
max_den(centers_t *centers, const RefSymmSCMatrix& dens, int s1, int s2)
{
  double t,max=0.0;

  int s1i = centers->func_num[s1];
  int s1u = s1i + INT_SH_NFUNC((centers),s1);

  int s2i = centers->func_num[s2];
  int s2u = s2i + INT_SH_NFUNC((centers),s1);

  for (int i=s1i; i < s1u; i++) {
    for (int j=s2i; j < s2u; j++) {
      if ((t=fabs(dens.get_element(i,j))) > max)
        max=t;
    }
  }

  if (max < 1.0e-15) max=1.0e-15;
  int ret = (int) (log(max)/log(2.0));
  if (ret < -126) ret=-126;

  return (signed char) ret;
}

void
OSSSCF::form_density(const RefSCMatrix& vec,
                      const RefSymmSCMatrix& density,
                      const RefSymmSCMatrix& density_diff,
                      const RefSymmSCMatrix& opena_density,
                      const RefSymmSCMatrix& opena_density_diff,
                      const RefSymmSCMatrix& openb_density,
                      const RefSymmSCMatrix& openb_density_diff)
{
  int nbasis = basis()->nbasis();

  // find out what type of matrices we're dealing with
  if (LocalSCMatrix::castdown(vec.pointer())) {
    LocalSCMatrix *lvec = LocalSCMatrix::require_castdown(
      vec.pointer(), "OSSSCF::form_density");
    LocalSymmSCMatrix *ldens = LocalSymmSCMatrix::require_castdown(
      density.pointer(), "OSSSCF::form_density");
    LocalSymmSCMatrix *ldensd = LocalSymmSCMatrix::require_castdown(
      density_diff.pointer(), "OSSSCF::form_density");
    LocalSymmSCMatrix *loadens = LocalSymmSCMatrix::require_castdown(
      opena_density.pointer(), "OSSSCF::form_density");
    LocalSymmSCMatrix *loadensd = LocalSymmSCMatrix::require_castdown(
      opena_density_diff.pointer(), "OSSSCF::form_density");
    LocalSymmSCMatrix *lobdens = LocalSymmSCMatrix::require_castdown(
      openb_density.pointer(), "OSSSCF::form_density");
    LocalSymmSCMatrix *lobdensd = LocalSymmSCMatrix::require_castdown(
      openb_density_diff.pointer(), "OSSSCF::form_density");

    for (int i=0; i < nbasis; i++) {
      for (int j=0; j <= i; j++) {
        double pt=0;
        for (int k=0; k < _ndocc; k++)
          pt += 2.0*lvec->get_element(i,k)*lvec->get_element(j,k);

        double ptoa=0;
        for (int k=_ndocc; k < _ndocc+1; k++)
          ptoa += lvec->get_element(i,k)*lvec->get_element(j,k);

        if (loadensd)
          loadensd->set_element(i,j,ptoa-loadens->get_element(i,j));
        loadens->set_element(i,j,ptoa);

        double ptob=0;
        for (int k=_ndocc+1; k < _ndocc+2; k++)
          ptob += lvec->get_element(i,k)*lvec->get_element(j,k);

        if (lobdensd)
          lobdensd->set_element(i,j,ptob-lobdens->get_element(i,j));
        lobdens->set_element(i,j,ptob);

        if (ldensd)
          ldensd->set_element(i,j,pt+ptoa+ptob-ldens->get_element(i,j));
        ldens->set_element(i,j,pt+ptoa+ptob);
      }
    }
  }
}

void
OSSSCF::form_ao_fock(centers_t *centers, double *intbuf)
{
  int inttol = int_bound_log(_energy.desired_accuracy()/100.0);

  char *shnfunc = new char[centers->nshell];
  for (int i=0; i < centers->nshell; i++)
    shnfunc[i] = INT_SH_NFUNC((centers),i);

  signed char *pmax = new signed char[ioff(centers->nshell)];
  for (int i=0; i < centers->nshell; i++) {
    int ij=ioff(i);
    for (int j=0; j <= i; j++,ij++) {
      pmax[ij] = max_den(centers,_gr_dens_diff,i,j);
    }
  }
  
  double tnint=0;
  
  for (int i=0; i < centers->nshell; i++) {
    for (int j=0; j <= i; j++) {
      int ij = ioff(i)+j;
      int ijbnd = int_Qvec[ij];
      int ijpmx = pmax[ij];

      for (int k=0; k <= i; k++) {
        int ijkpmx=ijpmx;
        if (pmax[(ioff(i)+k)]-2 > ijkpmx)
          ijkpmx = pmax[(ioff(i)+k)]-2;

        if (pmax[IOFF(j,k)]-2 > ijkpmx)
          ijkpmx = pmax[IOFF(j,k)]-2;
        
        for (int l=0; l <= ((k==i)?j:k); l++) {
          int kl = ioff(k)+l;
          int klbnd = int_Qvec[kl];
          int ijklpmx = (pmax[kl]>ijkpmx) ? pmax[kl] : ijkpmx;

          if (pmax[(ioff(i)+l)]-2 > ijklpmx)
            ijklpmx = pmax[(ioff(i)+l)]-2;

          if (pmax[IOFF(j,l)]-2 > ijklpmx)
            ijklpmx = pmax[IOFF(j,k)]-2;

          if (ijklpmx+ijbnd+klbnd < inttol)
            continue;

          int s1=i, s2=j, s3=k, s4=l;

          int_erep(INT_EREP|INT_NOBCHK|INT_NOPERM,&s1,&s2,&s3,&s4);

          int n1 = shnfunc[s1];
          int n2 = shnfunc[s2];
          int n3 = shnfunc[s3];
          int n4 = shnfunc[s4];

          int e12 = (s2==s1);
          int e13e24 = (s3==s1) && (s4==s2);
          int e34 = (s4 == s3);
          int e_any = (e12||e13e24||e34);
          
          int index=0;
          
          if (e_any) {
            for (int bf1=0; bf1<=INT_MAX1(n1) ; bf1++) {
              int i1 = centers->func_num[s1] + bf1;

              for (int bf2=0; bf2<=INT_MAX2(e12,bf1,n2) ; bf2++) {
                int j1 = centers->func_num[s2] + bf2;
                int ij1=ioff(i1)+j1;

                for (int bf3=0; bf3<=INT_MAX3(e13e24,bf1,n3) ; bf3++) {
                  int k1 = centers->func_num[s3] + bf3;

                  for (int bf4=0;
                       bf4<=INT_MAX4(e13e24,e34,bf1,bf2,bf3,n4);bf4++){
                    if (INT_NONZERO(intbuf[index])) {
                      int l1 = centers->func_num[s4] + bf4;

                      int ii,jj,kk,ll;
                      if (ij1 >= ioff(k1)+l1) {
                        ii = i1; jj = j1; kk = k1; ll = l1;
                      } else {
                        ii = k1; jj = l1; kk = i1; ll = j1;
                      }

                      double pki_int = intbuf[index];
                      double pkval;

                      int lij,lkl;
                      
                      if (jj == kk) {
                        /*
                         * if i=j=k or j=k=l, then this integral contributes
                         * to J, K1, and K2 of G(ij), so
                         * pkval = (ijkl) - 0.25 * ((ikjl)-(ilkj))
                         *       = 0.5 * (ijkl)
                         */
                        if (ii == jj || kk == ll) {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          pkval = (lij==lkl) ? 0.25*pki_int: 0.5*pki_int;

                          _gr_gmat.accumulate_element(ii,jj,pkval*
                                             _gr_dens_diff.get_element(kk,ll));
                          _gr_gmat.accumulate_element(kk,ll,pkval*
                                             _gr_dens_diff.get_element(ii,jj));

                          _gr_opa_gmat.accumulate_element(ii,jj,pkval*
                                         _gr_opa_dens_diff.get_element(kk,ll));
                          _gr_opa_gmat.accumulate_element(kk,ll,pkval*
                                         _gr_opa_dens_diff.get_element(ii,jj));

                          _gr_opa_gmat.accumulate_element(ii,jj,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(kk,ll));
                          _gr_opa_gmat.accumulate_element(kk,ll,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(ii,jj));

                          _gr_opb_gmat.accumulate_element(ii,jj,pkval*
                                         _gr_opb_dens_diff.get_element(kk,ll));
                          _gr_opb_gmat.accumulate_element(kk,ll,pkval*
                                         _gr_opb_dens_diff.get_element(ii,jj));

                          _gr_opb_gmat.accumulate_element(ii,jj,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(kk,ll));
                          _gr_opb_gmat.accumulate_element(kk,ll,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(ii,jj));

                        } else {
                          /*
                           * if j=k, then this integral contributes
                           * to J and K1 of G(ij)
                           *
                           * pkval = (ijkl) - 0.25 * (ikjl)
                           *       = 0.75 * (ijkl)
                           */
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          pkval = (lij==lkl) ? 0.375*pki_int: 0.75*pki_int;

                          _gr_gmat.accumulate_element(ii,jj,pkval*
                                             _gr_dens_diff.get_element(kk,ll));
                          _gr_gmat.accumulate_element(kk,ll,pkval*
                                             _gr_dens_diff.get_element(ii,jj));

                          pkval *= 0.33333333333333333;
                          _gr_opa_gmat.accumulate_element(ii,jj,pkval*
                                         _gr_opa_dens_diff.get_element(kk,ll));
                          _gr_opa_gmat.accumulate_element(kk,ll,pkval*
                                         _gr_opa_dens_diff.get_element(ii,jj));

                          _gr_opa_gmat.accumulate_element(ii,jj,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(kk,ll));
                          _gr_opa_gmat.accumulate_element(kk,ll,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(ii,jj));

                          _gr_opb_gmat.accumulate_element(ii,jj,pkval*
                                         _gr_opb_dens_diff.get_element(kk,ll));
                          _gr_opb_gmat.accumulate_element(kk,ll,pkval*
                                         _gr_opb_dens_diff.get_element(ii,jj));

                          _gr_opb_gmat.accumulate_element(ii,jj,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(kk,ll));
                          _gr_opb_gmat.accumulate_element(kk,ll,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(ii,jj));

                          /*
                           * this integral also contributes to K1 and K2 of
                           * G(il)
                           *
                           * pkval = -0.25 * ((ijkl)+(ikjl))
                           *       = -0.5 * (ijkl)
                           */
                          lij=ioff(ii)+ll;
                          lkl=IOFF(kk,jj);
                          pkval = (lij==lkl)? 0.25*pki_int: 0.5*pki_int;

                          _gr_gmat.accumulate_element(ii,ll,-pkval*
                                             _gr_dens_diff.get_element(kk,jj));
                          _gr_gmat.accumulate_element(kk,jj,-pkval*
                                             _gr_dens_diff.get_element(ii,ll));

                          _gr_opa_gmat.accumulate_element(ii,ll,pkval*
                                         _gr_opa_dens_diff.get_element(kk,jj));
                          _gr_opa_gmat.accumulate_element(kk,jj,pkval*
                                         _gr_opa_dens_diff.get_element(ii,ll));

                          _gr_opa_gmat.accumulate_element(ii,ll,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(kk,jj));
                          _gr_opa_gmat.accumulate_element(kk,jj,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(ii,ll));

                          _gr_opb_gmat.accumulate_element(ii,ll,pkval*
                                         _gr_opb_dens_diff.get_element(kk,jj));
                          _gr_opb_gmat.accumulate_element(kk,jj,pkval*
                                         _gr_opb_dens_diff.get_element(ii,ll));

                          _gr_opb_gmat.accumulate_element(ii,ll,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(kk,jj));
                          _gr_opb_gmat.accumulate_element(kk,jj,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(ii,ll));
                        }
                      } else if (ii == kk || jj == ll) {
                        /*
                         * if i=k or j=l, then this integral contributes
                         * to J and K2 of G(ij)
                         *
                         * pkval = (ijkl) - 0.25 * (ilkj)
                         *       = 0.75 * (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        pkval = (lij==lkl) ? 0.375*pki_int: 0.75*pki_int;
                        _gr_gmat.accumulate_element(ii,jj,pkval*
                                             _gr_dens_diff.get_element(kk,ll));
                        _gr_gmat.accumulate_element(kk,ll,pkval*
                                             _gr_dens_diff.get_element(ii,jj));

                        pkval *= 0.3333333333333333;
                        _gr_opa_gmat.accumulate_element(ii,jj,pkval*
                                         _gr_opa_dens_diff.get_element(kk,ll));
                        _gr_opa_gmat.accumulate_element(kk,ll,pkval*
                                         _gr_opa_dens_diff.get_element(ii,jj));

                        _gr_opa_gmat.accumulate_element(ii,jj,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(kk,ll));
                        _gr_opa_gmat.accumulate_element(kk,ll,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(ii,jj));

                        _gr_opb_gmat.accumulate_element(ii,jj,pkval*
                                         _gr_opb_dens_diff.get_element(kk,ll));
                        _gr_opb_gmat.accumulate_element(kk,ll,pkval*
                                         _gr_opb_dens_diff.get_element(ii,jj));

                        _gr_opb_gmat.accumulate_element(ii,jj,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(kk,ll));
                        _gr_opb_gmat.accumulate_element(kk,ll,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(ii,jj));

                        /*
                         * this integral also contributes to K1 and K2 of
                         * G(ik)
                         *
                         * pkval = -0.25 * ((ijkl)+(ilkj))
                         *       = -0.5 * (ijkl)
                         */
                        lij=ioff(ii)+kk;
                        lkl=IOFF(jj,ll);
                        pkval = (lij==lkl) ? 0.25*pki_int : 0.5*pki_int;

                        _gr_gmat.accumulate_element(ii,kk,-pkval*
                                             _gr_dens_diff.get_element(jj,ll));
                        _gr_gmat.accumulate_element(jj,ll,-pkval*
                                             _gr_dens_diff.get_element(ii,kk));

                        _gr_opa_gmat.accumulate_element(ii,kk,pkval*
                                         _gr_opa_dens_diff.get_element(jj,ll));
                        _gr_opa_gmat.accumulate_element(jj,ll,pkval*
                                         _gr_opa_dens_diff.get_element(ii,kk));

                        _gr_opa_gmat.accumulate_element(ii,kk,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(jj,ll));
                        _gr_opa_gmat.accumulate_element(jj,ll,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(ii,kk));

                        _gr_opb_gmat.accumulate_element(ii,kk,pkval*
                                         _gr_opb_dens_diff.get_element(jj,ll));
                        _gr_opb_gmat.accumulate_element(jj,ll,pkval*
                                         _gr_opb_dens_diff.get_element(ii,kk));

                        _gr_opb_gmat.accumulate_element(ii,kk,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(jj,ll));
                        _gr_opb_gmat.accumulate_element(jj,ll,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(ii,kk));
                      } else {
                        /*
                         * This integral contributes to J of G(ij)
                         *
                         * pkval = (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        pkval = (lij==lkl)? 0.5*pki_int : pki_int;

                        _gr_gmat.accumulate_element(ii,jj,pkval*
                                             _gr_dens_diff.get_element(kk,ll));
                        _gr_gmat.accumulate_element(kk,ll,pkval*
                                             _gr_dens_diff.get_element(ii,jj));

                        /*
                         * and to K1 of G(ik)
                         *
                         * pkval = -0.25 * (ijkl)
                         */
                        lij=ioff(ii)+kk;
                        lkl=IOFF(jj,ll);
                        pkval = (lij==lkl) ? 0.125*pki_int : 0.25*pki_int;

                        _gr_gmat.accumulate_element(ii,kk,-pkval*
                                             _gr_dens_diff.get_element(jj,ll));
                        _gr_gmat.accumulate_element(jj,ll,-pkval*
                                             _gr_dens_diff.get_element(ii,kk));

                        _gr_opa_gmat.accumulate_element(ii,kk,pkval*
                                         _gr_opa_dens_diff.get_element(jj,ll));
                        _gr_opa_gmat.accumulate_element(jj,ll,pkval*
                                         _gr_opa_dens_diff.get_element(ii,kk));

                        _gr_opa_gmat.accumulate_element(ii,kk,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(jj,ll));
                        _gr_opa_gmat.accumulate_element(jj,ll,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(ii,kk));

                        _gr_opb_gmat.accumulate_element(ii,kk,pkval*
                                         _gr_opb_dens_diff.get_element(jj,ll));
                        _gr_opb_gmat.accumulate_element(jj,ll,pkval*
                                         _gr_opb_dens_diff.get_element(ii,kk));

                        _gr_opb_gmat.accumulate_element(ii,kk,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(jj,ll));
                        _gr_opb_gmat.accumulate_element(jj,ll,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(ii,kk));

                        if ((ii != jj) && (kk != ll)) {
                          /*
                           * if i!=j and k!=l, then this integral also
                           * contributes to K2 of G(il)
                           *
                           * pkval = -0.25 * (ijkl)
                           *
                           * note: if we get here, then ik can't equal jl,
                           * so pkval wasn't multiplied by 0.5 above.
                           */
                          lij=ioff(ii)+ll;
                          lkl=IOFF(kk,jj);

                          _gr_gmat.accumulate_element(ii,ll,-pkval*
                                             _gr_dens_diff.get_element(kk,jj));
                          _gr_gmat.accumulate_element(kk,jj,-pkval*
                                             _gr_dens_diff.get_element(ii,ll));

                          _gr_opa_gmat.accumulate_element(ii,ll,pkval*
                                         _gr_opa_dens_diff.get_element(kk,jj));
                          _gr_opa_gmat.accumulate_element(kk,jj,pkval*
                                         _gr_opa_dens_diff.get_element(ii,ll));

                          _gr_opa_gmat.accumulate_element(ii,ll,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(kk,jj));
                          _gr_opa_gmat.accumulate_element(kk,jj,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(ii,ll));

                          _gr_opb_gmat.accumulate_element(ii,ll,pkval*
                                         _gr_opb_dens_diff.get_element(kk,jj));
                          _gr_opb_gmat.accumulate_element(kk,jj,pkval*
                                         _gr_opb_dens_diff.get_element(ii,ll));

                          _gr_opb_gmat.accumulate_element(ii,ll,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(kk,jj));
                          _gr_opb_gmat.accumulate_element(kk,jj,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(ii,ll));
                        }
                      }
                    }
                    index++;
                  }
                }
              }
            }
          } else {
            for (int i1=centers->func_num[s1], bf1=0; bf1<n1 ; bf1++, i1++) {
              for (int j1=centers->func_num[s2], bf2=0; bf2<n2 ; bf2++, j1++) {
                int ij1=ioff(i1)+j1;

                for (int k1=centers->func_num[s3],bf3=0; bf3<n3; bf3++,k1++) {
                  for (int l1=centers->func_num[s4],bf4=0;bf4<n4;bf4++,l1++) {
                    if (INT_NONZERO(intbuf[index])) {

                      int ii,jj,kk,ll;
                      if (ij1 >= ioff(k1)+l1) {
                        ii = i1; jj = j1; kk = k1; ll = l1;
                      } else {
                        ii = k1; jj = l1; kk = i1; ll = j1;
                      }

                      double pki_int = intbuf[index];
                      double pkval;
                      int lij,lkl;

                      if (jj == kk) {
                        /*
                         * if j=k, then this integral contributes
                         * to J and K1 of G(ij)
                         *
                         * pkval = (ijkl) - 0.25 * (ikjl)
                         *       = 0.75 * (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        pkval = 0.75*pki_int;

                        _gr_gmat.accumulate_element(ii,jj,pkval*
                                             _gr_dens_diff.get_element(kk,ll));
                        _gr_gmat.accumulate_element(kk,ll,pkval*
                                             _gr_dens_diff.get_element(ii,jj));

                        pkval *= 0.3333333333333333;
                        _gr_opa_gmat.accumulate_element(ii,jj,pkval*
                                         _gr_opa_dens_diff.get_element(kk,ll));
                        _gr_opa_gmat.accumulate_element(kk,ll,pkval*
                                         _gr_opa_dens_diff.get_element(ii,jj));

                        _gr_opa_gmat.accumulate_element(ii,jj,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(kk,ll));
                        _gr_opa_gmat.accumulate_element(kk,ll,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(ii,jj));

                        _gr_opb_gmat.accumulate_element(ii,jj,pkval*
                                         _gr_opb_dens_diff.get_element(kk,ll));
                        _gr_opb_gmat.accumulate_element(kk,ll,pkval*
                                         _gr_opb_dens_diff.get_element(ii,jj));

                        _gr_opb_gmat.accumulate_element(ii,jj,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(kk,ll));
                        _gr_opb_gmat.accumulate_element(kk,ll,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(ii,jj));

                        /*
                         * this integral also contributes to K1 and K2 of
                         * G(il)
                         *
                         * pkval = -0.25 * ((ijkl)+(ikjl))
                         *       = -0.5 * (ijkl)
                         */
                        lij=ioff(ii)+ll;
                        lkl=IOFF(kk,jj);
                        pkval = 0.5 * pki_int;

                        _gr_gmat.accumulate_element(ii,ll,-pkval*
                                             _gr_dens_diff.get_element(kk,jj));
                        _gr_gmat.accumulate_element(kk,jj,-pkval*
                                             _gr_dens_diff.get_element(ii,ll));

                        _gr_opa_gmat.accumulate_element(ii,ll,pkval*
                                         _gr_opa_dens_diff.get_element(kk,jj));
                        _gr_opa_gmat.accumulate_element(kk,jj,pkval*
                                         _gr_opa_dens_diff.get_element(ii,ll));

                        _gr_opa_gmat.accumulate_element(ii,ll,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(kk,jj));
                        _gr_opa_gmat.accumulate_element(kk,jj,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(ii,ll));

                        _gr_opb_gmat.accumulate_element(ii,ll,pkval*
                                         _gr_opb_dens_diff.get_element(kk,jj));
                        _gr_opb_gmat.accumulate_element(kk,jj,pkval*
                                         _gr_opb_dens_diff.get_element(ii,ll));

                        _gr_opb_gmat.accumulate_element(ii,ll,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(kk,jj));
                        _gr_opb_gmat.accumulate_element(kk,jj,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(ii,ll));

                      } else if (ii == kk || jj == ll) {
                        /*
                         * if i=k or j=l, then this integral contributes
                         * to J and K2 of G(ij)
                         *
                         * pkval = (ijkl) - 0.25 * (ilkj)
                         *       = 0.75 * (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        pkval = 0.75*pki_int;

                        _gr_gmat.accumulate_element(ii,jj,pkval*
                                             _gr_dens_diff.get_element(kk,ll));
                        _gr_gmat.accumulate_element(kk,ll,pkval*
                                             _gr_dens_diff.get_element(ii,jj));

                        pkval *= 0.3333333333333333;
                        _gr_opa_gmat.accumulate_element(ii,jj,pkval*
                                         _gr_opa_dens_diff.get_element(kk,ll));
                        _gr_opa_gmat.accumulate_element(kk,ll,pkval*
                                         _gr_opa_dens_diff.get_element(ii,jj));

                        _gr_opa_gmat.accumulate_element(ii,jj,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(kk,ll));
                        _gr_opa_gmat.accumulate_element(kk,ll,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(ii,jj));

                        _gr_opb_gmat.accumulate_element(ii,jj,pkval*
                                         _gr_opb_dens_diff.get_element(kk,ll));
                        _gr_opb_gmat.accumulate_element(kk,ll,pkval*
                                         _gr_opb_dens_diff.get_element(ii,jj));

                        _gr_opb_gmat.accumulate_element(ii,jj,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(kk,ll));
                        _gr_opb_gmat.accumulate_element(kk,ll,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(ii,jj));

                        /*
                         * this integral also contributes to K1 and K2 of
                         * G(ik)
                         *
                         * pkval = -0.25 * ((ijkl)+(ilkj))
                         *       = -0.5 * (ijkl)
                         */
                        lij=ioff(ii)+kk;
                        lkl=IOFF(jj,ll);
                        pkval = 0.5 * pki_int;

                        _gr_gmat.accumulate_element(ii,kk,-pkval*
                                             _gr_dens_diff.get_element(jj,ll));
                        _gr_gmat.accumulate_element(jj,ll,-pkval*
                                             _gr_dens_diff.get_element(ii,kk));

                        _gr_opa_gmat.accumulate_element(ii,kk,pkval*
                                         _gr_opa_dens_diff.get_element(jj,ll));
                        _gr_opa_gmat.accumulate_element(jj,ll,pkval*
                                         _gr_opa_dens_diff.get_element(ii,kk));

                        _gr_opa_gmat.accumulate_element(ii,kk,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(jj,ll));
                        _gr_opa_gmat.accumulate_element(jj,ll,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(ii,kk));

                        _gr_opb_gmat.accumulate_element(ii,kk,pkval*
                                         _gr_opb_dens_diff.get_element(jj,ll));
                        _gr_opb_gmat.accumulate_element(jj,ll,pkval*
                                         _gr_opb_dens_diff.get_element(ii,kk));

                        _gr_opb_gmat.accumulate_element(ii,kk,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(jj,ll));
                        _gr_opb_gmat.accumulate_element(jj,ll,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(ii,kk));

                      } else {
                        /*
                         * This integral contributes to J of G(ij)
                         *
                         * pkval = (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        pkval = pki_int;

                        _gr_gmat.accumulate_element(ii,jj,pkval*
                                             _gr_dens_diff.get_element(kk,ll));
                        _gr_gmat.accumulate_element(kk,ll,pkval*
                                             _gr_dens_diff.get_element(ii,jj));

                        /*
                         * and to K1 of G(ik)
                         *
                         * pkval = -0.25 * (ijkl)
                         */
                        lij=ioff(ii)+kk;
                        lkl=IOFF(jj,ll);
                        pkval *= 0.25;

                        _gr_gmat.accumulate_element(ii,kk,-pkval*
                                             _gr_dens_diff.get_element(jj,ll));
                        _gr_gmat.accumulate_element(jj,ll,-pkval*
                                             _gr_dens_diff.get_element(ii,kk));

                        _gr_opa_gmat.accumulate_element(ii,kk,pkval*
                                         _gr_opa_dens_diff.get_element(jj,ll));
                        _gr_opa_gmat.accumulate_element(jj,ll,pkval*
                                         _gr_opa_dens_diff.get_element(ii,kk));

                        _gr_opa_gmat.accumulate_element(ii,kk,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(jj,ll));
                        _gr_opa_gmat.accumulate_element(jj,ll,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(ii,kk));

                        _gr_opb_gmat.accumulate_element(ii,kk,pkval*
                                         _gr_opb_dens_diff.get_element(jj,ll));
                        _gr_opb_gmat.accumulate_element(jj,ll,pkval*
                                         _gr_opb_dens_diff.get_element(ii,kk));

                        _gr_opb_gmat.accumulate_element(ii,kk,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(jj,ll));
                        _gr_opb_gmat.accumulate_element(jj,ll,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(ii,kk));

                        /*
                         * and to K2 of G(il)
                         *
                         * pkval = -0.25 * (ijkl)
                         */
                        lij=ioff(ii)+ll;
                        lkl=IOFF(kk,jj);

                        _gr_gmat.accumulate_element(ii,ll,-pkval*
                                             _gr_dens_diff.get_element(kk,jj));
                        _gr_gmat.accumulate_element(kk,jj,-pkval*
                                             _gr_dens_diff.get_element(ii,ll));

                        _gr_opa_gmat.accumulate_element(ii,ll,pkval*
                                         _gr_opa_dens_diff.get_element(kk,jj));
                        _gr_opa_gmat.accumulate_element(kk,jj,pkval*
                                         _gr_opa_dens_diff.get_element(ii,ll));

                        _gr_opa_gmat.accumulate_element(ii,ll,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(kk,jj));
                        _gr_opa_gmat.accumulate_element(kk,jj,-3.0*pkval*
                                         _gr_opb_dens_diff.get_element(ii,ll));

                        _gr_opb_gmat.accumulate_element(ii,ll,pkval*
                                         _gr_opb_dens_diff.get_element(kk,jj));
                        _gr_opb_gmat.accumulate_element(kk,jj,pkval*
                                         _gr_opb_dens_diff.get_element(ii,ll));

                        _gr_opb_gmat.accumulate_element(ii,ll,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(kk,jj));
                        _gr_opb_gmat.accumulate_element(kk,jj,-3.0*pkval*
                                         _gr_opa_dens_diff.get_element(ii,ll));
                      }
                    }
                    index++;
                  }
                }
              }
            }
          }
          tnint += (double) (n1*n2*n3*n4);
        }
      }
    }
  }
  
  int_reduce_storage_threshold();

  printf("%20.0f integrals\n",tnint);

  delete[] shnfunc;
  delete[] pmax;

  _fock.assign(_gr_gmat);
  _fock.accumulate(_gr_hcore);

  _gr_opa_gmat.scale(-1.0);
  _opa_fock.assign(_fock);
  _opa_fock.accumulate(_gr_opa_gmat);
  _gr_opa_gmat.scale(-1.0);

  _gr_opb_gmat.scale(-1.0);
  _opb_fock.assign(_fock);
  _opb_fock.accumulate(_gr_opb_gmat);
  _gr_opb_gmat.scale(-1.0);

}
