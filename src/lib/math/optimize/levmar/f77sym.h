#ifndef FC_HEADER_INCLUDED
#define FC_HEADER_INCLUDED

/* Mangling for Fortran global symbols without underscores. */
#define FC_GLOBAL(name,NAME) name##_

/* Mangling for Fortran global symbols with underscores. */
#define FC_GLOBAL_(name,NAME) name##_

/* Mangling for Fortran module symbols without underscores. */
#define FC_MODULE(mod_name,name, mod_NAME,NAME) __##mod_name##__##name

/* Mangling for Fortran module symbols with underscores. */
#define FC_MODULE_(mod_name,name, mod_NAME,NAME) __##mod_name##__##name

/*--------------------------------------------------------------------------*/
/* Mangle some symbols automatically.                                       */
#define F77_PDSTEQR FC_GLOBAL(pdsteqr, PDSTEQR)
#define F77_DCOPY FC_GLOBAL(dcopy, DCOPY)
#define F77_DNRM2 FC_GLOBAL(dnrm2, DNRM2)
#define F77_DSCAL FC_GLOBAL(dscal, DSCAL)
#define F77_DGEMM FC_GLOBAL(dgemm, DGEMM)
#define F77_DGEMV FC_GLOBAL(dgemv, DGEMV)
#define F77_DAXPY FC_GLOBAL(daxpy, DAXPY)
#define F77_DDOT FC_GLOBAL(ddot, DDOT)
#define F77_DSPMV FC_GLOBAL(dspmv, DSPMV)
#define F77_DGESVD FC_GLOBAL(dgesvd, DGESVD)
#define F77_DSPSVX FC_GLOBAL(dspsvx, DSPSVX)
#define F77_DSYEVD FC_GLOBAL(dsyevd, DSYEVD)
#define F77_DSPTRF FC_GLOBAL(dsptrf, DSPTRF)
#define F77_DPPTRF FC_GLOBAL(dpptrf, DPPTRF)
#define F77_DSPTRI FC_GLOBAL(dsptri, DSPTRI)
#define F77_DPPTRI FC_GLOBAL(dpptri, DPPTRI)
#define F77_DLANSP FC_GLOBAL(dlansp, DLANSP)
#define F77_DSPCON FC_GLOBAL(dspcon, DSPCON)
#define F77_DPPCON FC_GLOBAL(dppcon, DPPCON)
#define F77_DLAMCH FC_GLOBAL(dlamch, DLAMCH)
#define F77_DLACPY FC_GLOBAL(dlacpy, DLACPY)
#define F77_DSPTRS FC_GLOBAL(dsptrs, DSPTRS)
#define F77_DPPTRS FC_GLOBAL(dpptrs, DPPTRS)
#define F77_DSPRFS FC_GLOBAL(dsprfs, DSPRFS)
#define F77_DPPRFS FC_GLOBAL(dpprfs, DPPRFS)
#define F77_DSYGV FC_GLOBAL(dsygv, DSYGV)
#define F77_DGEQRF FC_GLOBAL(dgeqrf, DGEQRF)
#define F77_DGETRF FC_GLOBAL(dgetrf, DGETRF)
#define F77_DGETRS FC_GLOBAL(dgetrs, DGETRS)
#define F77_DORGQR FC_GLOBAL(dorgqr, DORGQR)
#define F77_DPOTF2 FC_GLOBAL(dpotf2, DPOTF2)
#define F77_DTRTRS FC_GLOBAL(dtrtrs, DTRTRS)

#endif
