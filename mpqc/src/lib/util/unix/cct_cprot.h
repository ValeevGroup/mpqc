
/*
 * this code is from the the book
 * Advanced Programming in the UNIX Environment
 * by W. Richard Stevens
 *
 * it was obtained via anonymous ftp from ftp.uu.net, filename
 * /published/books/stevens.advprog.tar.Z
 *
 * I assume that it is public domain -- ETS
 *
 * this file was originally ourhdr.h, much was removed
 */

#ifndef  _libQC_cct_cprot_h
#define  _libQC_cct_cprot_h

#ifndef MAXLINE
#define MAXLINE 4096
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern char *strerror(int);

extern void close_open_files(void);

extern void  err_dump(const char *, ...);  /* {App misc_source} */
extern void  err_msg(const char *, ...);
extern void  err_quit(const char *, ...);
extern void  err_ret(const char *, ...);
extern void  err_sys(const char *, ...);

extern void  log_msg(const char *, ...);    /* {App misc_source} */
extern void  log_open(const char *, int, int);
extern void  log_quit(const char *, ...);
extern void  log_ret(const char *, ...);
extern void  log_sys(const char *, ...);

extern void pr_exit(int);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
