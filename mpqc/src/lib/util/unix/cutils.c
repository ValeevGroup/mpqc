
/* cutils.c -- some utility functions for libCCT
 *
 *      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
 *      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
 *      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
 *      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
 *      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
 *      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
 *
 *  Author:
 *      E. T. Seidl
 *      Bldg. 12A, Rm. 2033
 *      Computer Systems Laboratory
 *      Division of Computer Research and Technology
 *      National Institutes of Health
 *      Bethesda, Maryland 20892
 *      Internet: seidl@alw.nih.gov
 *      November, 1992
 */

/*
 * this code is from the the book
 * Advanced Programming in the UNIX Environment
 * by W. Richard Stevens
 *
 * it was obtained via anonymous ftp from ftp.uu.net, filename
 * /published/books/stevens.advprog.tar.Z
 *
 * I assume that it is public domain -- ETS
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include  <errno.h>    /* for definition of errno */
#include  <stdarg.h>    /* ANSI C header file */
#include <sys/param.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <unistd.h>

#include  <util/unix/cct_cprot.h>

#ifdef HAVE_SYSLOG_H
#include <syslog.h>
#endif

static void  err_doit(int, const char *, va_list);

char  *pname = NULL;    /* caller can set this from argv[0] */

/* Nonfatal error related to a system call.
 * Print a message and return. */

void
err_ret(const char *fmt, ...)
{
  va_list    ap;

  va_start(ap, fmt);
  err_doit(1, fmt, ap);
  va_end(ap);
  return;
}

/* Fatal error related to a system call.
 * Print a message and terminate. */

void
err_sys(const char *fmt, ...)
{
  va_list    ap;

  va_start(ap, fmt);
  err_doit(1, fmt, ap);
  va_end(ap);
  exit(1);
}

/* Fatal error related to a system call.
 * Print a message, dump core, and terminate. */

void
err_dump(const char *fmt, ...)
{
  va_list    ap;

  va_start(ap, fmt);
  err_doit(1, fmt, ap);
  va_end(ap);
  abort();    /* dump core and terminate */
  exit(1);    /* shouldn't get here */
}

/* Nonfatal error unrelated to a system call.
 * Print a message and return. */

void
err_msg(const char *fmt, ...)
{
  va_list    ap;

  va_start(ap, fmt);
  err_doit(0, fmt, ap);
  va_end(ap);
  return;
}

/* Fatal error unrelated to a system call.
 * Print a message and terminate. */

void
err_quit(const char *fmt, ...)
{
  va_list    ap;

  va_start(ap, fmt);
  err_doit(0, fmt, ap);
  va_end(ap);
  exit(1);
}

/* Print a message and return to caller.
 * Caller specifies "errnoflag". */

static void
err_doit(int errnoflag, const char *fmt, va_list ap)
{
  int    errno_save;
  char  buf[MAXLINE];

  errno_save = errno;    /* value caller might want printed */
  vsprintf(buf, fmt, ap);
  if (errnoflag)
    sprintf(buf+strlen(buf), ": %s", strerror(errno_save));
  strcat(buf, "\n");
  fflush(stdout);    /* in case stdout and stderr are the same */
  fputs(buf, stderr);
  fflush(stderr);    /* SunOS 4.1.* doesn't grok NULL argument */
  return;
}

/**************************************************************/


/* Error routines for programs that can run as a daemon. */

static void log_doit(int, int, const char *, va_list ap);

/* Initialize syslog(), if running as daemon. */

void
log_open(const char *ident, int option, int facility)
{
#ifdef HAVE_SYSLOG_H
  openlog(ident, option, facility);
#endif
}

/* Nonfatal error related to a system call.
 * Print a message with the system's errno value and return. */

void
log_ret(const char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
#ifdef HAVE_SYSLOG_H
  log_doit(1, LOG_ERR, fmt, ap);
#endif
  va_end(ap);
  return;
  }

/* Fatal error related to a system call.
 * Print a message and terminate. */

void
log_sys(const char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
#ifdef HAVE_SYSLOG_H
  log_doit(1, LOG_ERR, fmt, ap);
#endif
  va_end(ap);
  exit(2);
  }

/* Nonfatal error unrelated to a system call.
 * Print a message and return. */

void
log_msg(const char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
#ifdef HAVE_SYSLOG_H
  log_doit(0, LOG_ERR, fmt, ap);
#endif
  va_end(ap);
  return;
  }

/* Fatal error unrelated to a system call.
 * Print a message and terminate. */

void
log_quit(const char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
#ifdef HAVE_SYSLOG_H
  log_doit(0, LOG_ERR, fmt, ap);
#endif
  va_end(ap);
  exit(2);
  }

/* Print a message and return to caller.
 * Caller specifies "errnoflag" and "priority". */

static void
log_doit(int errnoflag, int priority, const char *fmt, va_list ap)
{
  int		errno_save;
  char	buf[MAXLINE];

  errno_save = errno;		/* value caller might want printed */
  vsprintf(buf, fmt, ap);
  if (errnoflag)
    sprintf(buf+strlen(buf), ": %s", strerror(errno_save));

  strcat(buf, "\n");

#ifdef HAVE_SYSLOG_H
  syslog(priority, buf);
#endif

  return;
  }

/**************************************************************/

#if !defined(HAVE_STRERROR)

static struct errordesc {
  int val;
  char *msg;
  } errors[] = {
  0, "Not an error",
#ifdef	EPERM
  1, "Not owner",
#endif
#ifdef	ENOENT
  2, "No such file or directory",
#endif
#ifdef	ESRCH
  3, "No such process",
#endif
#ifdef	EINTR
  4, "Interrupted system call",
#endif
#ifdef	EIO
  5, "I/O error",
#endif
#ifdef	ENXIO
  6, "No such device or address",
#endif
#ifdef	E2BIG
  7, "Arg list too long",
#endif
#ifdef	ENOEXEC
  8, "Exec format error",
#endif
#ifdef	EBADF
  9, "Bad file number",
#endif
#ifdef	ECHILD
 10, "No children",
#endif
#ifdef	EAGAIN
 11, "No more processes",
#endif
#ifdef	ENOMEM
 12, "Not enough core",
#endif
#ifdef	EACCES
 13, "Permission denied",
#endif
#ifdef	EFAULT
 14, "Bad address",
#endif
#ifdef	ENOTBLK
 15, "Block device required",
#endif
#ifdef	EBUSY
 16, "Mount device busy",
#endif
#ifdef	EEXIST
 17, "File exists",
#endif
#ifdef	EXDEV
 18, "Cross-device link",
#endif
#ifdef	ENODE
 19, "No such device",
#endif
#ifdef	ENOTDIR
 20, "Not a director",
#endif
#ifdef	EISDIR
 21, "Is a directory",
#endif
#ifdef	EINVAL
 22, "Invalid argument",
#endif
#ifdef	ENFILE
 23, "File table overflow",
#endif
#ifdef	EMFILE
 24, "Too many open files",
#endif
#ifdef	ENOTTY
 25, "Not a typewriter",
#endif
#ifdef	ETXTBSY
 26, "Text file busy",
#endif
#ifdef	EFBIG
 27, "File too large",
#endif
#ifdef	ENOSPC
 28, "No space left on device",
#endif
#ifdef	ESPIPE
 29, "Illegal seek",
#endif
#ifdef	EROFS
 30, "Read-only file system",
#endif
#ifdef	EMLINK
 31, "Too many links",
#endif
#ifdef	EPIPE
 32, "Broken pipe",
#endif
#ifdef	EDOM
 33, "Argument too large",
#endif
#ifdef	ERANGE
 34, "Result too large",
#endif
#ifdef	EWOULDBLOCK
 35, "Operation would block",
#endif
#ifdef	EINPROGRESS
 36, "Operation now in progress",
#endif
#ifdef	EALREADY
 37, "Operation already in progress",
#endif
#ifdef	ENOTSOCK
 38, "Socket operation on non-socket",
#endif
#ifdef	EDESTADDRREQ
 39, "Destination address required",
#endif
#ifdef	EMSGSIZE
 40, "Message too long",
#endif
#ifdef	EPROTOTYPE
 41, "Protocol wrong type for socket",
#endif
#ifdef	ENOPROTOOPT
 42, "Protocol not available",
#endif
#ifdef	EPROTONOSUPPORT
 43, "Protocol not supported",
#endif
#ifdef	ESOCKTNOSUPPORT
 44, "Socket type not supported",
#endif
#ifdef	EOPNOTSUPP
 45, "Operation not supported on socket",
#endif
#ifdef	EPFNOSUPPORT
 46, "Protocol family not supported",
#endif
#ifdef	EAFNOSUPPORT
 47, "Address family not supported by protocol family",
#endif
#ifdef	EADDRINUSE
 48, "Address already in use",
#endif
#ifdef	EADDRNOTAVAIL
 49, "Can't assign requested address",
#endif
#ifdef	ENETDOWN
 50, "Network is down",
#endif
#ifdef	ENETUNREACH
 51, "Network is unreachable",
#endif
#ifdef	ENETRESET
 52, "Network dropped connection on reset",
#endif
#ifdef	ECONNABORTED
 53, "Software caused connection abort",
#endif
#ifdef	ECONNRESET
 54, "Connection reset by peer",
#endif
#ifdef	ENOBUFS
 55, "No buffer space available",
#endif
#ifdef	EISCONN
 56, "Socket is already connected",
#endif
#ifdef	ENOTCONN
 57, "Socket is not connected",
#endif
#ifdef	ESHUTDOWN
 58, "Can't send after socket shutdown",
#endif
#ifdef	ETOOMANYREFS
 59, "Too many references: can't splice",
#endif
#ifdef	ETIMEDOUT
 60, "Connection timed out",
#endif
#ifdef	ECONNREFUSED
 61, "Connection refused",
#endif
#ifdef	ELOOP
 62, "Too many levels of symbolic links",
#endif
#ifdef	ENAMETOOLONG
 63, "File name too long",
#endif
#ifdef	EHOSTDOWN
 64, "Host is down",
#endif
#ifdef	EHOSTUNREACH
 65, "No route to host",
#endif
#ifdef	ENOTEMPTY
 66, "Directory not empty",
#endif
#ifdef	EPROCLIM
 67, "Too many processes",
#endif
#ifdef	EUSERS
 68, "Too many users",
#endif
#ifdef	EDQUOT
 69, "Disc quota exceeded",
#endif
#ifdef	ESTALE
 70, "Stale NFS file handle",
#endif
#ifdef	EREMOTE
 71, "Too many levels of remote in path",
#endif
#ifdef	ENOSTR
 72, "Device is not a stream",
#endif
#ifdef	ETIME
 73, "Timer expired",
#endif
#ifdef	ENOSR
 74, "Out of streams resources",
#endif
#ifdef	ENOMSG
 75, "No message of desired type",
#endif
#ifdef	EBADMSG
 76, "Trying to read unreadable message",
#endif
#ifdef EIDRM
 77, "Identifier removed",
#endif
#ifdef EDEADLK
 78, "Deadlock condition.",
#endif
#ifdef ENOLCK
 79, "No record locks available.",
#endif
#ifdef ENONET
 80, "Machine is not on the network",
#endif
#ifdef ERREMOTE
 81, "Object is remote",
#endif
#ifdef ENOLINK
 82, "the link has been severed",
#endif
#ifdef EADV
 83, "advertise error",
#endif
#ifdef ESRMNT
 84, "srmount error",
#endif
#ifdef ECOMM
 85, "Communication error on send",
#endif
#ifdef EPROTO
 86, "Protocol error",
#endif
#ifdef EMULTIHOP
 87, "multihop attempted",
#endif
#ifdef EDOTDOT
 88, "Cross mount point (not an error)",
#endif
#ifdef EREMCHG
 89, "Remote address changed",
#endif
#ifdef ENOSYS
 90, "function not implemented",
#endif
  -1, NULL };


char *
strerror(int error)
{
  int i;
  static char	mesg[30];

  for(i=0; errors[i].msg!=NULL ; i++) 
    if(error==errors[i].val) return errors[i].msg;

  sprintf(mesg, "Unknown error (%d)", error);
  return(mesg);
  }

#endif
/**************************************************************/

/*
 * NOFILE should be defined in <sys/param.h>
 */

#ifndef NOFILE
#define NOFILE 256
#endif

/* 
 * close all open file descriptors, and then reset errno to get rid of
 * any EBADFs
 */
void
close_open_files() 
{
  int i;
  for(i=0; i < NOFILE; i++) close(i);
  errno=0;
  }


/*
 * print a message about the exit status of a child process
 */

void
pr_exit(int status)
{
  int stat=status;

  if(WIFEXITED(stat))
   fprintf(stderr,"normal termination, exit status = %d\n",WEXITSTATUS(stat));

  else if(WIFSIGNALED(stat))
    fprintf(stderr,"abnormal termination, signal = %d%s\n",WTERMSIG(stat),
#ifdef WCOREDUMP
       WCOREDUMP(stat) ? " (core file generated)" : ""
#else
       ""
#endif
    );

  else if(WIFSTOPPED(stat))
    fprintf(stderr,"child stopped, signal = %d\n",WSTOPSIG(stat));

  }
