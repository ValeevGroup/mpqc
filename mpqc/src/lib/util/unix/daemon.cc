
/* daemon.C -- implementation of the Daemon class
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <fcntl.h>
#ifdef HAVE_SYSLOG_H
#include <syslog.h>
#endif
#ifdef HAVE_TERMIOS_H
#include <termios.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>  // for umask()
#endif

#include <iostream.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <util/unix/cct_cprot.h>
#include <util/unix/daemon.h>

#ifdef HAVE_SYS_IOCTL_H
#include <sys/ioctl.h>
#endif

#ifndef SIGCLD
#define SIGCLD SIGCHLD
#endif

////////////////////////////////////////////////////////

/*
 * set the default action for various nasty signals to the function "func".
 * ok if func is null.
 */

void Daemon::WhenKilled(Sigfunc_t func)
{
  unix_signal[SIGTERM].SigHandler(func);
  unix_signal[SIGHUP].SigHandler(func);
  unix_signal[SIGINT].SigHandler(func);
  unix_signal[SIGQUIT].SigHandler(func);
  unix_signal[SIGABRT].SigHandler(func);
#ifdef SIGILL
  unix_signal[SIGILL].SigHandler(func);
#endif
#ifdef SIGFPE
  unix_signal[SIGFPE].SigHandler(func);
#endif
#ifdef SIGBUS
  unix_signal[SIGBUS].SigHandler(func);
#endif
#ifdef SIGSEGV
  unix_signal[SIGSEGV].SigHandler(func);
#endif
#ifdef SIGSYS
  unix_signal[SIGSYS].SigHandler(func);
#endif
#ifdef SIGPIPE
  unix_signal[SIGPIPE].SigHandler(func);
#endif
  }

////////////////////////////////////////////////////////

/*
 * the following ignores the signals SIGTTIN, SIGTTOU, and SIGTSTP as all
 * good daemons should
 */

void Daemon::IgnoreTTY(void)
{
#ifdef SIGTTIN
  unix_signal[SIGTTIN].Ignore();
#endif
#ifdef SIGTTOU
  unix_signal[SIGTTOU].Ignore();
#endif
#ifdef SIGTSTP
  unix_signal[SIGTSTP].Ignore();
#endif
  }

////////////////////////////////////////////////////////


/*
 * this function will wait on dead children so they don't become zombies.
 * be aware that the calling process may get an interrupted system call
 * when we return, so they had better handle that.
 */

static SIG_HANDLER(daemon_sig_child__,sig)
{
#ifdef BSD
  union __wait status;
  while(wait3(&status,WNOHANG,(struct rusage *)0)>0) ;
#endif
  }

void Daemon::IgnoreChildren(void)
{
#ifdef BSD
  unix_signal[SIGCHLD].SigHandler(daemon_sig_child__);
#else
  unix_signal[SIGCLD].Ignore();
#endif
  }

////////////////////////////////////////////////////////

/*
 * disassociate ourselves from the parent's process group so that signals
 * that kill my parent won't kill me, then lose the controlling terminal,
 * and finally make things such that I can't reaquire one.
 */

void Daemon::LoseTerminal()
{
#ifdef BSD
  int fd;

  /* lose old process group */
  if(setpgrp(0,getpid())<0)
    err_sys("lose_terminal:  could not change process group");

  /* lose controlling terminal, after this we can't reaquire one */
  if((fd=open("/dev/tty",O_RDWR)) >= 0) {
    ioctl(fd, TIOCNOTTY,0);
    close(fd);
    }
#else
  int childpid;

  /*
   * for sysv, first let's fork, then, since we are guaranteed to not
   * be a process group leader, setpgrp will get rid of the controlling
   * terminal
   */

  if((childpid=fork())<0)
    err_sys("lose_terminal:  could not fork");
  else if(childpid)
    exit(0);  /* parent dies */

  if(setpgrp()<0)
    err_sys("lose_terminal:  could not change process group");

  /* now fork again to make sure we can't reaquire controlling terminal */
  if((childpid=fork())<0)
    err_sys("Daemon::LoseTerm:  could not fork a second time");
  else if(childpid)
    exit(0);  /* parent dies */
#endif
  }

////////////////////////////////////////////////////////


// start_daemon is a way to start a generic daemon.  It follows the
// recommendations in Chapter 2 of:
//   UNIX Network Programming
//   WR Stevens
//
//  see p. 83 for function daemon_start()

void Daemon::Start(int options)
{
  /*
   * if started by init (process 1) then we can skip a bit brother
   */

  if(getppid()==1) goto out;

  /*
   * ignore terminal stop signals (BSD)
   */

  if (!(options&NO_IGN_TTY)) IgnoreTTY();

  /*
   * run in the background
   */

  if(!(options&NO_BKGRND)) {
    int childpid;
    if((childpid=fork()) < 0)
      err_sys("StartDaemon:  could not fork off child");
    else if(childpid>0) {           // parent
      cerr << "[" << childpid << "]\n";
      exit(0);
      }
    }

  /*
   * disassociate from controlling terminal and process group.
   * ensure the process cannot reacquire a new controlling terminal
   */

  if(!(options&NO_LOSE_TERM)) LoseTerminal();

out:

  /*
   * close any open file descriptors
   */

  if(!(options&NO_CLOSE)) close_open_files();

  /*
   * chdir to root
   */

  chdir("/");

  /*
   * clear any inherited file mode creation mask
   */

  umask(0);

  /*
   * see if the caller cares about the exit status of any children.
   * for sys V, if we ignore SIGCLD, then dead children will not become
   * zombies if we don't wait() them.  For BSD, we have to catch SIGCLD
   * and perform a wait()
   */

  if(options&IGN_CHLD) IgnoreChildren();

  }

/////////////////////////////////////////////////////////////////////////

/*
 * daemon_open_logfile will open a file and send stdout and stderr to
 * that file.
 */
int Daemon::OpenLogfile(const char *path)
{
  int fd;
  if((fd=open(path, O_WRONLY|O_CREAT|O_TRUNC, 0644)) <0) {
    syslog(LOG_ERR,"daemon_open_logfile: open(%s) failed",path);
    err_ret("daemon_open_logfile: open(%s) failed",path);
    return -1;
    }

  /* don't buffer output to fd */
  setbuf(fdopen(fd,"a"),NULL);

  if(dup2(fd,1) == -1) { /* stdout */
    syslog(LOG_ERR, "daemon_open_logfile: dup2(fd,1) failed");
    err_ret("daemon_open_logfile: dup2(fd,1) failed");
    return -1;
    }
  if(dup2(1,2) == -1) {  /* stderr */
    syslog(LOG_ERR, "daemon_open_logfile: dup2(1,2) failed");
    err_ret("daemon_open_logfile: dup2(1,2) failed");
    return -1;
    }

  close(fd);
  return 0;
  }

int Daemon::OpenStdOut(void)
{
  int fd;
  if((fd=open("/dev/tty", O_WRONLY|O_CREAT|O_TRUNC, 0644)) <0) {
    syslog(LOG_ERR,"daemon_open_stdout_and_err: open() failed");
    err_ret("daemon_open_stdout_and_err: open() failed");
    return -1;
    }

#ifdef BSD
  ioctl(fd,TIOCSCTTY,0);
#endif

  /* don't buffer output to fd */
  setbuf(fdopen(fd,"a"),NULL);

  if(dup2(fd,1) == -1) { /* stdout */
    syslog(LOG_ERR, "daemon_open_stdout_and_err: dup2(fd,1) failed");
    err_ret("daemon_open_logfile: dup2(fd,1) failed");
    return -1;
    }
  if(dup2(1,2) == -1) {  /* stderr */
    syslog(LOG_ERR, "daemon_open_stdout_and_err: dup2(1,2) failed");
    err_ret("daemon_open_logfile: dup2(1,2) failed");
    return -1;
    }

  close(fd);
  return 0;
  }
