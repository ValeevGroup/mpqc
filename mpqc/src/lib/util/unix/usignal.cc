
/* usignal.C -- implementation of the unix signal classes
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

#include <signal.h>
#include "usignal.h"
#include "cct_cprot.h"

////////////////////////////////////////////////////////////////////////

UnixSignalSet::UnixSignalSet()
{
  sigemptyset(&sigset);
  sigemptyset(&save_sigset);
}

UnixSignalSet::UnixSignalSet(const sigset_t& s) :
  sigset(s)
{
  sigprocmask(SIG_SETMASK,&sigset,&save_sigset);
  }

UnixSignalSet::~UnixSignalSet()
{
}

int UnixSignalSet::BlockSet()
{
  if(sigprocmask(SIG_BLOCK,&sigset,&save_sigset)<0) {
    err_ret("UnixSignalSet::BlockSet(): sigprocmask call failed");
    return -1;
    }
  return 0;
  }

int UnixSignalSet::UnblockSet()
{
  if(sigprocmask(SIG_UNBLOCK,&sigset,&save_sigset)<0) {
    err_ret("UnixSignalSet::UnblockSet(): sigprocmask call failed");
    return -1;
    }
  return 0;
  }

int UnixSignalSet::SetMask(const sigset_t& s)
{
  sigset=s;
  if(sigprocmask(SIG_SETMASK,&sigset,&save_sigset)<0) {
    err_ret("UnixSignalSet::SetMask(): sigprocmask call failed");
    return -1;
    }
  return 0;
  }

int UnixSignalSet::IsPending(int sig)
{
  sigset_t pend;

  if(sigpending(&pend)<0)
    err_sys("UnixSignalSet::IsPending(%d): sigpending call failed",sig);

  return sigismember(&pend,sig);
  }

////////////////////////////////////////////////////////////////////////

UnixSignal::UnixSignal()
{
}

UnixSignal::~UnixSignal()
{
}

UnixSignal::UnixSignal(int sig, Sigfunc_t hand) :
  UnixSignalSet(), sig_typ(sig)
{
  if(!hand) hand=(Sigfunc_t)SIG_DFL;

  act.sa_handler=hand;
  act.sa_mask=sigset;
  act.sa_flags=0;

  if(sig_typ==SIGALRM) {
#ifdef SA_INTERRUPT // SunOS
    act.sa_flags |= SA_INTERRUPT;
#endif
    }
  else {
#ifdef SA_RESTART // SVR4 & 4.3+BSD
    act.sa_flags |= SA_RESTART;
#endif
    }

  if(sigaction(sig_typ,&act,&oact) < 0)
    err_sys("UnixSignal::UnixSignal(%d): sigaction call failed",sig_typ);
  }

Sigfunc_t UnixSignal::SigHandler(Sigfunc_t hand)
{
  if(!hand) hand=(Sigfunc_t)SIG_DFL;
  act.sa_handler=hand;
  act.sa_mask=sigset;
  if(sigaction(sig_typ,&act,&oact) < 0) {
    err_ret("UnixSignal::SigHandler: sigaction() failed");
    return((Sigfunc_t)SIG_ERR);
    }
  return((Sigfunc_t)oact.sa_handler);
  }

Sigfunc_t UnixSignal::SigHandlerInt(Sigfunc_t hand)
{
  if(!hand) hand=(Sigfunc_t)SIG_DFL;
  act.sa_handler=hand;
#ifdef SA_INTERRUPT // SunOS
  act.sa_flags |= SA_INTERRUPT;
#endif
#ifdef SA_RESTART
  act.sa_flags &= (~SA_RESTART);
#endif
  act.sa_mask=sigset;
  if(sigaction(sig_typ,&act,&oact) < 0) {
    err_ret("UnixSignal::SigHandlerInt: sigaction() failed");
    return((Sigfunc_t)SIG_ERR);
    }
  return((Sigfunc_t)oact.sa_handler);
  }

Sigfunc_t UnixSignal::Ignore()
{
  act.sa_handler=(Sigfunc_t)SIG_IGN;
  act.sa_mask=sigset;
  if(sigaction(sig_typ,&act,&oact) < 0) {
    err_ret("UnixSignal::Ignore: sigaction() failed");
    return((Sigfunc_t)SIG_ERR);
    }
  return((Sigfunc_t)oact.sa_handler);
  }
  
Sigfunc_t UnixSignal::Reset()
{
  act=oact;
  if(sigaction(sig_typ,&act,&oact) < 0) {
    err_ret("UnixSignal::Reset: sigaction() failed");
    return((Sigfunc_t)SIG_ERR);
    }
  return((Sigfunc_t)oact.sa_handler);
  }

////////////////////////////////////////////////////////////////////////

UnixSignal * Signals::sig[32];

void Signals::init()
{
  for(int i=1; i < 32; i++) 
    if(i!=SIGKILL && i!=SIGSTOP) 
      sig[i]=new UnixSignal(i);
  }

Signals unix_signal;
