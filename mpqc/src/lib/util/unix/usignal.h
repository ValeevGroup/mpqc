
/* usignal.h -- definition of the unix signal classes
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

#ifndef _libQC_usignal_h
#define _libQC_usignal_h

#ifdef __GNUC__
#pragma interface
#endif

#include <signal.h>

#if defined(SGI)
#ifdef __GNUC__
typedef  void  Sigfunc(int,...);  /* for signal handlers */
#else
typedef  void  Sigfunc(...);  /* for signal handlers */
#endif
#define SIG_HANDLER(a,b) void a(int b,...)
#else
typedef  void  Sigfunc(int);  /* for signal handlers */
#define SIG_HANDLER(a,b) void a(int b)
#endif
typedef Sigfunc* Sigfunc_t;


#if !defined(I860) || defined(PARAGON)
/////////////////////////////////////////////////////////////////////////

class UnixSignalSet {
    friend class UnixSignal;
  private:
    sigset_t sigset;
    sigset_t save_sigset;
  public:
    UnixSignalSet();
    UnixSignalSet(const sigset_t&);
    virtual ~UnixSignalSet();
    
    inline virtual int EmptySet() { return sigemptyset(&sigset); }
    inline virtual int FillSet() { return sigfillset(&sigset); }
    inline virtual int AddToSet(int s) { return sigaddset(&sigset,s); }
    inline virtual int DelFromSet(int s) { return sigdelset(&sigset,s); }
    inline virtual int IsMember(int s) { return sigismember(&sigset,s); }

    virtual int BlockSet();
    virtual int UnblockSet();
    virtual int SetMask(const sigset_t&);
    virtual int IsPending(int);

  };

/////////////////////////////////////////////////////////////////////////

class UnixSignal : public UnixSignalSet {
  private:
    int sig_typ;
    struct sigaction act;
    struct sigaction oact;

    UnixSignal();
   
  public:
    UnixSignal(int,Sigfunc_t =0);
    virtual ~UnixSignal();

    Sigfunc_t SigHandler(Sigfunc_t);
    Sigfunc_t SigHandlerInt(Sigfunc_t);
    Sigfunc_t Ignore();
    Sigfunc_t Reset();

    inline Sigfunc_t Handler() { return (Sigfunc_t) act.sa_handler; }
    inline void AddFlag(int f) { act.sa_flags |= f; }
    inline void DelFlag(int f) { act.sa_flags &= (~f); }
  };

/////////////////////////////////////////////////////////////////////////

class Signals {
    static UnixSignal *sig[32];
  public:
    static void init();
    UnixSignal& operator[](int i) { return *sig[i]; }
  };
 
extern Signals unix_signal;

#endif /* I860 */
#endif
