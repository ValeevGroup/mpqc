
/* daemon.h -- definition of the Daemon class
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

#ifndef _libQC_daemon_h
#define _libQC_daemon_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/unix/cct_cprot.h>
#include <util/unix/usignal.h>

class Daemon {
  public:
    enum { NO_IGN_TTY=1, NO_BKGRND=2, NO_LOSE_TERM=4, NO_CLOSE=8, IGN_CHLD=16 };
    static void WhenKilled(Sigfunc_t);
    static void IgnoreTTY(void);
    static void IgnoreChildren(void);
    static void Start(int);
    static void LoseTerminal(void);
    static int OpenLogfile(const char*);
    static int OpenStdOut(void);
  };

#endif
