#ifdef __GNUG__
#pragma interface
#endif

#ifndef _util_misc_bug_h
#define _util_misc_bug_h

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/ref/ref.h>

SavableState_REF_fwddec(Debugger);

//. The \clsnm{Debugger} class assists programmers in finding
// bugs in applications.  It can (sometimes) fork a debugger
// and print a traceback when a signal is received or at the
// request of the programmer.
class Debugger: public SavableState {
#define CLASSNAME Debugger
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/state/stated.h>
#include <util/class/classd.h>
  protected:
    char *prefix_;
    char *exec_;
    char *cmd_;
    int debugger_ready_;

    int debug_;
    int traceback_;
    int exit_on_signal_;
    int sleep_;
    int wait_for_debugger_;
    int *mysigs_;

    void init();

    static Debugger *default_debugger_;
  public:
    Debugger(const char *exec = 0);
    Debugger(const RefKeyVal&);
    Debugger(StateIn&);
    ~Debugger();

    //. The \srccd{debug} member attempts to start a debugger
    // running on the current process.
    virtual void debug(const char *reason = 0);
    //. The \srccd{traceback} member attempts a stack traceback
    // for the current process.  A symbol table must be saved for
    // the executable if any sense is to be made of the traceback.
    // Tracebacks are currently available only for a limited number
    // of architectures.
    virtual void traceback(const char *reason = 0);
    //. These set up the desired behavior upon receipt of a signal.
    // The default is that debug and traceback are both on.  Pass
    // these routines a zero to turn off these features.
    virtual void set_debug_on_signal(int);
    virtual void set_traceback_on_signal(int);
    //. After a signal is processed the default action is to
    // exit the program.  This can be changed by calling
    // \srccd{set\_exit\_on\_signal} with a nonzero argument.
    virtual void set_exit_on_signal(int);
    //. After the debugger is started the default is to wait
    // in a infinite loop.  Call this if with a nonzero argument
    // if you don't want to enter this loop.
    virtual void set_wait_for_debugger(int);

    //. The debugger should start a debugger when signal @var{sig}
    // is caught.
    virtual void handle(int sig);
    //. This calls \srccd{handle} with all of the major signals.
    virtual void handle_defaults();

    //. This sets a prefix which preceeds all messages printing
    // by \clsnm{Debugger}.
    virtual void set_prefix(const char *p);
    //. Set the prefix to the decimal represention of @var{p}
    // followed by a ``: ''.
    virtual void set_prefix(int p);

    //. Sets the command to be exectuted with \srccd{debug} is called.
    // The character sequence ``\$(EXEC)'' is replaced by the executable
    // name (see \srccd{set\_exec}), ``\$(PID)'' is replaced by the
    // current process id, and ``\$(PREFIX)'' is replaced by the
    // prefix.
    virtual void set_cmd(const char *);
    //. Calls \srccd{set\_cmd} with a hopefully suitable default.
    virtual void default_cmd();
    //. Set the name of the exectuble for the current process.
    // It is up to the programmer to set this, even if the \clsnm{Debugger}
    // is initialized with the \clsnmref{KeyVal} constructor.
    virtual void set_exec(const char *);

    //. Called with signal @var{sig} is received.  This will call
    // \srccd{debug}.  This is mainly for internal use.
    virtual void got_signal(int sig);

    //. Set the global default debugger.  The initial value is null.
    static void set_default_debugger(const RefDebugger &);
    //. Return the global default debugger.
    static Debugger *default_debugger();

    void save_data_state(StateOut&);
};
SavableState_REF_dec(Debugger);

#endif
