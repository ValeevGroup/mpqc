#ifdef __GNUG__
#pragma implementation
#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <fcntl.h>
#ifndef F_SETFD
#  define F_SETFD 2
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <iostream.h>
#include <signal.h>

#include <util/keyval/keyval.h>
#include <util/misc/bug.h>

//////////////////////////////////////////////////////////////////////
// static variables

static Debugger *signals[NSIG];

//////////////////////////////////////////////////////////////////////
// static routines

#ifdef SIGHASELLIP
// required for CC -64 on IRIX 6.0.1 and for gcc on IRIX 5.3
typedef RETSIGTYPE (*handler_type)(...);
#else
typedef RETSIGTYPE (*handler_type)(int);
#endif

static void
handler(int sig)
{
  if (signals[sig]) signals[sig]->got_signal(sig);
}

static char *
append(char *cmd, const char *a)
{
  int len = strlen(a)+1;
  if (cmd) len += strlen(cmd);
  char *n = new char[len];
  n[0] = '\0';
  if (cmd) strcat(n,cmd);
  strcat(n,a);
  return n;
}

static char *
append(char *cmd, char a)
{
  char aa[2];
  aa[0] = a;
  aa[1] = '\0';
  return append(cmd, aa);
}

static char *
append(char *cmd, int i)
{
  char a[128];
  sprintf(a,"%d",i);
  return append(cmd, a);
}

//////////////////////////////////////////////////////////////////////
// Debugger class definition

Debugger *Debugger::default_debugger_ = 0;

#define CLASSNAME Debugger
#define PARENTS public SavableState
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
Debugger::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

Debugger::Debugger(const char *exec)
{
  init();

  prefix_ = new char[1];
  prefix_[0] = '\0';

  set_exec(exec);

  default_cmd();
}

Debugger::Debugger(const RefKeyVal &keyval)
{
  init();

  debug_ = keyval->booleanvalue("debug");
  if (keyval->error() != KeyVal::OK) debug_ = 1;

  traceback_ = keyval->booleanvalue("traceback");
  if (keyval->error() != KeyVal::OK) traceback_ = 1;

  exit_on_signal_ = keyval->booleanvalue("exit");
  if (keyval->error() != KeyVal::OK) traceback_ = 1;

  wait_for_debugger_ = keyval->booleanvalue("wait_for_debugger");
  if (keyval->error() != KeyVal::OK) wait_for_debugger_ = 1;

  cmd_ = keyval->pcharvalue("cmd");

  prefix_ = keyval->pcharvalue("prefix");

  if (keyval->booleanvalue("handle_defaults")) handle_defaults();
  if (keyval->error() != KeyVal::OK) handle_defaults();

  if (cmd_ == 0) default_cmd();

  if (prefix_ == 0) {
      prefix_ = new char[1];
      prefix_[0] = '\0';
    }
}

Debugger::~Debugger()
{
  delete[] prefix_;
  delete[] exec_;
  delete[] cmd_;
  for (int i=0; i<NSIG; i++) {
      if (mysigs_[i]) signals[i] = 0;
    }
  delete[] mysigs_;
}

Debugger::Debugger(StateIn&s):
  SavableState(s)
{
  init();

  s.getstring(prefix_);
  s.getstring(exec_);
  s.getstring(cmd_);
  s.get(debug_);
  s.get(traceback_);
  s.get(exit_on_signal_);
  s.get(wait_for_debugger_);

  int i, nsig, tmp;
  s.get(nsig);
  for (i=0; i<nsig; i++) {
      s.get(tmp);
      handle(tmp);
    }
}

void
Debugger::save_data_state(StateOut&s)
{
  s.putstring(prefix_);
  s.putstring(exec_);
  s.putstring(cmd_);
  s.put(debug_);
  s.put(traceback_);
  s.put(exit_on_signal_);
  s.put(wait_for_debugger_);

  int i, nsig = 0;
  for (i=0; i<NSIG; i++) if (mysigs_[i]) nsig++;
  s.put(nsig);
  for (i=0; i<NSIG; i++) if (mysigs_[i]) s.put(i);
}

void
Debugger::init()
{
  exec_ = 0;
  prefix_ = 0;
  cmd_ = 0;

  exit_on_signal_ = 1;
  traceback_ = 1;
  debug_ = 1;
  wait_for_debugger_ = 1;

  mysigs_ = new int[NSIG];
  for (int i=0; i<NSIG; i++) {
      mysigs_[i] = 0;
    }
}

void
Debugger::handle(int sig)
{
  if (sig >= NSIG) return;
#ifdef HAVE_SIGNAL
  signal(sig, (handler_type)handler);
#endif
  signals[sig] = this;
  mysigs_[sig] = 1;
}

void
Debugger::handle_defaults()
{
  handle(SIGSEGV);
  handle(SIGFPE);
  handle(SIGQUIT);
  handle(SIGIOT);
  handle(SIGINT);
  handle(SIGHUP);
#ifdef SIGBUS
  handle(SIGBUS);
#endif
}

void
Debugger::set_exec(const char *exec)
{
  delete[] exec_;
  if (exec) {
      exec_ = new char[strlen(exec)+1];
      strcpy(exec_, exec);
    }
  else {
      exec_ = 0;
    }
}

void
Debugger::set_prefix(const char *p)
{
  delete[] prefix_;
  if (p) {
      prefix_ = new char[strlen(p)+1];
      strcpy(prefix_, p);
    }
  else {
      prefix_ = 0;
    }
}

void
Debugger::set_prefix(int i)
{
  char p[128];
  sprintf(p, "%3d: ", i);
  set_prefix(p);
}

void
Debugger::default_cmd()
{
#ifdef __GNUG__
  int gcc = 1;
#else
  int gcc = 0;
#endif
  int has_x11_display = (getenv("DISPLAY") != 0);

  if (!gcc && sizeof(void*) == 8 && has_x11_display) {
      set_cmd("xterm -title \"$(PREFIX)$(EXEC)\" -e dbx -p $(PID) $(EXEC) &");
    }
  else if (has_x11_display) {
      set_cmd("xterm -title \"$(PREFIX)$(EXEC)\" -e gdb $(EXEC) $(PID) &");
    }
  else {
      set_cmd(0);
    }
}

void
Debugger::set_cmd(const char *cmd)
{
  delete[] cmd_;
  if (cmd) {
      cmd_ = new char[strlen(cmd)+1];
      strcpy(cmd_, cmd);
    }
  else {
      cmd_ = 0;
    }
}

void
Debugger::debug(const char *reason)
{
  cout << prefix_ << "Debugger::debug: ";
  if (reason) cout << reason;
  else cout << "no reason given";
  cout << endl;

#ifndef HAVE_SYSTEM
  abort();
#else
  if (cmd_) {
      int pid = getpid();
      char cpid[128];
      sscanf(cpid,"%d",pid);
      // contruct the command name
      char *cmd = 0;
      for (char *c=cmd_; *c;) {
          if (!strncmp("$(PID)",c,6)) {
              cmd = append(cmd,pid);
              c += 6;
            }
          else if (!strncmp("$(EXEC)",c,7)) {
              if (exec_) cmd = append(cmd,exec_);
              c += 7;
            }
          else if (!strncmp("$(PREFIX)",c,9)) {
              if (prefix_) cmd = append(cmd,prefix_);
              c += 9;
            }
          else {
              cmd = append(cmd,*c);
              c++;
            }
        }
      // start the debugger
      cout << prefix_ << "Debugger: starting \"" << cmd << "\"" << endl;
      debugger_ready_ = 0;
      system(cmd);
      delete[] cmd;
      // wait until the debugger is ready
      if (wait_for_debugger_) while(!debugger_ready_);
    }
#endif
}

void
Debugger::got_signal(int sig)
{
  const char *signame;
  if (sig == SIGSEGV) signame = "SIGSEGV";
  else if (sig == SIGFPE) signame = "SIGFPE";
  else if (sig == SIGHUP) signame = "SIGHUP";
  else if (sig == SIGINT) signame = "SIGINT"; 
#ifdef SIGBUS
  else if (sig == SIGBUS) signame = "SIGBUS";
#endif
  else signame = "UNKNOWN SIGNAL";

  if (traceback_) {
      traceback(signame);
    }
  if (debug_) {
      debug(signame);
    }

  if (exit_on_signal_) {
      cout << prefix_ << "Debugger: exiting" << endl;
      exit(1);
    }
  else {
      cout << prefix_ << "Debugger: continuing" << endl;
    }

  //handle(sig);
}

void
Debugger::set_debug_on_signal(int v)
{
  debug_ = v;
}

void
Debugger::set_traceback_on_signal(int v)
{
  traceback_ = v;
}

void
Debugger::set_wait_for_debugger(int v)
{
  wait_for_debugger_ = v;
}

void
Debugger::set_exit_on_signal(int v)
{
  exit_on_signal_ = v;
}

void
Debugger::set_default_debugger(const RefDebugger &d)
{
  if (default_debugger_) {
      default_debugger_->dereference();
      // let a smart pointer figure out what to do with the old debugger
      RefDebugger old(default_debugger_);
    }
  if (d.pointer()) d.pointer()->reference();
  default_debugger_ = d.pointer();
}

Debugger *
Debugger::default_debugger()
{
  return default_debugger_;
}

#define SIMPLE_STACK (defined(linux) && defined(i386)) \
                     || (defined(__OSF1__) && defined(i860))

void
Debugger::traceback(const char *reason)
{
#if SIMPLE_STACK
  int bottom = 0x1234;
  void **topstack = (void**)0xffffffffL;
  void **botstack = (void**)0x70000000L;
  // signal handlers can put weird things in the return address slot,
  // so it is usually best to keep toptext large.
  void **toptext = (void**)0xffffffffL;
  void **bottext = (void**)0x00010000L;
#endif // SIMPLE_STACK

  cout << prefix_ << "Debugger::traceback:";
  if (reason) cout << reason;
  else cout << "no reason given";
  cout << endl;
#if (defined(linux) && defined(i386))
  topstack = (void**)0xc0000000;
  botstack = (void**)0xb0000000;
#endif
#if (defined(__OSF1__) && defined(i860))
  topstack = (void**)0x80000000;
  botstack = (void**)0x70000000;
#endif

#if SIMPLE_STACK
  // This will go through the stack assuming a simple linked list
  // of pointers to the previous frame followed by the return address.
  // It trys to be careful and avoid creating new execptions, but there
  // are no guarantees.
  void **stack = (void**) &bottom;

  void **frame_pointer = (void**) stack[1];
  while(frame_pointer >= botstack
        && frame_pointer < topstack
        && frame_pointer[1] >= bottext
        && frame_pointer[1] < toptext) {
      cout << prefix_ << "frame: " << (void*)frame_pointer;
      cout.flush();
      cout << "  retaddr: " << frame_pointer[1] << endl;
      frame_pointer = (void**)*frame_pointer;
    }
#else
  cout << prefix_ << "traceback not available for this arch" << endl;
#endif // SIMPLE_STACK
}
