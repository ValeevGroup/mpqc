//
// bug.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUG__
#pragma implementation
#endif

#ifdef HAVE_CONFIG_H
#include <scconfig.h>
#endif

#include <fcntl.h>
#ifndef F_SETFD
#  define F_SETFD 2
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include <signal.h>

#ifdef HAVE_BACKTRACE
#  include <execinfo.h>
#endif

#include <util/keyval/keyval.h>
#include <util/misc/bug.h>
#include <util/state/stateio.h>

// usually in signal.h, but not always.
#ifndef NSIG
#  define NSIG 100
#endif

using namespace std;
using namespace sc;

//////////////////////////////////////////////////////////////////////
// static variables

static Debugger *signals[NSIG];

//////////////////////////////////////////////////////////////////////
// static routines

extern "C" {
#ifdef SIGHASELLIP
// required for CC -64 on IRIX 6.0.1 and for gcc on IRIX 5.3
typedef RETSIGTYPE (*handler_type)(...);
#else
typedef RETSIGTYPE (*handler_type)(int);
#endif
}

static void
handler(int sig)
{
  if (signals[sig]) signals[sig]->got_signal(sig);
}

static void
append(char *cmd, const char *a, int len)
{
  int l = strlen(cmd) + strlen(a)+1;
  if (l > len) {
      ExEnv::outn() << "Debugger: command string too long" << endl;
      abort();
    }
  strcat(cmd,a);
}

static void
append(char *cmd, char a, int len)
{
  char aa[2];
  aa[0] = a;
  aa[1] = '\0';
  append(cmd, aa, len);
}

static void
append(char *cmd, int i, int len)
{
  char a[128];
  sprintf(a,"%d",i);
  append(cmd, a, len);
}

//////////////////////////////////////////////////////////////////////
// Debugger class definition

Debugger *Debugger::default_debugger_ = 0;

static ClassDesc Debugger_cd(
    typeid(Debugger),"Debugger",1,"public SavableState",
    0, create<Debugger>, create<Debugger>);

Debugger::Debugger(const char *exec)
{
  init();

  prefix_ = new char[1];
  prefix_[0] = '\0';

  set_exec(exec);

  default_cmd();
}

Debugger::Debugger(const Ref<KeyVal> &keyval)
{
  init();

  debug_ = keyval->booleanvalue("debug");
  if (keyval->error() != KeyVal::OK) debug_ = 1;

  traceback_ = keyval->booleanvalue("traceback");
  if (keyval->error() != KeyVal::OK) traceback_ = 1;

  exit_on_signal_ = keyval->booleanvalue("exit");
  if (keyval->error() != KeyVal::OK) exit_on_signal_ = 1;

  sleep_ = keyval->intvalue("sleep");
  if (keyval->error() != KeyVal::OK) sleep_ = 0;

  wait_for_debugger_ = keyval->booleanvalue("wait_for_debugger");
  if (keyval->error() != KeyVal::OK) wait_for_debugger_ = 1;

  cmd_ = keyval->pcharvalue("cmd");

  prefix_ = keyval->pcharvalue("prefix");

  handle_sigint_ = keyval->booleanvalue("handle_sigint");
  if (keyval->error() != KeyVal::OK) handle_sigint_=1;
  
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
  s.get(sleep_);
  s.get(debug_);
  s.get(traceback_);
  s.get(exit_on_signal_);
  s.get(wait_for_debugger_);
  s.get(handle_sigint_);

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
  s.put(sleep_);
  s.put(debug_);
  s.put(traceback_);
  s.put(exit_on_signal_);
  s.put(wait_for_debugger_);
  s.put(handle_sigint_);

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
  sleep_ = 0;

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
#ifdef SIGSEGV
  handle(SIGSEGV);
#endif
#ifdef SIGFPE
  handle(SIGFPE);
#endif
#ifdef SIGQUIT
  handle(SIGQUIT);
#endif
#ifdef SIGIOT
  handle(SIGIOT);
#endif
#ifdef SIGINT
  if (handle_sigint_)
      handle(SIGINT);
#endif
#ifdef SIGHUP
  handle(SIGHUP);
#endif
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
  ExEnv::outn() << prefix_ << "Debugger::debug: ";
  if (reason) ExEnv::outn() << reason;
  else ExEnv::outn() << "no reason given";
  ExEnv::outn() << endl;

#ifndef HAVE_SYSTEM
  abort();
#else
  if (cmd_) {
      int pid = getpid();
      // contruct the command name
      const int cmdlen = 512;
      char cmd[cmdlen];
      cmd[0] = '\0';
      for (char *c=cmd_; *c;) {
          if (!strncmp("$(PID)",c,6)) {
              append(cmd,pid,cmdlen);
              c += 6;
            }
          else if (!strncmp("$(EXEC)",c,7)) {
              if (exec_) append(cmd,exec_,cmdlen);
              c += 7;
            }
          else if (!strncmp("$(PREFIX)",c,9)) {
              if (prefix_) append(cmd,prefix_,cmdlen);
              c += 9;
            }
          else {
              append(cmd,*c,cmdlen);
              c++;
            }
        }
      // start the debugger
      ExEnv::outn() << prefix_ << "Debugger: starting \"" << cmd << "\"" << endl;
      debugger_ready_ = 0;
      system(cmd);
      // wait until the debugger is ready
      if (sleep_) {
          ExEnv::outn() << prefix_ << "Sleeping " << sleep_
               << " seconds to wait for debugger ..." << endl;
          sleep(sleep_);
      }
      if (wait_for_debugger_) {
          ExEnv::outn() << prefix_
                        << ": Spinning until debugger_ready_ is set ..." << endl;
          while(!debugger_ready_);
        }
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
      ExEnv::outn() << prefix_ << "Debugger: exiting" << endl;
      exit(1);
    }
  else {
      ExEnv::outn() << prefix_ << "Debugger: continuing" << endl;
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
Debugger::set_default_debugger(const Ref<Debugger> &d)
{
  if (default_debugger_) {
      default_debugger_->dereference();
      // let a smart pointer figure out what to do with the old debugger
      Ref<Debugger> old(default_debugger_);
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
#ifdef HAVE_BACKTRACE
  ExEnv::outn() << prefix_ << "Debugger::traceback(using backtrace):";
  if (reason) ExEnv::outn() << reason;
  else ExEnv::outn() << "no reason given";
  ExEnv::outn() << endl;

  const int n = 100;
  void *p[n];
  int nret = backtrace(p,n);
  if (nret == 0) {
      ExEnv::outn() << prefix_ << "backtrace returned no state information" << std::endl;
    }
  for (int i=0; i<nret; i++) {
      ExEnv::outn() << prefix_
                    << "frame " << i
                    << ": return address = " << p[i]
                    << std::endl;
    }
#else // HAVE_BACKTRACE
#if SIMPLE_STACK
  int bottom = 0x1234;
  void **topstack = (void**)0xffffffffL;
  void **botstack = (void**)0x70000000L;
  // signal handlers can put weird things in the return address slot,
  // so it is usually best to keep toptext large.
  void **toptext = (void**)0xffffffffL;
  void **bottext = (void**)0x00010000L;
#endif // SIMPLE_STACK

  ExEnv::outn() << prefix_ << "Debugger::traceback:";
  if (reason) ExEnv::outn() << reason;
  else ExEnv::outn() << "no reason given";
  ExEnv::outn() << endl;
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

  void **frame_pointer = (void**) stack[3];
  while(frame_pointer >= botstack
        && frame_pointer < topstack
        && frame_pointer[1] >= bottext
        && frame_pointer[1] < toptext) {
      ExEnv::outn() << prefix_ << "frame: " << (void*)frame_pointer;
      ExEnv::outn().flush();
      ExEnv::outn() << "  retaddr: " << frame_pointer[1] << endl;
      frame_pointer = (void**)*frame_pointer;
    }
#else
  ExEnv::outn() << prefix_ << "traceback not available for this arch" << endl;
#endif // SIMPLE_STACK
#endif // HAVE_BACKTRACE
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
