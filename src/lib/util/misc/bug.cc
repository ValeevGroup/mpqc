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

#ifdef HAVE_CONFIG_H
#include <mpqc_config.h>
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
#include <sstream>
#include <iterator>

#if defined(HAVE_LIBUNWIND)
#  define UNW_LOCAL_ONLY
#  include <libunwind.h>
#elif defined(HAVE_BACKTRACE)
#  include <execinfo.h>
#endif

#ifdef HAVE_CXA_DEMANGLE
#include <cxxabi.h>
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
append(std::string &cmd, int i)
{
  char a[128];
  sprintf(a,"%d",i);
  cmd += a;
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

  cmd_ = keyval->stringvalue("cmd");

  prefix_ = keyval->stringvalue("prefix");

  handle_sigint_ = keyval->booleanvalue("handle_sigint");
  if (keyval->error() != KeyVal::OK) handle_sigint_=1;
  
  if (keyval->booleanvalue("handle_defaults")) handle_defaults();
  if (keyval->error() != KeyVal::OK) handle_defaults();

  if (cmd_.empty()) default_cmd();

}

Debugger::~Debugger()
{
  for (int i=0; i<NSIG; i++) {
      if (mysigs_[i]) signals[i] = 0;
    }
  delete[] mysigs_;
}

Debugger::Debugger(StateIn&s):
  SavableState(s)
{
  init();

  s.get(prefix_);
  s.get(exec_);
  s.get(cmd_);
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
  s.put(prefix_);
  s.put(exec_);
  s.put(cmd_);
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
  exec_.resize(0);
  prefix_.resize(0);
  cmd_.resize(0);
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
  if (exec) {
      exec_ = exec;
    }
  else {
      exec_.resize(0);
    }
}

void
Debugger::set_prefix(const char *p)
{
  if (p) {
      prefix_ = p;
    }
  else {
      prefix_.resize(0);
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
  if (cmd) {
      cmd_ = cmd;
    }
  else {
      cmd_.resize(0);
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
  if (!cmd_.empty()) {
      int pid = getpid();
      // contruct the command name
      std::string cmd = cmd_;
      std::string::size_type pos;
      std::string pidvar("$(PID)");
      while ((pos = cmd.find(pidvar))
             != std::string::npos) {
          std::string pidstr;
          append(pidstr,pid);
          cmd.replace(pos,pidvar.size(),pidstr);
        }
      std::string execvar("$(EXEC)");
      while ((pos = cmd.find(execvar))
             != std::string::npos) {
          cmd.replace(pos,execvar.size(),exec_);
        }
      std::string prefixvar("$(PREFIX)");
      while ((pos = cmd.find(prefixvar))
             != std::string::npos) {
          cmd.replace(pos,prefixvar.size(),prefix_);
        }
      // start the debugger
      ExEnv::outn() << prefix_ << "Debugger: starting \"" << cmd << "\"" << endl;
      debugger_ready_ = 0;
      system(cmd.c_str());
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
  Debugger::__traceback(prefix_,reason);
}

void
Debugger::__traceback(const std::string& prefix, const char *reason)
{
  Backtrace result(prefix);
  const size_t nframes_to_skip = 2;
#if defined(HAVE_LIBUNWIND)
  ExEnv::outn() << prefix << "Debugger::traceback(using libunwind):";
#elif defined(HAVE_BACKTRACE) // !HAVE_LIBUNWIND
  ExEnv::outn() << prefix << "Debugger::traceback(using backtrace):";
#else // !HAVE_LIBUNWIND && !HAVE_BACKTRACE
# if defined(SIMPLE_STACK)
  ExEnv::outn() << prefix << "Debugger::traceback:";
# else
  ExEnv::outn() << prefix << "traceback not available for this arch" << endl;
  return;
# endif // SIMPLE_STACK
#endif // HAVE_LIBUNWIND, HAVE_BACKTRACE

  if (reason) ExEnv::outn() << reason;
  else ExEnv::outn() << "no reason given";
  ExEnv::outn() << endl;

  if (result.empty())
    ExEnv::outn() << prefix << "backtrace returned no state information" << std::endl;
  else
    ExEnv::outn() << result.str(nframes_to_skip) << std::endl;
}

/////////////////////////////////////////////////////////////////////////////

Debugger::Backtrace::Backtrace(const std::string& prefix) : prefix_(prefix)
{
#ifdef HAVE_LIBUNWIND
  {
    unw_cursor_t cursor; unw_context_t uc;
    unw_word_t ip, sp, offp;
    int frame = 0;

    unw_getcontext(&uc);
    unw_init_local(&cursor, &uc);
    while (unw_step(&cursor) > 0) {
      unw_get_reg(&cursor, UNW_REG_IP, &ip);
      unw_get_reg(&cursor, UNW_REG_SP, &sp);
      char name[1024];
      unw_get_proc_name(&cursor, name, 1024, &offp);
      std::ostringstream oss;
      oss << prefix_
          << "frame " << frame
          << ": " << scprintf("ip = 0x%lx sp = 0x%lx ", (long) ip, (long) sp)
          << " symbol = " << __demangle(name);
      frames_.push_back(oss.str());
      ++frame;
    }
  }
#elif defined(HAVE_BACKTRACE) // !HAVE_LIBUNWIND
  void* stack_addrs[1024];
  const int naddrs = backtrace(stack_addrs, 1024);
  char** frame_symbols = backtrace_symbols(stack_addrs, naddrs);
  // starting @ 1 to skip this function
  for(int i=1; i<naddrs; ++i) {
    // extract (mangled) function name
    std::string mangled_function_name;
    {
      std::istringstream iss(std::string(frame_symbols[i]), std::istringstream::in);
      std::string frame, file, address;
      iss >> frame >> file >> address >> mangled_function_name;
    }

    std::ostringstream oss;
    oss << prefix_
        << "frame " << i
        << ": return address = " << stack_addrs[i] << std::endl
        << "  symbol = " << __demangle(mangled_function_name);
    frames_.push_back(oss.str());
  }
  free(frame_symbols);
#else // !HAVE_LIBUNWIND && !HAVE_BACKTRACE
#if defined(SIMPLE_STACK)
  int bottom = 0x1234;
  void **topstack = (void**)0xffffffffL;
  void **botstack = (void**)0x70000000L;
  // signal handlers can put weird things in the return address slot,
  // so it is usually best to keep toptext large.
  void **toptext = (void**)0xffffffffL;
  void **bottext = (void**)0x00010000L;
#endif // SIMPLE_STACK

#if (defined(linux) && defined(i386))
  topstack = (void**)0xc0000000;
  botstack = (void**)0xb0000000;
#endif
#if (defined(__OSF1__) && defined(i860))
  topstack = (void**)0x80000000;
  botstack = (void**)0x70000000;
#endif

#if defined(SIMPLE_STACK)
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

    std::ostringstream oss;
    oss << prefix_ << "frame: " << (void*)frame_pointer;
    oss << "  retaddr: " << frame_pointer[1];
    frames_.push_back(oss.str());

    frame_pointer = (void**)*frame_pointer;
  }
#endif // SIMPLE_STACK
#endif // HAVE_BACKTRACE
}

Debugger::Backtrace::Backtrace(const Backtrace& other) :
    frames_(other.frames_),
    prefix_(other.prefix_)
{
}

std::string
Debugger::Backtrace::str(size_t nframes_to_skip) const {
  std::ostringstream oss;
  std::copy(frames_.begin() + nframes_to_skip, frames_.end(),
            std::ostream_iterator<std::string>(oss, "\n"));
  return oss.str();
}

std::string
Debugger::Backtrace::__demangle(const std::string& symbol) {
  std::string dsymbol;
#ifdef HAVE_CXA_DEMANGLE
  {
    int status;
    char* dsymbol_char = abi::__cxa_demangle(symbol.c_str(), 0, 0, &status);
    if (status == 0) { // success
      dsymbol = dsymbol_char;
      free(dsymbol_char);
    }
    else // fail
      dsymbol = symbol;
  }
#else
  dsymbol = symbol;
#endif
  return dsymbol;
}

/////////////////////////////////////////////////////////////////////////////
// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
