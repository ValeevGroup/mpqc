//
// bug.cpp
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

#include "bug.h"

#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iterator>
#include <sstream>
#include <unistd.h>

#include "mpqc/util/core/backtrace.h"
#include "mpqc/util/core/exenv.h"
#include "mpqc/util/keyval/forcelink.h"

MPQC_CLASS_EXPORT2("Debugger", mpqc::Debugger)

// usually in signal.h, but not always.
#ifndef NSIG
#define NSIG 100
#endif

using namespace std;
using namespace mpqc;

//////////////////////////////////////////////////////////////////////
// static variables

static Debugger *signals[NSIG];

//////////////////////////////////////////////////////////////////////
// Debugger class definition

std::shared_ptr<Debugger> Debugger::default_debugger_(nullptr);

Debugger::Debugger(const char *exec) : Debugger(KeyVal()) {
  set_exec(exec);
}

Debugger::Debugger(const KeyVal &keyval) {
  init();

  debug_ = keyval.value<bool>("debug", true);
  traceback_ = keyval.value<bool>("traceback", true);
  exit_on_signal_ = keyval.value<bool>("exit", true);
  sleep_ = keyval.value<bool>("sleep", false);
  wait_for_debugger_ = keyval.value<bool>("wait_for_debugger", true);
  cmd_ = keyval.value<std::string>("cmd", std::string());
  if (cmd_.empty()) default_cmd();
  prefix_ = keyval.value<std::string>("prefix", std::string());
  handle_sigint_ = keyval.value<bool>("handle_sigint", true);
  if (keyval.value<bool>("handle_defaults", true)) handle_defaults();
}

Debugger::~Debugger() {
  for (int i = 0; i < NSIG; i++) {
    if (mysigs_[i]) signals[i] = 0;
  }
  delete[] mysigs_;
}

void Debugger::init() {
  exec_.resize(0);
  prefix_.resize(0);
  cmd_.resize(0);
  sleep_ = 0;

  exit_on_signal_ = 1;
  traceback_ = 1;
  debug_ = 1;
  wait_for_debugger_ = 1;

  mysigs_ = new int[NSIG];
  for (int i = 0; i < NSIG; i++) {
    mysigs_[i] = 0;
  }
}

namespace {
static void
handler(int sig)
{
  if (signals[sig]) signals[sig]->got_signal(sig);
}
}

void Debugger::handle(int sig) {
  if (sig >= NSIG) return;
#ifdef HAVE_SIGNAL
  typedef void (*handler_type)(int);
  signal(sig, (handler_type)handler);
#endif
  signals[sig] = this;
  mysigs_[sig] = 1;
}

void Debugger::handle_defaults() {
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
  if (handle_sigint_) handle(SIGINT);
#endif
#ifdef SIGHUP
  handle(SIGHUP);
#endif
#ifdef SIGBUS
  handle(SIGBUS);
#endif
#ifdef SIGABRT
  handle(SIGABRT);
#endif
}

void Debugger::set_exec(const char *exec) {
  if (exec) {
    exec_ = exec;
  } else {
    exec_.resize(0);
  }
}

void Debugger::set_prefix(const char *p) {
  if (p) {
    prefix_ = p;
  } else {
    prefix_.resize(0);
  }
}

void Debugger::set_prefix(int i) {
  char p[128];
  sprintf(p, "%3d: ", i);
  set_prefix(p);
}

void Debugger::default_cmd() {
  int has_x11_display = (getenv("DISPLAY") != 0);

  if (has_x11_display) {
    set_cmd("xterm -title \"$(PREFIX)$(EXEC)\" -e gdb $(EXEC) $(PID) &");
  } else {
    set_cmd(0);
  }
}

void Debugger::set_cmd(const char *cmd) {
  if (cmd) {
    cmd_ = cmd;
  } else {
    cmd_.resize(0);
  }
}

void Debugger::debug(const char *reason) {
  ExEnv::outn() << prefix_ << "Debugger::debug: ";
  if (reason)
    ExEnv::outn() << reason;
  else
    ExEnv::outn() << "no reason given";
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
    while ((pos = cmd.find(pidvar)) != std::string::npos) {
      std::string pidstr;
      pidstr += std::to_string(pid);
      cmd.replace(pos, pidvar.size(), pidstr);
    }
    std::string execvar("$(EXEC)");
    while ((pos = cmd.find(execvar)) != std::string::npos) {
      cmd.replace(pos, execvar.size(), exec_);
    }
    std::string prefixvar("$(PREFIX)");
    while ((pos = cmd.find(prefixvar)) != std::string::npos) {
      cmd.replace(pos, prefixvar.size(), prefix_);
    }
    // start the debugger
    ExEnv::outn() << prefix_ << "Debugger: starting \"" << cmd << "\"" << endl;
    debugger_ready_ = 0;
    const auto system_retvalue = system(cmd.c_str());
    if (system_retvalue != 0) { // call to system() failed
      ExEnv::outn() << prefix_ << "Failed debugger launch: system() did not succeed ..." << endl;
    }
    else {  // call to system() succeeded
      // wait until the debugger is ready
      if (sleep_) {
        ExEnv::outn() << prefix_ << "Sleeping " << sleep_
                      << " seconds to wait for debugger ..." << endl;
        sleep(sleep_);
      }
      if (wait_for_debugger_) {
        ExEnv::outn() << prefix_ << ": Spinning until debugger_ready_ is set ..."
                      << endl;
        while (!debugger_ready_)
          ;
      }
    }
  }
#endif
}

void Debugger::got_signal(int sig) {
  const char *signame;
  if (sig == SIGSEGV)
    signame = "SIGSEGV";
  else if (sig == SIGFPE)
    signame = "SIGFPE";
  else if (sig == SIGHUP)
    signame = "SIGHUP";
  else if (sig == SIGINT)
    signame = "SIGINT";
  else if (sig == SIGABRT)
    signame = "SIGABRT";
#ifdef SIGBUS
  else if (sig == SIGBUS)
    signame = "SIGBUS";
#endif
  else
    signame = "UNKNOWN SIGNAL";

  if (traceback_) {
    traceback(signame);
  }
  if (debug_) {
    debug(signame);
  }

  if (exit_on_signal_) {
    ExEnv::outn() << prefix_ << "Debugger: exiting" << endl;
    exit(1);
  } else {
    ExEnv::outn() << prefix_ << "Debugger: continuing" << endl;
  }

  // handle(sig);
}

void Debugger::set_debug_on_signal(int v) { debug_ = v; }

void Debugger::set_traceback_on_signal(int v) { traceback_ = v; }

void Debugger::set_wait_for_debugger(int v) { wait_for_debugger_ = v; }

void Debugger::set_exit_on_signal(int v) { exit_on_signal_ = v; }

void Debugger::set_default_debugger(const std::shared_ptr<Debugger> &d) {
  default_debugger_ = d;
}

std::shared_ptr<Debugger> Debugger::default_debugger() {
  return default_debugger_;
}

#define SIMPLE_STACK \
  (defined(linux) && defined(i386)) || (defined(__OSF1__) && defined(i860))

void Debugger::traceback(const char *reason) {
  Debugger::__traceback(prefix_, reason);
}

void Debugger::__traceback(const std::string &prefix, const char *reason) {
  detail::Backtrace result(prefix);
  const size_t nframes_to_skip = 2;
#if defined(HAVE_LIBUNWIND)
  ExEnv::outn() << prefix << "Debugger::traceback(using libunwind):";
#elif defined(HAVE_BACKTRACE)  // !HAVE_LIBUNWIND
  ExEnv::outn() << prefix << "Debugger::traceback(using backtrace):";
#else                          // !HAVE_LIBUNWIND && !HAVE_BACKTRACE
#if defined(SIMPLE_STACK)
  ExEnv::outn() << prefix << "Debugger::traceback:";
#else
  ExEnv::outn() << prefix << "traceback not available for this arch" << endl;
  return;
#endif  // SIMPLE_STACK
#endif  // HAVE_LIBUNWIND, HAVE_BACKTRACE

  if (reason)
    ExEnv::outn() << reason;
  else
    ExEnv::outn() << "no reason given";
  ExEnv::outn() << endl;

  if (result.empty())
    ExEnv::outn() << prefix << "backtrace returned no state information"
                  << std::endl;
  else
    ExEnv::outn() << result.str(nframes_to_skip) << std::endl;
}

/////////////////////////////////////////////////////////////////////////////

namespace mpqc {

void launch_gdb_xterm() {
  auto debugger = std::make_shared<mpqc::Debugger>();
  if (ExEnv::initialized()) debugger->set_exec(ExEnv::argv()[0]);
  debugger->debug("Starting gdb ...");
}

void launch_lldb_xterm() {
  auto debugger = std::make_shared<mpqc::Debugger>();
  if (ExEnv::initialized()) debugger->set_exec(ExEnv::argv()[0]);
  debugger->set_cmd("xterm -title \"$(PREFIX)$(EXEC)\" -e lldb -p $(PID) &");
  debugger->debug("Starting gdb ...");
}

}  // namespace mpqc

/////////////////////////////////////////////////////////////////////////////
// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
