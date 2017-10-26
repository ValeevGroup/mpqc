//
// bug.h
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

#ifndef MPQC4_SRC_MPQC_UTIL_MISC_BUG_H_
#define MPQC4_SRC_MPQC_UTIL_MISC_BUG_H_

#include <string>
#include <vector>
#include "mpqc/util/keyval/keyval.h"

namespace mpqc {

/**
 * The Debugger class describes what should be done when a catastrophic
 * error causes unexpected program termination.  It can try things such as
 * start a debugger running where the program died or it can attempt to
 * produce a stack traceback showing roughly where the program died.  These
 * attempts will not always succeed.
 */
class Debugger: virtual public DescribedClass {
  protected:
    std::string prefix_;
    std::string exec_;
    std::string cmd_;
    volatile int debugger_ready_;

    bool debug_;
    bool traceback_;
    bool exit_on_signal_;
    bool sleep_;
    bool wait_for_debugger_;
    bool handle_sigint_;
    int *mysigs_;

    void init();

    static std::shared_ptr<Debugger> default_debugger_;

    /** prints out a backtrace
     *
     * @param prefix this string will be prepended at the beginning of each line of Backtrace
     * @param reason optional string specifying the reason for traceback
     * @return backtrace
     */
    static void __traceback(const std::string& prefix, const char *reason = 0);
  public:
    /** \brief The KeyVal constructor.

        The KeyVal object will be queried for the following keywords:
        <dl>
        <dt><tt>debug</tt><dd> Try to start a debugger when an error occurs.
        Doesn't work on all machines. The default is true, if possible.

        <dt><tt>traceback</tt><dd> Try to print out a traceback extracting
        return addresses from the call stack.  Doesn't work on most
        machines.  The default is true, if possible.

        <dt><tt>exit</tt><dd> Exit on errors.  The default is true.

        <dt><tt>wait_for_debugger</tt><dd> When starting a debugger go into
        an infinite loop to give the debugger a chance to attach to the
        process.  The default is true.

        <dt><tt>sleep</tt><dd> When starting a debugger wait this many
        seconds to give the debugger a chance to attach to the process.
        The default is 0.

        <dt><tt>handle_defaults</tt><dd> Handle a standard set of signals
        such as SIGBUS, SIGSEGV, etc.  The default is true.

        <dt><tt>prefix</tt><dd> Gives a string that is printed before each
        line that is printed by Debugger. The default is nothing.

        <dt><tt>cmd</tt><dd> Gives a command to be executed to start the
        debugger.  If environment variable <tt>DISPLAY</tt> is set, the default
        command is <tt>"xterm -title \"$(PREFIX)$(EXEC)\" -e gdb $(EXEC) $(PID) &"</tt>
        where token <tt>$(PREFIX)</tt> is replaced by the value of keyword <tt>prefix</tt>,
        <tt>$(EXEC)</tt> is replaced by the name of the current executable,
        and <tt>$(PID)</tt> is the process ID of each process. This will launch
        one <tt>xterm</tt> per MPI rank. There is no default if <tt>DISPLAY</tt> is not set (hence
        no debugger will be launched, only an informational statement will be printed). N.B. To use
        LLDB instead of GDB set <tt>cmd</tt> to <tt>"xterm -title \"$(PREFIX)$(EXEC)\" -e lldb -p $(PID) &"</tt>.

        </dl> */
    Debugger(const KeyVal&);
    /** @brief Programmatic construction of Debugger
     *
     * This is equivalent to:
     * @code
     * Debugger d(KeyVal());
     * d.set_exec(exec);
     * @endcode
     *
     * @param exec the executable name
     */
    Debugger(const char *exec = 0);
    ~Debugger();

    /** The debug member attempts to start a debugger
        running on the current process. */
    virtual void debug(const char *reason = 0);
    /** The traceback member attempts to produce a Backtrace
     for the current process.  A symbol table must be saved for
     the executable if any sense is to be made of the traceback.
     This feature is available on platforms with (1) libunwind,
     (2) backtrace, or (3) certain platforms with hardwired unwinding.
     @param reason optional string specifying the reason for traceback
     */
    virtual void traceback(const char *reason = 0);
    /// Turn on or off debugging on a signel.  The default is on.
    virtual void set_debug_on_signal(int);
    /// Turn on or off traceback on a signel.  The default is on.
    virtual void set_traceback_on_signal(int);
    /// Turn on or off exit after a signel.  The default is on.
    virtual void set_exit_on_signal(int);
    /** Turn on or off running an infinite loop after the debugger is started.
        This loop gives the debugger a chance to attack to the process.
        The default is on. */
    virtual void set_wait_for_debugger(int);

    /// The Debugger will be actived when sig is caught.
    virtual void handle(int sig);
    /// This calls handle(int) with all of the major signals.
    virtual void handle_defaults();

    /// This sets a prefix which preceeds all messages printing by Debugger.
    virtual void set_prefix(const char *p);
    /// Set the prefix to the decimal represention of p followed by a ": ".
    virtual void set_prefix(int p);

    /** Sets the command to be exectuted when debug is called.
        The character sequence "$(EXEC)" is replaced by the executable
        name (see set_exec), "$(PID)" is replaced by the
        current process id, and "$(PREFIX)" is replaced by the
        prefix. */
    virtual void set_cmd(const char *);
    /// Calls set_cmd with a hopefully suitable default.
    virtual void default_cmd();
    /** Set the name of the executable for the current process.
        It is up to the programmer to set this, even if the Debugger
        is initialized with the KeyVal constructor. */
    virtual void set_exec(const char *);

    /// Called when signal sig is received.  This is mainly for internal use.
    virtual void got_signal(int sig);

    /// Set the global default debugger.  The initial value is null.
    static void set_default_debugger(const std::shared_ptr<Debugger> &);
    /// Return the global default debugger.
    static std::shared_ptr<Debugger> default_debugger();
};

}

#endif  // MPQC4_SRC_MPQC_UTIL_MISC_BUG_H_

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
