//
// psiexenv.h
//
// Copyright (C) 2002 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _chemistry_qc_psi_exenv_h
#define _chemistry_qc_psi_exenv_h

using namespace std;

#include <string>
#include <sstream>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <chemistry/qc/psi/psiinput.h>
#include <chemistry/qc/psi/psifile11.h>

namespace sc {

/// PsiExEnv specifies a Psi calculation

class PsiExEnv: public DescribedClass {

    // Static Psi info
    static string defaultinputname_;
    static string defaultoutputname_;
    static string file11name_;
    static int ckptfile_;

    // Defaults
    static string defaultpsiprefix_;
    static string defaultcwd_;
    static string defaultfileprefix_;
    static string defaultstdout_;
    static string defaultstderr_;

    // Calculation-specific info
    string psiprefix_;
    string cwd_;        // working directory where all files will be placed
    string inputname_;
    string outputname_;
    string fileprefix_;
    string stdout_;     // Standard output of psi modules
    string stderr_;     // Standard error of psi modules
    int nscratch_;
    string *scratch_;
    Ref<PsiInput> psiinput_;
    Ref<PsiFile11> psifile11_;

    void config_psio();
    psi::PSIO psio_;
    psi::Chkpt* chkpt_;

    // Add the following to the PATH environmental variable
    void add_to_path(const string &);

  public:
    PsiExEnv(const Ref<KeyVal>&);
    PsiExEnv(char *cwd, char *fileprefix, int nscratch, char **scratch);
    ~PsiExEnv();

    /// Returns the PsiInput object which PsiExEnv uses
    Ref<PsiInput> get_psi_input() const { return psiinput_;};
    /// Returns the PsiFile11 object which PsiExEnv uses
    Ref<PsiFile11> get_psi_file11() const { return psifile11_;};
    
    /// Executes Psi input+driver
    int run_psi();
    /// Executes a Psi module
    int run_psi_module(const char *);

    /// Returns current working directory
    string get_cwd() const { return cwd_;};
    /// Returns the Psi file prefix
    string get_fileprefix() const { return fileprefix_; };
    /// Returns the number of scratch locations
    int get_nscratch() const { return nscratch_; };
    /// Returns the ith scratch location
    string get_scratch(int i) const { return scratch_[i]; };

    /// Returns an instance of psi::PSIO
    psi::PSIO& psio() { return psio_; }
    /// Returns an instance of psi::Chkpt
    psi::Chkpt& chkpt();
    
    void print(std::ostream&o=ExEnv::out0()) const;
};

}

#endif
