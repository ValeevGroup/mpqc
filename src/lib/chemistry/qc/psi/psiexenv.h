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

#include <string>
#include <sstream>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <chemistry/qc/psi/psiinput.h>
#include <chemistry/qc/psi/psifile11.h>
#include <chemistry/qc/wfn/spin.h>

namespace sc {

/// PsiExEnv specifies a Psi execution environment.
class PsiExEnv: public DescribedClass {

    /// The default constructor is identical to the KeyVal constructor with all parameters set to their default values.
    PsiExEnv();

    static Ref<PsiExEnv> default_instance_;

    // Static Psi info
    static std::string defaultinputname_;
    static std::string defaultoutputname_;
    static std::string file11name_;
    static int ckptfile_;

    // Defaults
    static std::string defaultpsiprefix_;
    static std::string defaultcwd_;
    static std::string defaultfileprefix_;
    static std::string defaultstdout_;
    static std::string defaultstderr_;

    // Calculation-specific info
    std::string psiprefix_;
    std::string cwd_;        // working directory where all files will be placed
    std::string inputname_;
    std::string outputname_;
    std::string fileprefix_;
    std::string stdout_;     // Standard output of psi modules
    std::string stderr_;     // Standard error of psi modules
    int nscratch_;
    std::vector<std::string> scratch_;
    Ref<PsiInput> psiinput_;
    Ref<PsiFile11> psifile11_;

    bool keep_output_;

    void config_psio();
    psi::PSIO psio_;
    psi::Chkpt* chkpt_;

    // Add the following to the PATH environmental variable
    void add_to_path(const std::string &);

    // task id
    int me_;

  public:
    /** A KeyVal constructor is used to generate a PsiExEnv
         object from the input. The full list of keywords
         that are accepted is below.

         <table border="1">

         <tr><td>%Keyword<td>Type<td>Default<td>Description

         <tr><td><tt>cwd</tt><td>string<td>.<td>the current working directory (small Psi files will be written here)

         <tr><td><tt>fileprefix</tt><td>string<td>filename.psi<td>the name prefix used for scratch files produced by Psi

         <tr><td><tt>stdout</tt><td>string<td>fileprefix.stdout<td>the file to which Psi standard output will be written

         <tr><td><tt>stderr</tt><td>string<td>fileprefix.stderr<td>the file to which Psi standard error will be written

         <tr><td><tt>scratch</tt><td>array of strings<td>[ <ConsumableResources::disk_location()> ]<td> the location to which large scratch files will be written

         </table>

         In addition, the location of Psi executables can be overridden by setting environmental variable <tt>PSIBIN</tt>
         to the desired value. By default, Psi executables in the location specified with the <tt>--with-psi</tt> configure
         script option will be used.
      */
    PsiExEnv(const Ref<KeyVal>&);
    ~PsiExEnv();

    static const Ref<PsiExEnv>& get_default_instance();

    /// Returns the PsiInput object which PsiExEnv uses
    Ref<PsiInput> psi_input() const { return psiinput_; }
    /// Returns the PsiFile11 object which PsiExEnv uses
    Ref<PsiFile11> psi_file11() const { return psifile11_; }
    /// Creates the PsiInput object which PsiExEnv uses
    Ref<PsiInput> get_psi_input();
    /// Creates the PsiFile11 object which PsiExEnv uses
    Ref<PsiFile11> get_psi_file11();

    /// Executes Psi input+driver.
    /// \param skip_input whether to skip running "input" module. The default is false.
    /// \sa PsiExEnv::run_psi_module()
    void run_psi(bool skip_input = false);
    /// Executes a Psi module using a system call. Throws if psi fails.
    void run_psi_module(const char * module, const std::vector<std::string>& args = std::vector<std::string>());
    /// cleans Psi scratch files using the same approach as psiclean.
    /// \param fullclean if set to true, clean out all files including the checkpoint file
    void run_psiclean(bool fullclean = true);

    /// Returns current working directory
    const std::string& get_cwd() const { return cwd_;};
    /// Returns the Psi file prefix
    const std::string& get_fileprefix() const { return fileprefix_; };
    /// Returns the number of scratch locations
    int get_nscratch() const { return nscratch_; };
    /// Returns the ith scratch location
    const std::string& get_scratch(int i) const { return scratch_[i]; };

    /// Returns an instance of psi::PSIO
    psi::PSIO& psio() { return psio_; }
    /// Returns an instance of psi::Chkpt
    psi::Chkpt& chkpt();

    void print(std::ostream&o=ExEnv::out0()) const;

    /// this will cause psi driver to keep the output next time it's run. Must be called every time to prevent the output file from truncation
    void keep_output();
};

/// PsiChkpt know to read data from Psi checkpoint file and convert it to conform to the representations expected in MPQC
class PsiChkpt : public RefCount {
  public:
    /// Assume environment described by exenv; integrals specifies the basis set and the factory with conventions compatible with Psi
    PsiChkpt(const Ref<PsiExEnv>& exenv,
             const Ref<Integral>& integral,
             int debug);
    ~PsiChkpt();

    int debug() const { return debug_; }

    const Ref<PsiExEnv>& exenv() const { return exenv_; }
    const Ref<Integral>& integral() const { return integral_; }

    /// are there spin-specific orbitals?
    bool have_spin_unrestricted_mos() const;
    /// read the orbital energies for spincase s and return as a diagonal matrix
    RefDiagSCMatrix evals(SpinCase1 s,
                          bool spin_restricted = true) const;
    /// read the orbital coefficients for spincase s
    RefSCMatrix coefs(SpinCase1 s,
                      bool spin_restricted = true) const;
    /// Returns a map from shells in Psi3 basis to std::pair<shell,contraction> in MPQC basis (note that Psi3 does not handle general contractions)
    std::vector< std::pair<unsigned int,unsigned int> > shell_map() const;
    /// Returns a map from AO in Psi3 basis to AO in MPQC basis
    std::vector<unsigned int> ao_map() const;

  private:
    Ref<PsiExEnv> exenv_;
    Ref<Integral> integral_;
    int debug_;
};

}

#endif
