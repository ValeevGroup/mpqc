//
// molden.h
//
// Copyright (C) 2013 MPQC authors
//
// Author: David Hollman <david.s.hollman@gmail.com
// Maintainer: DSH
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

#ifndef MOLDEN_H_
#define MOLDEN_H_

#include <util/misc/runnable.h>
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/wfn/orbitalspace.h>


// TODO get this right for different types of integrals
#define USE_OLD_SOLIDHARM_ORDERING 0
#if not USE_OLD_SOLIDHARM_ORDERING // CCA ordering
static inline int ipure(int l, int m) { return l+m; }
static inline int ipure_molden(int l, int m) { return m<0?2*-m:(m==0?0:2*m-1); }
#endif

namespace sc {

class WriteMolden: public Runnable {
  private:
    std::string filename_;
    Ref<OneBodyWavefunction> obwfn_;
    Ref<OrbitalSpace> mospace_;

    void initialize();

    Ref<Molecule> molecule() { return obwfn_->molecule(); }

    void write_atoms_section(std::ostream &out);
    void write_gto_section(std::ostream &out);
    void write_mo_section(std::ostream &out);

#if not USE_OLD_SOLIDHARM_ORDERING
    std::vector<int> bmap_;
#endif

  public:
    /** The KeyVal constructor

        <dl>
        <dt><tt>obwfn</tt></dt><dd> The OneBodyWavefunction whose orbitals are printed in the Molden file.
        There is no default for this option.</dd>

        <dt><tt>filename</tt></dt><dd> Specifies the filename of the file to
        write the output to. If "-" is given, the output will be written to the
        standard output. The default is "-".</dd>
        </dl>

     */
    WriteMolden(const Ref<KeyVal> &);

    /// Writes the molden file
    void run();


};

}

#endif /* MOLDEN_H_ */
