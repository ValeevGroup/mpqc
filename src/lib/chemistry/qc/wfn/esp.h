//
// esp.h
//
// Copyright (C) 2006 Toon Verstraelen.
//
// Author: Toon Verstraelen
//
// This file is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// This file is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with the MPQC; see the file COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#ifndef _chemistry_qc_wfn_esp_h
#define _chemistry_qc_wfn_esp_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/molecule/molecule.h>
#include <math/mmisc/grid.h>

namespace sc {

/** The WriteElectrostaticPotential class writes the electrostatic potential at
    user defined grid points to the standard output or to a separate file.*/
class WriteElectrostaticPotential: public WriteGrid {
  protected:
    Ref<Wavefunction> wfn_;
    Ref<SymmSCMatrix> ao_density_;
    Ref<SymmSCMatrix> pc_mat_;
    bool electronic_;
    bool nuclear_;

    void initialize();
    void label(char* buffer);
    Ref<Molecule> get_molecule();
    double calculate_value(SCVector3 point);
  public:
    /** The KeyVal constructor
        
        <dl>

        <dt><tt>wfn</tt></dt><dd> The Wavefunction of which the electrostatic
        potential is calculated. There is no default for this option.</dd>

        <dt><tt>nuclear</tt></dt><dd> Wether the nuclear terms should be
        included in the electrostatic potential. The default is yes.</dd>

        <dt><tt>electronic</tt></dt><dd> Wether the electronic terms should be
        included in the electrostatic potential. The default is yes.</dd>

        </dl> */
    WriteElectrostaticPotential(const Ref<KeyVal> &);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
