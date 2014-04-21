//
// orbital.h
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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

#ifndef _chemistry_qc_wfn_orbital_h
#define _chemistry_qc_wfn_orbital_h

#include <math/isosurf/volume.h>
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/wfn/orbitalspace.h>
#include <math/mmisc/grid.h>

namespace sc {

class Orbital: public Volume {
  protected:
    Ref<OneBodyWavefunction> wfn_;
    int orbital_;

    virtual void compute();
  public:
    Orbital(const Ref<KeyVal>&);
    Orbital(const Ref<OneBodyWavefunction>&, int orbital);
    ~Orbital();
    virtual void boundingbox(double valuemin,
                             double valuemax,
                             SCVector3& p1, SCVector3& p2);
};

/** The WriteOrbital class writes an orbital at
    user defined grid points to the standard output or to a separate file.*/
class WriteOrbital: public WriteGrid {
  protected:
    Ref<OneBodyWavefunction> obwfn_;
    int orbital_;

    void label(char* buffer);
    Ref<Molecule> get_molecule();
    double calculate_value(SCVector3 point);
    void initialize();

  public:
    /** The KeyVal constructor uses all keywords of WriteGrid, plus the following additional keywords:

        <dl>

        <dt><tt>obwfn</tt></dt><dd> The OneBodyWavefunction whose orbitals are calculated.
        There is no default for this option.</dd>

        <dt><tt>orbital</tt></dt><dd> Index of the orbital to be plotted. There is no default.</dd>

        </dl>

        N.B. Although WriteGrid requires keyword <tt>grid</tt> to be specified, omitting it for WriteOrbital
        will automatically construct a grid appropriate for the Molecule object of <tt>obwfn</tt>.
      */
    WriteOrbital(const Ref<KeyVal> &);
    ~WriteOrbital();

    /// if needed, creates default grid for the base class using the molecule of obwfn
    static Ref<KeyVal> process_keyval_for_base_class(const Ref<KeyVal> kv);

    /** constructs a default grid for the given molecule, using VDWShape for \c mol. The algorithm is very basic and does not
     * do any axis rotation, etc., so may be suboptimal for non-spherical systems.
     *
     * @param mol the Molecule object
     * @param resolution the grid voxel size
     * @param margin how much extra space to add around the VDWShape
     * @return Grid object that encloses \c mol
     */
    static Ref<Grid> make_default_grid(const Ref<Molecule>& mol,
                                       double resolution = 0.2,
                                       double margin = 1.0);
};

/** The WriteOrbitals class writes orbitals at
 *  user defined grid points to the standard output or to a separate file.
 */
class WriteOrbitals: public WriteVectorGrid {
  private:
    struct OrbitalMap : public DimensionMap {
      std::vector<int> map;
      int operator()(int o) const { return map[o]; }
      std::size_t ndim() const { return map.size(); }
    };
  protected:
    // provide either obwfn+omap or orbs
    Ref<OneBodyWavefunction> obwfn_;
    OrbitalMap omap_;
    Ref<OrbitalSpace> orbs_;

    void label(char* buffer);
    Ref<Molecule> get_molecule();
    void  calculate_values(const std::vector<SCVector3>& Points, std::vector<double>& Values);
    std::size_t ndim() const { return omap_.ndim(); }
    const DimensionMap& dimension_map() const { return omap_; }
    void initialize();
  public:
    /** The KeyVal constructor accepts keywords of WriteVectorGrid and the following additional keywords.
        <dl>

        <dt><tt>obwfn</tt></dt><dd> The OneBodyWavefunction whose orbitals are calculated.
        There is no default for this option.</dd>

        <dt><tt>first</tt></dt><dd> The index of the first orbital to be plotted. MOs are indexed according to their energy,
        hence <tt>first=1</tt> refers to the lowest-energy orbital. The default value is 1.</dd>

        <dt><tt>last</tt></dt><dd> The index of the last orbital to be plotted. By default the highest energy orbital is selected.
        MOs are indexed according to their energy, hence <tt>last=<# of orbitals></tt> refers to the highest-energy orbital.</dd>
        </dl>


        N.B. Although WriteVectorGrid requires keyword <tt>grid</tt> to be specified, omitting it for WriteOrbitals
        will automatically construct a grid appropriate for the Molecule object of <tt>obwfn</tt>.

      */
    WriteOrbitals(const Ref<KeyVal> &);
    /**
     * Evaluates orbitals \c orbs on \c grid and writes them to \c gridfile in format \c gridformat.
     * @param orbs the OrbitalSpace object that specifies the AO coefficients of the orbitals
     * @param labels vector of int's that maps orbitals in \c orbs to their "absolute" index; if empty,
     *               assume that \c orbs contains all orbitals (i.e. first orbital = 1, second = 2).
     * @param grid   Grid object
     * @param gridformat  output format for the grid data, the only supported value is "gaussian_cube"
     * @param gridfile    the file name to which the data will be written
     */
    WriteOrbitals(const Ref<OrbitalSpace> & orbs, const std::vector<int>& labels, const Ref<sc::Grid> & grig,
                  std::string gridformat, std::string gridfile);
    ~WriteOrbitals();
};


}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
