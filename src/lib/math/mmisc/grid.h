//
// grid.h
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

#ifndef _util_misc_grid_h
#define _util_misc_grid_h

#include <chemistry/molecule/molecule.h>
#include <util/class/class.h>
#include <util/misc/runnable.h>
#include <util/misc/units.h>


using boost::property_tree::ptree;

namespace sc {

/** The Grid class defines a finite regular Carthesian grid.
    The grid can be 0, 1, 2 or 3 dimensional. */
class Grid: public DescribedClass {
  public:
    int numx;
    int numy;
    int numz;
    SCVector3 origin;
    SCVector3 axisx;
    SCVector3 axisy;
    SCVector3 axisz;
    Ref<Units> unit;
    /** The KeyVal constructor.
    
        <dl>

        <dt><tt>numx</tt></dt><dd> The number of voxels along axisx
        (defined below). The default value is 1.</dd>

        <dt><tt>numy</tt></dt><dd> The number of voxels along axisy
        (defined below). The default value is 1.</dd>

        <dt><tt>numz</tt></dt><dd> The number of voxels along axisz
        (defined below). The default value is 1.</dd>

        <dt><tt>origin</tt></dt><dd> The origin of the grid. The default is
        [0 0 0].</dd>

        <dt><tt>axisx</tt></dt><dd> Defines volume element length and direction for the first axis of the grid.
        The default is [1 0 0].</dd>

        <dt><tt>axisy</tt></dt><dd> Same as <tt>axisx</tt>, but for the second axis. The default is [0 1 0].</dd>

        <dt><tt>axisz</tt></dt><dd> Same as <tt>axisx</tt>, but for the third axis. The default is [0 0 1].</dd>

        <dt><tt>unit</tt></dt><dd> The unit in which the parameters are given.
        Notice that this does not determine the unit used in the output file.
        The default value is bohr.</dd>

        </dl> */
    Grid(const Ref<KeyVal> &);
    Grid(int numx = 1, int numy = 1, int numz = 1,
         SCVector3 origin = SCVector3(0.0, 0.0, 0.0),
         SCVector3 axisx = SCVector3(1.0, 0.0, 0.0),
         SCVector3 axisy = SCVector3(0.0, 1.0, 0.0),
         SCVector3 axisz = SCVector3(0.0, 0.0, 1.0),
         Ref<Units> unit = new Units("bohr")
         );
};

/** The abstract WriteGrid class provides an interface for writing the value
    of a scalar function evaluated at a given set of grid points to a file. */
class WriteGrid: public Runnable {
  private:
    void wf_mpqc(std::ostream &out);
    void wf_gaussian_cube(std::ostream &out);
    void wf_vtk2(std::ostream &out);
    void wf_mpqc_raw(std::ostream &out);
  protected:
    std::string filename_;
    Ref<Grid> grid_;
    std::string format_;
    void (WriteGrid::*write_format_)(std::ostream &out);
    /** Prepares some pre-caculated values before the repetitive grid
        calculations are perfomed. */
    virtual void initialize() = 0;
    /** A label that identifies the scalar function evaluated at the grid
        points, is written to the buffer argument. The classname, concatenated
        with some important properties should be sufficient. No whitespace
        allowed, length of the string is limited to 256 characters. */
    virtual void label(char* buffer) = 0;
    /// Returns the molecule around which the grid values are calculated
    virtual Ref<Molecule> get_molecule() = 0;
    /// Returns the value of the scalar function at the given coordinate.
    virtual double calculate_value(SCVector3 point) = 0;
  public:
    /** The KeyVal constructor.
        <dl>

        <dt><tt>grid</tt></dt><dd> A Grid that specifies the grid points at
        which the scalar function should be calculated.</dd>

        <dt><tt>filename</tt></dt><dd> Specifies the filename of the file to
        write the output to. If "-" is given, the output will be written to the
        standard output. The default is "-".</dd>
        
        <dt><tt>format</tt></dt><dd> The format in which the grid data is to be
        written. Currently four formats have been implemented:
          <ul>
            <li><tt>mpqc</tt>: A very comprehensive format.</li>
            <li><tt>gaussian_cube</tt>: The format used by Gaussian.</li>
            <li><tt>vtk2</tt>: This format is usefull for vizualizing the grid
            data with the vtk library. (http://www.vtk.org/)</li>
            <li><tt>mpqc_raw</tt>: A very simple format that contains both the
            coordinates of the grid points and the value of the scalar function
            at each point.</li>
          </ul>
        </dd>
        </dl> */
    WriteGrid(const Ref<KeyVal> &);
    /// Writes the grid data.
    void run();
};

/**
  WriteVectorGrid provides an interface for writing the value
  of a vector function evaluated at a given set of grid points to a file (compare to WriteGrid).
 */
class WriteVectorGrid: public Runnable, virtual public DescribedXMLWritable {
  public:
    // see wf_gaussian_cube
    struct DimensionMap {
      virtual int operator()(int d) const {
        return d;
      }
    };

  private:
    /** for now, only Gaussian cube file format is implemented; other formats can be done similarly when necessary.
     *
     *  @param DimensionMap vector's dimensions are assumed to be a noncontiguous subset of a larger (super)set of dimensions.
     *         DimensionMap is a functor that maps dimensions to the superset.
     */
    void wf_gaussian_cube(std::ostream &out, const DimensionMap& dmap);
  protected:
    std::string filename_;
    Ref<Grid> grid_;
    std::string format_;
    void (WriteVectorGrid::*write_format_)(std::ostream &out, const DimensionMap& dmap);

    /** Prepares some pre-calculated values before the repetitive grid calculations
    are perfomed.*/
    virtual void initialize() = 0;
    /** A label that identifies the scalar function evaluated at the grid
        points, is written to the buffer argument. The classname, concatenated
        with some important properties should be sufficient. No whitespace
        allowed, length of the string is limited to 256 characters. */
    virtual void label(char* buffer) = 0;
    /// Returns the molecule around which the grid values are calculated
    virtual Ref<Molecule> get_molecule() = 0;
    /// Returns the value of the vector function at the given coordinate.
    virtual void calculate_values(const std::vector<SCVector3>& Points, std::vector<double>& Vals)=0;
    /// number of dimensions of the vector
    virtual std::size_t ndim() const =0;
    /// dimension map \sa DimensionMap
    virtual const DimensionMap& dimension_map() const =0;
  public:
    /** The KeyVal constructor.
        <dl>

        <dt><tt>grid</tt></dt><dd> A Grid that specifies the grid points at
        which the scalar function should be calculated.</dd>

        <dt><tt>filename</tt></dt><dd> Specifies the filename of the file to
        write the output to. If it is not given, the output will be written to the
        standard output. The default is "Grid".</dd>

        <dt><tt>format</tt></dt><dd> The format in which the grid data is to be
        written. Currently only guassian_cube format has been implemented:
          <ul>
            <li><tt>gaussian_cube</tt>: The format used by Gaussian.</li>
          </ul>
        </dd>
        </dl> */
    WriteVectorGrid(const Ref<KeyVal> &);
    WriteVectorGrid(const Ref<sc::Grid> & grid,
                    std::string gridformat, std::string gridfile);
    /// Writes the grid data.
    void run();
    virtual ptree& write_xml(ptree& parent, const XMLWriter& writer);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
