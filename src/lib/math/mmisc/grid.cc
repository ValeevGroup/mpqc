//
// grid.cc
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

#include <stdexcept>

#include <math/scmat/vector3.h>
#include <math/scmat/matrix.h>
#include <math/mmisc/grid.h>
#include <util/misc/scexception.h>
#ifdef MPQC_NEW_RUNTIME
#  include <util/misc/xmlwriter.h>
#endif

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////////
// GridDefinition

static ClassDesc Grid_cd(
    typeid(Grid),"Grid",1,
    "public DescribedClass", 0, create<Grid>, 0);

Grid::Grid(const Ref<KeyVal> &keyval):
  origin(0.,0.,0.),
  axisx(1.,0.,0.),
  axisy(0.,1.,0.),
  axisz(0.,0.,1.)
{
  KeyValValueint default_numx(1);
  numx = keyval->intvalue("numx", default_numx);
  KeyValValueint default_numy(1);
  numy = keyval->intvalue("numy", default_numy);
  KeyValValueint default_numz(1);
  numz = keyval->intvalue("numz", default_numz);

  if (keyval->exists("origin")) {
      for (int i=0; i<3; i++) {
          origin[i] = keyval->doublevalue("origin",i);
        }
    }
  if (keyval->exists("axisx")) {
      for (int i=0; i<3; i++) {
          axisx[i] = keyval->doublevalue("axisx",i);
        }
    }
  if (keyval->exists("axisy")) {
      for (int i=0; i<3; i++) {
          axisy[i] = keyval->doublevalue("axisy",i);
        }
    }
  if (keyval->exists("axisz")) {
      for (int i=0; i<3; i++) {
          axisz[i] = keyval->doublevalue("axisz",i);
        }
    }

  if (keyval->exists("unit")) {
      std::string tmp = keyval->stringvalue("unit");
      unit = new Units(tmp.c_str());
    }
  else {
      unit = new Units("bohr");
    }
}

Grid::Grid(int nx, int ny, int nz,
           SCVector3 o,
           SCVector3 ax,
           SCVector3 ay,
           SCVector3 az,
           Ref<Units> u) :
           numx(nx), numy(ny), numz(nz),
           origin(o),
           axisx(ax),
           axisy(ay),
           axisz(az),
           unit(u)
{
}

/////////////////////////////////////////////////////////////////////////////
// WriteGrid

static ClassDesc WriteGrid_cd(
    typeid(WriteGrid),"WriteGrid",1,
    "public Runnable");

WriteGrid::WriteGrid(const Ref<KeyVal> &keyval)
{
  grid_ << keyval->describedclassvalue("grid");
  if (grid_.null()) {
      InputError ex("valid \"grid\" missing",
                    __FILE__, __LINE__, "grid", "(null)", class_desc());
      try {
          ex.elaborate()
              << "WriteGrid KeyVal ctor requires"
              << " that \"grid\" specifies an object"
              << " of type Grid" << std::endl;
        }
      catch (...) {}
      throw ex;
    }
    
  if (keyval->exists("filename")) {
      filename_ = keyval->stringvalue("filename");
    }
  else {
      filename_ = "-";
    }

  if (keyval->exists("format")) {
      format_ = keyval->stringvalue("format");
      if (format_ == "mpqc") {
          write_format_ = &WriteGrid::wf_mpqc;
        }
      else if (format_ == "gaussian_cube") {
          write_format_ = &WriteGrid::wf_gaussian_cube;
        }
      else if (format_ == "vtk2") {
          write_format_ = &WriteGrid::wf_vtk2;
        }
      else if (format_ == "mpqc_raw") {
          write_format_ = &WriteGrid::wf_mpqc_raw;
        }
      else {
          InputError ex("valid \"format\" missing",
                        __FILE__, __LINE__, "format", "(null)", class_desc());
          try {
              ex.elaborate()
                  << "WriteElectronGrid KeyVal ctor requires"
                  << " that \"format\" is one of \"mpqc\", "
                  << " \"gaussian_cube\", \"vtk2\" or \"mpqc_raw\""
                  << ". The requested format was \""
                  << format_ << "\"." << std::endl;
            }
          catch (...) {}
          throw ex;
        }
    }
  else {
      format_ = "mpqc";
      write_format_ = &WriteGrid::wf_mpqc;
    }
}

void
WriteGrid::run()
{
  initialize();
  
  std::ostream *out;
  if (filename_ == "-") {
      out = &(ExEnv::out0());
    }
  else {
      char buffer[256];
      label(buffer);
      ExEnv::out0() << incindent << indent << buffer
                    << " is writing its output to \"" << filename_
                    << "\" using the \"" << format_
                    << "\" format." << std::endl;
      ExEnv::out0() << decindent;
      out = new std::ofstream(filename_.c_str());
    }

  // this is where the writting operation is done, by calling one of the next few functions.
  (*this.*write_format_)(*out);

  if (filename_ == "-") {
      *out << decindent;
    }
  else {
      delete out;
    }
}

void
WriteGrid::wf_mpqc(std::ostream &out) {
  double conv = grid_->unit->to_atomic_units();
  
  char buffer[256];
  label(buffer);
  out << buffer << ":" << std::endl;
  out << incindent;
  out << indent << "unit = " << grid_->unit->string_rep() << std::endl;

  out << indent << "origin = [";
  for (int i=0; i<3; i++) out << " " << grid_->origin[i];
  out << "]" << std::endl;

  out << indent << "numx = " << grid_->numx << std::endl;
  out << indent << "axisx = [";
  for (int i=0; i<3; i++) out << " " << grid_->axisx[i];
  out << "]" << std::endl;

  out << indent << "numy = " << grid_->numy << std::endl;
  out << indent << "axisy = [";
  for (int i=0; i<3; i++) out << " " << grid_->axisy[i];
  out << "]" << std::endl;

  out << indent << "numz = " << grid_->numz << std::endl;
  out << indent << "axisz = [";
  for (int i=0; i<3; i++) out << " " << grid_->axisz[i];
  out << "]" << std::endl;

  SCVector3 pointx;
  SCVector3 pointy;
  SCVector3 pointz;
  
  for (int i=0; i<grid_->numx; i++) {
      pointx = grid_->origin + i * grid_->axisx;
      for (int j=0; j<grid_->numy; j++) {
          pointy = pointx + j * grid_->axisy;
          for (int k=0; k<grid_->numz; k++) {
              pointz = pointy + k * grid_->axisz;
              out << indent
                  << scprintf("%16.12f", calculate_value(pointz*conv))
                  << std::endl;
            }
        }
    }  
}

void
WriteGrid::wf_gaussian_cube(std::ostream &out) {
  double to_atomic = grid_->unit->to_atomic_units();
  Ref<Units> angstrom = new Units("angstrom");
  double to_angstrom = angstrom->from_atomic_units();
  Ref<Molecule> mol = get_molecule();
  
  out << "Gaussian cube file generated by MPQC." << std::endl;
  char buffer[80];
  label(buffer);
  out << buffer << std::endl;
  out << std::endl;
  
  out.fill(' ');
  out << std::fixed << std::setprecision(6);
  
  out << std::setw( 4) << mol->natom() 
      << std::setw(12) << grid_->origin[0]*to_atomic
      << std::setw(12) << grid_->origin[1]*to_atomic
      << std::setw(12) << grid_->origin[2]*to_atomic << std::endl;
  out << std::setw( 4) << grid_->numx
      << std::setw(12) << grid_->axisx[0]*to_atomic
      << std::setw(12) << grid_->axisx[1]*to_atomic
      << std::setw(12) << grid_->axisx[2]*to_atomic << std::endl;
  out << std::setw( 4) << grid_->numy
      << std::setw(12) << grid_->axisy[0]*to_atomic
      << std::setw(12) << grid_->axisy[1]*to_atomic
      << std::setw(12) << grid_->axisy[2]*to_atomic << std::endl;
  out << std::setw( 4) << grid_->numz
      << std::setw(12) << grid_->axisz[0]*to_atomic
      << std::setw(12) << grid_->axisz[1]*to_atomic
      << std::setw(12) << grid_->axisz[2]*to_atomic;
      
  for (int atom=0; atom<mol->natom(); atom++) {
      out << std::endl
          << std::setw( 4) << mol->Z(atom)
          << std::setw(12) << 0.0  // this value is expected, but has no use.
          << std::setw(12) << mol->r(atom, 0)
          << std::setw(12) << mol->r(atom, 1)
          << std::setw(12) << mol->r(atom, 2);
    }
    
  SCVector3 pointx;
  SCVector3 pointy;
  SCVector3 pointz;
  
  out << std::scientific << std::uppercase << std::setprecision(5);
  for (int i=0; i<grid_->numx; i++) {
      pointx = grid_->origin + i * grid_->axisx;
      for (int j=0; j<grid_->numy; j++) {
          pointy = pointx + j * grid_->axisy;
          for (int k=0; k<grid_->numz; k++) {
              if (k%6==0) out << std::endl;
              pointz = pointy + k * grid_->axisz;
              out << setw(13)
                  << calculate_value(pointz*to_atomic);
            }
        }
    }  
}

void
WriteGrid::wf_vtk2(std::ostream &out) {
  double to_atomic = grid_->unit->to_atomic_units();

  SCVector3 pointx;
  SCVector3 pointy;
  SCVector3 pointz;
  
  char buffer[256];
  label(buffer);
  int num = grid_-> numx * grid_-> numy * grid_-> numz;

  out << "# vtk DataFile Version 2.0" << std::endl;
  out << "cube file: " << buffer << std::endl;
  out << "ASCII" << std::endl;
  out << "DATASET STRUCTURED_GRID" << std::endl;
  out << "DIMENSIONS " << grid_-> numx << " " << grid_-> numy 
      << " " << grid_-> numz << std::endl;
  out << "POINTS " << num << " float" << std::endl;
  
  for (int i=0; i<grid_->numx; i++) {
      pointx = grid_->origin + i * grid_->axisx;
      for (int j=0; j<grid_->numy; j++) {
          pointy = pointx + j * grid_->axisy;
          for (int k=0; k<grid_->numz; k++) {
              pointz = pointy + k * grid_->axisz;
              out << pointz[0] << " " << pointz[1] << " "
                  << pointz[2] << std::endl;
            }
        }
    }
      
  out << "POINT_DATA " << num << std::endl;
  out << "SCALARS " << buffer << " float 1" << std::endl;
  out << "LOOKUP_TABLE default" << std::endl;
  
  out << std::scientific << std::uppercase << std::setprecision(5);
  for (int i=0; i<grid_->numx; i++) {
      pointx = grid_->origin + i * grid_->axisx;
      for (int j=0; j<grid_->numy; j++) {
          pointy = pointx + j * grid_->axisy;
          for (int k=0; k<grid_->numz; k++) {
              pointz = pointy + k * grid_->axisz;
              out << calculate_value(pointz*to_atomic) << std::endl;
            }
        }
    }  
}

void
WriteGrid::wf_mpqc_raw(std::ostream &out) {
  double to_atomic = grid_->unit->to_atomic_units();

  SCVector3 pointx;
  SCVector3 pointy;
  SCVector3 pointz;
  
  char buffer[256];
  label(buffer);
  int num = grid_-> numx * grid_-> numy * grid_-> numz;

  out << "# MPQC raw grid data (atomic units)" << std::endl;
  out << "# " << buffer << std::endl;
  out << "# Number of records: " << num << std::endl;
  
  for (int i=0; i<grid_->numx; i++) {
      pointx = grid_->origin + i * grid_->axisx;
      for (int j=0; j<grid_->numy; j++) {
          pointy = pointx + j * grid_->axisy;
          for (int k=0; k<grid_->numz; k++) {
              pointz = (pointy + k * grid_->axisz)*to_atomic;
              out << pointz[0] << " " 
                  << pointz[1] << " "
                  << pointz[2] << " "
                  << calculate_value(pointz) << std::endl;
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// WriteVectorGrid

static ClassDesc WriteVectorGrid_cd(
    typeid(WriteVectorGrid),"WriteVectorGrid",1,
    "public Runnable");

WriteVectorGrid::WriteVectorGrid(const Ref<KeyVal> &keyval)
{
  grid_ << keyval->describedclassvalue("grid");
  if (grid_.null()) {
      InputError ex("valid \"grid\" missing",
                    __FILE__, __LINE__, "grid", "(null)", class_desc());
      try {
          ex.elaborate()
              << "WriteVectorGrid KeyVal ctor requires"
              << " that \"grid\" specifies an object"
              << " of type Grid" << std::endl;
        }
      catch (...) {}
      throw ex;
    }

  filename_ = keyval->stringvalue("filename", KeyValValuestring(""));
  if (filename_ == "")
    filename_ = SCFormIO::fileext_to_filename_string(".grid");

  if (keyval->exists("format")) {
      format_ = keyval->stringvalue("format");
      if (format_ == "gaussian_cube") {
          write_format_ = &WriteVectorGrid::wf_gaussian_cube;
        }
      else {
          InputError ex("valid \"format\" missing",
                        __FILE__, __LINE__, "format", "(null)", class_desc());
          try {
              ex.elaborate()
                  << "WriteVectorGrid KeyVal ctor requires"
                  << " that \"format\" is \"gaussian_cube\". "
                  << "The requested format was \""
                  << format_ << "\"." << std::endl;
            }
          catch (...) {}
          throw ex;
        }
    }
  else {
      format_ = "gaussian_cube";
      write_format_ = &WriteVectorGrid::wf_gaussian_cube;
    }
}



WriteVectorGrid::WriteVectorGrid(const Ref<sc::Grid> & grid,
                                 std::string gridformat,
                                 std::string gridfile) :
    grid_(grid),
    format_(gridformat),
    filename_(gridfile)
{
  if (format_ == "gaussian_cube") {
      write_format_ = &WriteVectorGrid::wf_gaussian_cube;
    }
  else {
    ProgrammingError("WriteVectorGrid: unrecognized format", __FILE__, __LINE__);
  }
  if (filename_ == "")
    filename_ = SCFormIO::fileext_to_filename_string(".grid");
}

void
WriteVectorGrid::run()
{
  initialize();

  std::ostream *out;
  if (filename_ == "-") {
      out = &(ExEnv::out0());
    }
  else {
      char buffer[256];
      label(buffer);
      ExEnv::out0() << incindent << indent << buffer
                    << " is writing its output to \"" << filename_
                    << "\" using the \"" << format_
                    << "\" format." << std::endl;
      ExEnv::out0() << decindent;
      out = new std::ofstream(filename_.c_str());
    }

  // this is where the writting operation is done, by calling one of the next few functions.
  (*this.*write_format_)(*out, this->dimension_map());
 // this->wf_gaussian_cube(*out);

  if (filename_ == "-") {
      *out << decindent;
    }
  else {
      delete out;
    }
}

void
WriteVectorGrid::wf_gaussian_cube(std::ostream &out, const DimensionMap& dmap) {
  const double to_atomic = grid_->unit->to_atomic_units();
  Ref<Units> angstrom = new Units("angstrom");
  const double to_angstrom = angstrom->from_atomic_units();
  Ref<Molecule> mol = get_molecule();
  const int nd = this->ndim();

  out << "Gaussian cube file generated by MPQC." << std::endl;
  char buffer[80];
  label(buffer);
  out << buffer << std::endl;

  out.fill(' ');
  out << std::fixed << std::setprecision(6);

  // multiple orbitals; depending on the number of orbitals, the gaussian cube format is different
  out << std::setw( 4) << ( (nd != 1 ? -1 : 1) * (int)mol->natom())
      << std::setw(12) << grid_->origin[0]*to_atomic
      << std::setw(12) << grid_->origin[1]*to_atomic
      << std::setw(12) << grid_->origin[2]*to_atomic << std::endl;
  out << std::setw( 4) << grid_->numx
      << std::setw(12) << grid_->axisx[0]*to_atomic
      << std::setw(12) << grid_->axisx[1]*to_atomic
      << std::setw(12) << grid_->axisx[2]*to_atomic << std::endl;
  out << std::setw( 4) << grid_->numy
      << std::setw(12) << grid_->axisy[0]*to_atomic
      << std::setw(12) << grid_->axisy[1]*to_atomic
      << std::setw(12) << grid_->axisy[2]*to_atomic << std::endl;
  out << std::setw( 4) << grid_->numz
      << std::setw(12) << grid_->axisz[0]*to_atomic
      << std::setw(12) << grid_->axisz[1]*to_atomic
      << std::setw(12) << grid_->axisz[2]*to_atomic;

  for (int atom=0; atom<mol->natom(); atom++)
  {
      out << std::endl
          << std::setw( 4) << mol->Z(atom)
          << std::setw(12) << 0.0  // this value is expected, but has no use.
          << std::setw(12) << mol->r(atom, 0)
          << std::setw(12) << mol->r(atom, 1)
          << std::setw(12) << mol->r(atom, 2);
  }
  if(nd != 1)
  {
      out << std::endl;
      out << std::setw(4) << nd; // give the number of dimensions
      for (int var = 0; var < nd; ++var)
      {
          out << std::setw(6) << dmap(var);
      }
  }

  out << std::scientific << std::uppercase << std::setprecision(5);

  std::vector<SCVector3> Points;
  for (int i=0; i<grid_->numx; i++)
  {
      SCVector3 pointx = grid_->origin + i * grid_->axisx;
      for (int j=0; j<grid_->numy; j++)
      {
          SCVector3 pointy = pointx + j * grid_->axisy;
          for (int k=0; k<grid_->numz; k++)
          {
              SCVector3 pointz = pointy + k * grid_->axisz;
              Points.push_back(pointz * to_atomic);
          }
      }
  }

  std::vector<double> Vals;
  this->calculate_values(Points, Vals);
  {
      const int zperiod = grid_->numz * this->ndim();
      int npts_on_line = 0;
      for (int pt = 0; pt < Vals.size(); ++pt)
      {
          if(pt%zperiod == 0 || npts_on_line == 6) {
            out << std::endl;
            npts_on_line = 0;
          }
          out << setw(13) << Vals[pt]; npts_on_line += 1;
      }
  }
}

#ifdef MPQC_NEW_FEATURES
boost::property_tree::ptree&
WriteVectorGrid::write_xml(
    boost::property_tree::ptree& parent,
    const XMLWriter& writer
)
{
  using boost::property_tree::ptree;
  ptree& my_tree = this->get_my_ptree(parent);
  //----------------------------------------//
  initialize();
  //----------------------------------------//
  const double to_atomic = grid_->unit->to_atomic_units();
  Ref<Units> angstrom = new Units("angstrom");
  const double to_angstrom = angstrom->from_atomic_units();
  Ref<Molecule> mol = get_molecule();
  const int nd = this->ndim();
  const int npoints = grid_->numx * grid_->numy * grid_->numz;
  //----------------------------------------//
  my_tree.put("nfunctions", nd);
  writer.insert_child(my_tree, grid_);
  //----------------------------------------//
  std::vector<SCVector3> Points;
  for (int i=0; i<grid_->numx; i++)
  {
      SCVector3 pointx = grid_->origin + i * grid_->axisx;
      for (int j=0; j<grid_->numy; j++)
      {
          SCVector3 pointy = pointx + j * grid_->axisy;
          for (int k=0; k<grid_->numz; k++)
          {
              SCVector3 pointz = pointy + k * grid_->axisz;
              Points.push_back(pointz * to_atomic);
          }
      }
  }
  //----------------------------------------//
  std::vector<double> Vals;
  this->calculate_values(Points, Vals);
  //----------------------------------------//
  for(int idim = 0; idim < nd; ++idim){
    ptree& ftree = my_tree.add_child("function", ptree());
    ftree.put("<xmlattr>.index", idim);
    double* fdata = allocate<double>(npoints);
    for(int ipoint = 0; ipoint < npoints; ++ipoint){
      fdata[ipoint] = Vals[ipoint*nd + idim];
    }
    writer.put_binary_data(
        ftree.add_child("data", ptree()),
        fdata,
        npoints
    );
  }
  //----------------------------------------//
  return my_tree;
}
#endif // MPQC_NEW_FEATURES



/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
