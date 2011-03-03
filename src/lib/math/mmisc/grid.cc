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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdexcept>

#include <math/scmat/vector3.h>
#include <math/scmat/matrix.h>
#include <math/mmisc/grid.h>
#include <util/class/scexception.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////////
// GridDefinition

static ClassDesc Grid_cd(
    typeid(Grid),"Grid",1,
    "public DescribedClass", 0, create<Grid>, 0);

Grid::Grid(const Ref<KeyVal> &keyval):
  origin(0.,0.,0.),
  axisx(0.,0.,0.),
  axisy(0.,0.,0.),
  axisz(0.,0.,0.)
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
  char buffer[256];
  label(buffer);
  out << buffer << std::endl;
  out << std::endl;
  
  out.fill(' ');
  out << std::fixed << std::setprecision(6);
  
  out << std::setw( 4) << mol->natom() 
      << std::setw(12) << grid_->origin[0]*to_atomic*to_angstrom
      << std::setw(12) << grid_->origin[1]*to_atomic*to_angstrom
      << std::setw(12) << grid_->origin[2]*to_atomic*to_angstrom << std::endl;
  out << std::setw( 4) << grid_->numx
      << std::setw(12) << grid_->axisx[0]*to_atomic*to_angstrom
      << std::setw(12) << grid_->axisx[1]*to_atomic*to_angstrom
      << std::setw(12) << grid_->axisx[2]*to_atomic*to_angstrom << std::endl;
  out << std::setw( 4) << grid_->numy
      << std::setw(12) << grid_->axisy[0]*to_atomic*to_angstrom
      << std::setw(12) << grid_->axisy[1]*to_atomic*to_angstrom
      << std::setw(12) << grid_->axisy[2]*to_atomic*to_angstrom << std::endl;
  out << std::setw( 4) << grid_->numz
      << std::setw(12) << grid_->axisz[0]*to_atomic*to_angstrom
      << std::setw(12) << grid_->axisz[1]*to_atomic*to_angstrom
      << std::setw(12) << grid_->axisz[2]*to_atomic*to_angstrom;
      
  for (int atom=0; atom<mol->natom(); atom++) {
      out << std::endl
          << std::setw( 4) << mol->Z(atom)
          << std::setw(12) << 0.0  // this value is expected, but has no use.
          << std::setw(12) << mol->r(atom, 0)*to_angstrom
          << std::setw(12) << mol->r(atom, 1)*to_angstrom
          << std::setw(12) << mol->r(atom, 2)*to_angstrom;
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
// WriteGrids

static ClassDesc WriteGrids_cd(
    typeid(WriteGrids),"WriteGrids",1,
    "public Runnable");

WriteGrids::WriteGrids(const Ref<KeyVal> &keyval)
{
  grid_ << keyval->describedclassvalue("grid");
  if (grid_.null()) {
      InputError ex("valid \"grid\" missing",
                    __FILE__, __LINE__, "grid", "(null)", class_desc());
      try {
          ex.elaborate()
              << "WriteGrids KeyVal ctor requires"
              << " that \"grid\" specifies an object"
              << " of type Grid" << std::endl;
        }
      catch (...) {}
      throw ex;
    }

  filename_ = keyval->stringvalue("filename", KeyValValuestring("Grid"));
  /** we start numbering orbitals from '1', instead of 0 */
  first_ = keyval->intvalue("first", KeyValValueint(1));
  last_ = keyval->intvalue("last", KeyValValueint(0));
  //whether these values are appropriete is checked in the contructor of the inherited class WriteOrbitals
  // since in this class we do not have information of the wavefunction and thus do not know the number of orbitals,
  // we use last_ == 0 to indicate that all orbitals will be written (by modifying the value in the
   // WriteOrbitals class constructor).

  if (keyval->exists("format")) {
      format_ = keyval->stringvalue("format");
      if (format_ == "gaussian_cube") {
          write_format_ = &WriteGrids::wf_gaussian_cube;
        }
      else {
          InputError ex("valid \"format\" missing",
                        __FILE__, __LINE__, "format", "(null)", class_desc());
          try {
              ex.elaborate()
                  << "WriteGrids KeyVal ctor requires"
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
      write_format_ = &WriteGrids::wf_gaussian_cube;
    }
}



void
WriteGrids::run()
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
 // this->wf_gaussian_cube(*out);

  if (filename_ == "-") {
      *out << decindent;
    }
  else {
      delete out;
    }
}





void
WriteGrids::wf_gaussian_cube(std::ostream &out) {
  double to_atomic = grid_->unit->to_atomic_units();
  Ref<Units> angstrom = new Units("angstrom");
  double to_angstrom = angstrom->from_atomic_units();
  Ref<Molecule> mol = get_molecule();

  out << "Gaussian cube file generated by MPQC." << std::endl;
  char buffer[256];
  label(buffer);
  out << buffer << std::endl;
  out << std::endl;

  out.fill(' ');
  out << std::fixed << std::setprecision(6);

  if(first_ != last_) // multiple orbitals; depending on the number of orbitals, the gaussian cube format is different
  {
      out << std::setw( 4) << (-1 * mol->natom())
          << std::setw(12) << grid_->origin[0]*to_atomic*to_angstrom
          << std::setw(12) << grid_->origin[1]*to_atomic*to_angstrom
          << std::setw(12) << grid_->origin[2]*to_atomic*to_angstrom << std::endl;
  }
  else // only one orbital
  {
      out << std::setw( 4) << mol->natom()
          << std::setw(12) << grid_->origin[0]*to_atomic*to_angstrom
          << std::setw(12) << grid_->origin[1]*to_atomic*to_angstrom
          << std::setw(12) << grid_->origin[2]*to_atomic*to_angstrom << std::endl;
  }
  out << std::setw( 4) << grid_->numx
      << std::setw(12) << grid_->axisx[0]*to_atomic*to_angstrom
      << std::setw(12) << grid_->axisx[1]*to_atomic*to_angstrom
      << std::setw(12) << grid_->axisx[2]*to_atomic*to_angstrom << std::endl;
  out << std::setw( 4) << grid_->numy
      << std::setw(12) << grid_->axisy[0]*to_atomic*to_angstrom
      << std::setw(12) << grid_->axisy[1]*to_atomic*to_angstrom
      << std::setw(12) << grid_->axisy[2]*to_atomic*to_angstrom << std::endl;
  out << std::setw( 4) << grid_->numz
      << std::setw(12) << grid_->axisz[0]*to_atomic*to_angstrom
      << std::setw(12) << grid_->axisz[1]*to_atomic*to_angstrom
      << std::setw(12) << grid_->axisz[2]*to_atomic*to_angstrom;

  for (int atom=0; atom<mol->natom(); atom++)
  {
      out << std::endl
          << std::setw( 4) << mol->Z(atom)
          << std::setw(12) << 0.0  // this value is expected, but has no use.
          << std::setw(12) << mol->r(atom, 0)*to_angstrom
          << std::setw(12) << mol->r(atom, 1)*to_angstrom
          << std::setw(12) << mol->r(atom, 2)*to_angstrom;
  }
  if(first_ != last_)
  {
      out << std::endl;
      out << std::setw(4) << last_ - first_ + 1; // give the number of orbitals
      for (int var = first_; var <= last_; ++var)
      {
          out << std::setw(6) << var;
      }
  }
  SCVector3 pointx;
  SCVector3 pointy;
  SCVector3 pointz;

  out << std::scientific << std::uppercase << std::setprecision(5);
  int numpoints = grid_->numx * grid_->numy * grid_->numz;
  int numorbs = last_ - first_ + 1;
  int totalnum = numorbs * numpoints;

  std::vector<int> Orbs(numorbs);
  std::vector<SCVector3> Points(numpoints);
  double * Vals = new double [totalnum];
  for (int i=0; i<grid_->numx; i++)
  {
      pointx = grid_->origin + i * grid_->axisx;
      for (int j=0; j<grid_->numy; j++)
      {
          pointy = pointx + j * grid_->axisy;
          for (int k=0; k<grid_->numz; k++)
          {
              static int countpoint = 0;
              pointz = pointy + k * grid_->axisz;
              Points[countpoint] = pointz;
              countpoint++;
          }
      }
  }

  for(int orbitalnum = first_; orbitalnum <= last_; orbitalnum++)
  {
      static int countt = 0;
      Orbs[countt] = orbitalnum-1; // for convience, we number first_ and last_ starting with '1';
      countt++;                    // but to be consistent with c++, we have decrease it by 1 here if
  }                                // we don't want to change in other functions

  this->calculate_values(Orbs, Points, Vals);
  {
      int val = 0;
      for (val = 0; val < totalnum; ++val)
      {
          if(val%6 == 0) out << std::endl;
          out << setw(13) << Vals[val];
      }
  }

  delete[] Vals;
}






/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
