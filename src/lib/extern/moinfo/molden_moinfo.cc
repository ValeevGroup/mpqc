//
// moinfo.cc
//
// Copyright (C) 2011 Edward Valeev
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

#include <extern/moinfo/molden_moinfo.h>
#include <iostream>
#include <map>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/split.h>

using namespace sc;

MOLDEN_ExternReadMOInfo::MOLDEN_ExternReadMOInfo(const std::string & filename)
{
  std::ifstream in(filename.c_str());

  //////
  // parse molecular geometry
  //////
  Ref<Molecule> molecule = new Molecule;
  const size_t nline = 256;
  char skipline[nline];
  for (int i = 0; i < 2; ++i) // skip the first 2 lines
  {
    in.getline(skipline, nline);
  }

  int natoms;
  in >> natoms;

  { // the third line; make sure the unit is bohr
    std::string temp1, temp2;
    in >> temp1;
    in >> temp2;
    if(temp2 != "(AU)")
      throw SCException("Unexpected MOLDEN input");
  }

  for (int i = 0; i < (2+natoms); ++i) // skip lines for Mulliken charges
  {
    in.getline(skipline, nline);
  }

  // the format: H1 sequence# charge x xy z
  for (int i = 0; i < natoms; ++i)
  {
    std::string atom;
    unsigned int offset;
    double charge, x, y, z;
    in >> atom >> offset >> charge >> x >> y >> z;
    molecule->add_atom(charge, x ,y, z);
  }

  //molecule->set_point_group(new PointGroup("d2h"));
  molecule->print();

  //////
  // parse basis set
  //////
  {
    Ref<AssignedKeyVal> tmpkv = new AssignedKeyVal;
    // add molecule
    tmpkv->assign("molecule", molecule.pointer());
    const unsigned int num_atoms = molecule->natom();
    bool have_gencon = false;

    // make name for atomic basis sets
    for (unsigned int a = 0; a < num_atoms; ++a)
    {
      std::ostringstream oss1, oss2;
      oss1 << "basis:" << a;
      oss2 << "custom" << a;
      tmpkv->assign(oss1.str().c_str(), oss2.str().c_str());
    }

    // read in first line
    std::string buffer, buffer2;
    const size_t nline = 256;
    char linebuf[nline];
    in.getline(linebuf, nline);
    buffer = linebuf;

    unsigned int nao = 0, nmo = 0;
    std::vector<double> mo_ao_coef; // to receive the mo coefficients (then put in coef_)
    std::vector<int> shell_offset(1, 0);
    std::vector<char> shells; // ['S', 'P', ...]
    std::map<char, int> shell_nao;
    shell_nao['S'] = 1;
    shell_nao['P'] = 3;
    shell_nao['D'] = 5;
    shell_nao['F'] = 7;
    shell_nao['G'] = 9;
    std::map<char, int *> am_reorder;
    am_reorder['S'] = {1}; // the spherical harmonic functins in molcas are in different order than in mpqc
    am_reorder['P'] = {2,3,1};
    am_reorder['D'] = {5,3,1,2,4};
    am_reorder['F'] = {7,5,3,1,2,4,6};
    am_reorder['G'] = {9,7,5,3,1,2,4,6,8};
    std::vector<int> mpqcorder; // this tells how the AOs are reordered


    // read shells for each atom
    int current_shell_index = -1;
    int current_atom_index = -1;
    bool have_more_atoms;
    bool have_more_shells;
    bool have_more_prims;
    std::map<int, int> current_shell_index_on_atom;
    std::string shell_prefix;
    std::string am;
    int puream = 1;

    do
    { // atom loop
      do
      { // shell loop
        int current_primitive_index = -1;
        do
        { // primitive loop
          const size_t first_not_whitespace = buffer.find_first_not_of(" ");
          if(first_not_whitespace == std::string::npos) // empty line; read one more line
          {
            const size_t nline = 256;
            char linebuf[nline];
            in.getline(linebuf, nline);
            buffer = linebuf;
          }
          const size_t first_not_whitespace = buffer.find_first_not_of(" ");
          const char first_char = buffer[first_not_whitespace];
          if (first_char == '[') // done with the basis set specification
          {
            have_more_atoms = false;
            have_more_shells = false;
            have_more_prims = false;
          }
          else
          { // have new shell or primitive? parse further
            const size_t first_whitespace = buffer.find_first_of(
                " ", first_not_whitespace);
            if(first_whitespace == std::string::npos)
            {
              have_more_atoms = true;
              have_more_shells = false;
              have_more_prims = false;
//              const size_t nline = 256;
//              char linebuf[nline];
//              in.getline(linebuf, nline);
//              buffer2 = linebuf; // the second line tells the angular momentum
            }
            else
            {
              have_more_atoms = false;
              const char first_char = buffer[first_not_whitespace];
              if(isalpha(first_char))
              {
                have_more_shells = true;
                have_more_prims = false;
              }
              else
              {
                have_more_shells = false;
                have_more_prims = true;
                ++current_primitive_index;
              }
            }
          }

          std::istringstream iss(buffer);

          // have new atom?
          if (have_more_atoms)
          {
            iss >> current_atom_index;
            --current_atom_index;
          }

          // have new shell?
          if (have_more_shells)
          {
            int nprim;
            iss >> am >> prim;
            if(am.size() > 1)
              throw SCException("Expected: angular momentum spec string of lenghth > 1");
            char amchar = toupper(am[0]); // convert 's' to 'S' etc
            shells.push_back(amchar);
            nao += shell_nao[amchar];
            int * am_mpqc_order = am_reorder[amchar];
            for (int i = 0; i < shell_nao[amchar]; ++i) // maps to the reordered sequence
            {
              mpqcorder.push_back(shell_offset.back() + am_mpqc_order[i]-1); // '-1': couting from 1
            }
            shell_offset.push_back(nao);


            if (current_shell_index_on_atom.find(current_atom_index)
                == current_shell_index_on_atom.end())
              current_shell_index_on_atom[current_atom_index] = -1;
            current_shell_index =
                ++current_shell_index_on_atom[current_atom_index];
            std::ostringstream oss;
            oss << ":basis:" << molecule->atom_name(current_atom_index)
                << ":custom" << current_atom_index << ":" << current_shell_index;
            shell_prefix = oss.str();

            tmpkv->assign((shell_prefix + ":type:0:am").c_str(),
                          std::string(1, amchar));
            tmpkv->assign((shell_prefix + ":type:0:puream").c_str(),
                            std::string("true"));
          }

          // have new primitive?
          if (have_more_prims)
          {
            double exponent, coef0, coef1;
            iss >> exponent >> coef0;
            { // exponent
              std::ostringstream oss;
              oss << shell_prefix << ":exp:" << current_primitive_index;
              tmpkv->assign(oss.str().c_str(), exponent);
            }
            { // coef 0
              std::ostringstream oss;
              oss << shell_prefix << ":coef:0:" << current_primitive_index;
              tmpkv->assign(oss.str().c_str(), coef0);
            }
          } // done with the current primitive

          // read next line
          if (have_more_prims || have_more_shells)
          {
            const size_t nline = 256;
            char linebuf[nline];
            in.getline(linebuf, nline);
            buffer = linebuf;
          }
        } while (have_more_prims); // loop over prims in this shell

      } while (have_more_shells); // loop over shells on this atom
    } while(have_more_atoms);
    Ref<KeyVal> kv = tmpkv;
    Ref<GaussianBasisSet> basis = new GaussianBasisSet(kv);

    basis_ = basis;
  }
  basis_->print();
  int nshells = shells.size();

  orbsym_.resize(nmo); // assume C1 symmetry for now.
  for(int o=0; o<nmo; ++o)
    orbsym_[o] = 0;

  bool done_reading = false
  while((not in.eof()) and (not done_reading))
  {
    { // not end of line yet;
       const size_t nline = 256;
       char linebuf[nline];
       in.getline(linebuf, nline);
       buffer = linebuf;
    }
    const size_t fnw = buffer.find_first_not_of(" ");
    if(fnw == std::string::npos) // finish reading MO coefficients
      done_reading = true;
    else
    {
      for (int j = 0; j < 4; ++j) // skip 3 lines: "Sym = ..."; "Ene = ..."; "Spin =="
      {
        const size_t nline = 256;
        char linebuf[nline];
        in.getline(linebuf, nline);
        buffer = linebuf;
      }
      for (int jj = 0; jj < nao; ++jj) // read AO coefficients for one MO
      {
        int xx, yy;
        in >> xx >> yy;
        mo_ao_coef.push_back(yy)
      }
      nmo++;
    }
  }


  RefSCDimension aodim = new SCDimension(nao, 1);
  aodim->blocks()->set_subdim(0, new SCDimension(nao));
  RefSCDimension modim = new SCDimension(nmo, 1);
  modim->blocks()->set_subdim(0, new SCDimension(nmo));
  coefs_ = basis_->so_matrixkit()->matrix(aodim, modim);
  coefs_.assign(0.0);
  for (int i = 0; i < nmo; ++i)
  {
    for (int j = 0; j < nao; ++j)
    {
      int ij = i*nao +j;
      coefs_.set_element(mpqcorder[j], i, mo_ao_coef[ij]);
    }
  }

  coefs_.print("MOLDEN_ExternReadMOInfo:: MO coefficients");
  in.close();
}


Ref<GaussianBasisSet> MOLDEN_ExternReadMOInfo::basis() const
{
    return basis_;
}

RefSCMatrix MOLDEN_ExternReadMOInfo::coefs() const
{
    return coefs_;
}

std::vector<unsigned int> MOLDEN_ExternReadMOInfo::orbsym() const
{
    return orbsym_;
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
