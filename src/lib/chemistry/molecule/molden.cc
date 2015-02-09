//
// molden.cc
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


#include <chemistry/molecule/molden.h>
#include <chemistry/molecule/atominfo.h>
#include <chemistry/qc/basis/petite.h>

using namespace std;
using namespace sc;


static ClassDesc WriteMolden_cd(
    typeid(WriteMolden),"WriteMolden",1,
    "public Runnable", 0, create<WriteMolden>, 0);

WriteMolden::WriteMolden(const Ref<KeyVal> &keyval){
  // Get the filename
  if (keyval->exists("filename")) {
    filename_ = keyval->stringvalue("filename");
  }
  else {
    filename_ = "-";
  }

  //----------------------------------------//
  // Get the OneBodyWavefunction

  obwfn_ << keyval->describedclassvalue("obwfn");
  if (obwfn_.null()) {
      InputError ex("valid \"obwfn\" missing",
                    __FILE__, __LINE__, "obwfn", "(null)", class_desc());
      try {
          ex.elaborate()
              << "WriteMolden KeyVal constructor requires"
              << " that \"obwfn\" specifies an object"
              << " of type OneBodyWavefunction" << std::endl;
        }
      catch (...) {}
      throw ex;
    }

}

void
WriteMolden::initialize(){
  // Get the mospace_ object from the OneBodyWavefunction
  Ref<PetiteList> pl = obwfn_->integral()->petite_list();
  RefSCMatrix aocoefs_full = pl->evecs_to_AO_basis(obwfn_->so_to_mo().t());
  mospace_ = new OrbitalSpace("p", "energy-ordered MOs to evaluate",
                                               aocoefs_full, obwfn_->basis(), obwfn_->integral(),
                                               obwfn_->eigenvalues(),
                                               0, 0,
                                               OrbitalSpace::energy);
}

void
WriteMolden::run(){

  initialize();

  //----------------------------------------//

  std::ostream *out;
  if (filename_ == "-") {
      out = &(ExEnv::out0());
  }
  else {
      ExEnv::out0() << incindent << indent << "WriteMolden"
                    << " is writing its output to \"" << filename_ << "\""
                    << std::endl;
      ExEnv::out0() << decindent;
      out = new std::ofstream(filename_.c_str());
  }

  //----------------------------------------//

  // Always write the [Molden Format] tag

  *out << "[Molden Format]" << endl;

  //----------------------------------------//

  // Write various sections

  write_atoms_section(*out);
  write_gto_section(*out);
  write_mo_section(*out);

  //----------------------------------------//

  if (filename_ == "-") {
      *out << decindent;
    }
  else {
      delete out;
    }

}

static inline int eq(const char* a, const char* b)
{
  return !strcmp(a,b);
}

void
WriteMolden::write_atoms_section(std::ostream &out)
{
  Ref<Molecule> mol = molecule();
  out << "[Atoms]";
  //if(eq(mol->geometry_units()->string_rep(), "angstrom")){
  //  out << " Angs";
  //}
  //else{
  //  out << " AU";
  //}
  out << " AU";
  out << endl;

  out.fill(' ');
  out << std::fixed << std::setprecision(8);
  for(int iatom = 0; iatom < mol->natom(); ++iatom){
    out << std::setw(4) << mol->atom_symbol(iatom)
        << std::setw(6) << (iatom + 1)
        << std::setw(6) << mol->Z(iatom)
        << setw(14) << mol->r(iatom, 0)
        << setw(14) << mol->r(iatom, 1)
        << setw(14) << mol->r(iatom, 2) << endl;
  }

}

void
WriteMolden::write_gto_section(std::ostream &out)
{
  Ref<GaussianBasisSet> basis = obwfn_->basis();
#if not USE_OLD_SOLIDHARM_ORDERING
      bmap_.resize(basis->nbasis());
#endif
  assert(basis->has_pure());

  out << "[GTO]" << endl;
  for(int center = 0; center < basis->ncenter(); ++center){
    out << setw(4) << (center + 1)
        << setw(4) << 0
        << endl;
    for(int ish = 0; ish < basis->nshell_on_center(center); ++ish){
      GaussianShell sh = basis->shell(center, ish);
      if(sh.ncontraction() == 1){
        out << setw(4) << sh.amchar(0)
            << setw(6) << sh.nprimitive()
            << setw(6) << "1.00"
            << endl;
        for(int iprim = 0; iprim < sh.nprimitive(); ++iprim){
          out << std::scientific << std::uppercase << std::setprecision(10);
          out << setw(20) << sh.exponent(iprim)
              << setw(20) << sh.coefficient_norm(0, iprim)
              << endl;
        } // end loop over primatives
#if not USE_OLD_SOLIDHARM_ORDERING
        // build the map
        const int ishoff = basis->shell_to_function(basis->shell_on_center(center, ish));
        const int l = sh.am(0);
        for(int m = -l; m <= l; ++m){
          bmap_[ishoff + ipure_molden(l, m)] = ishoff + ipure(l, m);
        }
#endif
      } // end if ncontraction == 1
      else{
        assert(false); // Not implemented
      }
    } //end loop over shells on center
    // format calls for a blank line here
    out << endl;
  } // end loop over centers

}

void
WriteMolden::write_mo_section(std::ostream &out)
{
  // TODO Open shell cases
  RefSCMatrix aocoeffs = mospace_->coefs_nb();
  RefDiagSCMatrix evals = mospace_->evals();
  out << "[MO]" << endl;
  const int nmo = aocoeffs.ncol();
  const int nao = obwfn_->basis()->nbasis();
  out << fixed << setprecision(10);
  for(int imo = 0; imo < nmo; ++imo){
    out << setw(8) << "Sym="
        << molecule()->point_group()->char_table().gamma(mospace_->orbsym()[imo]).symbol()
        << endl
        << setw(8) << "Ene="
        << setw(16) << evals->get_element(imo)
        << endl
        << setw(8) << "Spin="
        << "Alpha"
        << endl
        << setw(8) << "Occup="
        << obwfn_->occupation(imo)
        << endl;
    for(int iao = 0; iao < nao; ++iao){
      out << setw(6) << (iao+1)
#if not USE_OLD_SOLIDHARM_ORDERING
          << setw(16) << aocoeffs->get_element(iao, imo)
          //<< setw(16) << aocoeffs->get_element(bmap_[iao], imo)
#else
          << setw(16) << aocoeffs->get_element(iao, imo)
#endif
          << endl;
    } // end loop over AOs
  } // end loop over MOs
}
