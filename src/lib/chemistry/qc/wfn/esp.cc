//
// esp.cc
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
#include <chemistry/qc/wfn/esp.h>
#include <util/misc/scexception.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////////
// WriteElectrostaticPotential

static ClassDesc WriteElectrostaticPotential_cd(
    typeid(WriteElectrostaticPotential),"WriteElectrostaticPotential",1,
    "public WriteGrid", 0, create<WriteElectrostaticPotential>, 0);

WriteElectrostaticPotential::WriteElectrostaticPotential(const Ref<KeyVal> &keyval):
  WriteGrid(keyval)
{
  wfn_ << keyval->describedclassvalue("wfn");
  if (wfn_ == 0) {
      InputError ex("valid \"wfn\" missing",
                    __FILE__, __LINE__, "wfn", "(null)", class_desc());
      try {
          ex.elaborate()
              << "WriteElectrostaticPotential KeyVal ctor requires"
              << " that \"wfn\" specifies an object"
              << " of type Wavefunction" << std::endl;
        }
      catch (...) {}
      throw ex;
    }

  KeyValValueboolean default_electronic(1);
  electronic_ = keyval->booleanvalue("electronic", default_electronic);

  KeyValValueboolean default_nuclear(1);
  nuclear_ = keyval->booleanvalue("nuclear", default_nuclear);
}

void
WriteElectrostaticPotential::initialize()
{
	pc_mat_ = wfn_->ao_density()->copy();

  ao_density_ = wfn_->ao_density()->copy();
  ao_density_->scale(2.0);          // SCElementScalarProduct computes with unique matrix elements only -> in
  ao_density_->scale_diagonal(0.5); // symmetric matrices the off-diagonal elements are doubled to compensate this.
}

void
WriteElectrostaticPotential::label(char* buffer)
{
  sprintf(buffer, "WriteElectrostaticPotential");
}

Ref<Molecule>
WriteElectrostaticPotential::get_molecule()
{
  return wfn_->molecule();
}

double
WriteElectrostaticPotential::calculate_value(SCVector3 point)
{
  double result = 0.0;
  double charge = 1.0;
  
  const double* positions = point.data();
	Ref<PointChargeData> pcdata = new PointChargeData(1, &positions, &charge);
	
	Ref<Molecule> molecule = wfn_->molecule();
  Ref<Integral> integral = wfn_->integral();
  
  if (nuclear_) {
      Ref<GaussianBasisSet> atom_basis = wfn_->atom_basis();
    	if (atom_basis != NULL) {
        	integral->set_basis(atom_basis);
          Ref<OneBodyOneCenterInt> ob_atom = integral->point_charge1(pcdata);
          
          const double *atom_buffer = ob_atom->buffer();
          const double *atom_basis_coef = wfn_->atom_basis_coef();
          
          for (int i=0,icoef=0; i<atom_basis->ncenter(); i++) {
              if (atom_basis->nshell_on_center(i) > 0) {
                  int joff = atom_basis->shell_on_center(i,0);
                  for (int j=0; j<atom_basis->nshell_on_center(i); j++) {
                      int jsh = j + joff;
                      ob_atom->compute_shell(jsh);
                      int nfunc = atom_basis->shell(jsh).nfunction();
                      for (int k=0; k<nfunc; k++,icoef++) {
                          result -= atom_basis_coef[icoef] * atom_buffer[k];
                        }
                    }
                }
              else {
                  SCVector3 a(molecule->r(i));
                  result += molecule->charge(i) / a.dist(point);    
                }
            }
        }
      else {
          for (int i=0; i<molecule->natom(); i++) {
              SCVector3 a(molecule->r(i));
              result += molecule->charge(i) / a.dist(point);    
            }
        }
    }
    
  if (electronic_) {
      integral->set_basis(wfn_->basis());

      Ref<SCElementOp> pc_op = new OneBodyIntOp(integral->point_charge(pcdata));
      pc_mat_->assign(0.0);
      pc_mat_->element_op(pc_op);

      Ref<SCElementScalarProduct> sp = new SCElementScalarProduct;
      sp->init();
      pc_mat_->element_op(sp, ao_density_);

      result += sp->result();
    }
    
  return result;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
