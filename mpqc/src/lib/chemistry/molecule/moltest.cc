//
// moltest.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
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

#include <math.h>
#include <strstream.h>

#include <util/state/stateio.h>
#include <math/scmat/local.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/hess.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/molecule/coor.h>
#include <util/state/state_bin.h>
#include <util/render/object.h>
#include <util/render/oogl.h>
#include <util/misc/formio.h>
#include <chemistry/molecule/formula.h>

// force linkage 
#include <chemistry/molecule/linkage.h>

__builtin_delete(void*ptr)
{
  if (ptr>(void*)0 && ptr<(void*)0x100) abort();
  if (ptr) free(ptr);
}

void do_displacement(RefMolecularCoor&mc,int i);

void
print_atominfo(const RefAtomInfo &atominfo,
               const RefAtomInfo &refatominfo)
{
  cout << "Rvdw(H) = " << refatominfo->vdw_radius(1)
       << " " << atominfo->vdw_radius(1)
       << endl;
  cout << "Rvdw(C) = " << refatominfo->vdw_radius(6)
       << " " << atominfo->vdw_radius(6)
       << endl;
  cout << "Rb(H) = " << refatominfo->bragg_radius(1)
       << " " << atominfo->bragg_radius(1)
       << endl;
  cout << "Ra(H) = " << refatominfo->atomic_radius(1)
       << " " << atominfo->atomic_radius(1)
       << endl;
  cout << "mass(H) = " << refatominfo->mass(1)
       << " " << atominfo->mass(1)
       << endl;
  cout << "rgb(H) = "
       << "["
       << refatominfo->rgb(1,0) << " "
       << refatominfo->rgb(1,1) << " "
       << refatominfo->rgb(1,2) << " "
       << "] ["
       << atominfo->rgb(1,0) << " "
       << atominfo->rgb(1,1) << " "
       << atominfo->rgb(1,2) << " "
       << "]"
       << endl;
}

int 
main(int argc, char **argv)
{
  int i;

  // get the message group.  first try the commandline and environment
  RefMessageGrp grp = MessageGrp::initial_messagegrp(argc, argv);
  if (grp.nonnull())
    MessageGrp::set_default_messagegrp(grp);
  else
    grp = MessageGrp::get_default_messagegrp();

  RefKeyVal kv;
  if (argc == 2) {
      kv = new ParsedKeyVal(argv[1]);
    }
  else {
      kv = new ParsedKeyVal(SRCDIR "/moltest.in");
    }

  RefAtomInfo atominfo = kv->describedclassvalue("atominfo");
  if (atominfo.nonnull()) {
      RefAtomInfo refatominfo = new AtomInfo;
      cout << node0<< "-------------- testing atominfo --------------" << endl;
      if (grp->me() == 0) print_atominfo(atominfo, refatominfo);
      cout << node0 << "saving/restoring atominfo" << endl;
      ostrstream ostrs;
      StateOutBin so(ostrs);
      atominfo.save_state(so);
      atominfo = 0;
      so.close();
      ostrs.flush();
      istream istrs(ostrs.rdbuf());
      StateInBin si(istrs);
      atominfo.restore_state(si);
      if (grp->me() == 0) print_atominfo(atominfo, refatominfo);
      if (grp->n() > 1) {
          BcastState b(grp, 0);
          b.bcast(atominfo);
        }
      if (grp->me() == 1) {
          print_atominfo(atominfo, refatominfo);
        }
    }

  RefMolecule mol = kv->describedclassvalue("molecule");
  if (mol.nonnull()) {
      cout << "-------------- testing molecule --------------" << endl;

      MolecularFormula formula(mol);
      cout << "Molecular Formula" << endl << formula.formula() << endl;
      cout << "Number of Atomtypes" << endl << formula.natomtypes() << endl;
      cout << "Atomtype, Number of Atoms of This Type" << endl;
      for(i=0; i<formula.natomtypes(); i++) {
        cout << formula.Z(i) << "," << formula.nZ(i) << endl;
        }

      mol->cleanup_molecule();
      cout << "Clean Molecule:\n";
      mol->print();

      mol->transform_to_principal_axes();
      cout << "Clean Molecule wrt principal axes:\n";
      mol->print();

      int nunique = mol->nunique();

      cout << "nunique=" << nunique << ":";
      for (i=0; i < nunique; i++) cout << " " << mol->unique(i)+1;
      cout << endl;

      mol->point_group()->char_table().print();

      cout << "---------- testing molecule save/restore ----------" << endl;

      ostrstream ostrs;
      StateOutBin so(ostrs);
      cout << "saveing ..." << endl;
      mol.save_state(so);
      mol = 0;
      so.close();
      ostrs.flush();
      istream istrs(ostrs.rdbuf());
      StateInBin si(istrs);
      cout << "restoring ..." << endl;
      mol.restore_state(si);
      cout << "printing restored molecule:" << endl;
      mol->print();
    }

  cout << "-------------- initializing render tests --------------" << endl;
  RefRender ren = kv->describedclassvalue("renderer");
  RefRenderedObject renmol = kv->describedclassvalue("renderedmolecule");
  if (ren.nonnull() && renmol.nonnull()) {
      cout << "-------------- testing renderer --------------" << endl;
      ren->render(renmol);
    }

  //exit(0);

  RefSetIntCoor simp = kv->describedclassvalue("simp");
  if (simp.nonnull()) {
      cout << "-------------- testing simp  --------------" << endl;
      RefIntCoorGen gen = kv->describedclassvalue("generator");
      if (gen.nonnull()) {
          gen->print();
        }
      cout << "simp before update:\n";
      simp->print_details(mol);
      simp->update_values(mol);
      cout << "simp:\n";
      simp->print_details(mol);
    }

  // compare the analytic bmatrix to the finite displacement bmatrix
  RefSetIntCoor bmat_test = kv->describedclassvalue("bmat_test");
  if (bmat_test.nonnull()) {
      cout << "-------------- bmat_test  --------------" << endl;
      RefSCMatrixKit kit = SCMatrixKit::default_matrixkit();
      RefSCDimension dnc(new SCDimension(bmat_test->n()));
      RefSCDimension dn3(new SCDimension(mol->natom()*3));
      RefSCMatrix bmatrix(dnc,dn3,kit);
      RefSCMatrix fd_bmatrix(dnc,dn3,kit);
      cout << "testing bmat with:\n";
      bmat_test->update_values(mol);
      bmat_test->print();
      bmat_test->bmat(mol,bmatrix);
      bmat_test->fd_bmat(mol,fd_bmatrix);
      cout << "test bmatrix:\n";
      bmatrix.print();
      cout << "fd bmatrix:\n";
      fd_bmatrix.print();
      RefSCMatrix diff = fd_bmatrix - bmatrix;
      cout << "difference between test and finite displacement bmatrix:\n";
      diff.print();
      cout << "% difference between test and finite displacement bmatrix:\n";
      for (i=0; i<diff.nrow(); i++) {
          for (int j=0; j<diff.ncol(); j++) {
              double denom = fabs(fd_bmatrix(i,j));
              double num = fabs(diff(i,j));
              if (denom < 0.000001) denom = 0.000001;
              if (num < 0.00001) diff(i,j) = 0.0;
              else diff(i,j) = 100.0 * fabs(diff(i,j))/denom;
            }
        }
      diff.print();

      cout << "testing for translational invariance of each coordinate:"
           << endl;
      for (i=0; i<bmat_test->n(); i++) {
          cout << "  coor " << scprintf("%2d",i) << ":";
          for (int j=0; j<3; j++) {
              double sum = 0.0;
              for (int k=0; k<mol->natom(); k++) {
                  sum += bmatrix(i,k*3+j);
                }
              cout << scprintf(" % 16.12f",sum);
            }
          cout << endl;
        }
      bmatrix.gi().print("The inverse bmatrix");
    }

  cout.flush();
  cerr.flush();
  
  // now we get ambitious
  RefMolecularCoor mc = kv->describedclassvalue("molcoor");
  cout.flush();
  cerr.flush();

  if (mc.nonnull()) {
      cout << "-------------- testing molcoor  --------------" << endl;
      mc->print();

      cout.flush();
      cerr.flush();

      // do_displacement(mc,0);
      // do_displacement(mc,1);
      // do_displacement(mc,2);
      // do_displacement(mc,3);

      RefSCMatrixKit kit = SCMatrixKit::default_matrixkit();
      RefSymmSCMatrix hessian(mc->dim(),kit);
      mc->guess_hessian(hessian);

      // cout << "The guess hessian:\n";
      // hessian.print();
    }

  RefMolecularEnergy me = kv->describedclassvalue("energy");
  if (me.nonnull()) {
      cout << "-------------- testing energy  --------------" << endl;
      me->print();
    }

  RefMolecularHessian molhess = kv->describedclassvalue("hess");
  RefSymmSCMatrix xhessian;
  if (molhess.nonnull()) {
      xhessian = molhess->cartesian_hessian();
    }

  RefMolecularFrequencies molfreq = kv->describedclassvalue("freq");
  if (molfreq.nonnull() && xhessian.nonnull()) {
      cout << "-------------- testing freq  --------------" << endl;
      molfreq->compute_frequencies(xhessian);
    }

  return 0;
}


void
do_displacement(RefMolecularCoor&mc,int i)
{
  if (i>=mc->dim().n()) return;
  // now try to displace the geometry
  RefSCVector internal(mc->dim(),mc->matrixkit());
  mc->to_internal(internal);
  cout << "The initial internal coordinates:\n";
  internal.print();
  internal(i) = internal(i) + 0.2;
  cout << "The new internal coordinates:\n";
  internal.print();
  mc->to_cartesian(internal);
  mc->to_internal(internal);
  cout << "The actual new internal coordinates:\n";
  internal.print();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
