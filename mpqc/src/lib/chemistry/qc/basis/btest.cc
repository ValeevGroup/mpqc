//
// btest.cc --- test program
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

#include <util/keyval/keyval.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/symmint.h>
#include <chemistry/qc/intv3/intv3.h>

static void
test_overlap(const RefGaussianBasisSet& gbs, const RefGaussianBasisSet& gbs2,
             const RefIntegral& intgrl)
{
  intgrl->set_basis(gbs);

  // first form AO basis overlap
  RefSymmSCMatrix s(gbs->basisdim(), gbs->matrixkit());
  RefSCElementOp ov = new OneBodyIntOp(new OneBodyIntIter(intgrl->overlap()));
  s.assign(0.0);
  s.element_op(ov);
  ov=0;
  s.print("overlap");
      
  // now transform s to SO basis
  RefPetiteList pl = intgrl->petite_list();
  RefSymmSCMatrix sb = pl->to_SO_basis(s);
  sb.print("blocked s");
      
  // and back to AO basis
  s = pl->to_AO_basis(sb);
  s.print("reconstituted s");

  // form skeleton overlap
  ov = new OneBodyIntOp(new SymmOneBodyIntIter(intgrl->overlap(),pl));
  s.assign(0.0);
  s.element_op(ov);
  ov=0;
  s.print("overlap");

  // and symmetrize to get blocked overlap again
  sb.assign(0.0);
  pl->symmetrize(s,sb);
  sb.print("blocked again");

  s=0; sb=0;

  // now try overlap between two basis sets
  RefSCMatrix ssq(gbs2->basisdim(),gbs->basisdim(),gbs2->matrixkit());
  intgrl->set_basis(gbs2,gbs);

  ov = new OneBodyIntOp(new OneBodyIntIter(intgrl->overlap()));
  ssq.assign(0.0);
  ssq.element_op(ov);
  ssq.print("overlap sq");
  ov=0;

  RefPetiteList pl2 = intgrl->petite_list(gbs2);
  RefSCMatrix ssqb(pl2->AO_basisdim(), pl->AO_basisdim(), gbs->so_matrixkit());
  ssqb->convert(ssq);

  RefSCMatrix syms2 = pl2->aotoso().t() * ssqb * pl->aotoso();
  syms2.print("symm S2");
}

static void
test_eigvals(const RefGaussianBasisSet& gbs, const RefIntegral& intgrl)
{
  intgrl->set_basis(gbs);
  RefPetiteList pl = intgrl->petite_list();

  // form AO Hcore and evecs
  RefSymmSCMatrix hcore_ao(gbs->basisdim(), gbs->matrixkit());
  RefSCMatrix ao_evecs(gbs->basisdim(), gbs->basisdim(), gbs->matrixkit());
  RefDiagSCMatrix ao_evals(gbs->basisdim(), gbs->matrixkit());
  
  hcore_ao.assign(0.0);

  RefSCElementOp op = new OneBodyIntOp(new OneBodyIntIter(intgrl->kinetic()));
  hcore_ao.element_op(op);
  op=0;

  RefOneBodyInt nuc = intgrl->nuclear();
  nuc->reinitialize();
  op = new OneBodyIntOp(nuc);
  hcore_ao.element_op(op);
  op=0;
  
  hcore_ao.print("Hcore (AO)");
  
  hcore_ao.diagonalize(ao_evals, ao_evecs);
  ao_evecs.print("AO Evecs");
  ao_evals.print("AO Evals");

  // form SO Hcore and evecs
  RefSymmSCMatrix hcore_so(pl->SO_basisdim(), gbs->so_matrixkit());
  RefSCMatrix so_evecs(pl->SO_basisdim(), pl->SO_basisdim(),
                       gbs->so_matrixkit());
  RefDiagSCMatrix so_evals(pl->SO_basisdim(), gbs->so_matrixkit());
  
  // reuse hcore_ao to get skeleton Hcore
  hcore_ao.assign(0.0);

  op = new OneBodyIntOp(new SymmOneBodyIntIter(intgrl->kinetic(),pl));
  hcore_ao.element_op(op);
  op=0;

  nuc = intgrl->nuclear();
  nuc->reinitialize();
  op = new OneBodyIntOp(new SymmOneBodyIntIter(nuc,pl));
  hcore_ao.element_op(op);
  op=0;
  
  pl->symmetrize(hcore_ao, hcore_so);

  hcore_so.print("Hcore (SO)");
  
  hcore_so.diagonalize(so_evals, so_evecs);
  so_evecs.print("SO Evecs");
  so_evals.print("SO Evals");

  RefSCMatrix new_ao_evecs = pl->evecs_to_AO_basis(so_evecs);
  new_ao_evecs.print("AO Evecs again");

  //RefSCMatrix new_so_evecs = pl->evecs_to_SO_basis(ao_evecs);
  //new_so_evecs.print("SO Evecs again");

  pl->to_AO_basis(hcore_so).print("Hcore (AO) again");
}

int
main(int, char *argv[])
{
  char *filename = (argv[1]) ? argv[1] : SRCDIR "/btest.kv";
  
  RefKeyVal keyval = new ParsedKeyVal(filename);
  
  RefIntegral intgrl = new IntegralV3;

  for (int i=0; i<keyval->count("test"); i++) {
      RefGaussianBasisSet gbs = keyval->describedclassvalue("test", i);
      RefGaussianBasisSet gbs2 = keyval->describedclassvalue("test2", i);

      test_overlap(gbs,gbs2,intgrl);

      test_eigvals(gbs,intgrl);

      StateOutText out("btest.out");
      gbs.save_state(out);
      StateInText in("btest.out");
      gbs.restore_state(in);
      gbs->print();
      intgrl->petite_list()->print();
    }

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
