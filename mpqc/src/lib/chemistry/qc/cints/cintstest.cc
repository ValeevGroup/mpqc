//
// cintstest.cc
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#include <stdlib.h>
#include <string.h>

#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <util/keyval/keyval.h>
#include <util/group/message.h>
#include <util/group/pregtime.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/intv3/int1e.h>
#include <chemistry/qc/intv3/int2e.h>
#include <chemistry/qc/intv3/intv3.h>
#include <chemistry/qc/intv3/cartitv3.h>
#include <chemistry/qc/cints/cints.h>
#include <chemistry/qc/cints/int2e.h>
#include <chemistry/qc/cints/cartit.h>


using namespace std;
using namespace sc;

#define CINTS

void compare_1e_cints_vs_v3(Ref<OneBodyInt>& obcints, Ref<OneBodyInt>& obv3);
void compare_2e_cints_vs_v3(Ref<TwoBodyInt>& tbcints, Ref<TwoBodyInt>& tbv3);
void compare_2e_puream_cints_vs_v3(Ref<TwoBodyInt>& tbcints, Ref<TwoBodyInt>& tbv3);
void compare_2e_bufsum_cints_vs_v3(Ref<TwoBodyInt>& tbcints, Ref<TwoBodyInt>& tbv3);
void compare_2e_unique_bufsum_cints_vs_v3(Ref<TwoBodyInt>& tbcints, Ref<TwoBodyInt>& tbv3);
void print_grt_ints(Ref<TwoBodyInt>& tbcints);
void compare_2e_permute(Ref<Integral>& cints);
void test_int_shell_1e(const Ref<KeyVal>&, const Ref<Int1eV3> &int1ev3,
                       void (Int1eV3::*int_shell_1e)(int,int),
                       int permute);
void test_3_center(const Ref<KeyVal>&, const Ref<Int2eV3> &);
void test_4_center(const Ref<KeyVal>& keyval, const Ref<Int2eV3> &int2ev3);
void test_4der_center(const Ref<KeyVal>&, const Ref<Int2eV3> &int2ev3);

#define maxint 9

void
testint(const Ref<OneBodyInt>& in)
{
  if (in.null()) {
      cout << "null integral generator" << endl;
      abort();
    }
}

void
testint(const Ref<OneBodyDerivInt>& in)
{
  if (in.null()) {
      cout << "null integral generator" << endl;
      abort();
    }
}

void
testint(const Ref<TwoBodyInt>& in)
{
  if (in.null()) {
      cout << "null integral generator" << endl;
      abort();
    }
}

void
testint(const Ref<TwoBodyDerivInt>& in)
{
  if (in.null()) {
      cout << "null integral generator" << endl;
      abort();
    }
}

/*void
do_bounds_stats(const Ref<KeyVal>& keyval,
                const Ref<Int2eV3> &int2ev3)
{
  int i,j;
  int nshell = int2ev3->basis()->nshell();
  int eps = -10;
  int *nonzero = new int[nshell];
  for (i=0; i<nshell; i++) {
      if (i==0) nonzero[i] = 0;
      else nonzero[i] = nonzero[i-1];
      for (j=0; j<=i; j++) {
          if (int2ev3->erep_4bound(i,j,-1,-1) > eps) {
              nonzero[i]++;
            }
        }
    }
  for (i=0; i<nshell; i++) {
      int natom = 1 + int2ev3->basis()->shell_to_center(i);
      int npq = (i*(i+1))/2;
      cout<<scprintf("nsh=%2d nat=%2d npq=%4d npq>eps=%4d npq>eps/nsh=%9.4f /nat=%9.4f",
                     i, natom, npq, nonzero[i], double(nonzero[i])/i,
                     double(nonzero[i])/natom)
          << endl;
    }
  delete[] nonzero;
}
*/

int main(int argc, char **argv)
{
  int ii, i,j,k,l,m,n;

  Ref<MessageGrp> msg = MessageGrp::initial_messagegrp(argc,argv);
  if (msg.null()) msg = new ProcMessageGrp();
  MessageGrp::set_default_messagegrp(msg);

  Ref<RegionTimer> tim = new ParallelRegionTimer(msg,"cintstest", 1, 1);

  char *infile = new char[strlen(SRCDIR)+strlen("/cintstest.in")+1];
  sprintf(infile,SRCDIR "/cintstest.in");
  if (argc == 2) {
    infile = argv[1];
    }

  Ref<KeyVal> pkv(new ParsedKeyVal(infile));
  Ref<KeyVal> tkeyval(new PrefixKeyVal(":test", pkv));

  Ref<GaussianBasisSet> basis = require_dynamic_cast<GaussianBasisSet*>(
    tkeyval->describedclassvalue("basis").pointer(),"main\n");
  Ref<Molecule> mol = basis->molecule();

  int tproc = tkeyval->intvalue("test_processor");
  if (tproc >= msg->n()) tproc = 0;
  int me = msg->me();

  if (me == tproc) cout << "testing on processor " << tproc << endl;

  int storage = tkeyval->intvalue("storage");
  cout << "storage = " << storage << endl;
  /*  Ref<Integral> intgrlv3 = new IntegralV3(basis,basis,basis,basis);
  Ref<Int1eV3> int1ev3 = new Int1eV3(intgrlv3.pointer(),basis,basis,1);
  Ref<Int2eV3> int2ev3 = new Int2eV3(intgrlv3.pointer(),basis,basis,basis,basis,
                                   1, storage);
  

  int permute = tkeyval->booleanvalue("permute");
  tim->enter("overlap");
  if (me == tproc && tkeyval->booleanvalue("overlap")) {
      cout << scprintf("testing overlap:\n");
      test_int_shell_1e(tkeyval, int1ev3, &Int1eV3::overlap, permute);
    }
  tim->change("kinetic");
  if (me == tproc && tkeyval->booleanvalue("kinetic")) {
      cout << scprintf("testing kinetic:\n");
      test_int_shell_1e(tkeyval, int1ev3, &Int1eV3::kinetic, permute);
    }
  tim->change("hcore");
  if (me == tproc && tkeyval->booleanvalue("hcore")) {
      cout << scprintf("testing hcore:\n");
      test_int_shell_1e(tkeyval, int1ev3, &Int1eV3::hcore, permute);
    }
  tim->change("nuclear");
  if (me == tproc && tkeyval->booleanvalue("nuclear")) {
      cout << scprintf("testing nuclear:\n");
      test_int_shell_1e(tkeyval, int1ev3, &Int1eV3::nuclear, permute);
    }
  tim->change("3 center");
  if (me == tproc && tkeyval->booleanvalue("3")) {
      test_3_center(tkeyval, int2ev3);
    }
  tim->change("4 center");
  if (me == tproc && tkeyval->booleanvalue("4")) {
      test_4_center(tkeyval, int2ev3);
    }
  tim->change("4 center der");
  if (me == tproc && tkeyval->booleanvalue("4der")) {
      test_4der_center(tkeyval, int2ev3);
    }
  tim->change("bound stats");
  if (me == tproc && tkeyval->booleanvalue("boundstats")) {
      do_bounds_stats(tkeyval, int2ev3);
    }

    tim->change("IntegralV3");*/

  tim->enter("Integral");
  Ref<Integral> integral = new IntegralV3(basis);
#ifdef CINTS
  Ref<Integral> integralcints = new IntegralCints(basis);
#endif

  Ref<OneBodyInt> overlapv3 = integral->overlap();
  Ref<OneBodyInt> kineticv3 = integral->kinetic();
  Ref<OneBodyInt> nuclearv3 = integral->nuclear();
  Ref<OneBodyInt> hcorev3 = integral->hcore();

#ifdef CINTS
  Ref<OneBodyInt> overlapcints = integralcints->overlap();
  testint(overlapcints);
  Ref<OneBodyInt> kineticcints = integralcints->kinetic();
  testint(kineticcints);
  Ref<OneBodyInt> nuclearcints = integralcints->nuclear();
  testint(nuclearcints);
  Ref<OneBodyInt> hcorecints = integralcints->hcore();
  testint(hcorecints);
#endif

  Ref<TwoBodyInt> erepv3 = integral->electron_repulsion();

#ifdef CINTS
  int storage_needed = integralcints->storage_required_eri(basis);
  cout << scprintf("Need %d bytes to create EriCints\n",storage_needed);
  Ref<TwoBodyInt> erepcints = integralcints->electron_repulsion();
  testint(erepcints);
  storage_needed = integralcints->storage_required_grt(basis);
  cout << scprintf("Need %d bytes to create GRTCints\n",storage_needed);
  Ref<TwoBodyInt> grtcints = integralcints->grt();
  testint(grtcints);
#endif
  tim->exit();

  // Test iterators
  /*  CartesianIterCints citer(3);
  cout << "Cartesian f-shell:" << endl;
  for ( citer.start(); int(citer) ; citer.next() )
    cout << "nx = " << citer.a() << " ny = " << citer.b() << " nz = " << citer.c() << endl;
  RedundantCartesianIterCints rciter(3);
  cout << "Redundant Cartesian f-shell:" << endl;
  for ( rciter.start(); int(rciter) ; rciter.next() )
    cout << "nx = " << rciter.a() << " ny = " << rciter.b() << " nz = " << rciter.c() << endl;
  */

  //cout << "Testing Cints' overlap integrals against IntV3's" << endl;
  //  compare_1e_cints_vs_v3(overlapcints,overlapv3);
  //cout << "Testing Cints' kinetic energy integrals against IntV3's" << endl;
  //compare_1e_cints_vs_v3(kineticcints,kineticv3);
  //cout << "Testing Cints' nuclear attraction integrals against IntV3's" << endl;
  //compare_1e_cints_vs_v3(nuclearcints,nuclearv3);
  //cout << "Testing Cints' core hamiltonian integrals against IntV3's" << endl;
  //compare_1e_cints_vs_v3(hcorecints,hcorev3);

  //  compare_2e_permute(integralcints);

  cout << "Testing Cints' ERIs against IntV3's" << endl;
  compare_2e_cints_vs_v3(erepcints,erepv3);
  //compare_2e_puream_cints_vs_v3(erepcints,erepv3);
  cout << "Testing Cints' ERIs (from GRTCints) against IntV3's" << endl;
  compare_2e_cints_vs_v3(grtcints,erepv3);

#ifdef CINTS
  cout << "Testing sums of Cints' ERIs against IntV3's" << endl;
  compare_2e_bufsum_cints_vs_v3(erepcints,erepv3);
  cout << "Testing sums of Cints' ERIs (from GRTCints) against IntV3's" << endl;
  compare_2e_bufsum_cints_vs_v3(grtcints,erepv3);

  cout << "Testing sums of unique Cints' ERIs against IntV3's" << endl;
  erepcints->set_redundant(0);
  erepv3->set_redundant(0);
  grtcints->set_redundant(0);
  compare_2e_unique_bufsum_cints_vs_v3(erepcints,erepv3);
  cout << "Testing sums of unique Cints' ERIs (from GRTCints) against IntV3's" << endl;
  compare_2e_unique_bufsum_cints_vs_v3(grtcints,erepv3);

  cout << "Printing GRT integrals" << endl;
  print_grt_ints(grtcints);
#endif

  //  tim->print();
  return 0;
}

void
compare_1e_cints_vs_v3(Ref<OneBodyInt>& obcints, Ref<OneBodyInt>& obv3)
{
  Ref<GaussianBasisSet> basis = obcints->basis();
  for (int sh1=4; sh1<basis->nshell(); sh1++)
    for (int sh2=0; sh2<basis->nshell(); sh2++) {
      int nbf2 = basis->shell(sh2).nfunction();
      obv3->compute_shell(sh1,sh2);
      obcints->compute_shell(sh1,sh2);
      const double *buffercints = obcints->buffer();
      const double *bufferv3 = obv3->buffer();

      int bf1_offset = 0;
      for (int gc1=0; gc1<basis->shell(sh1).ncontraction(); gc1++) {
	int am1 = basis->shell(sh1).am(gc1);
	CartesianIterCints citer1(am1);
	CartesianIterV3 iter1(am1);
	for ( citer1.start(); int(citer1) ; citer1.next() ) {
	  int bf1cints = bf1_offset + citer1.bfn();
	  int bf1v3;
	  for( iter1.start(); int(iter1) ; iter1.next() ) {
	    if (iter1.a() == citer1.a() &&
		iter1.b() == citer1.b() &&
		iter1.c() == citer1.c()) {
	      bf1v3 = bf1_offset + iter1.bfn();
	      break;
	    }
	  }

	  int bf2_offset = 0;
	  for (int gc2=0; gc2<basis->shell(sh2).ncontraction(); gc2++) {
	    int am2 = basis->shell(sh2).am(gc2);
	    CartesianIterCints citer2(am2);
	    CartesianIterV3 iter2(am2);
	    
	    for ( citer2.start(); int(citer2) ; citer2.next() ) {
	      int bf2cints = bf2_offset + citer2.bfn();
	      int bf2v3;
	      for( iter2.start(); int(iter2) ; iter2.next() ) {
		if (iter2.a() == citer2.a() &&
		    iter2.b() == citer2.b() &&
		    iter2.c() == citer2.c()) {
		  bf2v3 = bf2_offset + iter2.bfn();
		  break;
		}
	      }
	      
	      double valuecints = buffercints[bf1cints*nbf2 + bf2cints];
	      double valuev3 = bufferv3[bf1v3*nbf2 + bf2v3];
	      if (fabs(valuecints-valuev3) > 1E-13) {
		cout << scprintf("Discrepancy in OEInt(sh1 = %d, sh2 = %d)\n",sh1,sh2);
		cout << scprintf("bf1 = %d   bf2 = %d  OEIntegral(cints) = %20.15lf\n",bf1cints,bf2cints,valuecints);
		cout << scprintf("bf1 = %d   bf2 = %d  OEIntegral(V3)    = %20.15lf\n\n",bf1v3,bf2v3,valuev3);
	      }
	    }
	    bf2_offset += basis->shell(sh2).nfunction(gc2);
	  }
	}
	bf1_offset += basis->shell(sh1).nfunction(gc1);
      }
    }
}

void
compare_2e_cints_vs_v3(Ref<TwoBodyInt>& tbcints, Ref<TwoBodyInt>& tbv3)
{
  const double *buffercints = tbcints->buffer();
  const double *bufferv3 = tbv3->buffer();

  Ref<GaussianBasisSet> basis = tbcints->basis();
  for (int sh1=0; sh1<basis->nshell(); sh1++)
    for (int sh2=0; sh2<basis->nshell(); sh2++)
      for (int sh3=0; sh3<basis->nshell(); sh3++)
	for (int sh4=0; sh4<basis->nshell(); sh4++)
	  {
	    //sh1=0;sh2=0;sh3=8;sh4=3;
	    //cout << scprintf("Computing TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",sh1,sh2,sh3,sh4);
	    tbv3->compute_shell(sh1,sh2,sh3,sh4);
	    tbcints->compute_shell(sh1,sh2,sh3,sh4);

	    int nbf2 = basis->shell(sh2).nfunction();
	    int nbf3 = basis->shell(sh3).nfunction();
	    int nbf4 = basis->shell(sh4).nfunction();

	    int bf1_offset = 0;
	    for(int gc1=0;gc1<basis->shell(sh1).ncontraction(); gc1++) {
	      int am1 = basis->shell(sh1).am(gc1);
	      CartesianIterCints citer1(am1);
	      CartesianIterV3 iter1(am1);
	      for ( citer1.start(); int(citer1) ; citer1.next() ) {
		int bf1cints = citer1.bfn();
		int bf1v3;
		for( iter1.start(); int(iter1) ; iter1.next() ) {
		  if (iter1.a() == citer1.a() &&
		      iter1.b() == citer1.b() &&
		      iter1.c() == citer1.c()) {
		    bf1v3 = iter1.bfn();
		    break;
		  }
		}
		bf1cints += bf1_offset;
		bf1v3 += bf1_offset;
	      
		int bf2_offset = 0;
		for(int gc2=0;gc2<basis->shell(sh2).ncontraction(); gc2++) {
		  int am2 = basis->shell(sh2).am(gc2);
		  CartesianIterCints citer2(am2);
		  CartesianIterV3 iter2(am2);
		  for ( citer2.start(); int(citer2) ; citer2.next() ) {
		    int bf2cints = citer2.bfn();
		    int bf2v3;
		    for( iter2.start(); int(iter2) ; iter2.next() ) {
		      if (iter2.a() == citer2.a() &&
			  iter2.b() == citer2.b() &&
			  iter2.c() == citer2.c()) {
			bf2v3 = iter2.bfn();
			break;
		      }
		    }
		    bf2cints += bf2_offset;
		    bf2v3 += bf2_offset;
		    
		    int bf3_offset = 0;
		    for(int gc3=0;gc3<basis->shell(sh3).ncontraction(); gc3++) {
		      int am3 = basis->shell(sh3).am(gc3);
		      CartesianIterCints citer3(am3);
		      CartesianIterV3 iter3(am3);
		      for ( citer3.start(); int(citer3) ; citer3.next() ) {
			int bf3cints = citer3.bfn();
			int bf3v3;
			for( iter3.start(); int(iter3) ; iter3.next() ) {
			  if (iter3.a() == citer3.a() &&
			      iter3.b() == citer3.b() &&
			      iter3.c() == citer3.c()) {
			    bf3v3 = iter3.bfn();
			    break;
			  }
			}
			bf3cints += bf3_offset;
			bf3v3 += bf3_offset;
			    
			int bf4_offset = 0;
			for(int gc4=0;gc4<basis->shell(sh4).ncontraction(); gc4++) {
			  int am4 = basis->shell(sh4).am(gc4);
			  CartesianIterCints citer4(am4);
			  CartesianIterV3 iter4(am4);
			  for ( citer4.start(); int(citer4) ; citer4.next() ) {
			    int bf4cints = citer4.bfn();
			    int bf4v3;
			    for( iter4.start(); int(iter4) ; iter4.next() ) {
			      if (iter4.a() == citer4.a() &&
				  iter4.b() == citer4.b() &&
				  iter4.c() == citer4.c()) {
				bf4v3 = iter4.bfn();
				break;
			      }
			    }
			    bf4cints += bf4_offset;
			    bf4v3 += bf4_offset;
		
			    double valuecints = buffercints[((bf1cints*nbf2 + bf2cints)*nbf3 + bf3cints)*nbf4 + bf4cints];
			    double valuev3 = bufferv3[((bf1v3*nbf2 + bf2v3)*nbf3 + bf3v3)*nbf4 + bf4v3];
			    if (fabs(valuecints-valuev3) > 1E-12) {
			      cout << scprintf("Discrepancy in TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",sh1,sh2,sh3,sh4);
			      cout << scprintf("bf1 = %d  bf2 = %d  bf3 = %d  bf4 = %d  TEIntegral(cints) = %20.15lf\n",
					       bf1cints,bf2cints,bf3cints,bf4cints,valuecints);
			      cout << scprintf("bf1 = %d  bf2 = %d  bf3 = %d  bf4 = %d  TEIntegral(V3)    = %20.15lf\n\n",
					       bf1v3,bf2v3,bf3v3,bf4v3,valuev3);
			    }
			  }
			  bf4_offset+=basis->shell(sh4).nfunction(gc4);
			}
		      }
		      bf3_offset+=basis->shell(sh3).nfunction(gc3);
		    }
		  }
		  bf2_offset+=basis->shell(sh2).nfunction(gc2);
		}
	      }
	      bf1_offset+=basis->shell(sh1).nfunction(gc1);
	    }
	    //return;
	  }
}


void
compare_2e_puream_cints_vs_v3(Ref<TwoBodyInt>& tbcints, Ref<TwoBodyInt>& tbv3)
{
  const double *buffercints = tbcints->buffer();
  const double *bufferv3 = tbv3->buffer();

  Ref<GaussianBasisSet> basis = tbcints->basis();
  for (int sh1=0; sh1<basis->nshell(); sh1++)
    for (int sh2=0; sh2<basis->nshell(); sh2++)
      for (int sh3=0; sh3<basis->nshell(); sh3++)
	for (int sh4=0; sh4<basis->nshell(); sh4++)
	  {
	    //	    	    	    sh1 = 0; sh2 = 0; sh3 = 6; sh4 = 13;
			    /*	    if ( !((basis->shell(sh1).has_pure() || basis->shell(sh1).max_am()==0) &&
		   (basis->shell(sh2).has_pure() || basis->shell(sh2).max_am()==0) &&
		   (basis->shell(sh3).has_pure() || basis->shell(sh3).max_am()==0) &&
		   (basis->shell(sh4).has_pure() || basis->shell(sh4).max_am()==0)) )
		   continue;*/
	    tbcints->compute_shell(sh1,sh2,sh3,sh4);
	    tbv3->compute_shell(sh1,sh2,sh3,sh4);

	    int nbf1 = basis->shell(sh1).nfunction();
	    int nbf2 = basis->shell(sh2).nfunction();
	    int nbf3 = basis->shell(sh3).nfunction();
	    int nbf4 = basis->shell(sh4).nfunction();

	    	    cout << scprintf("Computing TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",sh1,sh2,sh3,sh4);
	    	    cout << scprintf("size = %d\n",nbf1*nbf2*nbf3*nbf4);
	    for(int ijkl=0;ijkl<nbf1*nbf2*nbf3*nbf4;ijkl++) {
	      double valuecints = buffercints[ijkl];
	      double valuev3 = bufferv3[ijkl];
	      if (fabs(valuecints-valuev3) > 1E-13) {
		cout << scprintf("Discrepancy in TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",sh1,sh2,sh3,sh4);
		cout << scprintf("1234 = %d TEIntegral(cints) = %20.15lf\n",
				 ijkl,valuecints);
		cout << scprintf("TEIntegral(V3)    = %20.15lf\n\n",
				 valuev3);
	      }
	    }
	    //  	    	    return;
	  }
}


void
compare_2e_bufsum_cints_vs_v3(Ref<TwoBodyInt>& tbcints, Ref<TwoBodyInt>& tbv3)
{
  Ref<GaussianBasisSet> basis = tbcints->basis();
  const double *buffercints = tbcints->buffer();
  const double *bufferv3 = tbv3->buffer();

  for (int sh1=0; sh1<basis->nshell(); sh1++)
    for (int sh2=0; sh2<basis->nshell(); sh2++)
      for (int sh3=0; sh3<basis->nshell(); sh3++)
	for (int sh4=0; sh4<basis->nshell(); sh4++)
	  {
	    //	    sh1=12;sh2=12;sh3=12;sh4=12;
	    tbcints->compute_shell(sh1,sh2,sh3,sh4);
	    tbv3->compute_shell(sh1,sh2,sh3,sh4);

	    int nbf1 = basis->shell(sh1).nfunction();
	    int nbf2 = basis->shell(sh2).nfunction();
	    int nbf3 = basis->shell(sh3).nfunction();
	    int nbf4 = basis->shell(sh4).nfunction();
	    
	    double sum_cints = 0.0;
	    double sum_v3 = 0.0;

	    int index = 0;
	    for (int i=0; i<nbf1; i++) {
	      for (int j=0; j<nbf2; j++) {
		for (int k=0; k<nbf3; k++) {
		  for (int l=0; l<nbf4; l++) {
		    sum_cints += buffercints[index];
		    sum_v3 += bufferv3[index];
		    /*		    cout << scprintf("index = %d TEIntegral(cints) = %20.15lf\n",
				     index,buffercints[index]);
		    cout << scprintf("index = %d TEIntegral(V3) = %20.15lf\n\n",
				     index,bufferv3[index]);
		    */
		    index++;
		  }
		}
	      }
	    }
	    
	    if (fabs(sum_cints-sum_v3) > 1E-12) {
	      cout << scprintf("Discrepancy in TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",sh1,sh2,sh3,sh4);
	      cout << scprintf("TEIntegralSum(cints) = %20.15lf\n",
			       sum_cints);
	      cout << scprintf("TEIntegralSum(V3)    = %20.15lf\n\n",
			       sum_v3);
	    }
	    //return;
	  }
}

void
compare_2e_unique_bufsum_cints_vs_v3(Ref<TwoBodyInt>& tbcints, Ref<TwoBodyInt>& tbv3)
{
  Ref<GaussianBasisSet> basis = tbcints->basis();
  const double *buffercints = tbcints->buffer();
  const double *bufferv3 = tbv3->buffer();

  for (int sh1=0; sh1<basis->nshell(); sh1++)
    for (int sh2=0; sh2<=sh1; sh2++)
      for (int sh3=0; sh3<=sh1; sh3++)
	for (int sh4=0; sh4 <= ((sh1==sh3) ? sh2 : sh3) ; sh4++)
	  {
	    tbcints->compute_shell(sh1,sh2,sh3,sh4);
	    tbv3->compute_shell(sh1,sh2,sh3,sh4);

	    int nbf1 = basis->shell(sh1).nfunction();
	    int nbf2 = basis->shell(sh2).nfunction();
	    int nbf3 = basis->shell(sh3).nfunction();
	    int nbf4 = basis->shell(sh4).nfunction();
	    
	    double sum_cints = 0.0;
	    double sum_v3 = 0.0;

	    int e12 = (sh1 == sh2) ? 1 : 0;
	    int e34 = (sh3 == sh4) ? 1 : 0;
	    int e13e24 = (((sh1 == sh3)&&(sh2==sh4))||((sh1==sh4)&&(sh2==sh3))) ? 1 : 0;

	    int index = 0;
	    for (int i=0; i<nbf1; i++) {
	      int jmax = e12 ? i : nbf2-1;
	      for (int j=0; j<=jmax; j++) {
		int kmax = e13e24 ? i : nbf3-1;
		for (int k=0; k<=kmax; k++) {
		  int lmax = e34 ? ( (e13e24&&(i==k)) ? j : k) : ( (e13e24&&(i==k)) ? j : nbf4-1);
		  for (int l=0; l<=lmax; l++) {
		    sum_cints += buffercints[index];
		    sum_v3 += bufferv3[index];
		    /*		    cout << scprintf("index = %d TEIntegral(cints) = %20.15lf\n",
				     index,buffercints[index]);
		    cout << scprintf("index = %d TEIntegral(V3) = %20.15lf\n\n",
				     index,bufferv3[index]);
		    */
		    index++;
		  }
		}
	      }
	    }
	    
	    if (fabs(sum_cints-sum_v3) > 1E-12) {
	      cout << scprintf("Discrepancy in TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",sh1,sh2,sh3,sh4);
	      cout << scprintf("TEIntegralSum(cints) = %20.15lf\n",
			       sum_cints);
	      cout << scprintf("TEIntegralSum(V3)    = %20.15lf\n\n",
			       sum_v3);
	    }
	  }
}

void
print_grt_ints(Ref<TwoBodyInt>& tbcints)
{
  Ref<GaussianBasisSet> basis = tbcints->basis();
  const double *buffer[4];
  buffer[0] = tbcints->buffer(TwoBodyInt::eri);
  buffer[1] = tbcints->buffer(TwoBodyInt::r12);
  buffer[2] = tbcints->buffer(TwoBodyInt::r12t1);
  buffer[3] = tbcints->buffer(TwoBodyInt::r12t2);
  char teout_filename[] = "teout0.dat";
  FILE *teout[4];

  for(int te_type=0;te_type<4;te_type++) {
    teout_filename[5] = te_type + '0';
    teout[te_type] = fopen(teout_filename,"w");
  }

  for (int ush1=0; ush1<basis->nshell(); ush1++)
    for (int ush2=0; ush2<=ush1; ush2++)
      for (int ush3=0; ush3<=ush2; ush3++)
	for (int ush4=0; ush4 <=ush3 ; ush4++)
	  {
	    int S1[3], S2[3], S3[3], S4[4];
	    int num = 1;
	    S1[0] = ush1;
	    S2[0] = ush2;
	    S3[0] = ush3;
	    S4[0] = ush4;
	    if (ush1==ush2 && ush1==ush3 || ush2==ush3 && ush2==ush4)
	      num=1;
	    else if (ush1==ush3 || ush2==ush4) {
	      num = 2;
	      S1[1] = ush1;
	      S2[1] = ush3;
	      S3[1] = ush2;
	      S4[1] = ush4;
	    }
	    else if (ush2==ush3) {
	      num = 2;
	      S1[1] = ush1;
	      S2[1] = ush4;
	      S3[1] = ush2;
	      S4[1] = ush3;
	    }
	    else if (ush1==ush2 || ush3==ush4) {
	      num = 2;
	      S1[1] = ush1;
	      S2[1] = ush3;
	      S3[1] = ush2;
	      S4[1] = ush4;
	    }
	    else {
	      num = 3;
	      S1[1] = ush1;
	      S2[1] = ush3;
	      S3[1] = ush2;
	      S4[1] = ush4;
	      S1[2] = ush1;
	      S2[2] = ush4;
	      S3[2] = ush2;
	      S4[2] = ush3;
	    }
	      
	    for(int uq=0;uq<num;uq++) {
	      int sh1 = S1[uq];
	      int sh2 = S2[uq];
	      int sh3 = S3[uq];
	      int sh4 = S4[uq];

	      tbcints->compute_shell(sh1,sh2,sh3,sh4);

	      int nbf1 = basis->shell(sh1).nfunction();
	      int nbf2 = basis->shell(sh2).nfunction();
	      int nbf3 = basis->shell(sh3).nfunction();
	      int nbf4 = basis->shell(sh4).nfunction();
	    
	      int e12 = (sh1 == sh2) ? 1 : 0;
	      int e34 = (sh3 == sh4) ? 1 : 0;
	      int e13e24 = (((sh1 == sh3)&&(sh2==sh4))||((sh1==sh4)&&(sh2==sh3))) ? 1 : 0;

	      int index = 0;
	      for (int i=0; i<nbf1; i++) {
		int jmax = e12 ? i : nbf2-1;
		for (int j=0; j<=jmax; j++) {
		  int kmax = e13e24 ? i : nbf3-1;
		  for (int k=0; k<=kmax; k++) {
		    int lmax = e34 ? ( (e13e24&&(i==k)) ? j : k) : ( (e13e24&&(i==k)) ? j : nbf4-1);
		    for (int l=0; l<=lmax; l++) {
		      
		      double integral = buffer[0][index];
		      if (fabs(integral) > 1E-15) {
			fprintf(teout[0], "%5d%5d%5d%5d%20.10lf\n",
				basis->shell_to_function(sh1) + i+1, 
				basis->shell_to_function(sh2) + j+1, 
				basis->shell_to_function(sh3) + k+1, 
				basis->shell_to_function(sh4) + l+1, 
				integral);
		      }
		      integral = buffer[1][index];
		      if (fabs(integral) > 1E-15) {
			fprintf(teout[1], "%5d%5d%5d%5d%20.10lf\n",
				basis->shell_to_function(sh1) + i+1, 
				basis->shell_to_function(sh2) + j+1, 
				basis->shell_to_function(sh3) + k+1, 
				basis->shell_to_function(sh4) + l+1, 
				integral);
		      }
		      integral = buffer[2][index];
		      if (fabs(integral) > 1E-15) {
			fprintf(teout[2], "%5d%5d%5d%5d%20.10lf\n",
				basis->shell_to_function(sh1) + i+1, 
				basis->shell_to_function(sh2) + j+1, 
				basis->shell_to_function(sh3) + k+1, 
				basis->shell_to_function(sh4) + l+1, 
				integral);
		      }
		      integral = buffer[3][index];
		      if (fabs(integral) > 1E-15) {
			fprintf(teout[3], "%5d%5d%5d%5d%20.10lf\n",
				basis->shell_to_function(sh1) + i+1, 
				basis->shell_to_function(sh2) + j+1, 
				basis->shell_to_function(sh3) + k+1, 
				basis->shell_to_function(sh4) + l+1, 
				integral);
		      }
		      index++;
		    }
		  }
		}
	      }
	    }
	  }
  for(int te_type=0;te_type<4;te_type++)
    fclose(teout[te_type]);
}

void
compare_2e_permute(Ref<Integral>& cints)
{
  Ref<TwoBodyInt> tb1 = cints->electron_repulsion();
  Ref<TwoBodyInt> tb2 = cints->electron_repulsion();
  Ref<GaussianBasisSet> basis = tb1->basis();
  const double *buffer1 = tb1->buffer();
  const double *buffer2 = tb2->buffer();

  int sh1 = 0;
  int sh2 = 0;
  int sh3 = 4;
  int sh4 = 0;

  tb1->compute_shell(sh1,sh2,sh3,sh4);
  tb2->compute_shell(sh1,sh2,sh4,sh3);

  int nbf1 = basis->shell(sh1).nfunction();
  int nbf2 = basis->shell(sh2).nfunction();
  int nbf3 = basis->shell(sh3).nfunction();
  int nbf4 = basis->shell(sh4).nfunction();
	    
  for(int index = 0; index<nbf1*nbf2*nbf3*nbf4; index++)
    if (fabs(buffer1[index]-buffer2[index]) > 1E-13)
    {
      cout << scprintf("Discrepancy in TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",sh1,sh2,sh3,sh4);
      cout << scprintf("TEIntegral(cints1)    = %20.15lf\n",buffer1[index]);
      cout << scprintf("TEIntegral(cints2)    = %20.15lf\n\n",buffer2[index]);
    }
}


void
do_shell_test_1e(const Ref<Int1eV3> &int1ev3,
                 void (Int1eV3::*int_shell_1e)(int,int),
                 int permute, int i, int j, int na, int nb,
                 double *buf, double *pbuf)
{
  int ii = 0;
  int a;
  double *buffer = int1ev3->buffer();
  (int1ev3->*int_shell_1e)(i, j);
  for (a=0; a<na*nb; a++) {
      buf[a] = buffer[a];
    }
  (int1ev3->*int_shell_1e)(j, i);
  for (a=0; a<na*nb; a++) {
      pbuf[a] = buffer[a];
    }
  for (a=0; a<na; a++) {
      for (int b=0; b<nb; b++) {
          if (fabs(buf[ii] - pbuf[a + na*b]) > 1.0e-13) {
              cout << scprintf("----- 1e perm failed:"
                               "<%d %d|%d %d>:"
                               " %18.14f != %18.14f "
                               "<%d %d|%d %d>\n",
                               i, a, j, b,
                               buf[ii],
                               pbuf[a + na*b],
                               j, b, i, a);
            }
          if (fabs(buf[ii]) > 1.0e-15) {
              cout << scprintf(" <(%d %d)|(%d %d)> = %15.11f\n",
                               i,a,j,b, buf[ii]);
            }
          ii++;
        }
    }
}

void
test_int_shell_1e(const Ref<KeyVal>& keyval, const Ref<Int1eV3> &int1ev3,
                  void (Int1eV3::*int_shell_1e)(int,int),
                  int permute)
{
  int flags = 0;
  Ref<GaussianBasisSet> basis = int1ev3->basis();
  int maxfunc = basis->max_nfunction_in_shell();
  int size = maxfunc * maxfunc;
  double *buf = new double[size];
  double *pbuf = new double[size];
  int nshell = int1ev3->basis()->nshell();

  for (int i=0; i<nshell; i++) {
      int na = basis->shell(i).nfunction();
      for (int j=0; j<nshell; j++) {
          int nb = basis->shell(j).nfunction();
          do_shell_test_1e(int1ev3, int_shell_1e, permute,
                           i, j, na, nb, buf, pbuf);

        }
    }

  delete[] buf;
  delete[] pbuf;
}

void
test_3_center(const Ref<KeyVal>& keyval, const Ref<Int2eV3> &int2ev3)
{
  int ii, i,j,k,l,m,n;

  int2ev3->set_redundant(1);
  int2ev3->set_permute(0);
  double *buffer = int2ev3->buffer();
  int nshell = int2ev3->basis()->nshell();

  for (i=0; i<nshell; i++) {
      for (j=0; j<nshell; j++) {
          int sh[2], sizes[2];
          sh[0] = i;
          sh[1] = j;
          int2ev3->erep_2center(sh,sizes);
          ii = 0;
          for (k=0; k<sizes[0]; k++) {
              for (l=0; l<sizes[1]; l++) {
                  if (fabs(buffer[ii])>1.0e-15)
                      cout << scprintf(" ((%d %d)|(%d %d)) = %15.11f\n",
                                       sh[0],k,sh[1],l, buffer[ii]);
                  ii++;
                }
            }
        }
    }

  for (i=0; i<nshell; i++) {
      for (j=0; j<nshell; j++) {
          for (m=0; m<nshell; m++) {
              int sh[3], sizes[3];
              sh[0] = i;
              sh[1] = j;
              sh[2] = m;
              int2ev3->erep_3center(sh,sizes);
              ii = 0;
              for (k=0; k<sizes[0]; k++) {
                  for (l=0; l<sizes[1]; l++) {
                      for (n=0; n<sizes[2]; n++) {
                          if (fabs(buffer[ii])>1.0e-15)
                              cout << scprintf(
                                  " ((%d %d)|(%d %d)(%d %d)) = %15.11f\n",
                                  sh[0],k,sh[1],l,sh[2],n, buffer[ii]);
                          ii++;
                        }
                    }
                }
            }
        }
    }

}

void
init_shell_perm(const Ref<Int2eV3> &int2ev3, double *integrals,
                double buff[maxint][maxint][maxint][maxint],
                int sh[4], int sizes[4])
{
  int i, j, k, l;
  int oldp = int2ev3->permute();
  int2ev3->set_permute(0);
  int2ev3->erep(sh, sizes);
  int2ev3->set_permute(oldp);
  for (i=0; i<sizes[0]; i++) {
      for (j=0; j<sizes[1]; j++) {
          for (k=0; k<sizes[2]; k++) {
              for (l=0; l<sizes[3]; l++) {
                  buff[i][j][k][l] = *integrals++;
                }
            }
        }
    }
}

void
check_shell_perm(const Ref<Int2eV3> &int2ev3, double *integrals,
                 double buff[maxint][maxint][maxint][maxint],
                 int sh[4], int sizes[4], int p0, int p1, int p2, int p3)
{
  int ip[4], p[4];
  int psizes[4];
  int psh[4];
  int index = 0;
  int i[4];
  p[0] = p0;
  p[1] = p1;
  p[2] = p2;
  p[3] = p3;
  ip[p0] = 0;
  ip[p1] = 1;
  ip[p2] = 2;
  ip[p3] = 3;
  psh[0] = sh[p0];
  psh[1] = sh[p1];
  psh[2] = sh[p2];
  psh[3] = sh[p3];
  int oldp = int2ev3->permute();
  int2ev3->set_permute(0);
  int2ev3->erep(psh, psizes);
  int2ev3->set_permute(oldp);
  for (i[0]=0; i[0]<psizes[0]; i[0]++) {
      for (i[1]=0; i[1]<psizes[1]; i[1]++) {
          for (i[2]=0; i[2]<psizes[2]; i[2]++) {
              for (i[3]=0; i[3]<psizes[3]; i[3]++) {
                  if (fabs(buff[i[ip[0]]][i[ip[1]]][i[ip[2]]][i[ip[3]]]
                           - integrals[index]) > 1.0e-13) {
                      cout << scprintf("perm %d %d %d %d failed:"
                             "((%d %d)(%d %d)|(%d %d)(%d %d)):"
                             " %18.14f != %18.14f "
                             "((%d %d)(%d %d)|(%d %d)(%d %d))\n",
                             p0, p1, p2, p3,
                             sh[0],i[0], sh[1],i[1], sh[2],i[2], sh[3],i[3],
                             buff[i[ip[0]]][i[ip[1]]][i[ip[2]]][i[ip[3]]],
                             integrals[index],
                             psh[0],i[p[0]], psh[1],i[p[1]],
                             psh[2],i[p[2]], psh[3],i[p[3]]);
                    }
                  index++;
                }
            }
        }
    }
}

void
do_shell_quartet_test(const Ref<Int2eV3> &int2ev3,
                      int print, int printbounds, int bounds, int permute,
                      const Ref<KeyVal>& keyval,
                      int i, int j, int k, int l)
{
  int sh[4], sizes[4];
  int ibuf;
  int ii, jj, kk, ll;
  sh[0] = i;
  sh[1] = j;
  sh[2] = k;
  sh[3] = l;
  double maxintegral, integralbound;
  int boundijkl;
  if (bounds) {
      integralbound
          = int2ev3->logbound_to_bound(
              (boundijkl = int2ev3->erep_4bound(i,j,k,l))
              );
    }
  double *buffer = int2ev3->buffer();
  int2ev3->erep(sh,sizes);
  ibuf = 0;
  maxintegral = 0.0;
  for (ii=0; ii<sizes[0]; ii++) {
      for (jj=0; jj<sizes[1]; jj++) {
          for (kk=0; kk<sizes[2]; kk++) {
              for (ll=0; ll<sizes[3]; ll++) {
                  double absint = fabs(buffer[ibuf]);
                  if (absint > maxintegral) {
                      maxintegral = absint;
                    }
                  if (bounds &&  absint > integralbound) {
                      cout << scprintf("((%d %d)(%d %d)|(%d %d)(%d %d)) = %15.11f, "
                             "bound = %15.11f\n",
                             sh[0], ii, sh[1], jj, sh[2], kk, sh[3], ll,
                             buffer[ibuf], integralbound);
                      abort();
                    }
                  if (print && (absint > 1.0e-9
                                || (bounds && integralbound > 1.0e-9))) {
                      cout << scprintf(" ((%d %d)(%d %d)|(%d %d)(%d %d))"
                             " = %15.11f",
                             sh[0],ii,
                             sh[1],jj,
                             sh[2],kk,
                             sh[3],ll,
                             buffer[ibuf]);
                      if (bounds) {
                          cout << scprintf(" (%2d%% of bound)",
                                 (int)(100*(absint/integralbound)));
                        }
                      cout << scprintf("\n");
                    }
                  ibuf++;
                }
            }
        }
    }

  if (permute) {
      double buff1[maxint][maxint][maxint][maxint];
      sh[0] = i;
      sh[1] = j;
      sh[2] = k;
      sh[3] = l;
      init_shell_perm(int2ev3, buffer, buff1, sh, sizes);
      check_shell_perm(int2ev3, buffer, buff1, sh, sizes, 0, 1, 2, 3);
      check_shell_perm(int2ev3, buffer, buff1, sh, sizes, 1, 0, 2, 3);
      check_shell_perm(int2ev3, buffer, buff1, sh, sizes, 0, 1, 3, 2);
      check_shell_perm(int2ev3, buffer, buff1, sh, sizes, 1, 0, 3, 2);
      check_shell_perm(int2ev3, buffer, buff1, sh, sizes, 2, 3, 0, 1);
      check_shell_perm(int2ev3, buffer, buff1, sh, sizes, 2, 3, 1, 0);
      check_shell_perm(int2ev3, buffer, buff1, sh, sizes, 3, 2, 0, 1);
      check_shell_perm(int2ev3, buffer, buff1, sh, sizes, 3, 2, 1, 0);
    }

  if (bounds) {
      int boundij = int2ev3->erep_4bound(i,j,-1,-1);
      int boundkl = int2ev3->erep_4bound(-1,-1,k,l);
      int badbound = 0;
      if (boundij < boundijkl || boundkl < boundijkl) {
          badbound = 1;
        }
      if (badbound || printbounds) {
          cout << scprintf("max(%d,%d,%d,%d)=%7.4f, bnd=%7.4f, "
                 "bnd(%d,%d,*,*)=%7.4f, bnd(*,*,%d,%d)=%7.4f\n",
                 i, j, k, l, maxintegral, integralbound,
                 i,j, int2ev3->logbound_to_bound(boundij),
                 k,l, int2ev3->logbound_to_bound(boundkl));
        }
      if (badbound) {
          cout << scprintf("ERROR: bad bound\n");
          abort();
        }
    }
}

void
do_4_center_test(const Ref<Int2eV3> &int2ev3, int print, int printbounds,
                 int bounds, int permute,
                 const Ref<KeyVal>& keyval)
{
  int ii,jj,kk,ll, i,j,k,l, ibuf;
  int nshell = int2ev3->basis()->nshell();
  int unique = keyval->booleanvalue("unique");
  int timestats = keyval->booleanvalue("timestats");
  Ref<RegionTimer> timer = new RegionTimer();

  if (!timestats) {
      for (i=0; i<nshell; i++) {
          int jmax = nshell - 1;
          if (unique) jmax = i;
          for (j=0; j<=jmax; j++) {
              int kmax = nshell - 1;
              if (unique) kmax = i;
              for (k=0; k<=kmax; k++) {
                  int lmax = nshell - 1;
                  if (unique) {
                      if (k==i) lmax = j;
                      else lmax = k;
                    }
                  for (l=0; l<=lmax; l++) {
                      do_shell_quartet_test(int2ev3, print, printbounds,
                                            bounds, permute,
                                            keyval, i, j, k, l);
                    }
                }
            }
        }
    }
  if (timestats && nshell) {
      unsigned short seed = 1234;
      seed48(&seed);
      const int nsample = 5000;
      const int ntrials = 50;
      double times[ntrials];
      for (i=0; i<ntrials; i++) {
          double t1 = timer->get_cpu_time();
          for (j=0; j<nsample; j++) {
              // pick an integral at random
              int ish = int(drand48()*nshell);
              int jsh = int(drand48()*ish);
              int ksh = int(drand48()*ish);
              int lsh;
              if (ish==ksh) lsh = int(drand48()*jsh);
              else lsh = int(drand48()*ksh);
              int sh[4], sizes[4];
              if (ish >= nshell) ish = nshell-1;
              if (jsh >= nshell) jsh = nshell-1;
              if (ksh >= nshell) ksh = nshell-1;
              if (lsh >= nshell) lsh = nshell-1;
              sh[0] = ish;
              sh[1] = jsh;
              sh[2] = ksh;
              sh[3] = lsh;
              int2ev3->erep(sh,sizes);
            }
          double t2 = timer->get_cpu_time();
          times[i] = t2-t1;
        }
      double ave = 0.0;
      for (i=0; i<ntrials; i++) {
          ave += times[i];
        }
      ave /= ntrials;
      double sigma2 = 0.0;
      for (i=0; i<ntrials; i++) {
          double diff = times[i] - ave;
          sigma2 += diff*diff;
        }
      double sigma = sqrt(sigma2/ntrials);
      // adjust sigma and ave from the trial average results to
      // the integral results
      ave /= nsample;
      sigma /= sqrt(double(nsample));
      cout << scprintf(" ave = %10.8f sigma = %10.8f (microsecs)\n"
                       " sigma/ave = %10.4f",
                       ave*1e6, sigma*1e6,
                       sigma/ave)
           << endl;
    }
}

void
test_4_center(const Ref<KeyVal>& keyval, const Ref<Int2eV3> &int2ev3)
{
  int i;

  cout << scprintf("4 center test:\n");
  cout << scprintf("  on entry int2ev3 used %d bytes\n", int2ev3->used_storage());

  int2ev3->set_permute(0);
  int2ev3->set_redundant(1);

  int storage = keyval->intvalue("storage") - int2ev3->used_storage();
  if (storage < 0) storage = 0;
  if (keyval->booleanvalue("store_integrals")) storage = 0;
  int niter = keyval->intvalue("niter");
  int print = keyval->booleanvalue("print");
  int bounds = keyval->booleanvalue("bounds");
  int permute = keyval->booleanvalue("permute");
  int printbounds = keyval->booleanvalue("printbounds");

  cout << scprintf("  storage   = %d\n", storage);
  cout << scprintf("  niter     = %d\n", niter);
  cout << scprintf("  print     = %d\n", print);
  cout << scprintf("  bounds    = %d\n", bounds);
  cout << scprintf("  permute   = %d\n", permute);
  cout << scprintf("printbounds = %d\n", printbounds);

  if (bounds) int2ev3->init_bounds();

  int2ev3->init_storage(storage);

  for (i=0; i<niter; i++) {
      do_4_center_test(int2ev3, print, printbounds, bounds, permute, keyval);
    }

  if (keyval->count("quartet") == 4) {
      do_shell_quartet_test(int2ev3, print, printbounds, bounds, permute,
                            keyval,
                            keyval->intvalue("quartet", 0),
                            keyval->intvalue("quartet", 1),
                            keyval->intvalue("quartet", 2),
                            keyval->intvalue("quartet", 3));
    }

  int2ev3->done_storage();
  int2ev3->done_bounds();
}

void
do_shell_quartet_der_test(const Ref<Int2eV3> &int2ev3,
                          double* buffer, int print, int printbounds,
                          int bounds, int permute,
                          const Ref<KeyVal>& keyval,
                          int i, int j, int k, int l)
{
  int ii,jj,kk,ll, ibuf, ider, xyz;
  der_centersv3_t dercenters;

  int sh[4], sizes[4];
  sh[0] = i;
  sh[1] = j;
  sh[2] = k;
  sh[3] = l;
  double maxintegral = 0.0, integralbound;
  int boundijkl;
  if (bounds) {
      integralbound
          = int2ev3->logbound_to_bound(
              (boundijkl = int2ev3->erep_4bound_1der(i,j,k,l))
              );
    }
  int2ev3->erep_all1der(sh,sizes,&dercenters);
  ibuf = 0;
  for (ider=0; ider<dercenters.n; ider++) {
      for (xyz=0; xyz<3; xyz++) {
          for (ii=0; ii<sizes[0]; ii++) {
              for (jj=0; jj<sizes[1]; jj++) {
                  for (kk=0; kk<sizes[2]; kk++) {
                      for (ll=0; ll<sizes[3]; ll++) {
                          double absint = fabs(buffer[ibuf]);
                          if (absint > maxintegral) {
                              maxintegral = absint;
                            }
                          if (bounds &&  absint > integralbound) {
                              cout << scprintf("((%d %d)(%d %d)|(%d %d)(%d %d))"
                                     " = %15.11f, bound = %15.11f\n",
                                     sh[0], ii, sh[1], jj,
                                     sh[2], kk, sh[3], ll,
                                     buffer[ibuf], integralbound);
                              abort();
                            }
                          if (print && absint > 1.0e-15) {
                              cout << scprintf(" ((%d %d)(%d %d)"
                                     "|(%d %d)(%d %d))(%d %d)"
                                     " = %15.11f\n",
                                     sh[0],ii,
                                     sh[1],jj,
                                     sh[2],kk,
                                     sh[3],ll,
                                     dercenters.num[ider], xyz,
                                     buffer[ibuf]
                                  );
                            }
                          ibuf++;
                        }
                    }
                }
            }
        }
    }

  if (bounds) {
      int boundij = int2ev3->erep_4bound_1der(i,j,-1,-1);
      int boundkl = int2ev3->erep_4bound_1der(-1,-1,k,l);
      int badbound = 0;
      if (boundij < boundijkl || boundkl < boundijkl) {
          badbound = 1;
        }
      if (badbound || printbounds) {
          cout << scprintf("max(%d,%d,%d,%d)=%7.4f, bnd=%7.4f, "
                 "bnd(%d,%d,*,*)=%8.4f, bnd(*,*,%d,%d)=%8.4f\n",
                 i, j, k, l, maxintegral, integralbound,
                 i,j, int2ev3->logbound_to_bound(boundij),
                 k,l, int2ev3->logbound_to_bound(boundkl));
        }
      if (badbound) {
          cout << scprintf("ERROR: bad bound\n");
          abort();
        }
    }
}

void
do_test_4der_center(const Ref<Int2eV3> &int2ev3,
                    double* buffer, int print, int printbounds,
                    int bounds, int permute,
                    const Ref<KeyVal>& keyval)
{
  int i,j,k,l;
  int nshell = int2ev3->basis()->nshell();
  for (i=0; i<nshell; i++) {
      for (j=0; j<nshell; j++) {
          for (k=0; k<nshell; k++) {
              for (l=0; l<nshell; l++) {
                  do_shell_quartet_der_test(int2ev3, buffer,
                                            print, printbounds,
                                            bounds, permute,
                                            keyval,
                                            i, j, k, l);
                }
            }
        }
    }
}

void
test_4der_center(const Ref<KeyVal>& keyval, const Ref<Int2eV3> &int2ev3)
{
  int i;

  int2ev3->set_permute(0);
  int2ev3->set_redundant(1);
  double *buffer = int2ev3->buffer();

  int niter = keyval->intvalue("niter");
  int print = keyval->booleanvalue("print");
  int bounds = keyval->booleanvalue("bounds");
  int printbounds = keyval->booleanvalue("printbounds");
  int permute = keyval->booleanvalue("permute");

  cout << scprintf("4 center derivative test:\n");
  cout << scprintf("  niter      = %d\n", niter);
  cout << scprintf("  print      = %d\n", print);
  cout << scprintf("  bounds     = %d\n", bounds);
  cout << scprintf("printbounds  = %d\n", printbounds);
  cout << scprintf("  permute    = %d\n", permute);

  if (bounds) int2ev3->init_bounds_1der();

  for (i=0; i<niter; i++) {
      do_test_4der_center(int2ev3, buffer,
                          print, printbounds, bounds, permute, keyval);
    }

  if (keyval->count("quartet") == 4) {
      do_shell_quartet_der_test(int2ev3, buffer, print, printbounds,
                                bounds, permute,
                                keyval,
                                keyval->intvalue("quartet", 0),
                                keyval->intvalue("quartet", 1),
                                keyval->intvalue("quartet", 2),
                                keyval->intvalue("quartet", 3));
    }

  if (bounds) int2ev3->done_bounds_1der();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
