//
// tascf_test.cpp
//
// Copyright (C) 2013 Drew Lewis
//
// Authors: Drew Lewis
// Maintainer: Drew Lewis and Edward Valeev
//
// This file is part of the MPQC Toolkit.
//
// The MPQC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The MPQC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the MPQC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <util/madness/init.h>
#include <util/misc/regtime.h>
#include <chemistry/qc/scf/tascf.hpp>
#define BOOST_TEST_MODULE test_tascf
#include <boost/test/included/unit_test.hpp>

#include <chemistry/qc/libint2/linkage.h>

using namespace boost::unit_test;
using namespace sc;
using namespace mpqc;

struct Mock_SCF : public TA::SCF {
    Mock_SCF(const Ref<KeyVal> &kval) : TA::SCF(kval) {}
    virtual double scf_energy() { return 2.0; }
    virtual double iter_energy() {return 2.1; }

    static sc::ClassDesc class_desc_;
};

sc::ClassDesc Mock_SCF::class_desc_(typeid(Mock_SCF), "Mock_SCF",
                      1, "public TA.SCF", 0, sc::create<Mock_SCF>, 0);

struct MADConfig {
    MADConfig() {
      int argc = boost::unit_test::framework::master_test_suite().argc;
      char** argv = boost::unit_test::framework::master_test_suite().argv;
      ExEnv::init(argc,argv);
      mpqc::MADNESSRuntime::initialize();
    }
    ~MADConfig() {
      mpqc::MADNESSRuntime::finalize();
    }
};

BOOST_GLOBAL_FIXTURE( MADConfig );

BOOST_AUTO_TEST_CASE( construct_scf_programmatically ){

    // Make a molecule H2
    Ref<Molecule> mol = new Molecule;
    mol->add_atom(1,0,1,-1);
    mol->add_atom(1,0,1,1);

    // Make keyval
    Ref<AssignedKeyVal> akv = new AssignedKeyVal;
    akv->assign("name", "STO-3G");
    akv->assign("molecule", mol.pointer());
    Ref<GaussianBasisSet> basis =
                    new GaussianBasisSet(Ref<KeyVal>(akv));
    akv->assign("basis", basis.pointer());

    //Construct object
    Ref<Mock_SCF> tscf = new Mock_SCF(akv);
    tscf->print();
}

BOOST_AUTO_TEST_CASE( construct_scf_txtkeyval ){

  const char *input = "./tascf_test.kv";
  Ref<KeyVal> kv = new ParsedKeyVal(input);
  Ref<Mock_SCF> tscf; tscf << kv->describedclassvalue("rhf");

  tscf->print();
  std::cout << "Overlap matrix:" << std::endl;
  std::cout << tscf->overlap() << std::endl;
  std::cout << "Hcore matrix:" << std::endl;
  std::cout << tscf->hcore() << std::endl;
}

BOOST_AUTO_TEST_CASE( compute_tawfn_overlap_txtkeyval ){

  const char *input = "./tascf_test.kv";
  Ref<KeyVal> kv = new ParsedKeyVal(input);
  Ref<Mock_SCF> tscf; tscf << kv->describedclassvalue("rhf_large");

  tscf->print();

  // compute overlap, time it
  sc::Ref<sc::RegionTimer> rtim = new sc::RegionTimer;
  std::cout << "Computing overlap matrix .............. ";
  sc::Timer tim(rtim, "rhf_large overlap");
  TiledArray::Array<double,2> S = tscf->overlap();
  tim.exit("rhf_large overlap");
  std::cout << "done (" << tim.wall_time("rhf_large overlap") << " sec)" << std::endl;

#if 0 // keeping this here till I think of somewhere to put it.
  auto it = S.begin();
  auto end = S.end();
  std::size_t global_sparse_elems = 0;
  double global_elems_saved = 0;
  double svd_cut = 1e-12;
  double sparse_cut = 1e-12;
  std::size_t total_elems = 0;
  while(it != end){
    decltype(S)::value_type tile = *it++;
    Eigen::Map<TiledArray::EigenMatrixXd>
              tile_map(tile.data(), tile.range().size()[0],
                       tile.range().size()[1]);
    Eigen::MatrixXd temp_mat = tile_map;
    auto tile_ptr = tile.data();
    decltype(tile_ptr) tile_end = tile.data()+tile.size();
    std::size_t sparse_elems = 0;
    while(tile_ptr != tile_end){
      if(*tile_ptr++ < sparse_cut){
        ++global_sparse_elems;
      }
      ++total_elems;
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(temp_mat);
    std::size_t counter = 0;
    for(auto i = 0; i < svd.singularValues().size(); ++i){
      if(svd.singularValues()[i] < svd_cut){
        ++counter;
      }
    }
    double elements_saved =  0.5 * double(std::min(tile.range().size()[0], tile.range().size()[1]) * counter);
    global_elems_saved += elements_saved;
  }
  std::cout << "\nThe total number of elements reduced by svd(" << svd_cut << ") is " << global_elems_saved << std::endl;
  std::cout << "\tFraction saved by SVD = " << global_elems_saved/double(total_elems) << std::endl;
  std::cout << "The total number of elements smaller than " << sparse_cut << " is " << global_sparse_elems << std::endl;
  std::cout << "\tFraction saved by elem truncation = " <<  global_sparse_elems/double(total_elems) << std::endl;
#endif

}










