//
// molcas_pt2r12.cc
//
// Copyright (C) 2014 Chong Peng
//
// Authors: Chong Peng
// Maintainer: Chong Peng and Edward Valeev
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

#include <mpqc_config.h>
#include <chemistry/qc/mbptr12/molcas_pt2r12.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/lcao/wfnworld.h>
#include <chemistry/qc/basis/union.h>
#include <extern/moinfo/moinfo.h>

using namespace sc;

ClassDesc MolcasPT2R12::class_desc_(typeid(MolcasPT2R12),
                                "MolcasPT2R12",
                                1,
                                "public ExternPT2R12",
                                0,
                                create<MolcasPT2R12>,
                                0
                                );

MolcasPT2R12::MolcasPT2R12 (const Ref<KeyVal>& kv) :
  ExternPT2R12(construct_extern_pt2r12(kv))
{
  prefix_ =  kv->stringvalue("prefix", KeyValValuestring(std::string()));

};


Ref<KeyVal> MolcasPT2R12::construct_extern_pt2r12(const Ref<KeyVal>& kv){

  std::string filename_prefix = kv->stringvalue("prefix", KeyValValuestring(std::string()) );
  std::string obs_name = kv->stringvalue("obs", KeyValValuestring(std::string()) );
  std::string cabs_name = kv->stringvalue("cabs", KeyValValuestring(std::string()) );
  if (cabs_name.empty() && not obs_name.empty()) {
    cabs_name = R12Technology::default_cabs_name(obs_name);
  }

  std::string cabs_contraction = kv->stringvalue("cabs_contraction", KeyValValuestring(std::string()) );
  std::string dfbs_name = kv->stringvalue("dfbs", KeyValValuestring(std::string()) );
  if (dfbs_name.empty() && not obs_name.empty())
    dfbs_name = DensityFittingRuntime::default_dfbs_name(obs_name, 1); // for OBS with cardinal number X use DFBS with cardinal number X+1
  if (dfbs_name == "none") dfbs_name = "";

  std::string f12exp = kv->stringvalue("f12exp", KeyValValuestring(std::string()) );
  // if OBS given but F12 exponent is not, look up a default value
  if (f12exp.empty() && not obs_name.empty()) {
    const double f12exp_default = R12Technology::default_stg_exponent(obs_name);
    if (f12exp_default != 0.0) {
      std::ostringstream oss;
      oss << f12exp_default;
      f12exp = oss.str();
    }
  }

  std::string r12 = kv->stringvalue("r12", KeyValValuestring(std::string()) );
#if defined(HAVE_MPQC3_RUNTIME)
  std::string cabs_singles = kv->stringvalue("cabs_singles", KeyValValuestring(std::string()) );
  std::string cabs_singles_basis = kv->stringvalue("singles_basis", KeyValValuestring(std::string()) );
  std::string partition = kv->stringvalue("patitionH", KeyValValuestring(std::string()) );
#endif

  Ref<Integral> integral;
  integral << kv->describedclassvalue("integrals");
  //
  // Read molecule, basis, and orbitals
  //
  Ref<ExternMOInfo> rdorbs = new ExternMOInfo(filename_prefix + ".pt2r12.dat",
          integral,
          obs_name); // all MO info is contained in rdorbs
  Ref<OrbitalSpace> orbs = rdorbs->orbs();
  Ref<GaussianBasisSet> basis = orbs->basis();
  RefSCMatrix C_ao = orbs->coefs();
  const std::vector<unsigned int>& fzcpi = rdorbs->fzcpi();
  const std::vector<unsigned int>& inactpi = rdorbs->inactpi();
  const std::vector<unsigned int>& actpi = rdorbs->actpi();
  const std::vector<unsigned int>& fzvpi = rdorbs->fzvpi();
  const unsigned int nfzc = std::accumulate(fzcpi.begin(), fzcpi.end(), 0.0);
  const unsigned int ninact = std::accumulate(inactpi.begin(), inactpi.end(), 0.0);
  const unsigned int nact = std::accumulate(actpi.begin(), actpi.end(), 0.0);
  const unsigned int nfzv = std::accumulate(fzvpi.begin(), fzvpi.end(), 0.0);
  const unsigned int nmo = orbs->rank();
  const unsigned int nuocc = nmo - nfzc - ninact - nact - nfzv;

  basis = orbs->basis();
  C_ao = orbs->coefs();

  /////////////////////////////////////////////
  // Read 2-RDM
  /////////////////////////////////////////////

  // molcas reports 2-RDM in terms of active occupied orbitals only, indexed occording to molcas convention
  // thus use the map from molcas active occupied orbitals to MPQC occupied range
  // first make an OrbitalSpace for MPQC occupied orbitals
  Ref<OrbitalSpace> occ_orbs = new OrbitalSpace(std::string("z(sym)"),
          std::string("symmetry-ordered occupied MOInfo orbitals"),
          orbs->coefs(),
          orbs->basis(),
          orbs->integral(), orbs->evals(),
          0, nuocc + nfzv,
          OrbitalSpace::symmetry);

  Ref<ExternSpinFreeRDMTwo> rdrdm2;
  rdrdm2 = new ExternSpinFreeRDMTwo(filename_prefix + ".pt2r12.rdm2.dat",
          rdorbs->actindexmap_occ(),
          occ_orbs);
  // create World in which we will compute
  // use defaults for all params
  Ref<WavefunctionWorld> world;
  {
    Ref<AssignedKeyVal> kva = new AssignedKeyVal;
    if (dfbs_name.empty() == false) {
      Ref<AssignedKeyVal> tmpkv = new AssignedKeyVal;
      tmpkv->assign("name", dfbs_name.c_str());
      if (dfbs_name.find("aug-cc-pV") != std::string::npos &&
              dfbs_name.find("Z-RI") != std::string::npos) { // if aug-cc-pVXZ-RI, make one as a union of
        // cc-pVXZ-RI and augmentation-cc-pVXZ-RI
        std::string ccpvxzri_name(dfbs_name, 4, dfbs_name.size()-4);

        Ref<AssignedKeyVal> tmpkv1 = new AssignedKeyVal;
        tmpkv1->assign("name", ccpvxzri_name);
        tmpkv1->assign("molecule", basis->molecule().pointer());

        Ref<GaussianBasisSet> ccpvxzri = new GaussianBasisSet(tmpkv1);

        Ref<AssignedKeyVal> tmpkv2 = new AssignedKeyVal;
        tmpkv2->assign("name", std::string("augmentation-") + ccpvxzri_name);
        tmpkv2->assign("molecule", basis->molecule().pointer());
        Ref<GaussianBasisSet> augmentationccpvxzri = new GaussianBasisSet(tmpkv2);

        Ref<GaussianBasisSet> df_basis = new UnionBasisSet(ccpvxzri, augmentationccpvxzri);
        kva->assign("df_basis", df_basis.pointer());
      }
      else { // otherwise assume the basis exists in the library
        tmpkv->assign("molecule", basis->molecule().pointer());
        Ref<KeyVal> kv = tmpkv;
        Ref<GaussianBasisSet> df_basis = new GaussianBasisSet(kv);
        kva->assign("df_basis", df_basis.pointer());
      }
    }
    Ref<KeyVal> kv = kva;
    world = new WavefunctionWorld(kv);
  }

  {
    Ref<AssignedKeyVal> kva = new AssignedKeyVal;
    kva->assign("orbs_info", rdorbs.pointer());
    kva->assign("rdm2", rdrdm2.pointer());
    kva->assign("world", world.pointer());
    kva->assign("f12exp", f12exp);
    kva->assign("cabs", cabs_name);
    kva->assign("basis", orbs->basis().pointer());
    kva->assign("molecule", orbs->basis()->molecule().pointer());
    if(not r12.empty())
      kva->assign("pt2_correction", r12);
    if(not cabs_contraction.empty())
      kva->assign("cabs_contraction", cabs_contraction);

#if defined(HAVE_MPQC3_RUNTIME)
    if(not cabs_singles.empty())
      kva->assign("cabs_singles", cabs_singles);
    if(not cabs_singles_basis.empty())
      kva->assign("cabs_singles_basis", cabs_singles_basis);
    if(not partition.empty())
      kva->assign("cabs_singles_h0", partition);
#endif
    Ref<KeyVal> kv = kva;
    return kv;
  }
}