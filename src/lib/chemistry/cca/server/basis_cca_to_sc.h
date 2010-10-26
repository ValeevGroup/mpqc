#ifndef _chemistry_cca_server_basisccatosc_h
#define _chemistry_cca_server_basisccatosc_h

namespace MPQC {
  sc::Ref<sc::GaussianBasisSet>
  basis_cca_to_sc( Chemistry::QC::GaussianBasis::MolecularInterface& );
}

#endif
