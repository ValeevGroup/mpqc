//
// Created by Chong Peng on 7/11/16.
//

#ifndef MPQC_DBMP2F12_H
#define MPQC_DBMP2F12_H

#include <mpqc/chemistry/qc/f12/mp2f12.h>


namespace mpqc{
namespace f12{


template <typename Tile>
class DBMP2F12 : protected MP2F12<Tile>{
public:

  using Policy = TA::SparsePolicy;
  using TArray = TA::DistArray<Tile, Policy>;
  using MolecularIntegralClass = integrals::MolecularIntegral<Tile,Policy>;

  DBMP2F12() = default;

  DBMP2F12(MolecularIntegralClass& mo_int) : MP2F12<Tile>(mo_int) {}

  virtual double compute(const rapidjson::Document& in){
    std::cout << "Not Implemented Yet!" << std::endl;
  }

private:


};

}
}



#endif //MPQC_DBMP2F12_H
