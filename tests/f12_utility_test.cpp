//
// Created by Chong Peng on 10/30/15.
//

#include <catch.hpp>
#include <mpqc/chemistry/qc/f12/f12_utility.h>

using namespace mpqc;

TEST_CASE("F12 Utility", "[f12_utility]") {
  SECTION("basis to exponant") {
    std::string b1 = "cc-pVDZ-F12";
    REQUIRE(f12::basis_to_f12exponent(b1) == 0.9);

    std::string b2 = "cc-pVDZ";
    REQUIRE(f12::basis_to_f12exponent(b2) == 1.2);

    std::string b3 = "  aug-cc-pVQZ";
    REQUIRE(f12::basis_to_f12exponent(b3) == 1.4);
  }

  SECTION("STG_NG_FIT") {
    auto result = f12::stg_ng_fit(6, 0.9);
    auto result00 = std::to_string(result[0].first);
    auto result01 = std::to_string(result[0].second);
    REQUIRE(result01 == "-0.383776");
    REQUIRE(result00 == "0.194846");

    auto result50 = std::to_string(result[5].first);
    auto result51 = std::to_string(result[5].second);
    REQUIRE(result51 == "-0.037943");
    REQUIRE(result50 == "245.372177");
  }
}