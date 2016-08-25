//
// Created by Chong Peng on 8/25/16.
//

#include <catch.hpp>
#include <mpqc/util/misc/string.h>

TEST_CASE("Util String", "[util_string]") {
  SECTION("to_string") {
    std::wstring wstr = L"αβ12";
    std::string str = mpqc::utility::to_string(wstr);

    REQUIRE(str == "αβ12");
  }

  SECTION("to_wstring") {
    std::string str = "αβ12";
    std::wstring wstr = mpqc::utility::to_wstring(str);
    REQUIRE(wstr == L"αβ12");

    wstr = mpqc::utility::to_wstring(str.c_str());
    REQUIRE(wstr == L"αβ12");
  }
}
