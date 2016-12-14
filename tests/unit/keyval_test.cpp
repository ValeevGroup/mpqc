#include <vector>
#include <sstream>

#include "mpqc/util/keyval/keyval.h"
#include "mpqc/util/external/madworld/parallel_file.h"
#include <madness/world/world.h>
#include "catch.hpp"
#include "mpqc/chemistry/qc/basis/basis.h"
#include "mpqc/chemistry/molecule/linkage.h"

using std::vector;
using std::array;
using std::string;
using std::stringstream;
using mpqc::KeyVal;
using mpqc::DescribedClass;

struct Base : public DescribedClass {
  Base(const KeyVal& kv) : DescribedClass(), value_(kv.value<int>("value")) {}
  Base(int v) : value_(v) {}
  virtual ~Base() {}

  int value() const { return value_; }

 private:
  int value_;
};
MPQC_CLASS_EXPORT_KEY(Base);

template <size_t tag>
struct Derived : public Base {
  Derived(const KeyVal& kv) : Base(kv), value_(kv.value<double>("dvalue")) {}
  ~Derived() {}

  double value() const { return value_; }

 private:
  double value_;
};

// only register Derived<0> (you could in principle register Derived generically
// not recommended due to complications with the static data initialization,etc)
MPQC_CLASS_EXPORT_KEY(Derived<0>);

struct Nested : public DescribedClass{

  Nested(const KeyVal& kv) {
    auto base = kv.keyval("base").class_ptr<Base>();
    base_ = base;
  }

  std::shared_ptr<Base> base_;
};

MPQC_CLASS_EXPORT_KEY(Nested)

TEST_CASE("KeyVal", "[keyval]") {
  // first, test basic programmatic construction

  KeyVal kv;

  kv.assign("x", 0);  // equiv JSON: "x":"0"
  REQUIRE(kv.value<int>("x") == 0);

  kv.assign(":x", "0");  // equiv JSON: "": { "x":"0" }
  REQUIRE(kv.value<string>(":x") == "0");
  REQUIRE(kv.value<int>(":x") == 0);

  // can chain multiple assign calls
  kv.assign(":z:0", true).assign(":z:1", -1.75);

  REQUIRE(kv.value<string>(":z:0") == "true");
  REQUIRE_THROWS(kv.value<int>(":z:0"));  // cannot obtain bool as int
  REQUIRE(kv.value<bool>(":z:0") == true);
  // reals can be read in any format, but use the one with enough precision
  REQUIRE(kv.value<float>(":z:1") == -1.75);
  REQUIRE(kv.value<double>(":z:1") == -1.75);

  // can overwrite
  kv.assign(":z:0", false).assign(":z:1", +2.35);  // overwrite old values
  REQUIRE(kv.value<bool>(":z:0") == false);
  REQUIRE(kv.value<double>(":z:1") == +2.35);

  // sequences are written as Arrays, types are lost, hence can write a vector
  // and read as an array
  kv.assign(":z:a:0", vector<int>{{0, 1, 2}});
  typedef array<int, 3> iarray3;
  kv.assign(":z:a:1", iarray3{{1, 2, 3}});
  REQUIRE(kv.value<vector<int>>(":z:a:1") == vector<int>({1, 2, 3}));
  {
    auto ref_array = iarray3{{0, 1, 2}};
    REQUIRE(kv.value<iarray3>(":z:a:0") == ref_array);
  }

  // can write arrays using 0-based keys instead of empty keys in JSON case
  kv.assign(":z:a:2", vector<int>{{7, 6, 5, 4}}, false);
  REQUIRE(kv.value<vector<int>>(":z:a:2") == vector<int>({7, 6, 5, 4}));

  // can count children
  REQUIRE(kv.count(":z") == 3);
  REQUIRE(kv.count(":z:0") == 0);
  REQUIRE(kv.count(":z:1") == 0);
  REQUIRE(kv.count(":z:a") == 3);
  REQUIRE(kv.count(":z:a:0") == 3);
  REQUIRE(kv.count(":z:a:1") == 3);
  REQUIRE(kv.count(":z:a:2") == 4);

  // can write pointers
  {
    std::unique_ptr<double> x(new double);
    kv.assign("double*", x.get());
    *x = 0.0;
    auto x_copy = kv.value<double*>("double*");
    REQUIRE(x_copy == x.get());
  }

  SECTION("JSON read/write") {
    stringstream oss;
    REQUIRE_NOTHROW(kv.write_json(oss));
    // std::cout << oss.str();
  }

  SECTION("JSON read/write") {
    stringstream oss;
    REQUIRE_NOTHROW(kv.write_xml(oss));
    // std::cout << oss.str();
  }

  SECTION("making subtree KeyVal") {
    auto kv_z = kv.keyval(":z");
    REQUIRE(kv_z.value<bool>("0") == false);
    REQUIRE(kv_z.value<double>("1") == +2.35);
  }

  SECTION("make classes") {
    {  // construct Base
      KeyVal kv;
      kv.assign("value", 1).assign("type", "Base");
      auto x1 = kv.class_ptr<Base>();
      REQUIRE(x1->value() == 1);
      auto x2 = kv.class_ptr<Base>();  // this returns the cached ptr
      REQUIRE(x1 == x2);
    }
    {  // construct Derived<0>
      KeyVal kv;
      kv.assign("value", 2).assign("dvalue", 2.0).assign("type", "Derived<0>");
      auto x1 = kv.class_ptr<Derived<0>>();
      REQUIRE(x1->value() == 2.0);
      auto x2 = kv.class_ptr<Derived<0>>();  // this returns the cached ptr
      REQUIRE(x1 == x2);
      auto x3 =
          kv.class_ptr<Base>();  // this returns the cached ptr, cast to Base
      REQUIRE(std::dynamic_pointer_cast<Base>(x1) == x3);
    }
    {  // programmatically add DescribedClasses to KeyVal
      std::shared_ptr<DescribedClass> bptr = std::make_shared<Base>(2);
      KeyVal kv;
      kv.assign("base", bptr);
      auto bptr1 = kv.class_ptr<Base>("base");
      REQUIRE(bptr == bptr1);
    }
  }

  SECTION("check references") {
    KeyVal kv;
    kv.assign("i1", 1);
    kv.assign("i2:a", true);
    kv.assign("i2:b", 2);
    kv.assign("i3", "$:i2:a");
    kv.assign("i4", "$:i2:..:i1");
    kv.assign("i5", "$i2:b");
    kv.assign("i2:c", "$..:i3");
    kv.assign("i 6", 1.33);  // spaces in keys are allowed

    kv.assign("c1:type", "Base").assign("c1:value", 1);
    kv.assign("i2:c2", "$:c1");

    stringstream oss;
    REQUIRE_NOTHROW(kv.write_json(oss));
    //    std::cout << oss.str();

    REQUIRE(kv.value<int>("i4") == 1);
    REQUIRE(kv.value<int>("i5") == 2);
    REQUIRE(kv.value<bool>("i2:c") == true);
    REQUIRE(kv.value<double>("i 6") == 1.33);

    auto kv_c1 = kv.keyval("c1");
    auto c1 = kv_c1.class_ptr<Base>();
    REQUIRE(c1->value() == 1);
    auto c2 = kv.class_ptr<Base>(
        "c1");  // another way, without making keyval for the object
    REQUIRE(c2 == c1);
    auto c3 = kv.class_ptr<Base>("i2:c2");  // returns cached ptr
    REQUIRE(c3 == c1);
  }

  SECTION("read complex JSON") {
    KeyVal kv;

    // clang-format off
    const char input[] =
"{                             \
  \"a\": 0,                    \
  \"b\": 1.25,                 \
  \"base\": {                  \
     \"type\":\"Base\",        \
     \"value\":\"$:a\"         \
  },                           \
  \"nested\" : {               \
      \"type\" : \"Nested\",   \
      \"base\" : \"$:base\"    \
  },                           \
  \"deriv0\": {                \
     \"type\":\"Derived<0>\",  \
     \"value\":\"$:a\",        \
     \"dvalue\":\"$:b\"        \
  },                           \
  \"mpqc\": {                  \
     \"base\":\"$..:base\",    \
     \"deriv\":\"$:deriv0\"    \
  },                           \
  \"a\": 1                     \
}";
    // clang-format on

    stringstream iss(input);
    REQUIRE_NOTHROW(kv.read_json(iss));

    stringstream oss;
    REQUIRE_NOTHROW(kv.write_json(oss));
    // std::cout << oss.str();

    REQUIRE(kv.value<int>("a") == 0);  // "a" specified twice, make sure
                                       // KeyVal::value gets the first
                                       // specification

    auto nested = kv.keyval("nested").class_ptr<Nested>();

    auto b1 = kv.keyval("mpqc:base").class_ptr<Base>();
    auto b2 = kv.keyval("base").class_ptr<Base>();
    REQUIRE(b1 == b2);
    auto b3 = kv.keyval("mpqc:base").class_ptr<Base>();
    REQUIRE(b1 == b3);
    Derived<0> x(kv.keyval("mpqc:deriv"));
    auto d1 = kv.keyval("mpqc:deriv").class_ptr<Derived<0>>();
    auto d2 = kv.keyval("deriv0").class_ptr<Derived<0>>();
    Derived<0>* d3 = &x;
    REQUIRE(d1 == d2);
    REQUIRE(d1.get() != d3);
    REQUIRE(d1->value() == kv.value<double>("b"));

    REQUIRE(kv.exists("c")==false);

    kv.keyval("mpqc:deriv").assign("d",0);
    REQUIRE(kv.keyval("deriv0").exists("d"));
    REQUIRE(kv.value<int>("deriv0:d")==0);
  }


  SECTION("Basis Test"){

    auto& world = TiledArray::get_default_world();
    std::string filename = "keyval_test.json";

    std::stringstream ss;
    mpqc::utility::parallel_read_file(world, filename ,ss);

    KeyVal kv;
    kv.read_json(ss);

    kv.assign("world", &world);

    REQUIRE_NOTHROW(kv.keyval("basis").class_ptr<::mpqc::lcao::gaussian::Basis>());

  }

}  // end of test case
