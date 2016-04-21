#include <vector>

#include <mpqc/util/keyval/keyval.hpp>
#include <catch.hpp>

using std::cout;
using std::endl;
using std::vector;
using std::string;
using mpqc::KeyVal;
using mpqc::DescribedClass;

struct Base : public DescribedClass {
    Base(const KeyVal& kv) : DescribedClass(), value_(kv.value<int>("value")) {}
    virtual ~Base() {}

    int value_;
    int value() const { return value_; }
};

template <size_t tag>
struct Derived : public Base {
    Derived(const KeyVal& kv) : Base(kv), value_(kv.value<double>("value")) {}
    ~Derived() {}

    double value_;
    double value() const { return value_; }
};

TEST_CASE("KeyVal", "[keyval]"){

  DescribedClass::register_keyval_ctor<Base>("Base");
  DescribedClass::register_keyval_ctor<Derived<0>>("Derived<0>");

  // setup tests programmatic construction

  KeyVal kv;

  kv.assign ("x", 0);
  REQUIRE(kv.value<int>("x") == 0);

  kv.assign (":x", "0");
  REQUIRE(kv.value<string>(":x") == "0");
  REQUIRE(kv.value<int>(":x") == 0);

  kv.assign (":z:0", true).assign (":z:1", -1.75);

  REQUIRE(kv.value<string>(":z:0") == "true");
  REQUIRE_THROWS(kv.value<int>(":z:0")); // cannot obtain bool as int
  REQUIRE(kv.value<bool>(":z:0") == true);
  REQUIRE(kv.value<float>(":z:1") == -1.75);
  REQUIRE(kv.value<double>(":z:1") == -1.75);

  kv.assign (":z:a:0", vector<int> ( {0, 1, 2}));
  //REQUIRE_THROWS(kv.value<vector<int>>(":z:a:0")); // not yet implemented

  SECTION("JSON read/write"){
    kv.write_json("keyvaltest0.json");
  }

  SECTION("JSON read/write"){
    kv.write_xml("keyvaltest0.xml");
  }

  SECTION("making subtree KeyVal"){
    auto kv_z = kv.keyval (":z");
    kv_z.write_json ("keyvaltest1.json");
    REQUIRE(kv_z.value<bool> ("0") == true);
    REQUIRE(kv_z.value<double> ("1") == -1.75);
  }

  SECTION("make classes"){

    {
      KeyVal kv;
      kv.assign("value", 1).assign("type", "Base");
      auto x1 = kv.class_ptr<Base>();
      REQUIRE(x1->value() == 1);
      auto x2 = kv.class_ptr<Base>(); // this returns the cached ptr
      REQUIRE(x1 == x2);
    }
    {
      KeyVal kv;
      kv.assign("value", 2.0).assign("type", "Derived<0>");
      auto x1 = kv.class_ptr<Derived<0>>();
      REQUIRE(x1->value() == 2.0);
      auto x2 = kv.class_ptr<Derived<0>>(); // this returns the cached ptr
      REQUIRE(x1 == x2);
      auto x3 = kv.class_ptr<Base>(); // this returns the cached ptr, cast to Base
      REQUIRE(std::dynamic_pointer_cast<Base>(x1) == x3);
    }
  }

  SECTION("check references"){
    KeyVal kv;
    kv.assign ("i1", 1);
    kv.assign ("i2:a", true);
    kv.assign ("i2:b", 2);
    kv.assign ("i3", "$:i2:a");
    kv.assign ("i4", "$:i2:..:i1");
    kv.assign ("i5", "$i2:b");
    kv.assign ("i2:c", "$..:i3");

    kv.assign("c1:type", "Base").assign("c1:value", 1);
    kv.assign ("i2:c2", "$:c1");

    kv.write_json ("keyvaltest2_orig.json");

    REQUIRE(kv.value<int>("i4") == 1);
    REQUIRE(kv.value<int>("i5") == 2);
    REQUIRE(kv.value<bool>("i2:c") == true);

    auto kv_c1 = kv.keyval("c1");
    auto c1 = kv_c1.class_ptr<Base>();
    REQUIRE(c1->value() == 1);
    auto kv_c2 = kv.keyval("i2:c2");
    auto c2 = kv_c2.class_ptr<Base>(); // returns cached ptr
    REQUIRE(c1 == c2);
  }

} // end of test case

