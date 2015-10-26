//
// Created by Chong Peng on 10/19/15.
//

#include <catch.hpp>

#include "../expression/orbital_index.h"

TEST_CASE("Orbital Index", "[orbital_index]"){

    SECTION("single letter case"){
        mpqc::OrbitalIndex i(L"i");
        REQUIRE( i.index() == mpqc::OrbitalIndex::Index::occ);

        mpqc::OrbitalIndex l(L"i");
        REQUIRE( l.index() == mpqc::OrbitalIndex::Index::occ);

        mpqc::OrbitalIndex a(L"a");
        REQUIRE( a.index() == mpqc::OrbitalIndex::Index::virt);

        mpqc::OrbitalIndex d(L"d");
        REQUIRE( d.index() == mpqc::OrbitalIndex::Index::virt);

        mpqc::OrbitalIndex p(L"p");
        REQUIRE( p.index() == mpqc::OrbitalIndex::Index::any);

        mpqc::OrbitalIndex s(L"p");
        REQUIRE( s.index() == mpqc::OrbitalIndex::Index::any);
    }

    SECTION("letter with number"){

        mpqc::OrbitalIndex i(L"i1");
        REQUIRE( i.index() == mpqc::OrbitalIndex::Index::occ);

        mpqc::OrbitalIndex l(L"l2");
        REQUIRE( l.index() == mpqc::OrbitalIndex::Index::occ);

        mpqc::OrbitalIndex a(L"a4");
        REQUIRE( a.index() == mpqc::OrbitalIndex::Index::virt);

        mpqc::OrbitalIndex d(L"d5");
        REQUIRE( d.index() == mpqc::OrbitalIndex::Index::virt);

        mpqc::OrbitalIndex p(L"p6");
        REQUIRE( p.index() == mpqc::OrbitalIndex::Index::any);

        mpqc::OrbitalIndex s(L"s7");
        REQUIRE( s.index() == mpqc::OrbitalIndex::Index::any);
    }

    SECTION("letter with prime"){
        mpqc::OrbitalIndex P(L"P'");
        REQUIRE( P.index() == mpqc::OrbitalIndex::Index::allany);

        mpqc::OrbitalIndex S(L"S'");
        REQUIRE( S.index() == mpqc::OrbitalIndex::Index::allany);

        mpqc::OrbitalIndex a(L"a'");
        REQUIRE( a.index() == mpqc::OrbitalIndex::Index::othervirt);

        mpqc::OrbitalIndex d(L"d'");
        REQUIRE( d.index() == mpqc::OrbitalIndex::Index::othervirt);

        mpqc::OrbitalIndex A(L"A'");
        REQUIRE( A.index() == mpqc::OrbitalIndex::Index::allvirt);

        mpqc::OrbitalIndex D(L"D'");
        REQUIRE( D.index() == mpqc::OrbitalIndex::Index::allvirt);
    }


    SECTION("letter with number and prime"){
        mpqc::OrbitalIndex P(L"P1'");
        REQUIRE( P.index() == mpqc::OrbitalIndex::Index::allany);

        mpqc::OrbitalIndex S(L"S2'");
        REQUIRE( S.index() == mpqc::OrbitalIndex::Index::allany);

        mpqc::OrbitalIndex a(L"a4'");
        REQUIRE( a.index() == mpqc::OrbitalIndex::Index::othervirt);

        mpqc::OrbitalIndex d(L"d5'");
        REQUIRE( d.index() == mpqc::OrbitalIndex::Index::othervirt);

        mpqc::OrbitalIndex A(L"A6'");
        REQUIRE( A.index() == mpqc::OrbitalIndex::Index::allvirt);

        mpqc::OrbitalIndex D(L"D7'");
        REQUIRE( D.index() == mpqc::OrbitalIndex::Index::allvirt);
    }

    SECTION("error handling"){
        // key not allowed
        REQUIRE_THROWS(mpqc::OrbitalIndex(L"m"));
        REQUIRE_THROWS(mpqc::OrbitalIndex(L"I"));
        REQUIRE_THROWS(mpqc::OrbitalIndex(L"I'"));
        REQUIRE_THROWS(mpqc::OrbitalIndex(L"I1"));

        //wrong length
        REQUIRE_THROWS(mpqc::OrbitalIndex(L"i123"));

        // wrong format
        REQUIRE_THROWS(mpqc::OrbitalIndex(L"a13"));
        REQUIRE_THROWS(mpqc::OrbitalIndex(L"a'4"));
    }

}