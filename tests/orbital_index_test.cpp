//
// Created by Chong Peng on 10/19/15.
//

#include <catch.hpp>

#include "../expression/orbital_index.h"

TEST_CASE("Orbital Index", "[orbital_index]"){

    SECTION("single letter case"){

        mpqc::OrbitalIndex m(L"m");
        REQUIRE( m.index() == mpqc::OrbitalIndex::Index::occ);

        mpqc::OrbitalIndex n(L"n");
        REQUIRE( n.index() == mpqc::OrbitalIndex::Index::occ);

        mpqc::OrbitalIndex i(L"i");
        REQUIRE( i.index() == mpqc::OrbitalIndex::Index::actocc);

        mpqc::OrbitalIndex l(L"l");
        REQUIRE( l.index() == mpqc::OrbitalIndex::Index::actocc);

        mpqc::OrbitalIndex x(L"x");
        REQUIRE( x.index() == mpqc::OrbitalIndex::Index::active);

        mpqc::OrbitalIndex y(L"y");
        REQUIRE( y.index() == mpqc::OrbitalIndex::Index::active);

        mpqc::OrbitalIndex a(L"a");
        REQUIRE( a.index() == mpqc::OrbitalIndex::Index::virt);

        mpqc::OrbitalIndex d(L"d");
        REQUIRE( d.index() == mpqc::OrbitalIndex::Index::virt);

        mpqc::OrbitalIndex p(L"p");
        REQUIRE( p.index() == mpqc::OrbitalIndex::Index::any);

        mpqc::OrbitalIndex s(L"p");
        REQUIRE( s.index() == mpqc::OrbitalIndex::Index::any);

        mpqc::OrbitalIndex k1(L"κ");
        REQUIRE( k1.index() == mpqc::OrbitalIndex::Index::obs);

        mpqc::OrbitalIndex v1(L"ν");
        REQUIRE( v1.index() == mpqc::OrbitalIndex::Index::obs);

        mpqc::OrbitalIndex K1(L"Κ");
        REQUIRE( K1.index() == mpqc::OrbitalIndex::Index::dfbs);

        mpqc::OrbitalIndex V1(L"Ν");
        REQUIRE( V1.index() == mpqc::OrbitalIndex::Index::dfbs);

        mpqc::OrbitalIndex a1(L"α");
        REQUIRE( a1.index() == mpqc::OrbitalIndex::Index::abs);

        mpqc::OrbitalIndex d1(L"δ");
        REQUIRE( d1.index() == mpqc::OrbitalIndex::Index::abs);
    }

    SECTION("letter with number"){

        mpqc::OrbitalIndex m(L"m2");
        REQUIRE( m.index() == mpqc::OrbitalIndex::Index::occ);

        mpqc::OrbitalIndex n(L"n6");
        REQUIRE( n.index() == mpqc::OrbitalIndex::Index::occ);

        mpqc::OrbitalIndex x(L"x7");
        REQUIRE( x.index() == mpqc::OrbitalIndex::Index::active);

        mpqc::OrbitalIndex y(L"y9");
        REQUIRE( y.index() == mpqc::OrbitalIndex::Index::active);

        mpqc::OrbitalIndex i(L"i1");
        REQUIRE( i.index() == mpqc::OrbitalIndex::Index::actocc);

        mpqc::OrbitalIndex l(L"l2");
        REQUIRE( l.index() == mpqc::OrbitalIndex::Index::actocc);

        mpqc::OrbitalIndex a(L"a4");
        REQUIRE( a.index() == mpqc::OrbitalIndex::Index::virt);

        mpqc::OrbitalIndex d(L"d5");
        REQUIRE( d.index() == mpqc::OrbitalIndex::Index::virt);

        mpqc::OrbitalIndex p(L"p6");
        REQUIRE( p.index() == mpqc::OrbitalIndex::Index::any);

        mpqc::OrbitalIndex s(L"s7");
        REQUIRE( s.index() == mpqc::OrbitalIndex::Index::any);

        mpqc::OrbitalIndex k1(L"κ1");
        REQUIRE( k1.index() == mpqc::OrbitalIndex::Index::obs);

        mpqc::OrbitalIndex v1(L"ν2");
        REQUIRE( v1.index() == mpqc::OrbitalIndex::Index::obs);

        mpqc::OrbitalIndex a1(L"α4");
        REQUIRE( a1.index() == mpqc::OrbitalIndex::Index::abs);

        mpqc::OrbitalIndex d1(L"δ5");
        REQUIRE( d1.index() == mpqc::OrbitalIndex::Index::abs);

        mpqc::OrbitalIndex K1(L"Κ1");
        REQUIRE( K1.index() == mpqc::OrbitalIndex::Index::dfbs);

        mpqc::OrbitalIndex V1(L"Ν2");
        REQUIRE( V1.index() == mpqc::OrbitalIndex::Index::dfbs);
    }

    SECTION("letter with prime"){

        mpqc::OrbitalIndex m(L"m'");
        REQUIRE( m.index() == mpqc::OrbitalIndex::Index::inactocc);

        mpqc::OrbitalIndex n(L"n'");
        REQUIRE( n.index() == mpqc::OrbitalIndex::Index::inactocc);

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

        mpqc::OrbitalIndex m(L"m3'");
        REQUIRE( m.index() == mpqc::OrbitalIndex::Index::inactocc);

        mpqc::OrbitalIndex n(L"n4'");
        REQUIRE( n.index() == mpqc::OrbitalIndex::Index::inactocc);

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
        REQUIRE_THROWS(mpqc::OrbitalIndex(L"e"));
        REQUIRE_THROWS(mpqc::OrbitalIndex(L"ε"));
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