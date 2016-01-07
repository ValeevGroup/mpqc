//
// Created by Chong Peng on 10/19/15.
//

#include <catch.hpp>

#include "../expression/orbital_index.h"

TEST_CASE("Orbital Index", "[orbital_index]"){

    SECTION("single letter case"){

        mpqc::OrbitalIndex m(L"m");
        REQUIRE( m.index() == mpqc::OrbitalIndex::Index::occ);
        REQUIRE( m.is_ao() == false);
        REQUIRE( m.is_mo() == true);
        REQUIRE( m.is_mo_in_obs()  == true);
        REQUIRE( m.is_mo_in_abs()  == false);
        REQUIRE( m.is_mo_in_ribs()  == false);

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
        REQUIRE( k1.is_mo() == false);
        REQUIRE( k1.is_ao() == true);

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

        mpqc::OrbitalIndex p1(L"ρ");
        REQUIRE( p1.index() == mpqc::OrbitalIndex::Index::ribs);

        mpqc::OrbitalIndex v2(L"υ");
        REQUIRE( v2.index() == mpqc::OrbitalIndex::Index::ribs);
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

        mpqc::OrbitalIndex a(L"a14");
        REQUIRE( a.index() == mpqc::OrbitalIndex::Index::virt);

        mpqc::OrbitalIndex d(L"d25");
        REQUIRE( d.index() == mpqc::OrbitalIndex::Index::virt);

        mpqc::OrbitalIndex p(L"p36");
        REQUIRE( p.index() == mpqc::OrbitalIndex::Index::any);

        mpqc::OrbitalIndex s(L"s7");
        REQUIRE( s.index() == mpqc::OrbitalIndex::Index::any);

        mpqc::OrbitalIndex k1(L"κ1");
        REQUIRE( k1.index() == mpqc::OrbitalIndex::Index::obs);

        mpqc::OrbitalIndex v1(L"ν52");
        REQUIRE( v1.index() == mpqc::OrbitalIndex::Index::obs);

        mpqc::OrbitalIndex a1(L"α14");
        REQUIRE( a1.index() == mpqc::OrbitalIndex::Index::abs);

        mpqc::OrbitalIndex d1(L"δ5");
        REQUIRE( d1.index() == mpqc::OrbitalIndex::Index::abs);

        mpqc::OrbitalIndex K1(L"Κ1");
        REQUIRE( K1.index() == mpqc::OrbitalIndex::Index::dfbs);

        mpqc::OrbitalIndex V1(L"Ν2");
        REQUIRE( V1.index() == mpqc::OrbitalIndex::Index::dfbs);

        mpqc::OrbitalIndex p1(L"ρ123");
        REQUIRE( p1.index() == mpqc::OrbitalIndex::Index::ribs);

        mpqc::OrbitalIndex v2(L"υ45");
        REQUIRE( v2.index() == mpqc::OrbitalIndex::Index::ribs);
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
        REQUIRE( a.is_mo_in_ribs() == false);
        REQUIRE( a.is_mo_in_abs() == true);
        REQUIRE( a.is_mo_in_obs() == false);

        mpqc::OrbitalIndex d(L"d'");
        REQUIRE( d.index() == mpqc::OrbitalIndex::Index::othervirt);

        mpqc::OrbitalIndex A(L"A'");
        REQUIRE( A.index() == mpqc::OrbitalIndex::Index::allvirt);
        REQUIRE( A.is_mo_in_ribs() == true);
        REQUIRE( A.is_mo_in_abs() == false);
        REQUIRE( A.is_mo_in_obs() == false);

        mpqc::OrbitalIndex D(L"D'");
        REQUIRE( D.index() == mpqc::OrbitalIndex::Index::allvirt);
    }


    SECTION("letter with number and prime"){

        mpqc::OrbitalIndex m(L"m'3");
        REQUIRE( m.index() == mpqc::OrbitalIndex::Index::inactocc);

        mpqc::OrbitalIndex n(L"n'4");
        REQUIRE( n.index() == mpqc::OrbitalIndex::Index::inactocc);

        mpqc::OrbitalIndex P(L"P'1");
        REQUIRE( P.index() == mpqc::OrbitalIndex::Index::allany);

        mpqc::OrbitalIndex S(L"S'2");
        REQUIRE( S.index() == mpqc::OrbitalIndex::Index::allany);

        mpqc::OrbitalIndex a(L"a'4");
        REQUIRE( a.index() == mpqc::OrbitalIndex::Index::othervirt);

        mpqc::OrbitalIndex d(L"d'5");
        REQUIRE( d.index() == mpqc::OrbitalIndex::Index::othervirt);

        mpqc::OrbitalIndex A(L"A'16");
        REQUIRE( A.index() == mpqc::OrbitalIndex::Index::allvirt);

        mpqc::OrbitalIndex D(L"D'178");
        REQUIRE( D.index() == mpqc::OrbitalIndex::Index::allvirt);
    }

    SECTION("error handling"){
        // key not allowed
        REQUIRE_THROWS(mpqc::OrbitalIndex(L"e"));
        REQUIRE_THROWS(mpqc::OrbitalIndex(L"ε"));
        REQUIRE_THROWS(mpqc::OrbitalIndex(L"I"));
        REQUIRE_THROWS(mpqc::OrbitalIndex(L"I'"));
        REQUIRE_THROWS(mpqc::OrbitalIndex(L"I1"));

        // wrong format
        REQUIRE_THROWS(mpqc::OrbitalIndex(L"1a3"));
        REQUIRE_THROWS(mpqc::OrbitalIndex(L"a4'"));
    }

}