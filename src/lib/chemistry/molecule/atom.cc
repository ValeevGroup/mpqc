/*
 * atom.cc
 *
 *  Created on: Jul 25, 2013
 *      Author: drewlewis
 */

#include <chemistry/molecule/atom.h>

void
sc::ToStateOut(const Atom &a, StateOut &so, int &count) {
    so.put_array_double(a.r(), 3);
    so.put(a.Z());
    so.put(a.have_charge());
    so.put(a.have_fragment());
    so.put(a.charge());
    so.put(a.fragment());
    so.put(a.mass());
    so.put(a.label());
}

void
sc::FromStateIn(Atom &a, StateIn &si, int &count){
    si.get_array_double(a.r_,3);
    si.get(a.Z_);
    si.get(a.have_charge_);
    si.get(a.have_fragment_);
    si.get(a.charge_);
    si.get(a.fragment_);
    si.get(a.mass_);
    si.get(a.label_);
}

