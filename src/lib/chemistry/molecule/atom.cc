/*
 * atom.cc
 *
 *  Created on: Jul 25, 2013
 *      Author: drewlewis
 */

#include <chemistry/molecule/atom.h>

void
sc::ToStateOut(const Atom &a, StateOut &so, int &count) {
    count += so.put_array_double(a.r(), 3);
    count += so.put(a.Z());
    count += so.put(a.have_charge());
    count += so.put(a.have_fragment());
    count += so.put(a.charge());
    count += so.put(a.fragment());
    count += so.put(a.mass());
    count += so.put(a.label());
}

void
sc::FromStateIn(Atom &a, StateIn &si, int &count){
    count += si.get_array_double(a.r_,3);
    count += si.get(a.Z_);
    count += si.get(a.have_charge_);
    count += si.get(a.have_fragment_);
    count += si.get(a.charge_);
    count += si.get(a.fragment_);
    count += si.get(a.mass_);
    count += si.get(a.label_);
}

