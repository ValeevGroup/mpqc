/*
 * atom_class.h
 *
 *  Created on: Jul 22, 2013
 *      Author: drewlewis
 */

#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <util/state/stateio.h>

namespace sc {
    class Atom {

    private:
        double r_[3];
        int Z_;
        bool have_charge_;
        bool have_fragment_;

        // optional prameters.
        double charge_;
        int fragment_;
        double mass_;
        std::string label_;

    public:
        Atom(int Z, double x, double y, double z, char *label = nullptr,
             double mass = 0, int have_charge = 0, double charge = 0,
             int have_fragment = 0, int fragment = 0)
             : Z_(Z), mass_(mass), label_(label ? label : ""), charge_(charge),
               have_charge_(have_charge), have_fragment_(have_fragment),
               fragment_(fragment)
        {
            r_[0] = x;
            r_[1] = y;
            r_[2] = z;
        }

        Atom(int Z, double x, double y, double z, const std::string &label,
             double mass = 0, int have_charge = 0, double charge = 0,
             int have_fragment = 0, int fragment = 0)
             : Z_(Z), mass_(mass), label_(label), charge_(charge),
               have_charge_(have_charge), have_fragment_(have_fragment),
               fragment_(fragment)
        {
            r_[0] = x;
            r_[1] = y;
            r_[2] = z;
        }

        //Don't use this guy.
        Atom() : Z_(-1), mass_(-1), label_(), have_charge_(-1), charge_(-1),
                        have_fragment_(-1), fragment_(-1)
        { r_[0] = -1; r_[1] = -1; r_[2] = -1; }

        double& xyz(int xyz){return r_[xyz];}
        const double& xyz(int xyz) const {return r_[xyz];}

        double* r(){return r_;}
        const double* r() const {return r_;}

        int Z() const {return Z_;}

        double mass() const {return mass_;}
        void set_mass(double mass){mass_ = mass;}

        bool have_charge() const {return have_charge_;}

        double charge() const {return charge_;}
        void set_charge(double charge){charge_ = charge;}

        bool have_fragment() const {return have_fragment_;}

        double fragment() const {return fragment_;}
        void set_fragment(double fragment){fragment_ = fragment;}

        const std::string label() const {return label_;}
        const char* label_c_str() const {return label_.c_str();}

        friend void FromStateIn(Atom &a, StateIn &so, int &count);

   };

    void ToStateOut(const Atom &a, StateOut &so, int &count);

    void FromStateIn(Atom &a, StateIn &so, int &count);
}


#endif /* ATOM_H */
