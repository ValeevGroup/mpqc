/*
 * atom_class.h
 *
 *  Created on: Jul 22, 2013
 *      Author: drewlewis
 */

#ifndef ATOM_CLASS_H_
#define ATOM_CLASS_H_

#include <string>

namespace sc {
    class Atom {

    private:
        double r_[3];
        int Z_;
        int have_charge_;
        int have_fragment_;

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

        Atom();

        double& xyz(int xyz){return r_[xyz];}
        const double& xyz(int xyz) const {return r_[xyz];}

        double* r(){return r_;}
        const double* r() const {return r_;}

        int Z() const {return Z_;}

        double mass() const {return mass_;}
        void set_mass(double mass){mass_ = mass;}

        int have_charge() const {return have_charge_;}

        double charge() const {return charge_;}
        void set_charge(double charge){charge_ = charge;}

        int have_fragment() const {return have_fragment_;}

        double fragment() const {return fragment_;}
        void set_fragment(double fragment){fragment_ = fragment;}

        const char* label() const {return label_.c_str();}
    };
}


#endif /* ATOM_CLASS_H_ */
