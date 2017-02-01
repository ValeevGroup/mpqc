
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_ATOM_MASSES_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_ATOM_MASSES_H_

#include <vector>

namespace mpqc {
namespace molecule {
namespace masses {

// Data comes from J. Meija et al Atomic weights of the elements 2013 Pure and
// Applied Chemistry 88, 265-291 (2016) unless otherwise commented.  If there
// was more than one option for the standard atomic weight the first list
// option was used. 
// Make zeroth atom be 0 so we can look up based on atomic number
static const std::vector<double> masses = {
    0.0,          1.00784,     4.002602,    6.938,       9.0121831,   
    10.806,       12.0096,     14.00643,    15.99903,    18.998403163,
    20.1797,      22.98976928, 24.304,      26.9815385,  28.084,      
    30.973761998, 32.059,      35.446,      39.948,      39.0983,     
    40.078,       44.955908,   47.867,      50.9415,     51.9961,     
    54.938044,    55.845,      58.933194,   58.6934,     63.546,      
    65.38,        69.723,      72.630,      74.921595,   78.971,      
    79.901,       83.798,      85.4678,     87.62,       88.90584,    
    91.224,       92.90637,    95.95,                                 
    98.9062547, // TC Is not in the table but this value is reported as the relative atomic mass for Tc99 in Atomic weights of the elements 2011 (IUPAC Technical Report) Pure Appl. Chem., Vol. 85, No. 5, pp. 1047–1078, 2013 
    101.07,       102.90550,   106.42,      107.8682,    112.414,     114.818, 
    118.710,      121.760,     127.60,      126.90447,   131.293,
    132.9054519,  137.327,     138.90547,   140.116,     140.90766, 144.242,      
    146.9151385, // Pm Is from isotope 147, but has more digits than in the IUPAC report from 2013
    150.36,      151.964,     157.25,
    158.92535,    162.500,     164.93033,   167.259,     168.93422,
    173.054,      174.9668,    178.49,      180.94788,   183.84,
    186.207,      190.23,      192.217,     195.084,     196.966569,
    200.592,      204.382,     207.2,       208.98040,   
    209,  // This is just the PO 209 with out any real values
    210.9874963,  // At isotope 211, but with more digits than the IUPAC report from 2013
    222.0175777, // Rn isotope 222, once again same, but with extra digits 
    223, 226, 227, // Just isotope labels in place of masses
    232.0377,     231.03588,   238.02891,   
    237,         244, // Just isotope labels in place of masses
    243.0613811, // I am unsure about the rest of the values, thankfully I don't think these have been used as of 02/01/2017.  
    243.0613891, 249.0749867, 249.0748535, 252.082980,
    257.095105,   258.098431,  259.10103,   262.10963,   265.11670,
    268.12545,    271.13347,   272.13803,   270.13465,   276.15116,
    281.16206,    280.16447,   285.17411,   284.17808,   289.18728,
    288.19249};

}  // namespace masses
}  // namespace molecule
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_ATOM_MASSES_H_
