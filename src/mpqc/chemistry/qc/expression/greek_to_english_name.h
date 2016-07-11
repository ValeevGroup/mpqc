//
// Created by Chong Peng on 1/25/16.
//

#ifndef MPQC_GREEK_TO_ENGLISH_NAME_H
#define MPQC_GREEK_TO_ENGLISH_NAME_H

#include <string>
#include <unordered_map>

namespace mpqc{
    const std::unordered_map<wchar_t, std::string> greek_to_english_name = {
            {L'Α',"ALPHA"},
            {L'Β',"BETA"},
            {L'Γ', "GAMMA"},
            {L'Δ', "DELTA"},
            {L'Ε', "EPSILON"},
            {L'Ζ', "ZETA"},
            {L'Η', "ETA"},
            {L'Θ',"THETA"},
            {L'Ι', "IOTA"},
            {L'Κ', "KAPPA"},
            {L'Λ', "LAMDA"},
            {L'Μ', "MU"},
            {L'Ν', "NU"},
            {L'Ξ',"XI"},
            {L'Ο', "OMIRCON"},
            {L'Π', "PI"},
            {L'Ρ', "RHO"},
            {L'Σ', "SIGMA"},
            {L'Τ', "TAU"},
            {L'Υ', "UPSILON"},
            {L'Φ', "PHI"},
            {L'Χ', "CHI"},
            {L'Ψ', "PSI"},
            {L'Ω', "OMEGA"},
            {L'α',"alpha"},
            {L'β',"beta"},
            {L'γ', "gamma"},
            {L'δ', "delta"},
            {L'ε', "epsilon"},
            {L'ζ', "zeta"},
            {L'η', "eta"},
            {L'θ',"theta"},
            {L'ι', "iota"},
            {L'κ', "kappa"},
            {L'λ', "lamda"},
            {L'μ', "mu"},
            {L'ν', "nu"},
            {L'ξ',"xi"},
            {L'ο', "omicron"},
            {L'π', "pi"},
            {L'ρ', "rho"},
            {L'σ', "sigma"},
            {L'τ', "tau"},
            {L'υ', "upsilon"},
            {L'φ', "phi"},
            {L'χ', "chi"},
            {L'ψ', "psi"},
            {L'ω', "omega"}
    };
}

#endif //MPQC_GREEK_TO_ENGLISH_NAME_H
