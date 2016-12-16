[![Build Status](https://travis-ci.com/ValeevGroup/mpqc4.svg?token=2pDpbi3swi4zsJxpapq7&branch=master)](https://travis-ci.com/ValeevGroup/mpqc4)

<img src=https://github.com/ValeevGroup/mpqc4/wiki/images/mpqc_logo_med.png width=115>

### Synopsis
The Massively Parallel Quantum Chemistry (MPQC) platform is a research package for *ab initio* simulation of the electronic structure of molecules and periodic solids. The current (4th) version of the package, or __MPQC4__, is a modern reenvisioning of the conceptual design of the original MPQC platform using the massively-parallel tensor framework [TiledArray](https://github.com/ValeevGroup/tiledarray), distributed task-based programming model and runtime [MADWorld](https://github.com/m-a-d-n-e-s-s/madness), and the Gaussian integrals library [Libint](https://github.com/evaleev/libint).

### Developers
MPQC4 is developed by the [Valeev Group](http://research.valeyev.net) at [Virginia Tech](http://www.vt.edu).

### License

MPQC4 is freely available under the terms of the GPL v3+ licence. See the the included LICENSE file for details. If you are interested in using MPQC4 under different licensing terms, please contact us.

### How to Cite

Cite MPQC4 as
> "MPQC4: Massively Parallel Quantum Chemistry", Edward F. Valeev, Cannada A. Lewis, Chong Peng, Justus A. Calvin, Jinmei Zhang, https://github.com/valeevgroup/mpqc4 .

The development of electronic structure methods in MPQC4 is partially described in the following publications:
* Cannada A. Lewis , Justus A. Calvin , and Edward F. Valeev, "Clustered Low-Rank Tensor Format: Introduction and Application to Fast Construction of Hartree-Fock Exchange.", *J. Chem. Theor. Comp.*, DOI 10.1021/acs.jctc.6b00884;
* Chong Peng, Justus A. Calvin, Fabijan Pavosevic, Jinmei Zhang, and Edward F. Valeev, "Massively Parallel Implementation of Explicitly Correlated Coupled-Cluster Singles and Doubles Using TiledArray Framework.", *J. Phys. Chem. A*, DOI 10.1021/acs.jpca.6b10150.

### Performance

Excellent strong scaling performance of the electronic structure methods in MPQC is demonstrated below for the coupled-cluster singles and doubles (CCSD) wave function solver. Parallel speed-up of 1 iteration of CCSD solver for uracil trimer in 6-31G* AO basis was measured on ["BlueRidge" cluster](https://secure.hosting.vt.edu/www.arc.vt.edu/computing/blueridge-sandy-bridge/) at Virginia Tech (wall time on 1 16-core node = 1290 sec):

![CCSD:UracilTrimer-speedup](https://github.com/ValeevGroup/tiledarray/wiki/images/uracil-trimer-ccsd-blueridge-speedup.png)

This figure was obtained with the help of an allocation from [Advanced Research Computing](https://secure.hosting.vt.edu/www.arc.vt.edu/) at Virginia Tech.

### Acknowledgements
Development of MPQC4 and its key components is made possible by past and present contributions from the National Science Foundation (awards CHE-0847295, CHE-0741927, OCI-1047696, CHE-1362655, ACI-1450262, and ACI-1550456), the Alfred P. Sloan Foundation, the Camille and Henry Dreyfus Foundation, the Department of Energy Exascale Computing Project (NWChemEx subproject), and the Department of Energy INCITE Program.
