
import numpy
from pympqc import *
from ase.units import Bohr, Hartree
from ase.visualize import view

def atoms_to_molecule(atoms):
    # ASE uses Angstrom throughout, but by default Molecule uses Bohr; apply conversion factors below
    m = Molecule()

    zs = atoms.get_atomic_numbers()
    ps = atoms.get_positions()

    # label='X' is needed to avoid segfault -- passing arguments by reference is not yet implemented properly
    for i in range(len(zs)):
        m.add_atom(int(zs[i]), ps[i,0] * Bohr, ps[i,1] * Bohr, ps[i,2] * Bohr, label='X');

    return m

class Calculator:
    """This is the ASE-calculator frontend for doing an MPQC calculation.
    """

    def __init__(self, mole=None):
        self.mole = mole
        # MPQC cannot compute stress
        self.stress = numpy.empty((3, 3))

    def set_x(self, atoms):
        # ASE uses angstroms throughout, hence convert to the internal units of Molecule
        m = self.mole.molecule()
        from_bohr = m.geometry_units().from_atomic_units()
        current_x = self.mole.get_x()
        positions = atoms.get_positions()
        for i in range(len(positions)):
            for j in range(3):
                current_x.set_element(i*3+j,(positions[i,j]/Bohr)* from_bohr)
        self.mole.set_x(current_x)

    def get_potential_energy(self, atoms=None, force_consistent=False):
        self.set_x(atoms)
        # ASE uses eV throughout, hence convert from Hartree
        return self.mole.energy() * Hartree
        
    def get_forces(self, atoms):
        self.set_x(atoms)
        positions = atoms.get_positions()
        natom = len(positions)
        grad = self.mole.gradient()
        ase_grad = numpy.empty((natom,3),dtype=float)
        for i in range(natom):
            for j in range(3):
                # ASE uses eV and Angstrom throughout, hence convert from Hartree /Bohr
                ase_grad[i,j] = -grad.get_element(i*3+j) * Hartree / Bohr
        return ase_grad

    def get_stress(self, atoms):
        return self.stress

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    
    return None

def ase_gui_view(atoms):
    if which("ase-gui") == None:
        print "ase-gui is not found in PATH, will skip visualization"
    else:
        view(atoms)
