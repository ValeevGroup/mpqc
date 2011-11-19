
#ifdef HAVE_CONFIG_H
#include <scconfig.h>
#endif

#include <fstream>

#include <util/keyval/keyval.h>
#include <math/isosurf/shape.h>
#include <chemistry/qc/dft/solvent.h>
#include <chemistry/molecule/formula.h>

#ifdef USING_NAMESPACE_STD
using namespace std;
#endif
using namespace sc;

static inline double
get_ki(int z)
{
  // The ki values (used in the computation of the dispersion coefficients)
  // for H, C, and N were taken from Vigne-Maeder and Claverie, JACS 1987, v109, pp24-28
  // and the value for O from Huron and Claverie, J. Phys. Chem. 1974, v78, p1862

  double ki;
  
  if (z <= 0) {
      ExEnv::errn() << "Non-positive nuclear charge encountered in computation of"
           << " dispersion coefficient" << endl;
      abort();
    }
  else if (z == 1) ki = 1.0;
  else if (z == 6) ki = 1.0;
  else if (z == 7) ki = 1.18;
  else if (z == 8) ki = 1.36;  // from Huron & Claverie, J.Phys.Chem v78, 1974, p1862
  else if (z > 1 && z < 6) {
      ki = 1.0;
      ExEnv::out0() << "Warning: No d6 dispersion coefficient available for atomic number " <<
              z << "; using value for carbon instead" << endl;
    }
  else {
      ki = 1.18;
      ExEnv::out0() << "Warning: No d6 dispersion coefficient available for atomic number " <<
              z << "; using value for nitrogen instead" << endl;
    }
  
  return ki;
}

static inline double
get_d6ii(int z, double r_vdw)
{
  // The dispersion coefficient d6 for a pair of atoms ij can be computed
  // from the dispersion coefficient d6ii for atom pair ii and d6jj for
  // atom pair jj by the formula: d6 = sqrt(d6ii*d6jj).
  // The dispersion coefficients d8 and d10 can be obtained from d6.
  // The d6ii values given below were taken from: Vigne-Maeder and Claverie
  // JACS 1987, v. 109, pp. 24-28.

  const double a6 = 0.143; // [kcal/mol]
  double d6ii;
  double ki;

  Ref<Units> unit = new Units("kcal/mol");
  
  ki = get_ki(z);
  d6ii = ki*ki*a6*pow(4*r_vdw*r_vdw,3.0);  // units of (kcal mol^-1)*bohr^6
  d6ii *= unit->to_atomic_units();  // convert to atomic units
  return d6ii;
}

static inline double
get_d8ii(double d6ii, double r_vdw)
{
  // The value of c8 was taken from Vigne-Maeder and Claverie, JACS 1987,
  // v. 109, pp 24-28 and is here obtained in atomic units by using
  // atomic units for d6ii and r_vdw

  double d8ii;
  const double c8 = 0.26626;

  d8ii = d6ii*c8*4*pow(r_vdw,2.0);
  
  return d8ii;
}

static inline double
get_d10ii(double d6ii, double r_vdw)
{
  // The value of c10 was taken from Vigne-Maeder and Claverie, JACS 1987,
  // v. 109, pp 24-28 and is here obtained in atomic units by using
  // atomic units for d6ii and r_vdw

  double d10ii;
  const double c10 = 0.095467;

  d10ii = d6ii*c10*16*pow(r_vdw,4.0);
  
  return d10ii;
}

// For debugging compute 6, 8, and 10 contributions separately
static inline double
disp6_contrib(double rasnorm, double d6)
{
  double edisp6_contrib;

  edisp6_contrib = d6/(3*pow(rasnorm,6.0)); // atomic units
  
  return edisp6_contrib;
}

static inline double
disp8_contrib(double rasnorm, double d8)
{
  double edisp8_contrib;

  edisp8_contrib = d8/(5*pow(rasnorm,8.0)); // atomic units
  
  return edisp8_contrib;
}

static inline double
disp10_contrib(double rasnorm, double d10)
{
  double edisp10_contrib;

  edisp10_contrib = d10/(7*pow(rasnorm,10.0)); // atomic units
  
  return edisp10_contrib;
}

static inline double
disp_contrib(double rasnorm, double d6, double d8, double d10)
{
  double edisp_contrib;

  edisp_contrib = d6/(3*pow(rasnorm,6.0)) + d8/(5*pow(rasnorm,8.0))
                + d10/(7*pow(rasnorm,10.0));
  
  return edisp_contrib;
}

static inline double
rep_contrib(double rasnorm, double ri_vdw, double rj_vdw, double ki, double kj,
            double kcalpermol_to_hartree)
{
  // The expression and the parameters used for the repulsion energy
  // were taken from Vigne-Maeder and Claverie, JACS 1987, v109, pp24-28
  // NB: We have omitted the factor Gij
  
  const double c = 90000; // [kcal/mol]
  const double gamma = 12.35;
  double erep_contrib;
  double tmp;
  
  tmp = gamma*rasnorm/(2.0*sqrt(ri_vdw*rj_vdw));

  erep_contrib = -ki*kj*c*(1.0/tmp + 2.0/(tmp*tmp) + 2.0/(tmp*tmp*tmp))*exp(-tmp);
  erep_contrib *= kcalpermol_to_hartree; // convert from kcal/mol to atomic units
  
  return erep_contrib;
}

double
BEMSolvent::disprep() 
{
  double edisprep = 0.0;
  double edisprep_contrib;
  double edisp6_contrib, edisp8_contrib, edisp10_contrib; // for debugging
  double erep_contrib;
  double edisp6 = 0.0; // for debugging
  double edisp8 = 0.0; // for debugging
  double edisp10 = 0.0; // for debugging
  double erep = 0.0;
  double proberadius;
  double radius;
  double rasnorm;
  double weight;
  double d6, d8, d10; // dispersion coefficients
  double d6aa, d8aa, d10aa; // dispersion coefficients for atom pair aa
  double d6ss, d8ss, d10ss; // dispersion coefficients for atom pair ss
  int i, iloop, isolute;
  int natomtypes;
  int z_solvent_atom;

  Ref<Units> unit = new Units("kcal/mol");
  double kcalpermol_to_hartree = unit->to_atomic_units();

  Ref<AtomInfo> atominfo = solute_->atominfo();
  Ref<AtomInfo> solventatominfo = solvent_->atominfo();
  MolecularFormula formula(solvent_);

  // Compute number of different atom types in solvent molecule
  natomtypes = formula.natomtypes();

  double *solute_d6ii  = new double[solute_->natom()];
  double *solute_d8ii  = new double[solute_->natom()];
  double *solute_d10ii = new double[solute_->natom()];
  double *solute_ki = new double[solute_->natom()];

  for (isolute=0; isolute<solute_->natom(); isolute++) {
      int Z_solute = solute_->Z(isolute);
      double radius = atominfo->vdw_radius(Z_solute);
      solute_d6ii[isolute] = get_d6ii(Z_solute,radius);
      solute_d8ii[isolute] = get_d8ii(solute_d6ii[isolute],radius);
      solute_d10ii[isolute] = get_d10ii(solute_d6ii[isolute],radius);
      solute_ki[isolute] = get_ki(Z_solute);
    }
  
  // Loop over atom types in solvent molecule
  for (iloop=0; iloop<natomtypes; iloop++) {

      // define the shape of the surface for current atom type
      Ref<UnionShape> us = new UnionShape;
      z_solvent_atom = formula.Z(iloop);
      proberadius = solventatominfo->vdw_radius(z_solvent_atom);
      for (i=0; i<solute_->natom(); i++) {
          us->add_shape(new SphereShape(solute_->r(i),
                                        atominfo->vdw_radius(solute_->Z(i))+proberadius));
        }
      
      // triangulate the surface
      Ref<AssignedKeyVal> keyval = new AssignedKeyVal;
      keyval->assign("volume", us.pointer());
      keyval->assign("order", 2);
      keyval->assign("remove_short_edges", 1);
      keyval->assign("remove_small_triangles", 1);
      keyval->assign("remove_slender_triangles", 1);
      keyval->assign("short_edge_factor", 0.8);
      keyval->assign("small_triangle_factor", 0.8);
      keyval->assign("slender_triangle_factor", 0.8);
      Ref<TriangulatedImplicitSurface> ts = new TriangulatedImplicitSurface(keyval.pointer());
      ts->init();

      // Debug print: check the triangulated surface
//      if (iloop == 0) {
//          ofstream geomviewfile("geomview.input");
//          ts->print_geomview_format(geomviewfile);
//        }
      
      ExEnv::out0().setf(ios::scientific,ios::floatfield); // use scientific format
      ExEnv::out0() << "Area of disp-rep surface generated with atom number "
           << setw(3) << setfill(' ') << z_solvent_atom
           << " as probe: " << setprecision(4) << ts->area()
           << " bohr^2" << endl;
      
      edisprep_contrib = 0.0;
      edisp6_contrib = 0.0;  // for debugging
      edisp8_contrib = 0.0;  // for debugging
      edisp10_contrib = 0.0; // for debugging
      erep_contrib = 0.0;
      TriangulatedSurfaceIntegrator triint(ts.pointer());

      double solvent_ki = get_ki(z_solvent_atom);
      d6ss = get_d6ii(z_solvent_atom,proberadius);
      d8ss = get_d8ii(d6ss, proberadius);
      d10ss = get_d10ii(d6ss, proberadius);
              
      // integrate the surface
      for (triint=0; triint.update(); triint++) {
          SCVector3 dA = triint.dA();
          SCVector3 location = triint.current()->point();
          weight = triint.weight();
          
          //Loop over atoms in solute
          for (isolute=0; isolute<solute_->natom(); isolute++) {

              SCVector3 atom(solute_->r(isolute)); 
              SCVector3 ras = location - atom;
              rasnorm = ras.norm();
              radius = atominfo->vdw_radius(solute_->Z(isolute));
              d6aa = solute_d6ii[isolute];
              d8aa = solute_d8ii[isolute];
              d10aa = solute_d10ii[isolute];
              d6 = sqrt(d6aa*d6ss);
              d8 = sqrt(d8aa*d8ss);
              d10 = sqrt(d10aa*d10ss);

              double f = ras.dot(dA)*weight;
              double tdisp6 = f*disp6_contrib(rasnorm,d6);
              double tdisp8 = f*disp8_contrib(rasnorm,d8);
              double tdisp10 = f*disp10_contrib(rasnorm,d10);
              double trep = f*rep_contrib(rasnorm,radius,proberadius,
                                          solute_ki[isolute],solvent_ki,
                                          kcalpermol_to_hartree);
              double tdisp = tdisp6+tdisp8+tdisp10;

              // add in contributions to various energies; the minus sign
              // is there to get the normal pointing into the cavity
              edisprep_contrib -= tdisp+trep;
              edisp6_contrib -= tdisp6;
              edisp8_contrib -= tdisp8;
              edisp10_contrib -= tdisp10;
              erep_contrib -= trep;
              
            }
        }
      
      edisprep += edisprep_contrib*formula.nZ(iloop);
      edisp6 += edisp6_contrib*formula.nZ(iloop);
      edisp8 += edisp8_contrib*formula.nZ(iloop);
      edisp10 += edisp10_contrib*formula.nZ(iloop);
      erep += erep_contrib*formula.nZ(iloop);
    }

  delete[] solute_d6ii;
  delete[] solute_d8ii;
  delete[] solute_d10ii;
  delete[] solute_ki;

  // Multiply energies by number density of solvent
  // Print out individual energy contributions in kcal/mol
  
  ExEnv::out0().setf(ios::scientific,ios::floatfield); // use scientific format
  ExEnv::out0().precision(5);
  ExEnv::out0() << "Edisp6:  " << edisp6*solvent_density_*unit->from_atomic_units()
       << " kcal/mol" << endl; 
  ExEnv::out0() << "Edisp8:  " << edisp8*solvent_density_*unit->from_atomic_units()
       << " kcal/mol" << endl;
  ExEnv::out0() << "Edisp10: " << edisp10*solvent_density_*unit->from_atomic_units()
       << " kcal/mol" << endl;


  ExEnv::out0() << "Total dispersion energy: "
       << (edisp6 + edisp8 + edisp10)*solvent_density_*unit->from_atomic_units()
       << " kcal/mol" << endl;
  ExEnv::out0() << "Repulsion energy:        " << setw(12) << setfill(' ') 
       << erep*solvent_density_*unit->from_atomic_units() << " kcal/mol" << endl;
  
  return edisprep*solvent_density_; // atomic units

}


