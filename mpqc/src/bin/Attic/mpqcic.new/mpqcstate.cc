
/* mpqcstate.cc -- implementation of mpqc state class
 *
 *      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
 *      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
 *      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
 *      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
 *      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
 *      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
 *
 *  Author:
 *      E. T. Seidl
 *      Bldg. 12A, Rm. 2033
 *      Computer Systems Laboratory
 *      Division of Computer Research and Technology
 *      National Institutes of Health
 *      Bethesda, Maryland 20892
 *      Internet: seidl@alw.nih.gov
 *      November, 1993
 */

#include <stdio.h>
#include <math.h>

extern "C" {
#include <util/misc/libmisc.h>
#include <comm/picl/picl.h>
#include <comm/picl/ext/piclext.h>
}

#include "mpqcstate.h"

#include <chemistry/molecule/simple.h>
#include <chemistry/molecule/symm.h>
#include <chemistry/molecule/simpleQCList.h>
#include <chemistry/molecule/symmQCList.h>
#include <math/nihmatrix/lmath.h>

#define CLASSNAME MPQC_state
#define PARENTS virtual public SavableState
#define HAVE_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
MPQC_state::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

MPQC_state::MPQC_state() :
  coords_(0), forces_(0), error_(0), nsave(0), state_(UNK)
{
  wfn_converged_=geom_converged_=0;
}
  
MPQC_state::MPQC_state(StateIn& si) :
  coords_(0), forces_(0), error_(0)
{
  si.get((int)state_); si.get(wfn_converged_); si.get(geom_converged_);
  hessian_ = DMatrix(si);
  coord_ = DVector(si);
  deltax_ = DVector(si);
  force_ = DVector(si);
  mol_ = Molecule(si);

  si.get(nsave);
  if (!nsave) return;

  coords_ = new DVector[nsave];
  forces_ = new DVector[nsave];
  error_ = new DVector[nsave];

  for (int i=0; i < nsave; i++) {
    coords_[i] = DVector(si);
    forces_[i] = DVector(si);
    error_[i] = DVector(si);
    }
}

MPQC_state::MPQC_state(const char *input, MPQC_input& mpqcin,
                       MPQC_geom_input& geomin, centers_t& centers) :
  coords_(0), forces_(0), error_(0)
{
  if (mynode0()==0) {
    if (!input)
      input = "mpqc.in";

    ParsedKeyVal rawin(input);
    PrefixKeyVal pref1(":mpqc",rawin);
    PrefixKeyVal pref2(":default",rawin);
    AggregateKeyVal keyval(pref1,pref2);

    read_input(keyval,mpqcin,geomin,centers);
    }
}

MPQC_state::MPQC_state(KeyVal& keyval, MPQC_input& mpqcin,
                       MPQC_geom_input& geomin, centers_t& centers) :
  coords_(0), forces_(0), error_(0)
{
  if (mynode0()==0)
    read_input(keyval,mpqcin,geomin,centers);
}

void
MPQC_state::save_data_state(StateOut& so)
{
  so.put(state_); so.put(wfn_converged_); so.put(geom_converged_);
  hessian_.save_object_state(so);
  coord_.save_object_state(so);
  deltax_.save_object_state(so);
  force_.save_object_state(so);
  mol_.save_object_state(so);

  so.put(nsave);
  for (int i=0; i < nsave; i++) {
    coords_[i].save_object_state(so);
    forces_[i].save_object_state(so);
    error_[i].save_object_state(so);
    }
}


// this is to be used when the scf wfn is not the one specified by the
// "wavefunction" keyword, eg. for an mp2
void 
MPQC_state::did_scf(MPQC_input& mpqcin)
{

  if (mpqcin.wavefunction() == MPQC_input::WFN_SCF) {
    did_wfn(mpqcin);
    return;
    }
  
  state_ = NEED_WFN;
}

void 
MPQC_state::did_wfn(MPQC_input& mpqcin)
{
  if (mpqcin.dertype()==1) 
    state_ = NEED_GRADIENT;
  else if (mpqcin.dertype()==2)
    state_ = NEED_HESSIAN;
  else if (mpqcin.optimize_geometry() && !geom_converged())
    state_ = NEED_NEXT_GEOM;
  else if (mpqcin.do_properties())
    state_ = NEED_PROPER;
  else 
    state_ = DONE;
}

void 
MPQC_state::did_gradient(MPQC_input& mpqcin)
{
  if (mpqcin.dertype()==2)
    state_ = NEED_HESSIAN;
  else if (mpqcin.optimize_geometry() && !geom_converged())
    state_ = NEED_NEXT_GEOM;
  else if (mpqcin.do_properties())
    state_ = NEED_PROPER;
  else 
    state_ = DONE;
}

void 
MPQC_state::did_hessian(MPQC_input& mpqcin)
{
  if (mpqcin.optimize_geometry() && !geom_converged())
    state_ = NEED_NEXT_GEOM;
  else if (mpqcin.do_properties())
    state_ = NEED_PROPER;
  else 
    state_ = DONE;
}

void 
MPQC_state::did_next_geom(MPQC_input& mpqcin)
{
  if (mpqcin.wavefunction() != MPQC_input::WFN_SCF)
    state_ = NEED_SCF_WFN;
  else
    state_ = NEED_WFN;
}

void 
MPQC_state::did_proper(MPQC_input& mpqcin)
{
  state_ = DONE;
}

void 
MPQC_state::broadcast()
{
  bcast0(&state_,sizeof(int),mtype_get(),0);
  bcast0(&wfn_converged_,sizeof(int),mtype_get(),0);
  bcast0(&geom_converged_,sizeof(int),mtype_get(),0);
}

centers_t * MPQC_state::centers;

void
MPQC_state::read_input(KeyVal& keyval, MPQC_input& mpqcin,
                       MPQC_geom_input& geomin, centers_t& cs)
{
  int i;

  MPQC_state::centers = &cs;

 // ok, first let's create the molecule object using the coords in centers
  for (i=0; i < centers->n; i++) {
    char* label = keyval.pcharvalue("atom_labels",i);
    AtomicCenter ac(centers->center[i].atom,
                    centers->center[i].r[0],
                    centers->center[i].r[1],
                    centers->center[i].r[2],
                    label
                    );
    delete[] label;
    mol_.add_atom(i,ac);
    }

  if (mpqcin.optimize_geometry()) {

   // now that we have a geometry, let's make the coordinates
    int ncoord=0;
    if(!geomin.cartesians()) {
      if(get_symmco(keyval,mpqcin,geomin) < 0) {
        err_msg("MPQC_state::read_input: yikes!  can't make symm coords");
        exit(1);
        }

      for(SymmCoListIter scp=symm_coords; scp; scp++,ncoord++) ;
      }
    else
      ncoord = mol_.natom()*3;

    coord_.resize(ncoord);
    coord_.zero();

    force_.resize(ncoord);
    force_.zero();

    deltax_.resize(ncoord);
    deltax_.zero();

    if (!geomin.cartesians()) {
      i=0;
      for(SymmCoListIter scp=symm_coords; scp; scp++,i++) {
        coord_[i] = scp->value();
        }
      }
    }
}

int
MPQC_state::get_symmco(KeyVal& keyval, MPQC_input& mpqcin,
                       MPQC_geom_input& geomin, centers_t& cs)
{
  SimpleCoList *list=0;

 // if the simples are defined in the input, read them, otherwise generate
 // them based on the geometry

  if(keyval.count("simp"))
    list = Geom_read_simples(&keyval);
  else
    list = Geom_form_simples(mol_);

  int nadd;
  if(nadd=keyval.count("add_simp")) {
    for(int i=0; i < nadd; i++) {
      char *val = keyval.pcharvalue("add_simp",i,0);

      if (!strcmp("stre",val))
        list->add(new Stre(&keyval,"add_simp",i));
      else if (!strcmp("bend",val))
        list->add(new Bend(&keyval,"add_simp",i));
      else if (!strcmp("tors",val))
        list->add(new Tors(&keyval,"add_simp",i));
      else if (!strcmp("out",val))
        list->add(new Out(&keyval,"add_simp",i));
      else if (!strcmp("linip",val))
        list->add(new LinIP(&keyval,"add_simp",i));
      else if (!strcmp("linop",val))
        list->add(new LinOP(&keyval,"add_simp",i));
      delete[] val;
      }
    }

  if(keyval.count("symm")) {
    symm_coords = Geom_read_symm(&keyval,"symm",list);
    if(!redundant) symm_coords = Geom_form_symm(mol,symm_coords,justa1);
    }
  else if(redundant)
    symm_coords = Geom_symm_from_simple(list,1);
  else
    symm_coords = Geom_form_symm(mol,list,justa1);

  if(!symm_coords) {
    err_msg("get_symmco: yikes!  no symmetrized internal coords");
    return -1;
    }

  guess_hessian(list);

  if(list) {
    printf("  Simple Internal Coordinates\n");
    Geom_calc_simples(list,mol);
    Geom_print_pretty(list);
    delete list;
    }

  printf("\n  in get_symmco:  %d symmetrized internal coordinates\n",
                hessian.nrow());

  return 0;
}
