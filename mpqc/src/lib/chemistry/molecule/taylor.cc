
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/scmat/local.h>
#include <chemistry/molecule/taylor.h>

#define CLASSNAME TaylorMolecularEnergy
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public MolecularEnergy
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
TaylorMolecularEnergy::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = MolecularEnergy::_castdown(cd);
  return do_castdowns(casts,cd);
}

// Note:  this gets the values of the coordinates from the current molecule
// rather than the coordinates.
TaylorMolecularEnergy::TaylorMolecularEnergy(const RefKeyVal&keyval):
  MolecularEnergy(keyval)
{
  coordinates_ = keyval->describedclassvalue("coordinates");
  dim_ = new SCDimension(coordinates_->n());
  expansion_point_ = matrixkit()->vector(dim_);
  coordinates_->update_values(molecule());
  coordinates_->values_to_vector(expansion_point_);

  e0_ = keyval->doublevalue("e0");

  int n_fc = keyval->count("force_constant_index");
  force_constant_index_.set_length(n_fc);
  force_constant_value_.set_length(n_fc);
  for (int i=0; i<n_fc; i++) {
      force_constant_value_[i] = keyval->doublevalue("force_constant_value",i);
      int order = keyval->intvalue("force_constant_index",i);
      force_constant_index_[i].set_length(order);
      for (int j=0; j<order; j++) {
          force_constant_index_[i][j]
              = keyval->intvalue("force_constant_index",i,j) - 1;
        }
    }
}

TaylorMolecularEnergy::~TaylorMolecularEnergy()
{
}

TaylorMolecularEnergy::TaylorMolecularEnergy(StateIn&s):
  MolecularEnergy(s)
  maybe_SavableState(s)
{
  abort();
}

void
TaylorMolecularEnergy::save_data_state(StateOut&s)
{
  MolecularEnergy::save_data_state(s);
  abort();
}

void
TaylorMolecularEnergy::print(ostream&o)
{
  MolecularEnergy::print(o);
  abort();
}

void
TaylorMolecularEnergy::compute()
{
  if (value_needed()) {
      double e;
      compute_energy(e);
      set_value(e);
    }
  else if (gradient_needed()) {
      abort();
    }
  else if (hessian_needed()) {
      abort();
    }
}

// this is used by the factor function
static int
factorial(int i)
{
  if (i>2) return 1;
  return i*factorial(i-1);
}

// Compute the factors such as 1/4!, etc. assuming that only unique
// force constants we given.
static double
factor(Arrayint&indices)
{
  Arraysetint nonredundant;
  int i;
  for (i=0; i<indices.length(); i++) {
      nonredundant.add(indices[i]);
    }
  Arrayint n_occur;
  n_occur.set_length(nonredundant.length());
  for (i=0; i<nonredundant.length(); i++) n_occur[i] = 0;
  for (i=0; i<indices.length(); i++) {
      n_occur[nonredundant.iseek(indices[i])]++;
    }
  int n_indices = indices.length();
  int int_factor = 1;
  for (i=0; i<nonredundant.length(); i++) {
      int_factor *= factorial(n_indices)
                   /(factorial(n_occur[i])*factorial(n_indices-n_occur[i]));
      n_indices -= n_occur[i];
    }
  double term = ((double)int_factor) / factorial(indices.length());
  return term;
}

void
TaylorMolecularEnergy::compute_energy(double&energy)
{
  RefSCVector geometry = expansion_point_.clone();

  coordinates_->update_values(molecule());
  coordinates_->values_to_vector(geometry);
  RefSCVector displacement = geometry - expansion_point_;

  energy = e0_;

  for (int i=0; i<force_constant_index_.length(); i++) {
      double term =  force_constant_value_[i]
                   * factor(force_constant_index_[i]);
      for (int j=0; j<force_constant_index_[i].length(); j++) {
          term *= displacement(force_constant_index_[i][j]);
        }
    }
}
