
#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/basis/symgaussbas.h>

SavableState_REF_def(SymmGaussianBasisSet);

#define CLASSNAME SymmGaussianBasisSet
#define PARENTS public GaussianBasisSet
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define VERSION 2
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
SymmGaussianBasisSet::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = GaussianBasisSet::_castdown(cd);
  return do_castdowns(casts,cd);
}

SymmGaussianBasisSet::SymmGaussianBasisSet(const GaussianBasisSet& gbs)
  : GaussianBasisSet(gbs), pl(*this)
{
}

SymmGaussianBasisSet::SymmGaussianBasisSet(const RefKeyVal&topkeyval)
  : GaussianBasisSet(topkeyval), pl(*this)
{
}

SymmGaussianBasisSet::SymmGaussianBasisSet(StateIn&s):
  GaussianBasisSet(s), pl(*this)
{
}

SymmGaussianBasisSet::~SymmGaussianBasisSet()
{
}

void
SymmGaussianBasisSet::save_data_state(StateOut&s)
{
  GaussianBasisSet::save_data_state(s);
}

PetiteList&
SymmGaussianBasisSet::petite_list()
{
  return pl;
}

////////////////////////////////////////////////////////////////////////////

SymmetryOrbitals::SymmetryOrbitals() :
  sos_(0), som_(0)
{
}

SymmetryOrbitals::SymmetryOrbitals(const RefSymmGaussianBasisSet& gbs) :
  gbs_(gbs)
{
  sos_ = gbs_->petite_list().aotoso();

  CharacterTable ct = gbs_->molecule()->point_group().char_table();
  
  som_ = new SO_block*[ct.nirrep()];

  SO_block *sp = sos_;
  for (int i=0; i < ct.nirrep(); i++) {
    som_[i] = sp;
    sp += ct.gamma(i).degeneracy();
  }
}

SymmetryOrbitals::~SymmetryOrbitals()
{
  gbs_=0;
  if (sos_) {
    delete[] sos_;
    sos_=0;
  }
  if (som_) {
    delete[] som_;
    som_=0;
  }
}

void
SymmetryOrbitals::print(FILE *out)
{
  char label[80];
  
  fprintf(out,"SymmetryOrbitals:\n");
  CharacterTable ct = gbs_->molecule()->point_group().char_table();

  for (int i=0; i < ct.nirrep(); i++) {
    fprintf(out,"  irrep %s:\n",ct.gamma(i).symbol());
    for (int j=0; j < ct.gamma(i).degeneracy(); j++) {
      sprintf(label,"%d",j+1);
      som_[i][j].print(label);
    }
    fprintf(out,"\n");
  }
}
    
RefBlockedSCDimension
SymmetryOrbitals::AO_basisdim()
{
  RefBlockedSCDimension ret =
    new BlockedSCDimension(gbs_->matrixkit(),gbs_->nbasis());

  return ret;
}

RefBlockedSCDimension
SymmetryOrbitals::SO_basisdim()
{
  return gbs_->petite_list().SO_basisdim();
}
  
////////////////////////////////////////////////////////////////////////////

AOSO_Transformation::AOSO_Transformation(
  const RefSymmGaussianBasisSet& gbs)
  : sos(gbs)
{
  ct = gbs->molecule()->point_group().char_table();
}

void
AOSO_Transformation::process(SCMatrixBlockIter&)
{
  fprintf(stderr,"AOSO_Transformation::process(SCMatrixBlockIter&):"
          " can only handle RectBlocks\n");
  abort();
}

void
AOSO_Transformation::process(SCMatrixRectBlock* blk)
{
  SO_block& sob = sos.sos_[current_block()];

  if (blk->jend > sob.len) {
    fprintf(stderr,"AOSO_Transformation::process(SCMatrixRectBlock*):"
            " blk is wrong size\n");
    abort();
  }
  
  int isize = blk->iend - blk->istart;
  int jsize = blk->jend - blk->jstart;
  
  memset(blk->data,0,isize*jsize*sizeof(double));

  for (int j=blk->jstart; j < blk->jend; j++) {
    contribution *tc = sob.so[j].cont;

    for (int i=0; i < sob.so[j].len; i++,tc++) {
      int bfn = (*tc).bfn;
      if (bfn >= blk->istart && bfn < blk->iend)
        blk->data[(bfn-blk->istart)*jsize+j] = (*tc).coef;
    }
  }
}

////////////////////////////////////////////////////////////////////////////

AOSO_Unit::AOSO_Unit(const RefBlockedSCDimension& rd1,
                     const RefBlockedSCDimension& rd2)
  : d1(rd1), d2(rd2)
{
}

void
AOSO_Unit::process(SCMatrixBlockIter&)
{
  fprintf(stderr,"AOSO_Unit::process(SCMatrixBlockIter&):"
          " can only handle RectBlocks\n");
  abort();
}

void
AOSO_Unit::process(SCMatrixRectBlock* blk)
{
  int fi = (d1->nblocks()==1) ? d1->first(0) : d1->first(current_block());
  int fj = (d2->nblocks()==1) ? d2->first(0) : d2->first(current_block());

  int isize = blk->iend - blk->istart;
  int jsize = blk->jend - blk->jstart;
  
  memset(blk->data,0,isize*jsize*sizeof(double));

  for (int i=blk->istart; i < blk->iend; i++)
    for (int j=blk->jstart; j < blk->jend; j++)
      if (i+fi==j+fj)
        blk->data[(i-blk->istart)*jsize+j-blk->jstart]=1;
}

