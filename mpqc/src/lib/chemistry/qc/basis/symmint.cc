
#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/basis/symmint.h>

////////////////////////////////////////////////////////////////////////////
// SymmOneBodyIntIter

SymmOneBodyIntIter::SymmOneBodyIntIter(const RefPetiteList& p) :
  pl(p)
{
}

SymmOneBodyIntIter::~SymmOneBodyIntIter()
{
}

void
SymmOneBodyIntIter::next()
{
  OneBodyIntIter::next();
  while (!pl->lambda(icur,jcur))
    OneBodyIntIter::next();
}

void
SymmOneBodyIntIter::next_ltri()
{
  OneBodyIntIter::next_ltri();
  while (!pl->lambda(icur,jcur))
    OneBodyIntIter::next_ltri();
}

void
SymmOneBodyIntIter::start()
{
  OneBodyIntIter::start();
  while (!pl->lambda(icur,jcur))
    OneBodyIntIter::next();
}

void
SymmOneBodyIntIter::start_ltri()
{
  OneBodyIntIter::start_ltri();
  while (!pl->lambda(icur,jcur))
    OneBodyIntIter::next_ltri();
}

double
SymmOneBodyIntIter::scale() const
{
  return (double) pl->lambda(icur,jcur) / (double) pl->order();
}

////////////////////////////////////////////////////////////////////////////

SymmetryOrbitals::SymmetryOrbitals() :
  sos_(0), som_(0)
{
}

SymmetryOrbitals::SymmetryOrbitals(const RefGaussianBasisSet& gbs,
                                   const RefPetiteList& pl) :
  gbs_(gbs),
  pl_(pl)
{
  sos_ = pl_->aotoso();

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
  return pl_->SO_basisdim();
}
  
////////////////////////////////////////////////////////////////////////////

AOSO_Transformation::AOSO_Transformation(
  const RefGaussianBasisSet& gbs, const RefPetiteList& pl)
  : sos(gbs,pl)
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

