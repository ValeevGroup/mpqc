
#include <math/scmat/elemop.h>
#include <math/scmat/blocked.h>
#include <math/optimize/scextrapmat.h>

#define CLASSNAME SymmSCMatrixSCExtrapData
#define PARENTS public SCExtrapData
#define HAVE_STATEIN_CTOR
#include <util/class/classi.h>
void *
SymmSCMatrixSCExtrapData::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCExtrapData::_castdown(cd);
  return do_castdowns(casts,cd);
}

SymmSCMatrixSCExtrapData::SymmSCMatrixSCExtrapData(StateIn& s) :
  SCExtrapData(s)
{
  RefSCMatrixKit k = SCMatrixKit::default_matrixkit();
  RefSCDimension dim;
  dim.restore_state(s);
  m = k->symmmatrix(dim);
  m.restore(s);
}

SymmSCMatrixSCExtrapData::SymmSCMatrixSCExtrapData(const RefSymmSCMatrix& mat)
{
  m = mat;
}

void
SymmSCMatrixSCExtrapData::save_data_state(StateOut& s)
{
  SCExtrapData::save_data_state(s);
  m.dim().save_state(s);
  m.save(s);
}

void
SymmSCMatrixSCExtrapData::zero()
{
  m.assign(0.0);
}

SCExtrapData*
SymmSCMatrixSCExtrapData::copy()
{
  return new SymmSCMatrixSCExtrapData(m.copy());
}

void
SymmSCMatrixSCExtrapData::accumulate_scaled(double scale,
                                            const RefSCExtrapData& data)
{
  SymmSCMatrixSCExtrapData* a
      = SymmSCMatrixSCExtrapData::require_castdown(
          data.pointer(), "SymmSCMatrixSCExtrapData::accumulate_scaled");

  RefSymmSCMatrix am = a->m.copy();
  am.scale(scale);
  m.accumulate(am);
}

///////////////////////////////////////////////////////////////////////////

#define CLASSNAME SymmSCMatrix2SCExtrapData
#define PARENTS public SCExtrapData
#define HAVE_STATEIN_CTOR
#include <util/class/classi.h>
void *
SymmSCMatrix2SCExtrapData::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCExtrapData::_castdown(cd);
  return do_castdowns(casts,cd);
}

SymmSCMatrix2SCExtrapData::SymmSCMatrix2SCExtrapData(StateIn&s) :
  SCExtrapData(s)
{
  RefSCMatrixKit k = SCMatrixKit::default_matrixkit();
  RefSCDimension dim;
  dim.restore_state(s);
  m1 = k->symmmatrix(dim);
  m2 = k->symmmatrix(dim);
  m1.restore(s);
  m2.restore(s);
}

SymmSCMatrix2SCExtrapData::SymmSCMatrix2SCExtrapData(
    const RefSymmSCMatrix& mat1,
    const RefSymmSCMatrix& mat2)
{
  m1 = mat1;
  m2 = mat2;
}

void
SymmSCMatrix2SCExtrapData::save_data_state(StateOut& s)
{
  SCExtrapData::save_data_state(s);
  m1.dim().save_state(s);
  m1.save(s);
  m2.save(s);
}

void
SymmSCMatrix2SCExtrapData::zero()
{
  m1.assign(0.0);
  m2.assign(0.0);
}

SCExtrapData*
SymmSCMatrix2SCExtrapData::copy()
{
  return new SymmSCMatrix2SCExtrapData(m1.copy(), m2.copy());
}

void
SymmSCMatrix2SCExtrapData::accumulate_scaled(double scale,
                                             const RefSCExtrapData& data)
{
  SymmSCMatrix2SCExtrapData* a
      = SymmSCMatrix2SCExtrapData::require_castdown(
          data.pointer(), "SymmSCMatrix2SCExtrapData::accumulate_scaled");

  RefSymmSCMatrix am = a->m1.copy();
  am.scale(scale);
  m1.accumulate(am);
  am = 0;

  am = a->m2.copy();
  am.scale(scale);
  m2.accumulate(am);
}

///////////////////////////////////////////////////////////////////////////

#define CLASSNAME SymmSCMatrixSCExtrapError
#define PARENTS public SCExtrapError
#define HAVE_STATEIN_CTOR
#include <util/class/classi.h>
void *
SymmSCMatrixSCExtrapError::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCExtrapError::_castdown(cd);
  return do_castdowns(casts,cd);
}

SymmSCMatrixSCExtrapError::SymmSCMatrixSCExtrapError(StateIn& s) :
  SCExtrapError(s)
{
  RefSCMatrixKit k = SCMatrixKit::default_matrixkit();
  RefSCDimension dim;
  dim.restore_state(s);

  int blocked;
  s.get(blocked);
  
  if (blocked)
    k = new BlockedSCMatrixKit(SCMatrixKit::default_matrixkit());
  
  m = k->symmmatrix(dim);
  m.restore(s);
}

SymmSCMatrixSCExtrapError::SymmSCMatrixSCExtrapError(
    const RefSymmSCMatrix& mat)
{
  m = mat;
}

void
SymmSCMatrixSCExtrapError::save_data_state(StateOut& s)
{
  SCExtrapError::save_data_state(s);
  m.dim().save_state(s);

  int blocked = (BlockedSymmSCMatrix::castdown(m)) ? 1 : 0;
  s.put(blocked);
  
  m.save(s);
}

double
SymmSCMatrixSCExtrapError::error()
{
  return m->maxabs();
}

double
SymmSCMatrixSCExtrapError::scalar_product(const RefSCExtrapError& arg)
{
  SymmSCMatrixSCExtrapError* a
      = SymmSCMatrixSCExtrapError::require_castdown(
          arg.pointer(), "SymmSCMatrixSCExtrapError::scalar_product");
  RefSCElementScalarProduct sp(new SCElementScalarProduct);
  m->element_op(sp.pointer(), a->m.pointer());
  return sp->result();
}
