//
// orthog.cc
//
// Copyright (C) 1998 Sandia National Laboratories
//
// Author: Curtis Janssen <cljanss@aros.ca.sandia.gov>
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <util/misc/formio.h>
#include <math/scmat/blkiter.h>
#include <chemistry/qc/scf/scf.h>

class SCElementMaxDiff: public SCElementOp {
#   define CLASSNAME SCElementMaxDiff
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int deferred_;
    double diff_1_diag_;
    double diff_0_diag_;
    double diff_0_offdiag_;
  public:
    SCElementMaxDiff();
    SCElementMaxDiff(StateIn&);
    ~SCElementMaxDiff();
    void save_data_state(StateOut&);
    void process(SCMatrixBlockIter&);
    int has_collect();
    void defer_collect(int);
    void collect(const RefMessageGrp&);
    double diff_1_diag() { return diff_1_diag_; }
    double diff_0_diag() { return diff_0_diag_; }
    double diff_0_offdiag() { return diff_0_offdiag_; }
};
DescribedClass_REF_dec(SCElementMaxDiff);

DescribedClass_REF_def(SCElementMaxDiff);
#define CLASSNAME SCElementMaxDiff
#define PARENTS   public SCElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

SCElementMaxDiff::SCElementMaxDiff():
  deferred_(0),
  diff_1_diag_(0.0),
  diff_0_diag_(0.0),
  diff_0_offdiag_(0.0)
{
}
SCElementMaxDiff::SCElementMaxDiff(StateIn&s):
  SCElementOp(s)
{
  s.get(diff_1_diag_);
  s.get(diff_0_diag_);
  s.get(diff_0_offdiag_);
  s.get(deferred_);
}
void
SCElementMaxDiff::save_data_state(StateOut&s)
{
  s.put(diff_1_diag_);
  s.put(diff_0_diag_);
  s.put(diff_0_offdiag_);
  s.put(deferred_);
}
void *
SCElementMaxDiff::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp::_castdown(cd);
  return do_castdowns(casts,cd);
}
SCElementMaxDiff::~SCElementMaxDiff() {}
void
SCElementMaxDiff::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      double v = i.get();
      double av = fabs(v);
      if (i.i() == i.j()) {
          double avm1 = fabs(v-1.0);
          if (v > 0.5) { if (avm1 > diff_1_diag_) diff_1_diag_ = avm1; }
          else if (av > diff_0_diag_) diff_0_diag_ = av;
        }
      else {
          if (av > diff_0_offdiag_) diff_0_offdiag_ = av;
        }
    }
}
int
SCElementMaxDiff::has_collect()
{
  return 1;
}
void
SCElementMaxDiff::defer_collect(int h)
{
  deferred_=h;
}
void
SCElementMaxDiff::collect(const RefMessageGrp&msg)
{
  if (!deferred_) {
    msg->max(diff_1_diag_);
    msg->max(diff_0_diag_);
    msg->max(diff_0_offdiag_);
    }
}

///////////////////////////////////////////////////////////////////////////


// C.t() * S * C should be a unit matrix.  There can be zeros
// on the diagonal if the basis set is linearly dependent and
// symmetric orthogonalization is used.
void
SCF::orthog_vector(RefSCMatrix &C, const RefSymmSCMatrix &S,
                   double diag_one_tolerance,
                   double offdiag_tolerance,
                   double diag_zero_tolerance)
{
  RefSymmSCMatrix CtSC(oso_dimension(), basis_matrixkit());
  CtSC.assign(0.0);
  CtSC.accumulate_transform(C, S, SCMatrix::TransposeTransform);
  if (debug_ > 1) CtSC.print("CtSC");

  RefSCElementMaxDiff maxdiff = new SCElementMaxDiff;
  CtSC.element_op(maxdiff.pointer());

  if (maxdiff->diff_1_diag() > diag_one_tolerance
      || maxdiff->diff_0_diag() > diag_zero_tolerance
      || maxdiff->diff_0_offdiag() > offdiag_tolerance) {
      cout << node0 << indent << "SCF: WARNING: nonorthogonal orbitals:"
           << " using GS orthogonalization" << endl;
      if (debug_ > 0) {
        cout << node0 << indent
             << "diff_1_diag    = " << maxdiff->diff_1_diag() << endl;
        cout << node0 << indent
             << "diff_0_diag    = " << maxdiff->diff_0_diag() << endl;
        cout << node0 << indent
             << "diff_0_offdiag = " << maxdiff->diff_0_offdiag() << endl;
        C.print("C before GS");
        CtSC.print("CtSC before GS");
      }
      C->schmidt_orthog(S.pointer(),basis()->nbasis());
      if (debug_ > 0) {
        C.print("C after GS");
        CtSC.print("CtSC after GS");
      }
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
