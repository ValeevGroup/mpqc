
#include <stdio.h>
#include <chemistry/qc/basis/basis.h>
#include "integral.h"

///////////////////////////////////////////////////////////////////////

OneBodyIntShellIter::OneBodyIntShellIter(OneBodyInt*i):
  obi(i)
{
  // allocate a buffer
  int biggest_shell = obi->basis().max_nfunction_in_shell();
  if (biggest_shell) {
      buffer_ = new double[biggest_shell * biggest_shell];
    }
  else {
      buffer_ = 0;
    }
  
  // set the initial range for the computation
  i_start_ = j_start_ = 0;
  i_end_ = j_end_ = obi->basis().nshell();
}

void OneBodyIntShellIter::set_range(
  int i_start,
  int i_length,
  int j_start,
  int j_length)
{
  i_start_ = i_start;
  j_start_ = j_start;
  i_end_ = i_start + i_length;
  j_end_ = j_start + j_length;
}

OneBodyIntShellIter::~OneBodyIntShellIter()
{
  if (buffer_) delete[] buffer_;
}

///////////////////////////////////////////////////////////////////////

OneBodyIntRedundantShellIter::OneBodyIntRedundantShellIter(OneBodyInt*obi):
  OneBodyIntShellIter(obi)
{
}
OneBodyIntRedundantShellIter::~OneBodyIntRedundantShellIter()
{
}

void OneBodyIntRedundantShellIter::start()
{
  i_ = i_start_;
  j_ = j_start_;

  i_len_ = obi->basis()[i_].nfunction();
  j_len_ = obi->basis()[j_].nfunction();
  i_function_ = obi->basis().shell_to_function(i_);
  j_function_ = obi->basis().shell_to_function(j_);
  obi->compute_shell(i_,j_,buffer_);
}

void OneBodyIntRedundantShellIter::next()
{
  j_++;
  if (j_ == j_end_) {
      i_++;
      j_ = j_start_;
      if (!*this) return;
    }
  
  if (*this) {
      i_len_ = obi->basis()[i_].nfunction();
      j_len_ = obi->basis()[j_].nfunction();
      i_function_ = obi->basis().shell_to_function(i_);
      j_function_ = obi->basis().shell_to_function(j_);
      obi->compute_shell(i_,j_,buffer_);
    }
}

void OneBodyIntRedundantShellIter::operator int()
{
  return i_ != i_end_;
}

///////////////////////////////////////////////////////////////////////

OneBodyIntNonredundantNonrepeatedShellIter::
  OneBodyIntNonredundantNonrepeatedShellIter(OneBodyInt*obi):
  OneBodyIntShellIter(obi),
  diagonal(1)
{
}
OneBodyIntNonredundantNonrepeatedShellIter::
  ~OneBodyIntNonredundantNonrepeatedShellIter()
{
}

void OneBodyIntNonredundantNonrepeatedShellIter::set_range(
  int i_start,
  int i_length,
  int j_start,
  int j_length)
{
  i_start_ = i_start;
  j_start_ = j_start;
  i_end_ = i_start + i_length;
  j_end_ = j_start + j_length;

  // test for a legal block:
  if (i_start_ == j_start_ && i_end_ == j_end_) {
      // diagonal block, OK
      diagonal = 1;
    }
  else if (j_end_ <= i_start_) {
      // part of the lower triangle, OK
      diagonal = 0;
    }
  else {
      fprintf(stderr,"OneBodyIntNonredundantNonrepeatedShellIter::set_range: "
	      "illegal range\n");
      exit(1);
    }
}

void OneBodyIntNonredundantNonrepeatedShellIter::start()
{
  if (diagonal) {
      i_ = i_start_ + 1;
      j_ = j_start_;
    }
  else {
      i_ = i_start_;
      j_ = j_start_;
    }

  if (*this) {
      i_len_ = obi->basis()[i_].nfunction();
      j_len_ = obi->basis()[j_].nfunction();
      i_function_ = obi->basis().shell_to_function(i_);
      j_function_ = obi->basis().shell_to_function(j_);
      obi->compute_shell(i_,j_,buffer_);
    }
}

void OneBodyIntNonredundantNonrepeatedShellIter::next()
{
  j_++;
  if (j_ == j_end_ || j_ == i_) {
      i_++;
      j_ = j_start_;
      if (!*this) return;
    }

  i_len_ = obi->basis()[i_].nfunction();
  j_len_ = obi->basis()[j_].nfunction();
  i_function_ = obi->basis().shell_to_function(i_);
  j_function_ = obi->basis().shell_to_function(j_);
  obi->compute_shell(i_,j_,buffer_);
}

void OneBodyIntNonredundantNonrepeatedShellIter::operator int()
{
  return i_ != i_end_;
}

///////////////////////////////////////////////////////////////////////

OneBodyIntNonredundantRestShellIter::
  OneBodyIntNonredundantRestShellIter(OneBodyInt*obi):
  OneBodyIntShellIter(obi),
  diagonal(1)
{
}
OneBodyIntNonredundantRestShellIter::~OneBodyIntNonredundantRestShellIter()
{
}

void OneBodyIntNonredundantRestShellIter::set_range(
  int i_start,
  int i_length,
  int j_start,
  int j_length)
{
  i_start_ = i_start;
  j_start_ = j_start;
  i_end_ = i_start + i_length;
  j_end_ = j_start + j_length;

  // test for a legal block:
  if (i_start_ == j_start_ && i_end_ == j_end_) {
      // diagonal block, OK
      diagonal = 1;
    }
  else if (j_end_ <= i_start_) {
      // part of the lower triangle, OK
      diagonal = 0;
    }
  else {
      fprintf(stderr,"OneBodyIntNonredundantRestShellIter::set_range: "
	      "illegal range\n");
      exit(1);
    }
}

void OneBodyIntNonredundantRestShellIter::start()
{
  if (diagonal) {
      i_ = i_start_;
      j_ = j_start_;
    }

  if (*this) {
      i_len_ = obi->basis()[i_].nfunction();
      j_len_ = obi->basis()[j_].nfunction();
      i_function_ = obi->basis().shell_to_function(i_);
      j_function_ = obi->basis().shell_to_function(j_);
      obi->compute_shell(i_,j_,buffer_);
    }
}

void OneBodyIntNonredundantRestShellIter::next()
{
  if (!diagonal) return;

  j_++; i_++;

  if (!*this) return;

  i_len_ = obi->basis()[i_].nfunction();
  j_len_ = obi->basis()[j_].nfunction();
  i_function_ = obi->basis().shell_to_function(i_);
  j_function_ = obi->basis().shell_to_function(j_);
  obi->compute_shell(i_,j_,buffer_);
}

void OneBodyIntNonredundantRestShellIter::operator int()
{
  return i_ != i_end_ && diagonal;
}

///////////////////////////////////////////////////////////////////////

OneBodyIntShellIter* OneBodyInt::get_nonredundant_nonrepeated_shell_iter()
{
  return new OneBodyIntNonredundantNonrepeatedShellIter(this);
}

OneBodyIntShellIter* OneBodyInt::get_nonredundant_rest_shell_iter()
{
  return new OneBodyIntNonredundantRestShellIter(this);
}

OneBodyIntShellIter* OneBodyInt::get_redundant_shell_iter()
{
  return new OneBodyIntRedundantShellIter(this);
}

OneBodyInt::OneBodyInt(const GaussianBasisSet*b):
  bs(b)
{
}
OneBodyInt::~OneBodyInt()
{
}
int OneBodyInt::nbasis()
{
  return bs->nbasis();
}
int OneBodyInt::nshell()
{
  return bs->nshell();
}
GaussianBasisSet& OneBodyInt::basis()
{
  return *bs;
}

void OneBodyInt::compute(SymmetricMatrix&matrix)
{
  // set up the matrix
  matrix.ReDimension(bs->nbasis());

  OneBodyIntShellIter& obinrepi = *get_nonredundant_nonrepeated_shell_iter();

  for (obinrepi.start(); obinrepi; obinrepi.next()) {
      for (int i=obinrepi.i_function(); i<obinrepi.i_function_fence(); i++) {
	  for (int j=obinrepi.j_function(); j<obinrepi.j_function_fence(); j++) {
	      matrix.element(i,j) = obinrepi.val_by_overall_bf(i,j);
	    }
	}
    }

  OneBodyIntShellIter& obiresti = *get_nonredundant_rest_shell_iter();

  for (obiresti.start(); obiresti; obiresti.next()) {
      for (int i=obiresti.i_function(); i<obiresti.i_function_fence(); i++) {
	  for (int j=obiresti.j_function(); j<obiresti.j_function_fence(); j++) {
	      matrix.element(i,j) = obiresti.val_by_overall_bf(i,j);
	    }
	}
    }
}

///////////////////////////////////////////////////////////////////////

TestOneBodyInt::TestOneBodyInt(const GaussianBasisSet*bs):
  OneBodyInt(bs)
{
}
TestOneBodyInt::~TestOneBodyInt()
{
}

void TestOneBodyInt::compute_shell(int i,int j,double*buf)
{
  int ii,jj,index=0;
  for (ii=0; ii<basis()[i].nfunction(); ii++) {
      for (jj=0; jj<basis()[j].nfunction(); jj++) {
	  buf[index] = 0.0;
	  index++;
	}
    }
}
