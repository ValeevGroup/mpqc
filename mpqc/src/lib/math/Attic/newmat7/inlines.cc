
#define WANT_MATH

#include "newmat.h"
#include "except.h"
#include "newmatap.h"
#include "newmatrm.h"
#include "newmatrc.h"
#include "precisio.h"

// These functions were orginally inlined in newmat.h.

MatrixDetails::~MatrixDetails() {}

LogAndSign::LogAndSign() { log_value=0.0; sign=1; }
void LogAndSign::ChangeSign() { sign = -sign; }
Real LogAndSign::LogValue() const { return log_value; }
int LogAndSign::Sign() const { return sign; }

void ExeCounter::operator++() { nexe++; }

MatrixConversionCheck::MatrixConversionCheck() { DoCheck=TRUE; } 
MatrixConversionCheck::~MatrixConversionCheck() { DoCheck=FALSE; }
void MatrixConversionCheck::CleanUp() { DoCheck=FALSE; }
Boolean MatrixConversionCheck::IsOn() { return DoCheck; }

MatrixType::nTypes() { return 9; }
MatrixType::MatrixType () {}
MatrixType::MatrixType (int i) : attribute(i) {}
int MatrixType::operator+() const { return attribute; }
MatrixType MatrixType::operator+(MatrixType mt) const
      { return MatrixType(attribute & mt.attribute); }
Boolean MatrixType::operator>=(MatrixType mt) const
      { return ( attribute & mt.attribute ) == attribute; }
Boolean MatrixType::operator==(MatrixType t) const
      { return (attribute == t.attribute); }
Boolean MatrixType::operator!=(MatrixType t) const
      { return (attribute != t.attribute); }
Boolean MatrixType::operator!() const { return (attribute & Valid) == 0; }
MatrixType MatrixType::i() const                        // type of inverse
      { return MatrixType(attribute & ~(Band+LUDeco)); }
MatrixType MatrixType::AddEqualEl() const          // Add constant to matrix
      { return MatrixType(attribute & (Valid + Symmetric)); }
MatrixType MatrixType::sub() const                      // type of submatrix
      { return MatrixType(attribute & Valid); }
MatrixType MatrixType::ssub() const            // type of sym submatrix
      { return MatrixType(attribute); }        // not for selection matrix
MatrixType::~MatrixType() {};

MatrixBandWidth::MatrixBandWidth(const int l, const int u) : lower(l), upper (u) {}
MatrixBandWidth::MatrixBandWidth(const int i) : lower(i), upper(i) {}
MatrixBandWidth MatrixBandWidth::t() const { return MatrixBandWidth(upper,lower); }
Boolean MatrixBandWidth::operator==(const MatrixBandWidth& bw) const
      { return (lower == bw.lower) && (upper == bw.upper); }

int ArrayLengthSpecifier::Value() const { return value; }
ArrayLengthSpecifier::ArrayLengthSpecifier(int l) : value(l) {}

#ifdef __GNUG__
GeneralMatrix* BaseMatrix::Evaluate() { return Evaluate(MatrixTypeUnSp); }
#endif
void BaseMatrix::CleanUp() {}                     // to clear store

int GeneralMatrix::Nrows() const { return nrows; }          // get dimensions
int GeneralMatrix::Ncols() const { return ncols; }
int GeneralMatrix::Storage() const { return storage; }
Real* GeneralMatrix::Store() const { return store; }
void GeneralMatrix::Protect() { tag=-1; }                   // can't delete or reuse
int GeneralMatrix::Tag() const { return tag; }
void GeneralMatrix::Release() { tag=1; }                    // del store after next use
void GeneralMatrix::Release(int t) { tag=t; }               // del store after t accesses
void GeneralMatrix::ReleaseAndDelete() { tag=0; }           // delete matrix after use
void GeneralMatrix::operator<<(const BaseMatrix& X) { Eq(X,this->Type()); }
void GeneralMatrix::Solver(MatrixRowCol&, const MatrixRowCol&) {}
void GeneralMatrix::RestoreRow(MatrixRowCol&) {}    // Restore matrix row
void GeneralMatrix::RestoreCol(MatrixRowCol&) {}    // Restore matrix col
void GeneralMatrix::SetParameters(const GeneralMatrix*) {}

Matrix::Matrix() {}
void Matrix::operator=(Real f) { GeneralMatrix::operator=(f); }
void Matrix::operator=(const Matrix& m) { operator=((const BaseMatrix&)m); }
#ifdef SETUP_C_SUBSCRIPTS
Real* Matrix::operator[](int m) { return store+m*ncols; }
const Real* Matrix::operator[](int m) const { return store+m*ncols; }
#endif
Matrix::Matrix(const Matrix& gm) { GetMatrix(&gm); }

nricMatrix::nricMatrix() {}
nricMatrix::nricMatrix(int m, int n)                     // standard declaration
   :  Matrix(m,n) { MakeRowPointer(); }
nricMatrix::nricMatrix(const BaseMatrix& bm)             // evaluate BaseMatrix
   :  Matrix(bm) { MakeRowPointer(); }
void nricMatrix::operator=(const BaseMatrix& bm)
   { DeleteRowPointer(); Matrix::operator=(bm); MakeRowPointer(); }
void nricMatrix::operator=(Real f) { GeneralMatrix::operator=(f); }
void nricMatrix::operator=(const nricMatrix& m) { operator=((const BaseMatrix&)m); }
void nricMatrix::operator<<(const BaseMatrix& X)
   { DeleteRowPointer(); Eq(X,this->Type()); MakeRowPointer(); }
nricMatrix::nricMatrix(const nricMatrix& gm) { GetMatrix(&gm); MakeRowPointer(); }
void nricMatrix::ReDimension(int m, int n)               // change dimensions
   { DeleteRowPointer(); Matrix::ReDimension(m,n); MakeRowPointer(); }
nricMatrix::~nricMatrix() { DeleteRowPointer(); }
#ifndef __ZTC__
Real** nricMatrix::nric() const { CheckStore(); return row_pointer-1; }
#endif

SymmetricMatrix::SymmetricMatrix() {}
void SymmetricMatrix::operator=(Real f) { GeneralMatrix::operator=(f); }
void SymmetricMatrix::operator=(const SymmetricMatrix& m) { operator=((const BaseMatrix&)m); }
#ifdef SETUP_C_SUBSCRIPTS
Real* SymmetricMatrix::operator[](int m) { return store+(m*(m+1))/2; }
const Real* SymmetricMatrix::operator[](int m) const { return store+(m*(m+1))/2; }
#endif
SymmetricMatrix::SymmetricMatrix(const SymmetricMatrix& gm) { GetMatrix(&gm); }

UpperTriangularMatrix::UpperTriangularMatrix() {}
void UpperTriangularMatrix::operator=(const UpperTriangularMatrix& m)
      { operator=((const BaseMatrix&)m); }
UpperTriangularMatrix::UpperTriangularMatrix(const UpperTriangularMatrix& gm) { GetMatrix(&gm); }
void UpperTriangularMatrix::operator=(Real f) { GeneralMatrix::operator=(f); }
#ifdef SETUP_C_SUBSCRIPTS
Real* UpperTriangularMatrix::operator[](int m) { return store+m*ncols-(m*(m+1))/2; }
const Real* UpperTriangularMatrix::operator[](int m) const { return store+m*ncols-(m*(m+1))/2; }
#endif
GeneralMatrix* UpperTriangularMatrix::MakeSolver() { return this; } // for solving

LowerTriangularMatrix::LowerTriangularMatrix() {}
LowerTriangularMatrix::LowerTriangularMatrix(const LowerTriangularMatrix& gm) { GetMatrix(&gm); }
void LowerTriangularMatrix::operator=(Real f) { GeneralMatrix::operator=(f); }
void LowerTriangularMatrix::operator=(const LowerTriangularMatrix& m)
     { operator=((const BaseMatrix&)m); }
#ifdef SETUP_C_SUBSCRIPTS
Real* LowerTriangularMatrix::operator[](int m) { return store+(m*(m+1))/2; }
const Real* LowerTriangularMatrix::operator[](int m) const { return store+(m*(m+1))/2; }
#endif
GeneralMatrix* LowerTriangularMatrix::MakeSolver() { return this; } // for solving

   DiagonalMatrix::DiagonalMatrix() {}
   DiagonalMatrix::DiagonalMatrix(const DiagonalMatrix& gm) { GetMatrix(&gm); }
   void DiagonalMatrix::operator=(Real f) { GeneralMatrix::operator=(f); }
   void DiagonalMatrix::operator=(const DiagonalMatrix& m) { operator=((const BaseMatrix&)m); }
#ifdef SETUP_C_SUBSCRIPTS
   Real& DiagonalMatrix::operator[](int m) { return store[m]; }
   const Real& DiagonalMatrix::operator[](int m) const { return store[m]; }
#endif
   GeneralMatrix* DiagonalMatrix::MakeSolver() { return this; } // for solving
#ifndef __ZTC__
   Real* DiagonalMatrix::nric() const
      { CheckStore(); return store-1; }         // for use by NRIC
#endif

   RowVector::RowVector() {}
   RowVector::RowVector(ArrayLengthSpecifier n) : Matrix(1,n.Value()) {}
   RowVector::RowVector(const RowVector& gm) { GetMatrix(&gm); }
   void RowVector::operator=(Real f) { GeneralMatrix::operator=(f); }
   void RowVector::operator=(const RowVector& m) { operator=((const BaseMatrix&)m); }
#ifdef SETUP_C_SUBSCRIPTS
   Real& RowVector::operator[](int m) { return store[m]; }
   const Real& RowVector::operator[](int m) const { return store[m]; }
#endif
#ifndef __ZTC__
   Real* RowVector::nric() const
      { CheckStore(); return store-1; }         // for use by NRIC
#endif

   ColumnVector::ColumnVector() {}
   ColumnVector::ColumnVector(ArrayLengthSpecifier n) : Matrix(n.Value(),1) {}
   ColumnVector::ColumnVector(const ColumnVector& gm) { GetMatrix(&gm); }
   void ColumnVector::operator=(Real f) { GeneralMatrix::operator=(f); }
   void ColumnVector::operator=(const ColumnVector& m) { operator=((const BaseMatrix&)m); }
#ifdef SETUP_C_SUBSCRIPTS
   Real& ColumnVector::operator[](int m) { return store[m]; }
   const Real& ColumnVector::operator[](int m) const { return store[m]; }
#endif
#ifndef __ZTC__
   Real* ColumnVector::nric() const
      { CheckStore(); return store-1; }         // for use by NRIC
#endif

   GeneralMatrix* CroutMatrix::MakeSolver() { return this; } // for solving
   void CroutMatrix::operator=(const CroutMatrix& m) { operator=((const BaseMatrix&)m); }

   BandMatrix::BandMatrix() { lower=0; upper=0; CornerClear(); }
   BandMatrix::BandMatrix(int n,int lb,int ub) { ReDimension(n,lb,ub); CornerClear(); }
   void BandMatrix::operator=(Real f) { GeneralMatrix::operator=(f); }
   void BandMatrix::operator=(const BandMatrix& m) { operator=((const BaseMatrix&)m); }
#ifdef SETUP_C_SUBSCRIPTS
   Real* BandMatrix::operator[](int m) { return store+(upper+lower)*m+lower; }
   const Real* BandMatrix::operator[](int m) const { return store+(upper+lower)*m+lower; }
#endif
   BandMatrix::BandMatrix(const BandMatrix& gm) { GetMatrix(&gm); }
   Real BandMatrix::SumSquare() const { CornerClear(); return GeneralMatrix::SumSquare(); }
   Real BandMatrix::SumAbsoluteValue() const
      { CornerClear(); return GeneralMatrix::SumAbsoluteValue(); }
   Real BandMatrix::MaximumAbsoluteValue() const
      { CornerClear(); return GeneralMatrix::MaximumAbsoluteValue(); }
   void BandMatrix::operator<<(const BaseMatrix& X) { GeneralMatrix::operator<<(X); }

   UpperBandMatrix::UpperBandMatrix() {}
   UpperBandMatrix::UpperBandMatrix(int n, int ubw)              // standard declaration
      : BandMatrix(n, 0, ubw) {}
   void UpperBandMatrix::operator=(Real f) { GeneralMatrix::operator=(f); }
   void UpperBandMatrix::operator=(const UpperBandMatrix& m)
      { operator=((const BaseMatrix&)m); }
    UpperBandMatrix::UpperBandMatrix(const UpperBandMatrix& gm) { GetMatrix(&gm); }
   GeneralMatrix* UpperBandMatrix::MakeSolver() { return this; }
   void UpperBandMatrix::ReDimension(int n,int ubw)              // change dimensions
      { BandMatrix::ReDimension(n,0,ubw); }
#ifdef SETUP_C_SUBSCRIPTS
   Real* UpperBandMatrix::operator[](int m) { return store+upper*m; }
   const Real* UpperBandMatrix::operator[](int m) const { return store+upper*m; }
#endif

   LowerBandMatrix::LowerBandMatrix() {}
   LowerBandMatrix::LowerBandMatrix(int n, int lbw)              // standard declaration
      : BandMatrix(n, lbw, 0) {}
   void LowerBandMatrix::operator=(Real f) { GeneralMatrix::operator=(f); }
   void LowerBandMatrix::operator=(const LowerBandMatrix& m)
      { operator=((const BaseMatrix&)m); }
   LowerBandMatrix::LowerBandMatrix(const LowerBandMatrix& gm) { GetMatrix(&gm); }
   GeneralMatrix* LowerBandMatrix::MakeSolver() { return this; }
   void LowerBandMatrix::ReDimension(int n,int lbw)             // change dimensions
      { BandMatrix::ReDimension(n,lbw,0); }
#ifdef SETUP_C_SUBSCRIPTS
   Real* LowerBandMatrix::operator[](int m) { return store+lower*(m+1); }
   const Real* LowerBandMatrix::operator[](int m) const { return store+lower*(m+1); }
#endif

   SymmetricBandMatrix::SymmetricBandMatrix() { lower=0; CornerClear(); }
   SymmetricBandMatrix::SymmetricBandMatrix(int n, int lb) { ReDimension(n,lb); CornerClear(); }
   void SymmetricBandMatrix::operator=(Real f) { GeneralMatrix::operator=(f); }
   void SymmetricBandMatrix::operator=(const SymmetricBandMatrix& m)
      { operator=((const BaseMatrix&)m); }
#ifdef SETUP_C_SUBSCRIPTS
   Real* SymmetricBandMatrix::operator[](int m) { return store+lower*(m+1); }
   const Real* SymmetricBandMatrix::operator[](int m) const { return store+lower*(m+1); }
#endif
   SymmetricBandMatrix::SymmetricBandMatrix(const SymmetricBandMatrix& gm) { GetMatrix(&gm); }
   Real SymmetricBandMatrix::MaximumAbsoluteValue() const
      { CornerClear(); return GeneralMatrix::MaximumAbsoluteValue(); }

   GeneralMatrix* BandLUMatrix::MakeSolver() { return this; } // for solving
   void BandLUMatrix::operator=(const BandLUMatrix& m) { operator=((const BaseMatrix&)m); }

   MultipliedMatrix::MultipliedMatrix(const BaseMatrix* bm1x, const BaseMatrix* bm2x)
      : bm1(bm1x),bm2(bm2x) {}

   AddedMatrix::AddedMatrix(const BaseMatrix* bm1x, const BaseMatrix* bm2x)
      : MultipliedMatrix(bm1x,bm2x) {}

   SolvedMatrix::SolvedMatrix(const BaseMatrix* bm1x, const BaseMatrix* bm2x)
      : MultipliedMatrix(bm1x,bm2x) {}

   SubtractedMatrix::SubtractedMatrix(const BaseMatrix* bm1x, const BaseMatrix* bm2x)
      : AddedMatrix(bm1x,bm2x) {}

   ShiftedMatrix::ShiftedMatrix(const BaseMatrix* bmx, Real fx) : bm(bmx),f(fx) {}

   ScaledMatrix::ScaledMatrix(const BaseMatrix* bmx, Real fx) : ShiftedMatrix(bmx,fx) {}

   NegatedMatrix::NegatedMatrix(const BaseMatrix* bmx) : bm(bmx) {}

   TransposedMatrix::TransposedMatrix(const BaseMatrix* bmx) : NegatedMatrix(bmx) {}

   InvertedMatrix::InvertedMatrix(const BaseMatrix* bmx) : NegatedMatrix(bmx) {}

   RowedMatrix::RowedMatrix(const BaseMatrix* bmx) : NegatedMatrix(bmx) {}

   ColedMatrix::ColedMatrix(const BaseMatrix* bmx) : NegatedMatrix(bmx) {}

   DiagedMatrix::DiagedMatrix(const BaseMatrix* bmx) : NegatedMatrix(bmx) {}

   MatedMatrix::MatedMatrix(const BaseMatrix* bmx, int nrx, int ncx)
      : NegatedMatrix(bmx), nr(nrx), nc(ncx) {}

   ConstMatrix::ConstMatrix(const GeneralMatrix* cgmx) : cgm(cgmx) {}

#ifndef TEMPS_DESTROYED_QUICKLY
   ReturnMatrixX::ReturnMatrixX(const ReturnMatrixX& tm) : gm(tm.gm) {}
#endif
   ReturnMatrixX::ReturnMatrixX(const GeneralMatrix* gmx) : gm((GeneralMatrix*&)gmx) {}

GetSubMatrix::GetSubMatrix
      (const BaseMatrix* bmx, int rs, int rn, int cs, int cn, Boolean is)
      : NegatedMatrix(bmx),
      row_skip(rs), row_number(rn), col_skip(cs), col_number(cn), IsSym(is) {}
   GetSubMatrix::GetSubMatrix(const GetSubMatrix& g)
      : NegatedMatrix(g.bm), row_skip(g.row_skip), row_number(g.row_number),
      col_skip(g.col_skip), col_number(g.col_number), IsSym(g.IsSym) {}
   void GetSubMatrix::operator=(const GetSubMatrix& m) { operator=((const BaseMatrix&)m); }

   MatrixInput::MatrixInput() { depth++; }
   MatrixInput::MatrixInput(const MatrixInput&) { depth++; }

   long SpaceException::st_type() { return 2; }
   long SpaceException::type() const { return 2; }
   void SpaceException::SetAction(int a) { action=a; }

   long MatrixException::st_type() { return 3; }
   long MatrixException::type() const { return 3; }

   long DataException::st_type() { return 3*53; }
   long DataException::type() const { return 3*53; }
   void DataException::SetAction(int a) { action=a; }

   long SingularException::st_type() { return 3*53*109; }
   long SingularException::type() const { return 3*53*109; }

   long NPDException::st_type() { return 3*53*113; }
   long NPDException::type() const { return 3*53*113; }

   long ConvergenceException::st_type() { return 3*59; }
   long ConvergenceException::type() const { return 3*59; }
   void ConvergenceException::SetAction(int a) { action=a; }

   long ProgramException::st_type() { return 3*61; }
   long ProgramException::type() const { return 3*61; }
   void ProgramException::SetAction(int a) { action=a; }

   long IndexException::st_type() { return 3*61*101; }
   long IndexException::type() const { return 3*61*101; }

   long VectorException::st_type() { return 3*61*107; }
   long VectorException::type() const { return 3*61*107; }

   long NotSquareException::st_type() { return 3*61*109; }
   long NotSquareException::type() const { return 3*61*109; }

   long SubMatrixDimensionException::st_type() { return 3*61*113; }
   long SubMatrixDimensionException::type() const { return 3*61*113; }

   long IncompatibleDimensionsException::st_type() { return 3*61*127; }
   long IncompatibleDimensionsException::type() const { return 3*61*127; }

   long NotDefinedException::st_type() { return 3*61*131; }
   long NotDefinedException::type() const { return 3*61*131; }

   long CannotBuildException::st_type() { return 3*61*137; }
   long CannotBuildException::type() const { return 3*61*137; }

   long InternalException::st_type() { return 3*67; }
   long InternalException::type() const { return 3*67; }
   void InternalException::SetAction(int a) { action=a; }

 LogAndSign LogDeterminant(const BaseMatrix& B)
   { return B.LogDeterminant(); }
 Real SumSquare(const BaseMatrix& B) { return B.SumSquare(); }
 Real Trace(const BaseMatrix& B) { return B.Trace(); }
 Real SumAbsoluteValue(const BaseMatrix& B)
   { return B.SumAbsoluteValue(); }
 Real MaximumAbsoluteValue(const BaseMatrix& B)
   { return B.MaximumAbsoluteValue(); }
 Real Norm1(const BaseMatrix& B) { return B.Norm1(); }
 Real Norm1(RowVector& RV) { return RV.MaximumAbsoluteValue(); }
 Real NormInfinity(const BaseMatrix& B) { return B.NormInfinity(); }
 Real NormInfinity(ColumnVector& CV)
   { return CV.MaximumAbsoluteValue(); } 

//JCM

 Real Norm2(const ColumnVector& CV) {return sqrt(CV.SumSquare());}
 Real Dot(ColumnVector& CV1, ColumnVector& CV2) {return CV1.Dot(CV2);}

   NEW_DELETE_IMPL(BaseMatrix)
   NEW_DELETE_IMPL(GeneralMatrix)
   NEW_DELETE_IMPL(Matrix)
   NEW_DELETE_IMPL(nricMatrix)
   NEW_DELETE_IMPL(SymmetricMatrix)
   NEW_DELETE_IMPL(UpperTriangularMatrix)
   NEW_DELETE_IMPL(LowerTriangularMatrix)
   NEW_DELETE_IMPL(DiagonalMatrix)
   NEW_DELETE_IMPL(RowVector)
   NEW_DELETE_IMPL(ColumnVector)
   NEW_DELETE_IMPL(CroutMatrix)
   NEW_DELETE_IMPL(BandMatrix)
   NEW_DELETE_IMPL(UpperBandMatrix)
   NEW_DELETE_IMPL(LowerBandMatrix)
   NEW_DELETE_IMPL(SymmetricBandMatrix)
   NEW_DELETE_IMPL(BandLUMatrix)
   NEW_DELETE_IMPL(MultipliedMatrix)
   NEW_DELETE_IMPL(AddedMatrix)
   NEW_DELETE_IMPL(SolvedMatrix)
   NEW_DELETE_IMPL(SubtractedMatrix)
   NEW_DELETE_IMPL(ShiftedMatrix)
   NEW_DELETE_IMPL(ScaledMatrix)
   NEW_DELETE_IMPL(NegatedMatrix)
   NEW_DELETE_IMPL(TransposedMatrix)
   NEW_DELETE_IMPL(InvertedMatrix)
   NEW_DELETE_IMPL(RowedMatrix)
   NEW_DELETE_IMPL(ColedMatrix)
   NEW_DELETE_IMPL(DiagedMatrix)
   NEW_DELETE_IMPL(MatedMatrix)
   NEW_DELETE_IMPL(ConstMatrix)
   NEW_DELETE_IMPL(ReturnMatrixX)
   NEW_DELETE_IMPL(GetSubMatrix)

// these functions are from except.h
   long Exception::st_type() { return 1; }
   long Exception::type() const { return 1; }

Tracer::Tracer(char* e)
   : entry(e), previous(Exception::last) { Exception::last = this; }

Tracer::~Tracer() { Exception::last = previous; }

void Tracer::ReName(char* e) { entry=e; }

   JumpItem::JumpItem() : trace(0), janitor(0), ji(JumpBase::jl)
      { JumpBase::jl = this; }
   JumpItem::~JumpItem() { JumpBase::jl = ji; }

   void Janitor::CleanUp() {}

// these are from newmatap.h

#ifndef __GNUG__
 void SVD(const Matrix& A, DiagonalMatrix& D, Matrix& U,
   Boolean withU) { SVD(A, D, U, U, withU, FALSE); }
#else
 void SVD(const Matrix& A, DiagonalMatrix& D, Matrix& U,
   Boolean withU) { SVD(A, D, U, U, withU, FALSE); }
#endif

   int LinearEquationSolver::search(const BaseMatrix*) const { return 0; }
   LinearEquationSolver::~LinearEquationSolver() { delete gm; }
   void LinearEquationSolver::CleanUp() { delete gm; } 
   GeneralMatrix* LinearEquationSolver::Evaluate(MatrixType) { return gm; }
   NEW_DELETE_IMPL(LinearEquationSolver)

// these are from newmatrm.h
   RectMatrixRowCol::RectMatrixRowCol(Real* st, int nx, int sp, int sh)
      : store(st), n(nx), spacing(sp), shift(sh) {}
   void RectMatrixRowCol::Reset(Real* st, int nx, int sp, int sh)
      { store=st; n=nx; spacing=sp; shift=sh; }
   Real& RectMatrixRowCol::operator[](int i) { return *(store+i*spacing); } // element
   Real& RectMatrixRowCol::First() { return *store; }                       // get first element
   void RectMatrixRowCol::DownDiag() { store += (shift+spacing); n--; }
   void RectMatrixRowCol::UpDiag() { store -= (shift+spacing); n++; }

   Real& RectMatrixRow::operator[](int i) { return *(store+i); }
   void RectMatrixRow::Down() { store += shift; }
   void RectMatrixRow::Right() { store++; n--; }
   void RectMatrixRow::Up() { store -= shift; }
   void RectMatrixRow::Left() { store--; n++; }

   void RectMatrixCol::Down() { store += spacing; n--; }
   void RectMatrixCol::Right() { store++; }
   void RectMatrixCol::Up() { store -= spacing; n++; }
   void RectMatrixCol::Left() { store--; }

   RectMatrixDiag::RectMatrixDiag(const DiagonalMatrix& D)
      : RectMatrixRowCol(D.Store(), D.Nrows(), 1, 1) {}
   Real& RectMatrixDiag::operator[](int i) { return *(store+i); }
   void RectMatrixDiag::DownDiag() { store++; n--; }
   void RectMatrixDiag::UpDiag() { store--; n++; }


 RectMatrixRow::RectMatrixRow
   (const Matrix& M, int row, int skip, int length)
   : RectMatrixRowCol( M.Store()+row*M.Ncols()+skip, length, 1, M.Ncols() ) {}

 RectMatrixRow::RectMatrixRow (const Matrix& M, int row)
   : RectMatrixRowCol( M.Store()+row*M.Ncols(), M.Ncols(), 1, M.Ncols() ) {}

 RectMatrixCol::RectMatrixCol
   (const Matrix& M, int skip, int col, int length)
   : RectMatrixRowCol( M.Store()+col+skip*M.Ncols(), length, M.Ncols(), 1 ) {}

 RectMatrixCol::RectMatrixCol (const Matrix& M, int col)
   : RectMatrixRowCol( M.Store()+col, M.Nrows(), M.Ncols(), 1 ) {}

 Real square(Real x) { return x*x; }
 Real sign(Real x, Real y)
   { return (y>=0) ? x : -x; }                    // assume x >=0

// these are from newmatrc.h

   LoadAndStoreFlag::LoadAndStoreFlag() {}
   LoadAndStoreFlag::LoadAndStoreFlag(int i) : ControlWord(i) {}
   LoadAndStoreFlag::LoadAndStoreFlag(LSF lsf) : ControlWord(lsf) {}
   LoadAndStoreFlag::LoadAndStoreFlag(const ControlWord& cwx) : ControlWord(cwx) {}

   void MatrixRowCol::IncrMat() { rowcol++; store += storage; }
					       // used by NextRow
   void MatrixRowCol::IncrDiag() { rowcol++; skip++; }
   void MatrixRowCol::IncrUT() { rowcol++; storage--; store += storage; skip++; }
   void MatrixRowCol::IncrLT() { rowcol++; store += storage; storage++; }
   Real* MatrixRowCol::operator()() { return store+skip; }   // pointer to first element
   Real* MatrixRowCol::Store() { return store; }
   int MatrixRowCol::Skip() { return skip; }                 // number of elements skipped
   int MatrixRowCol::Storage() { return storage; }           // number of elements stored
   void MatrixRowCol::Skip(int i) { skip=i; }
   void MatrixRowCol::Storage(int i) { storage=i; }
   MatrixRowCol::MatrixRowCol() {}                           // to stop warning messages


 MatrixRow::MatrixRow(GeneralMatrix* gmx, LoadAndStoreFlag cwx, int row)
{ gm=gmx; cw=cwx; rowcol=row; gm->GetRow(*this); } 

 void MatrixRow::Next() { gm->NextRow(*this); }

 MatrixCol::MatrixCol(GeneralMatrix* gmx, LoadAndStoreFlag cwx, int col)
{ gm=gmx; cw=cwx; rowcol=col; gm->GetCol(*this); } 

 MatrixCol::MatrixCol(GeneralMatrix* gmx, Real* r,
   LoadAndStoreFlag cwx, int col)
{ gm=gmx; store=r; cw=cwx+StoreHere; rowcol=col; gm->GetCol(*this); } 

 void MatrixCol::Next() { gm->NextCol(*this); }

// here are some dtors from newmat.h and newmatap.h

BaseMatrix::~BaseMatrix(){}
Matrix::~Matrix(){}
SymmetricMatrix::~SymmetricMatrix(){}
UpperTriangularMatrix::~UpperTriangularMatrix(){}
LowerTriangularMatrix::~LowerTriangularMatrix(){}
DiagonalMatrix::~DiagonalMatrix(){}
RowVector::~RowVector(){}
ColumnVector::~ColumnVector(){}
BandMatrix::~BandMatrix(){}
UpperBandMatrix::~UpperBandMatrix(){}
LowerBandMatrix::~LowerBandMatrix(){}
SymmetricBandMatrix::~SymmetricBandMatrix(){}
MultipliedMatrix::~MultipliedMatrix(){}
AddedMatrix::~AddedMatrix(){}
SolvedMatrix::~SolvedMatrix(){}
SubtractedMatrix::~SubtractedMatrix(){}
ShiftedMatrix::~ShiftedMatrix(){}
ScaledMatrix::~ScaledMatrix(){}
NegatedMatrix::~NegatedMatrix(){}
TransposedMatrix::~TransposedMatrix(){}
InvertedMatrix::~InvertedMatrix(){}
RowedMatrix::~RowedMatrix(){}
ColedMatrix::~ColedMatrix(){}
DiagedMatrix::~DiagedMatrix(){}
MatedMatrix::~MatedMatrix(){}
ConstMatrix::~ConstMatrix(){}
ReturnMatrixX::~ReturnMatrixX(){}
GetSubMatrix::~GetSubMatrix(){}
SymmetricEigenAnalysis::~SymmetricEigenAnalysis(){}

// inlines from boolean.h
   Boolean::Boolean(const int b) { value = b ? 1 : 0; }
   Boolean::Boolean(const void* b) { value = b ? 1 : 0; }
   Boolean::Boolean() {}
   Boolean::operator int() const { return value; }
   int Boolean::operator!() const { return !value; }

// inlines from controlw.h
   ControlWord::ControlWord() : cw(0) {}                     // do nothing
   ControlWord::ControlWord(int i) : cw(i) {}                // load an integer

      // select specific bits (for testing at least one set)
   ControlWord ControlWord::operator*(ControlWord i) const
      { return ControlWord(cw & i.cw); }
   void ControlWord::operator*=(ControlWord i)  { cw &= i.cw; }

      // set bits
   ControlWord ControlWord::operator+(ControlWord i) const
      { return ControlWord(cw | i.cw); }
   void ControlWord::operator+=(ControlWord i)  { cw |= i.cw; }

      // reset bits
   ControlWord ControlWord::operator-(ControlWord i) const
      { return ControlWord(cw - (cw & i.cw)); }
   void ControlWord::operator-=(ControlWord i) { cw -= (cw & i.cw); }

      // check if all of selected bits set or reset
   Boolean ControlWord::operator>=(ControlWord i) const { return (cw & i.cw) == i.cw; }
   Boolean ControlWord::operator<=(ControlWord i) const { return (cw & i.cw) == cw; }

      // flip selected bits
   ControlWord ControlWord::operator^(ControlWord i) const
      { return ControlWord(cw ^ i.cw); }
   ControlWord ControlWord::operator~() const { return ControlWord(~cw); }

      // convert to integer
   int ControlWord::operator+() const { return cw; }
   int ControlWord::operator!() const { return cw==0; }

// these are from precisio.h

#ifndef SystemV                    // if there is float.h

#ifdef USING_FLOAT
   int FloatingPointPrecision::Dig()
      { return FLT_DIG; }        // number of decimal digits or precision
   Real FloatingPointPrecision::Epsilon()
      { return FLT_EPSILON; }    // smallest number such that 1+Eps!=Eps
   int FloatingPointPrecision::Mantissa()
      { return FLT_MANT_DIG; }   // bits in mantisa
   Real FloatingPointPrecision::Maximum()
      { return FLT_MAX; }        // maximum value
   int FloatingPointPrecision::MaximumDecimalExponent()
      { return FLT_MAX_10_EXP; } // maximum decimal exponent
   int FloatingPointPrecision::MaximumExponent()
      { return FLT_MAX_EXP; }    // maximum binary exponent
   Real FloatingPointPrecision::Minimum()
      { return FLT_MIN; }        // minimum positive value
   int FloatingPointPrecision::MinimumDecimalExponent()
      { return FLT_MIN_10_EXP; } // minimum decimal exponent
   int FloatingPointPrecision::MinimumExponent()
      { return FLT_MIN_EXP; }    // minimum binary exponent
   int FloatingPointPrecision::Radix()
      { return FLT_RADIX; }      // exponent radix
   int FloatingPointPrecision::Rounds()
      { return FLT_ROUNDS; }     // addition rounding (1 = does round)
#endif

#ifdef USING_DOUBLE

   int FloatingPointPrecision::Dig()
      { return DBL_DIG; }        // number of decimal digits or precision
   Real FloatingPointPrecision::Epsilon()
      { return DBL_EPSILON; }    // smallest number such that 1+Eps!=Eps
   int FloatingPointPrecision::Mantissa()
      { return DBL_MANT_DIG; }   // bits in mantisa
   Real FloatingPointPrecision::Maximum()
      { return DBL_MAX; }        // maximum value
   int FloatingPointPrecision::MaximumDecimalExponent()
      { return DBL_MAX_10_EXP; } // maximum decimal exponent
   int FloatingPointPrecision::MaximumExponent()
      { return DBL_MAX_EXP; }    // maximum binary exponent
   Real FloatingPointPrecision::Minimum()
   {
#ifdef __BCPLUSPLUS__
       return 2.225074e-308;     // minimum positive value
#else
       return DBL_MIN;
#endif
   }
   int FloatingPointPrecision::MinimumDecimalExponent()
      { return DBL_MIN_10_EXP; } // minimum decimal exponent
   int FloatingPointPrecision::MinimumExponent()
      { return DBL_MIN_EXP; }    // minimum binary exponent
   int FloatingPointPrecision::Radix()
      { return FLT_RADIX; }      // exponent radix
   int FloatingPointPrecision::Rounds()
      { return FLT_ROUNDS; }     // addition rounding (1 = does round)

#endif

#endif


#ifdef SystemV                    // if there is no float.h

#ifdef USING_FLOAT

class FloatingPointPrecision
{
public:
   Real FloatingPointPrecision::Epsilon()
      { return pow(2.0,1-FSIGNIF); }  // smallest number such that 1+Eps!=Eps
   Real FloatingPointPrecision::Maximum()
      { return MAXFLOAT; }        // maximum value
   Real FloatingPointPrecision::Minimum()
      { return MINFLOAT; }        // minimum positive value
   FREE_CHECK(FloatingPointPrecision)
};

#endif


#ifdef USING_DOUBLE

   Real FloatingPointPrecision::Epsilon()
      { return pow(2.0,1-DSIGNIF); }  // smallest number such that 1+Eps!=Eps
   Real FloatingPointPrecision::Maximum()
      { return MAXDOUBLE; }          // maximum value
   Real FloatingPointPrecision::Minimum()
      { return MINDOUBLE; }
   FREE_CHECK(FloatingPointPrecision)

#endif

#endif

