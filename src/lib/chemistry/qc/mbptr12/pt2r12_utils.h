//
// pt2r12_utils.h
//
// Copyright (C) 2011 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_pt2r12_utils_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_pt2r12_utils_h

namespace {

  /// return the one-d offset for the index (i,j) of a triangle matrix
  int triang_half_INDEX_ordered(int i, int j) {
    return(i*(i+1)/2+j);
  }

  /// return the one-d offset for the index (i,j) of a triangle matrix
  int triang_half_INDEX(int i, int j) {
    return((i>j) ? triang_half_INDEX_ordered(i,j) : triang_half_INDEX_ordered(j,i));
  }

  /// return the one-d offset for the index (i,j) of a matrix
  int ordinary_INDEX(int i, int j, int coldim) {
    return(i*coldim+j);
  }

  /// tpdm_index and init_ioff: for indexing of density matrices.
  int tpdm_index(int i, int j, int k, int l,int coldim){
    //int ind_half1=triang_half_INDEX(i,j);
    //int ind_half2=triang_half_INDEX(k,l);
    int ind_half1=ordinary_INDEX(i,j,coldim);
    int ind_half2=ordinary_INDEX(k,l,coldim);
    return(triang_half_INDEX(ind_half1,ind_half2));
  }

  /// unpack a vector to a symmetric matrix; element-wise
  void vector_to_symmmatrix(RefSymmSCMatrix &matrix, const RefSCVector &vector) {
    int dim = matrix.dim().n();
    for(int i=0; i<dim; i++){
      for(int j=0; j<=i; j++) {
        matrix.set_element(i,j,vector.get_element(triang_half_INDEX(i,j)));
      }
    }
  }

  /// pack a symmetric matrix to a vector
  void symmmatrix_to_vector(RefSCVector &vector, const RefSymmSCMatrix &matrix) {
    int dim = matrix.dim().n();
    for(int i=0; i<dim; i++){
      for(int j=0; j<=i; j++) {
        vector.set_element(triang_half_INDEX(i,j),matrix.get_element(i,j));
      }
    }
  }

  /// unpack a vector to a matrix
  void vector_to_matrix(RefSCMatrix &matrix, const RefSCVector &vector) {
    int dim1 = matrix.rowdim().n();
    int dim2 = matrix.coldim().n();
    for(int i=0; i<dim1; i++) {
      for(int j=0; j<dim2; j++) {
        matrix.set_element(i,j,vector.get_element(ordinary_INDEX(i,j,dim2)));
      }
    }
  }

  int lowerupper_index(int p, int q);

  /// unpack a vector a matrix, depending on SpinCase2
  void vector_to_matrix(RefSCMatrix &matrix,const RefSCVector &vector,const SpinCase2 &pairspin) {
    int dim1 = matrix.rowdim().n();
    int dim2 = matrix.coldim().n();
    if(pairspin==AlphaBeta) {
      for(int i=0; i<dim1; i++) {
        for(int j=0; j<dim2; j++) {
          matrix.set_element(i,j,vector.get_element(ordinary_INDEX(i,j,dim2)));
        }
      }
    }
    else {  // pairspin==AlphaAlpha || pairspin==BetaBeta
      matrix->assign(0.0);
      for(int i=0; i<dim1; i++) {
        for(int j=0; j<i; j++) {
          const double value = vector.get_element(lowerupper_index(i,j));
          matrix.set_element(i,j,value);
          matrix.set_element(j,i,-value);
        }
      }
    }
  }

  /// pack a matrix to a vector
  void matrix_to_vector(RefSCVector &vector, const RefSCMatrix &matrix) {
    int dim1 = matrix.rowdim().n();
    int dim2 = matrix.coldim().n();
    for(int i=0; i<dim1; i++) {
      for(int j=0; j<dim2; j++) {
        vector.set_element(ordinary_INDEX(i,j,dim2),matrix.get_element(i,j));
      }
    }
  }

  /// pack a matrix to a vector
  void matrix_to_vector(RefSCVector &vector, const RefSymmSCMatrix &matrix) {
    int n = matrix.dim().n();
    for(int i=0; i<n; i++) {
      for(int j=0; j<=i; j++) {
        const double value = matrix.get_element(i,j);
        vector.set_element(ordinary_INDEX(i,j,n),value);
        vector.set_element(ordinary_INDEX(j,i,n),value);
      }
    }
  }

  /// pack a matrix to a vector
  void matrix_to_vector(RefSCVector &vector, const RefDiagSCMatrix &matrix) {
    int n = matrix.dim().n();
    for(int i=0; i<n; i++) {
      const double value = matrix.get_element(i);
      vector.set_element(ordinary_INDEX(i,i,n),value);
    }
  }

  /// pack a matrix to a vector depending on SpinCase2
  void matrix_to_vector(RefSCVector &vector, const RefSCMatrix &matrix,const SpinCase2 &pairspin) {
    int dim1 = matrix.rowdim().n();
    int dim2 = matrix.coldim().n();
    if(pairspin==AlphaBeta) {
      for(int i=0; i<dim1; i++) {
        for(int j=0; j<dim2; j++) {
          vector.set_element(ordinary_INDEX(i,j,dim2),matrix.get_element(i,j));
        }
      }
    }
    else {  // pairspin==AlphaAlpha || pairspin==BetaBeta
      for(int i=0; i<dim1; i++) {
        for(int j=0; j<i; j++) {
          vector.set_element(lowerupper_index(i,j),matrix.get_element(i,j));
        }
      }
    }
  }


  int lowertriang_index(int p,int q) {
    if(q>=p){
      throw ProgrammingError("lowertriang_index(p,q) -- q must be smaller than p.",__FILE__,__LINE__);
    }
    int index=p*(p+1)/2+q-p;
    return(index);
  }

  int lowerupper_index(int p, int q) {
    if(p>q) {
      return(lowertriang_index(p,q));
    }
    else if(q>p) {
      return(lowertriang_index(q,p));
    }
    else {
      throw ProgrammingError("lowerupper_index(p,q) -- p and q are not allowed to be equal.",__FILE__,__LINE__);
    }
  }
  double indexsizeorder_sign(int p,int q) {
    if(p>q) {
      return(1.0);
    }
    else if(q>p) {
      return(-1.0);
    }
    else {
      return(0.0);
    }
  }

  int antisym_pairindex(int i, int j) {
    int max_ij = std::max(i, j);
    int min_ij = std::min(i, j);
    return (max_ij -1)* max_ij/2 + min_ij;
  }

  /** this is a wrapper method so that we can get elements from 4-index antisymmetric
      matrix just as other matrices;
      as long as MM supports the method get_element */
  template <typename MatrixType>
  double get_4ind_antisym_matelement(MatrixType & MM, const int & U1, const int & U2,
                                                      const int & L1, const int & L2)
  {
      if(U1 == U2 || L1 == L2) return 0.0;
      else
      {
          int uppind = antisym_pairindex(U1, U2);
          int lowind = antisym_pairindex(L1, L2);
          const double totalsign = indexsizeorder_sign(U1,U2) * indexsizeorder_sign(L1, L2);
          return totalsign * MM.get_element(uppind, lowind);
      }
  }

  /** this is a wrapper method so that we can easily retrieve elements
      from 4-index matrix; as long as MM supports the method get_element */
  template <typename MatrixType>
  double get_4ind_matelement(MatrixType & MM, const int & U1, const int & U2,
                                              const int & L1, const int & L2,
                                              const int & UppDim, const int & LowDim)
  {
      int uppind = U1 * UppDim + U2;
      int lowind = L1 * LowDim + L2;
      return MM.get_element(uppind, lowind);
  }

  RefSCMatrix convert_RefSC_to_local_kit(const RefSCMatrix& A)
  {
    RefSCMatrix result;
    Ref<LocalSCMatrixKit> kit_cast_to_local; kit_cast_to_local << A.kit();
    if (kit_cast_to_local == 0) {
      Ref<LocalSCMatrixKit> local_kit = new LocalSCMatrixKit();
      RefSCMatrix A_local = local_kit->matrix(A.rowdim(), A.coldim());
      A_local->convert(A);
      result = A_local;
    }
    else
      result = A;
    return result;
  }


  RefSCMatrix transform_one_ind(RefSCMatrix AA, RefSCMatrix BB, int whichindex, int onedim) // transform one index. A is a two-index tensor;
  {                       // B is a 4-ind tensor; whichindex tells which to transform; onedim (1-4) tells the dimension of one ind; the second ind of A is the dummy
#if 0
    ExEnv::out0() << "test transform_one_ind:\n";
    ExEnv::out0() << (onedim*onedim)<< ", " << BB->nrow() << ", "<< BB->ncol() << "\n";
    AA.print(prepend_spincase(AlphaBeta, "transform_one_ind: AA").c_str());
    BB.print(prepend_spincase(AlphaBeta, "transform_one_ind: BB").c_str());
#endif
    MPQC_ASSERT(((onedim*onedim) == BB->nrow()) and (BB->nrow() == BB->ncol()));
    RefSCMatrix res = BB->clone();
    res.assign(0.0);
    int ext_ind, int_ind, row, col, Ap, Bp, Cp, Dp, A, B, C, D, a, b, c, d, f;
    ext_ind = int_ind = row = col = Ap = Bp = Cp = Dp = A = B = C = D = a = b = c = d = f = 0;
    double xx = 0;
    for (a = 0; a < onedim; ++a) // the external index of A
    {
      for (b = 0; b < onedim; ++b)
      {
        for (c = 0; c < onedim; ++c)
        {
          for (d = 0; d < onedim; ++d)
          {
            xx = 0;
            for (f = 0; f < onedim; ++f) // dummy index
            {
              switch (whichindex)
              {
                case 1: // Gamma^AB_CD * C_A^Ap
                  ext_ind = Ap = a; int_ind = A = f;  B = b; C = c; D = d;// by renaming, we have BB always of the same form, convenient
                  xx += AA->get_element(ext_ind, int_ind) * BB->get_element(A*onedim + B, C*onedim + D);
                  break;
                case 2: // Gamma^AB_CD * C_B^Bp
                  ext_ind = Bp = a; int_ind = B = f; A = b; C = c; D = d;
                  xx += AA->get_element(ext_ind, int_ind) *BB->get_element(A*onedim + B, C*onedim + D);
                  break;
                case 3: // C_Cp^C * Gamma^AB_CD
                  ext_ind = Cp = a; int_ind = C = f; A = b; B = c; D = d;
                  xx += BB->get_element(A*onedim + B, C*onedim + D)*AA->get_element(ext_ind, int_ind);//AA->(ext,int) == Transpose(AA)(int, ext)
                  break;
                case 4: // C_Dp^D * Gamma^AB_CD
                  ext_ind = Dp = a; int_ind = D = f; A = b; B = c; C = d;
                  xx += BB->get_element(A*onedim + B, C*onedim + D)*AA->get_element(ext_ind, int_ind);
                  break;
                default: abort();
              }
           }
           switch (whichindex)
           {
             case 1:// Gamma^AB_CD * C_A^Ap
               col = C*onedim + D; row = Ap*onedim +B;
               break;
             case 2:// Gamma^AB_CD * C_B^Bp
               col = C*onedim + D; row = A*onedim + Bp;
               break;
             case 3: // C_Cp^C * Gamma^AB_CD
               col = Cp *onedim + D; row = A*onedim + B;
               break;
             case 4: // C_Dp^D * Gamma^AB_CD
               col = C*onedim + Dp; row = A*onedim + B;
               break;
             default: abort();
           }
           res->set_element(row,col, xx);
//           ExEnv::out0() << "row, col, xx: " << row << ", " << col << ", " << xx << "\n";
          }
        }
      }
    }
    //res.print(prepend_spincase(AlphaBeta, "transform_one_ind: res").c_str());
    return res;
  }

  RefSymmSCMatrix convert_to_local_kit(const RefSymmSCMatrix& A) {   //forward declaration
     RefSymmSCMatrix result;
     Ref<LocalSCMatrixKit> kit_cast_to_local; kit_cast_to_local << A.kit();
     if (kit_cast_to_local == 0) {
       Ref<LocalSCMatrixKit> local_kit = new LocalSCMatrixKit();
       RefSymmSCMatrix A_local = local_kit->symmmatrix(A.dim());
       A_local->convert(A);
       result = A_local;
     }
     else
       result = A;
     return result;
   }

  /// mat corresponds to a 4 ind tensor (b1 b2, k1 k2). Return a reshaped matrix (b1, b2k1k2)
  RefSCMatrix RefSCMAT_combine234(RefSCMatrix mat, const int b1, const int b2,
                                  const int k1, const int k2)
  {
    const int ncol = b2 * k1 * k2;
    RefSCDimension rowdim = new SCDimension(b1);
    RefSCDimension coldim = new SCDimension(b2*k1*k2);
    RefSCMatrix res = mat->kit()->matrix(rowdim, coldim);
    res->assign(0.0);
    const int num_e = k2*k1;

    for (int I1 = 0; I1 < b1; ++I1)
    {
      for (int I2 = 0; I2 < b2; ++I2)
      {
        const int mat_row = I1 * b2 + I2;
        const int blockbegin = I2*num_e;
        const int blockend = blockbegin + num_e - 1;

//        RefSCVector k1k2row = mat->get_row(mat_row);
//        RefSCMatrix temp_mat = mat->kit()->matrix(new SCDimension(1), mat->coldim());
//        temp_mat->assign_row(k1k2row, 0);
//        res->assign_subblock(temp_mat, I1, I1, blockbegin, blockend);
        res->assign_subblock(mat, I1, I1, blockbegin, blockend, mat_row, 0);
      }
    }
    return res;
  }


} // end of anonymous namespace

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
