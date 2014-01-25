//
// nao.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
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
#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/transform.h>

using namespace std;
using namespace sc;

#undef DEBUG

// namespace sc {
// static RefSCMatrix
// operator *(const RefDiagSCMatrix &d, const RefSymmSCMatrix &s)
// {
//   RefSCMatrix ret(s.dim(), s.dim(), s.kit());
//   int n = s.dim()->n();
//   for (int i=0; i<n; i++) {
//       for (int j=0; j<n; j++) {
//           ret.set_element(i,j, d.get_element(i)*s.get_element(i,j));
//         }
//     }
//   return ret;
// }
// }

static RefSymmSCMatrix
weight_matrix(const RefDiagSCMatrix &d, const RefSymmSCMatrix &s)
{
  RefSymmSCMatrix ret = s.clone();
  int n = s.dim()->n();
  for (int i=0; i<n; i++) {
      for (int j=0; j<=i; j++) {
          ret.set_element(i,j, s.get_element(i,j)
                               *d.get_element(i)*d.get_element(j));
        }
    }
  return ret;
}

static int
nnmb_atom(int z, int l)
{
  if (l==0) {
      if (z <= 2) return 1;
      else if (z <= 10) return 2;
      else if (z <= 18) return 3;
    }
  else if (l==1) {
      if (z <= 4) return 0;
      else if (z <= 12) return 1;
      else if (z <= 20) return 2;
    }
  else if (l==2) {
      if (z <= 20) return 0;
    }
  else if (l==3) {
      if (z <= 56) return 0;
    }
  else {
      return 0;
    }
  ExEnv::errn() << "NAO: z too big" << endl;
  abort();
  return 0;
}

static int
nnmb_all_atom(int z, int maxl)
{
  int ret = 0;
  for (int i=0; i<=maxl; i++) {
      ret += nnmb_atom(z,i) * (2*i+1);
    }
  return ret;
}

static RefSymmSCMatrix
mhalf(const RefSymmSCMatrix &S)
{
  RefSCDimension tdim = S.dim();
  Ref<SCMatrixKit> kit = S.kit();

  // find a symmetric orthogonalization transform
  RefSCMatrix trans(tdim,tdim,kit);
  RefDiagSCMatrix eigval(tdim,kit);

  S.diagonalize(eigval,trans);

  Ref<SCElementOp> squareroot = new SCElementSquareRoot;
  eigval.element_op(squareroot);

  Ref<SCElementOp> invert = new SCElementInvert(1.0e-12);
  eigval.element_op(invert);

  RefSymmSCMatrix OL(tdim,kit);
  OL.assign(0.0);
  // OL = trans * eigval * trans.t();
  OL.accumulate_transform(trans, eigval);
  return OL;
}

static void
delete_partition_info(int natom, int *maxam_on_atom,
                      int **nam_on_atom, int ***amoff_on_atom)
{
  int i, j;
  for (i=0; i<natom; i++) {
      for (j=0; j<=maxam_on_atom[i]; j++) {
          delete[] amoff_on_atom[i][j];
        }
      delete[] nam_on_atom[i];
      delete[] amoff_on_atom[i];
    }
  delete[] maxam_on_atom;
  delete[] nam_on_atom;
  delete[] amoff_on_atom;
}

#ifdef DEBUG
static double
ttrace(const RefSCMatrix &N,
       const RefSymmSCMatrix &P,
       const RefSymmSCMatrix &S)
{
  RefSCMatrix Nt = N.t();
  RefSymmSCMatrix Pt = P.clone();
  Pt.assign(0.0);
  Pt.accumulate_transform(Nt, P);
  RefSymmSCMatrix St = S.clone();
  St.assign(0.0);
  St.accumulate_transform(Nt, S);
  return (mhalf(St)*Pt*mhalf(St)).trace();
}

// for N giving an orthonormal basis
static double
ttrace(const RefSCMatrix &N,
       const RefSymmSCMatrix &P)
{
  RefSCMatrix Nt = N.t();
  RefSymmSCMatrix Pt = P.clone();
  Pt.assign(0.0);
  Pt.accumulate_transform(Nt, P);
  return Pt.trace();
}
#endif

static RefSCMatrix
assemble(const RefSCDimension dim,
         const RefSCMatrix &Nm, int *Nm_map,
         const RefSCMatrix &Nr1, int *Nr1_map,
         const RefSCMatrix &Nr2 = 0,  int *Nr2_map = 0)
{
  int nnmb = Nm.ncol();
  int nr1 = Nr1.ncol();
  int nr2 = (Nr2.null()?0:Nr2.ncol());
  int nb = dim.n();
  if (nb != nnmb + nr1 + nr2) {
      ExEnv::errn() << "assemble: dim mismatch" << endl;
      abort();
    }
  RefSCMatrix N(Nm.rowdim(), Nm.rowdim(), Nm.kit());
  // collect Nm, Nr1, and Nr2 back into N
  int i;
  for (i=0; i<nnmb; i++) {
      if (Nm_map[i] < 0 || Nm_map[i] >= nb) {
          ExEnv::errn() << "assemble: bad Nm_map" << endl;
          abort();
        }
      N.assign_column(Nm.get_column(i), Nm_map[i]);
    }
  for (i=0; i<nr1; i++) {
      if (Nr1_map[i] < 0 || Nr1_map[i] >= nb) {
          ExEnv::errn() << "assemble: bad Nr1_map" << endl;
          abort();
        }
      N.assign_column(Nr1.get_column(i), Nr1_map[i]);
    }
  for (i=0; i<nr2; i++) {
      if (Nr2_map[i] < 0 || Nr2_map[i] >= nb) {
          ExEnv::errn() << "assemble: bad Nr2_map" << endl;
          abort();
        }
      N.assign_column(Nr2.get_column(i), Nr2_map[i]);
    }
  return N;
}

// form symmetry average NAO for each atom
static void
form_nao(const RefSymmSCMatrix &P, const RefSymmSCMatrix &S,
         const RefSCMatrix &N, const RefDiagSCMatrix &W, int natom,
         int *maxam_on_atom, int **nam_on_atom, int ***amoff_on_atom,
         const Ref<SCMatrixKit>& kit)
{
  int i,j,k,l,m;

  N.assign(0.0);
  W.assign(0.0);

  for (i=0; i<natom; i++) {
      for (j=0; j<=maxam_on_atom[i]; j++) {
          int nfunc = 2*j + 1;
          double oonfunc = 1.0/nfunc;
          int nt = nam_on_atom[i][j];
          RefSCDimension tdim(new SCDimension(nt));
          RefSymmSCMatrix Pt(tdim, kit);
          RefSymmSCMatrix St(tdim, kit);
          Pt.assign(0.0);
          St.assign(0.0);
          for (k=0; k<nt; k++) {
              for (l=0; l<nt; l++) {
                  double Stmp = 0.0;
                  double Ptmp = 0.0;
                  for (m=0; m<nfunc; m++) {
                      int ii = amoff_on_atom[i][j][k] + m;
                      int jj = amoff_on_atom[i][j][l] + m;
                      Stmp += S.get_element(ii,jj);
                      Ptmp += P.get_element(ii,jj);
                    }
                  St.set_element(k,l,Stmp*oonfunc);
                  Pt.set_element(k,l,Ptmp*oonfunc);
                }
            }
          // find a symmetric orthogonalization transform
          RefSymmSCMatrix OL = mhalf(St);

          // transform Pt to the orthogonal basis
          RefSymmSCMatrix PtL(tdim,kit);
          PtL.assign(0.0);
          PtL.accumulate_transform(OL, Pt);

          // diagonalize PtL
          RefSCMatrix trans(tdim,tdim,kit);
          RefDiagSCMatrix eigval(tdim,kit);
          PtL.diagonalize(eigval, trans);

          // transform trans to the nonortho basis
          trans = OL * trans;

#         ifdef DEBUG
          eigval.print("eigval");
#         endif
          // fill in the elements of W
          for (k=0; k<nt; k++) {
              // the eigenvalues come out backwards: reverse them
              int krev = nt-k-1;
              double elem = eigval.get_element(krev);
              for (m=0; m<nfunc; m++) {
                  int ii = amoff_on_atom[i][j][k] + m;
#                 ifdef DEBUG
                  ExEnv::outn().form("W(%2d) = %12.8f\n", ii, elem);
#                 endif
                  W.set_element(ii, elem);
                }
            }

          // fill in the elements of N
          for (k=0; k<nt; k++) {
              for (l=0; l<nt; l++) {
                  // the eigenvalues come out backwards: reverse them
                  int lrev = nt-l-1;
                  double elem = trans.get_element(k,lrev);
                  for (m=0; m<nfunc; m++) {
                      int ii = amoff_on_atom[i][j][k] + m;
                      int jj = amoff_on_atom[i][j][l] + m;
                      N.set_element(ii,jj, elem);
                    }
                }
            }
        }
    }
}

// From "Natural Population Analysis", Alan E. Reed, Robert B. Weinstock,
// Frank Weinhold, JCP, 83 (1985), p 735.
RefSCMatrix
Wavefunction::nao(double *atom_charges)
{

  Ref<GaussianBasisSet> b = basis();
  Ref<PetiteList> pl = integral()->petite_list();

  // compute S, the ao basis overlap
  RefSymmSCMatrix blockedS = pl->to_AO_basis(overlap());
  RefSymmSCMatrix S
      = dynamic_cast<BlockedSymmSCMatrix*>(blockedS.pointer())->block(0);
  blockedS = 0;
# ifdef DEBUG
  S.print("S");
# endif

  // compute P, the ao basis density
  RefSymmSCMatrix P
      = dynamic_cast<BlockedSymmSCMatrix*>(ao_density().pointer())->block(0);

  // why?  good question.
  RefSymmSCMatrix Ptmp = P->clone();
  Ptmp.assign(0.0);
  Ptmp->accumulate_transform(S, P);
# ifdef DEBUG
  P.print("P");
  ExEnv::out0() << "nelec = " << (mhalf(S) * Ptmp * mhalf(S)).trace() << endl;
  ExEnv::out0() << "nelec(2) = " << (P * S).trace() << endl;
# endif
  P = Ptmp;
  Ptmp = 0;

  int i,j,k,l;
  int nb = b->nbasis();
  int nsh = b->nshell();
  int natom = molecule()->natom();

# ifdef DEBUG
  ExEnv::out0() << "nb = " << nb << endl;
  ExEnv::out0() << "nsh = " << nsh << endl;
  ExEnv::out0() << "natom = " << natom << endl;
# endif

  // Step 2a. Transform to solid harmonics.
  // -- for now program will abort if basis does not use only S.H and cart d.
  RefSCDimension aodim = P.dim();
  RefSCMatrix Tdfg(aodim, aodim, matrixkit());
  Tdfg->unit();
  for (i=0; i<nsh; i++) {
      const GaussianShell &shell = b->shell(i);
      int off = b->shell_to_function(i);
      for (j=0; j<shell.ncontraction(); j++) {
          if (shell.am(j) == 2 && ! shell.is_pure(j)) {
              for (k=0; k<6; k++) {
                  for (l=0; l<6; l++) {
                      Tdfg.set_element(off+k,off+l,0.0);
                    }
                }
              // this will put the s function first and the d second
              // first grab the s function
              SphericalTransformIter *sti;
              sti = integral()->new_spherical_transform_iter(2,0,0);
              for (sti->begin(); sti->ready(); sti->next()) {
                  Tdfg->set_element(off + sti->pureindex(),
                                    off + sti->cartindex(),
                                    sti->coef());
                }
              delete sti;
              // now for the pure d part of the cartesian d shell
              sti = integral()->new_spherical_transform_iter(2,0,2);
              for (sti->begin(); sti->ready(); sti->next()) {
                  Tdfg->set_element(off + sti->pureindex() + 1,
                                    off + sti->cartindex(),
                                    sti->coef());
                }
              delete sti;
            }
          else if (shell.am(j) > 2 && ! shell.is_pure(j)) {
              ExEnv::errn() << "NAOs can only be computed for puream if am > 2" << endl;
              abort();
            }
          off += shell.nfunction(j);
        }
    }

  // Tdfg should already be orthogonal, normalize them
//   RefSCMatrix Tdfgo = Tdfg*Tdfg.t();
//   RefDiagSCMatrix Tdfg_norm(Tdfg.rowdim(), matrixkit());
//   for (i=0; i<nb; i++) {
//       double o = Tdfgo.get_element(i,i);
//       Tdfg_norm.set_element(i,1.0/sqrt(o));
//     }
//   Tdfgo = 0;
//   Tdfg = Tdfg_norm * Tdfg;

# ifdef DEBUG
  Tdfg.print("Tdfg");
  (Tdfg.t() * Tdfg).print("Tdfg.t() * Tdfg");
  (Tdfg * Tdfg.t()).print("Tdfg * Tdfg.t()");
# endif

  RefSymmSCMatrix Pdfg(aodim, matrixkit());
  // Pdfp = Tdfp.t() * P * Tdfp
  Pdfg.assign(0.0); Pdfg.accumulate_transform(Tdfg, P);
  RefSymmSCMatrix Sdfg(aodim, matrixkit());
  // Sdfp = Tdfp.t() * S * Tdfp
  Sdfg.assign(0.0); Sdfg.accumulate_transform(Tdfg, S);
# ifdef DEBUG
  ExEnv::out0() << "nelec = " << (mhalf(Sdfg) * Pdfg * mhalf(Sdfg)).trace() << endl;
# endif

  // Step 2b. Partitioning and symmetry averaging of P and S
  // Partitioning:
  int *maxam_on_atom = new int[natom];
  int **nam_on_atom = new int*[natom];
  int ***amoff_on_atom = new int**[natom];
  int maxam = -1;
  for (i=0; i<natom; i++) {
      maxam_on_atom[i] = -1;
      for (j=0; j<b->nshell_on_center(i); j++) {
          GaussianShell &shell = b->shell(i,j);
          int maxam_on_shell = shell.max_angular_momentum();
          if (maxam_on_atom[i] < maxam_on_shell)
              maxam_on_atom[i] = maxam_on_shell;
        }
      if (maxam_on_atom[i] > maxam) maxam = maxam_on_atom[i];
      nam_on_atom[i] = new int[maxam_on_atom[i]+1];
      for (j=0; j<=maxam_on_atom[i]; j++) {
          nam_on_atom[i][j] = 0;
          for (k=0; k<b->nshell_on_center(i); k++) {
              GaussianShell &shell = b->shell(i,k);
              for (l=0; l<shell.ncontraction(); l++) {
                  if (shell.am(l) == j) nam_on_atom[i][j]++;
                  if (shell.am(l) == 2 && !shell.is_pure(l) && j == 0) {
                      // the s component of a cartesian d
                      nam_on_atom[i][0]++;
                    }
                }
            }
        }
      amoff_on_atom[i] = new int*[maxam_on_atom[i]+1];
      for (j=0; j<=maxam_on_atom[i]; j++) {
          amoff_on_atom[i][j] = new int[nam_on_atom[i][j]];
          int nam = 0;
          for (k=0; k<b->nshell_on_center(i); k++) {
              GaussianShell &shell = b->shell(i,k);
              int function_offset
                  = b->shell_to_function(b->shell_on_center(i,k));
              int conoffset = 0;
              for (l=0; l<shell.ncontraction(); l++) {
                  if (shell.am(l) == j) {
                      amoff_on_atom[i][j][nam]
                          = function_offset + conoffset;
                      if (j == 2 && !shell.is_pure(l)) {
                          // the pure d part of a cartesian d shell is offset
                          amoff_on_atom[i][j][nam]++;
                        }
                      nam++;
                    }
                  if (shell.am(l) == 2 && (!shell.is_pure(l)) && j == 0) {
                      // the s component of a cartesian d
                      amoff_on_atom[i][j][nam]
                          = function_offset + conoffset;
                      nam++;
                    }
                  conoffset += shell.nfunction(l);
                }
            }
        }
    }

# ifdef DEBUG
  ExEnv::out0() << indent << "Basis set partitioning:" << endl;
  ExEnv::out0() << incindent;
  for (i=0; i<natom; i++) {
      ExEnv::out0() << indent <<  "atom " << i
           << " maxam = " << maxam_on_atom[i] << endl;
      ExEnv::out0() << incindent;
      for (j=0; j<=maxam_on_atom[i]; j++) {
          ExEnv::out0() << indent <<  "am = " << j
               << " n = " << nam_on_atom[i][j] << endl;
          ExEnv::out0() << incindent;
          ExEnv::out0() << indent << "offsets =";
          for (k=0; k<nam_on_atom[i][j]; k++) {
              ExEnv::out0() << " " << amoff_on_atom[i][j][k];
            }
          ExEnv::out0() << endl;
          ExEnv::out0() << decindent;
        }
      ExEnv::out0() << decindent;
    }
  ExEnv::out0() << decindent;
# endif

  // Symmetry averaging and Step 2c: Formation of pre-NAO's
  RefSCMatrix N(aodim, aodim, matrixkit());
  RefDiagSCMatrix W(aodim, matrixkit());
  form_nao(Pdfg, Sdfg, N, W, natom, maxam_on_atom, nam_on_atom, amoff_on_atom,
           matrixkit());
# ifdef DEBUG
  N.print("N");
  W.print("W");
  ExEnv::out0() << "nelec = " << ttrace(N, Pdfg, Sdfg) << endl;
# endif

  // Step 3a: selection of NMB orbitals

  // count the size of nmb
  int nnmb = 0;
  for (i=0; i<natom; i++) {
      nnmb += nnmb_all_atom(molecule()->Z(i),
                            maxam_on_atom[i]);
    }
  int nnrb = nb - nnmb;

# ifdef DEBUG
  ExEnv::out0() << "nnmb = " << nnmb << endl;
  ExEnv::out0() << "nnrb = " << nnrb << endl;
# endif

  RefSCDimension nmbdim = new SCDimension(nnmb);
  RefSCDimension nrbdim = new SCDimension(nnrb);

  // split N into the nmb and nrb parts
  RefSCMatrix Nm(aodim, nmbdim, matrixkit());
  RefSCMatrix Nr(aodim, nrbdim, matrixkit());
  RefDiagSCMatrix Wm(nmbdim, matrixkit());
  RefDiagSCMatrix Wr(nrbdim, matrixkit());
  int *Nm_map = new int[nnmb];
  int *Nr_map = new int[nnrb];

  int im = 0;
  int ir = 0;
  for (i=0; i<natom; i++) {
      int z = molecule()->Z(i);
      for (j=0; j<=maxam_on_atom[i]; j++) {
          int nnmb_zj = nnmb_atom(z,j);
          int nt = nam_on_atom[i][j];
          for (k=0; k<nt; k++) {
              int iN = amoff_on_atom[i][j][k];
              if (k<nnmb_zj) {
                  for (l=0; l<(2*j+1); l++) {
                      Nm_map[im] = iN;
                      Wm.set_element(im, W.get_element(iN));
                      Nm.assign_column(N.get_column(iN++),im++);
                    }
                }
              else {
                  for (l=0; l<(2*j+1); l++) {
                      Nr_map[ir] = iN;
                      Wr.set_element(ir, W.get_element(iN));
                      Nr.assign_column(N.get_column(iN++),ir++);
                    }
                }
            }
        }
    }
# ifdef DEBUG
  ExEnv::out0() << "Nmmap:"; for (i=0;i<nnmb;i++) ExEnv::out0()<<" "<<Nm_map[i]; ExEnv::out0()<<endl;
  ExEnv::out0() << "Nrmap:"; for (i=0;i<nnrb;i++) ExEnv::out0()<<" "<<Nr_map[i]; ExEnv::out0()<<endl;
  Wm.print("Wm");
  Wr.print("Wr");
  Nm.print("Nm");
  Nr.print("Nr");
  (Nm.t() * Sdfg * Nr).print("3a Smr");
  ExEnv::out0() << "nelec = "
       << ttrace(assemble(aodim,Nm,Nm_map,Nr,Nr_map), Pdfg, Sdfg) << endl;
# endif

  // Step 3b: Schmidt interatomic orthogonalization of NRB to NMB orbs

  // orthogonalize the NMB orbs (temporarily, to project them out of NRB)
  int ii=0;
  for (i=0; i<nnmb; i++,ii++) {
      N.assign_column(Nm.get_column(i),ii);
    }
  for (i=0; i<nnrb; i++,ii++) {
      N.assign_column(Nr.get_column(i),ii);
    }
  N->schmidt_orthog(Sdfg.pointer(),nnmb);

  RefSCMatrix Nmo = Nm.clone();
  for (i=0; i<nnmb; i++) {
      Nmo.assign_column(N.get_column(i),i);
    }
  RefSCMatrix OSmr = Nmo.t() * Sdfg * Nr;
  OSmr.scale(-1.0);
  Nr.accumulate(Nmo * OSmr);
# ifdef DEBUG
  OSmr.print("OSmr");
  Nmo.print("Nmo = Nm after temporay orthog");
  Nr.print("Nr after orthogonalization to NMB");
  (Nm.t() * Sdfg * Nr).print("3b Smr");
# endif
  Nmo = 0;

  // Step 3c: Restoration of natural character of the NRB
  // Partitioning:
  int *r_maxam_on_atom = new int[natom];
  int **r_nam_on_atom = new int*[natom];
  int ***r_amoff_on_atom = new int**[natom];
  int r_offset = 0;
  for (i=0; i<natom; i++) {
      int z = molecule()->Z(i);
      r_maxam_on_atom[i] = maxam_on_atom[i];
      r_nam_on_atom[i] = new int[r_maxam_on_atom[i]+1];
      for (j=0; j<=r_maxam_on_atom[i]; j++) {
          r_nam_on_atom[i][j] = nam_on_atom[i][j] - nnmb_atom(z,j);
          if (r_nam_on_atom[i][j] < 0) {
              ExEnv::errn() << "NAO: < 0 rydberg orbitals of a given type" << endl;
              abort();
            }
        }
      r_amoff_on_atom[i] = new int*[r_maxam_on_atom[i]+1];
      for (j=0; j<=r_maxam_on_atom[i]; j++) {
          r_amoff_on_atom[i][j] = new int[r_nam_on_atom[i][j]];
          for (k=0; k<r_nam_on_atom[i][j]; k++) {
              r_amoff_on_atom[i][j][k] = r_offset;
              r_offset += 2*j + 1;
            }
        }
    }
  RefSymmSCMatrix Pr(nrbdim, matrixkit());
  // Pr = Nr.t() * Tdfg.t() * P * Tdfg * Nr;
  Pr.assign(0.0); Pr.accumulate_transform(Nr.t(), Pdfg);
  RefSymmSCMatrix Sr(nrbdim, matrixkit());
  // Sr = Nr.t() * Tdfg.t() * S * Tdfg * Nr;
  Sr.assign(0.0); Sr.accumulate_transform(Nr.t(), Sdfg);

  // Symmetry averaging and restoration of natural character of NRB
  RefSCMatrix Nrr(nrbdim, nrbdim, matrixkit());
  form_nao(Pr, Sr, Nrr, Wr,
           natom, r_maxam_on_atom, r_nam_on_atom, r_amoff_on_atom,
           matrixkit());
  Nr = Nr * Nrr;
  // these are out-of-date
  Pr = 0; Sr = 0;
# ifdef DEBUG
  Wr.print("Wr after restoring natural character");
  Nr.print("Nr after restoring natural character");
  (Nm.t() * Sdfg * Nr).print("3c Smr");
# endif

  // Step 4a: Weighted interatomic orthogonalization
  // nmb part of OW
  RefSymmSCMatrix Sm(nmbdim, matrixkit());
  Sm.assign(0.0); Sm.accumulate_transform(Nm.t(), Sdfg);
  RefSymmSCMatrix SWm = weight_matrix(Wm, Sm);
  RefSCMatrix OWm = Wm * mhalf(SWm);
# ifdef DEBUG
  Sm.print("Sm before 4a");
  OWm.print("OWm");
  (OWm.t() * Sm * OWm).print("Sm after 4a");
  ExEnv::out0() << "nelec = "
       << ttrace(assemble(aodim,Nm,Nm_map,Nr,Nr_map), Pdfg, Sdfg) << endl;
# endif

  // put OWm into Nm
  Nm = Nm * OWm;

# ifdef DEBUG
  Nm.print("Nm after interatomic orthog");
  (Nm.t() * Sdfg * Nr).print("4a Smr before r orthog");
  ExEnv::out0() << "nelec = "
       << ttrace(assemble(aodim,Nm,Nm_map,Nr,Nr_map), Pdfg, Sdfg)
       << endl;
# endif

  // nrb part of OW
  // based on Wr, r is split into r1 and r2
  double tw = 1.0e-4; // the tolerance used for the split
  int nr1 = 0;
  int nr2 = 0;
  for (i=0; i<nnrb; i++) {
      if (fabs(Wr.get_element(i)) >= tw) nr1++;
      else nr2++;
    }
  RefSCDimension r1dim(new SCDimension(nr1));
  RefSCDimension r2dim(new SCDimension(nr2));
  RefSCMatrix Nr1(aodim, r1dim, matrixkit());
  RefSCMatrix Nr2(aodim, r2dim, matrixkit());
  RefDiagSCMatrix Wr1(r1dim, matrixkit());
  int *Nr1_map = new int[nr1];
  int *Nr2_map = new int[nr2];
  int ir1 = 0;
  int ir2 = 0;
  for (i=0; i<nnrb; i++) {
      if (fabs(Wr.get_element(i)) >= tw) {
          Nr1_map[ir1] = Nr_map[i];
          Wr1.set_element(ir1, Wr.get_element(i));
          Nr1.assign_column(Nr.get_column(i),ir1++);
        }
      else {
          Nr2_map[ir2] = Nr_map[i];
          Nr2.assign_column(Nr.get_column(i),ir2++);
        }
    }
# ifdef DEBUG
  ExEnv::out0() << "Nr1map:"; for (i=0;i<nr1;i++) ExEnv::out0()<<" "<<Nr1_map[i]; ExEnv::out0()<<endl;
  ExEnv::out0() << "Nr2map:"; for (i=0;i<nr2;i++) ExEnv::out0()<<" "<<Nr2_map[i]; ExEnv::out0()<<endl;
  Nr1.print("Nr1");
  Nr2.print("Nr2");
  (Nm.t() * Sdfg * Nr1).print("4a Smr1 before r orthog");
  (Nm.t() * Sdfg * Nr2).print("4a Smr2 before r orthog");
  ExEnv::out0() << "nelec = "
       << ttrace(assemble(aodim,Nm,Nm_map,Nr1,Nr1_map,Nr2,Nr2_map), Pdfg, Sdfg)
       << endl;
# endif

  // Schmidt orthogonalization of r2 to r1
  // Collect Nr together again (but in the order: r1, r2)
  ii=0;
  for (i=0; i<nr1; i++,ii++) {
      Nr.assign_column(Nr1.get_column(i),ii);
    }
  for (i=0; i<nr2; i++,ii++) {
      Nr.assign_column(Nr2.get_column(i),ii);
    }
  Nr->schmidt_orthog(Sdfg.pointer(),nr1);
  RefSCMatrix Nr1o = Nr1.copy();
  for (i=0; i<nr1; i++) {
      Nr1o.assign_column(Nr.get_column(i), i);
    }
  RefSCMatrix Or1r2 = Nr1o.t() * Sdfg * Nr2;
  Or1r2.scale(-1.0);
# ifdef DEBUG
  (Nm.t() * Sdfg * Nr2).print("4a Smr2 before orthog of r2 to r1");
  (Nr1.t() * Sdfg * Nr2).print("4a Sr1r2 before orthog of r2 to r1");
# endif
  Nr2.accumulate(Nr1o * Or1r2);
# ifdef DEBUG
  Nr2.print("Nr2 after orthogonalization to r1");
  (Nm.t() * Sdfg * Nr2).print("4a Smr2 after orthog of r2 to r1");
  (Nr1.t() * Sdfg * Nr2).print("4a Sr1r2 after orthog of r2 to r1");
  ExEnv::out0() << "nelec = "
       << ttrace(assemble(aodim,Nm,Nm_map,Nr1,Nr1_map,Nr2,Nr2_map), Pdfg, Sdfg)
       << endl;
# endif

  // weighted symmetric orthog of r1
  RefSymmSCMatrix Sr1(r1dim, matrixkit());
  Sr1.assign(0.0); Sr1.accumulate_transform(Nr1.t(), Sdfg);
  RefSymmSCMatrix SWr1 = weight_matrix(Wr1, Sr1);
  RefSCMatrix OWr1 = Wr1 * mhalf(SWr1);
# ifdef DEBUG
  OWr1.print("OWr1");
  (Nr1.t() * Sdfg * Nr1).print("Nr1.t() * Sdfg * Nr1");
# endif
  // Put OWr1 into Nr1
  Nr1 = Nr1 * OWr1;
# ifdef DEBUG
  Nr1.print("Nr1 after weighted symmetric orthogonalization");
  ExEnv::out0() << "nelec = "
       << ttrace(assemble(aodim,Nm,Nm_map,Nr1,Nr1_map,Nr2,Nr2_map), Pdfg, Sdfg)
       << endl;
# endif

  // symmetric orthog of r1
  RefSymmSCMatrix Sr2(r2dim, matrixkit());
  Sr2.assign(0.0); Sr2.accumulate_transform(Nr2.t(), Sdfg);
  RefSymmSCMatrix OWr2 = mhalf(Sr2);
# ifdef DEBUG
  OWr2.print("OWr2");
# endif

  // Put OWr2 into Nr2
  Nr2 = Nr2 * OWr2;
# ifdef DEBUG
  Nr2.print("Nr2 after weighted symmetric orthogonalization");
  ExEnv::out0() << "nelec = "
       << ttrace(assemble(aodim,Nm,Nm_map,Nr1,Nr1_map,Nr2,Nr2_map), Pdfg, Sdfg)
       << endl;
  ExEnv::out0() << "nelec(o) = "
       << ttrace(assemble(aodim,Nm,Nm_map,Nr1,Nr1_map,Nr2,Nr2_map), Pdfg)
       << endl;
# endif

  // Step 4b. restoration of the natural character of the naos

  // collect Nm, Nr1, and Nr2 back into N
  N = assemble(aodim, Nm,Nm_map, Nr1,Nr1_map, Nr2,Nr2_map);
# ifdef DEBUG
  N.print("N after 4a");
# endif

  // compute the density and overlap in the current basis
  // N currently has the entire transform, starting from the dfg basis
  P.assign(0.0);
  P.accumulate_transform(N.t(), Pdfg);
  S.assign(0.0);
  S.accumulate_transform(N.t(), Sdfg);
# ifdef DEBUG
  P.print("P after 4a");
  S.print("S after 4a");
  (Nm.t() * Sdfg * Nm).print("4a Sm");
  (Nr1.t() * Sdfg * Nr1).print("4a Sr1");
  (Nr2.t() * Sdfg * Nr2).print("4a Sr2");
  (Nm.t() * Sdfg * Nr1).print("4a Smr1");
  (Nm.t() * Sdfg * Nr2).print("4a Smr2");
  (Nr1.t() * Sdfg * Nr2).print("4a Sr1r2");
# endif

  RefSCMatrix Nred(aodim, aodim, matrixkit());
  form_nao(P, S, Nred, W, natom, maxam_on_atom, nam_on_atom, amoff_on_atom,
           matrixkit());
  N = N * Nred;

  RefSymmSCMatrix Pfinal(aodim, matrixkit());
  Pfinal.assign(0.0);
  Pfinal.accumulate_transform(N.t(), Pdfg);
# ifdef DEBUG
  Nred.print("Nred");
  N.print("N after 4b");
  ExEnv::out0() << "nelec = " << ttrace(N, Pdfg, Sdfg) << endl;
  ExEnv::out0() << "nelec(o) = " << ttrace(N, Pdfg) << endl;
  Pfinal.print("final P");
  (N.t() * Sdfg * N).print("final S");
  ExEnv::out0().form("nelec = trace(final P) = %14.8f", (N.t() * Pdfg * N).trace());

  (mhalf(Sdfg) * Pdfg * mhalf(Sdfg)).print("P in symm orth basis");
# endif

# ifdef DEBUG
  ExEnv::out0() << "nb   = " << nb << endl;
  ExEnv::out0() << "nnmb = " << nnmb << endl;
  ExEnv::out0() << "nnrb = " << nnrb << endl;
  ExEnv::out0() << "nr1  = " << nr1 << endl;
  ExEnv::out0() << "nr2  = " << nr2 << endl;
# endif

  ExEnv::out0() << indent << "Natural Population Analysis:" << endl;
  ExEnv::out0() << incindent;
  ExEnv::out0() << indent << " n   atom    charge ";
  for (i=0; i<=maxam; i++) {
      const char *am = "SPDFGH?";
      int index;
      if (i>6) index = 6;
      else index = i;
      ExEnv::out0() << "    ne(" << am[index] << ") ";
    }
  ExEnv::out0() << endl;
  for (i=0; i<natom; i++) {
      double e = 0.0;
      for (j=0; j<=maxam_on_atom[i]; j++) {
          for (k=0; k<nam_on_atom[i][j]; k++) {
              for (l=0; l<(2*j+1); l++) {
                  e += Pfinal.get_element(amoff_on_atom[i][j][k] + l,
                                          amoff_on_atom[i][j][k] + l);
                }
            }
        }
      std::string symbol(molecule()->atom_symbol(i));
      ExEnv::out0() << indent
           << scprintf("%3d   %2s   % 8.6f",i + 1,
                       symbol.c_str(),
                       double(molecule()->Z(i)) - e);
      if (atom_charges) {
          atom_charges[i] = molecule()->Z(i) - e;
        }
      for (j=0; j<=maxam_on_atom[i]; j++) {
          e = 0.0;
          for (k=0; k<nam_on_atom[i][j]; k++) {
              for (l=0; l<(2*j+1); l++) {
                  e += Pfinal.get_element(amoff_on_atom[i][j][k] + l,
                                          amoff_on_atom[i][j][k] + l);
                }
            }
          ExEnv::out0() << scprintf(" % 8.6f",e);
        }
      ExEnv::out0() << endl;
    }
  ExEnv::out0() << endl;
  ExEnv::out0() << decindent;

  delete[] Nm_map;
  delete[] Nr_map;
  delete[] Nr1_map;
  delete[] Nr2_map;
  delete_partition_info(natom,maxam_on_atom,nam_on_atom,amoff_on_atom);
  delete_partition_info(natom,r_maxam_on_atom,r_nam_on_atom,r_amoff_on_atom);

  return N;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
