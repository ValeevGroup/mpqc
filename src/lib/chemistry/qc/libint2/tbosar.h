//
// tbosar.h
//
// Copyright (C) 2001 Edward Valeev
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

#include <limits.h>
#include <stdexcept>

#include <util/ref/ref.h>
#include <util/misc/scexception.h>
#include <util/misc/consumableresources.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/libint2/shellpairs.h>
#include <chemistry/qc/basis/fjt.h>
#include <chemistry/qc/libint2/int2e.h>
#include <chemistry/qc/libint2/macros.h>
#include <chemistry/qc/libint2/libint2_utils.h>
#include <libint2/libint2.h>
#include <libint2/boys.h>

#if LIBINT2_SUPPORT_ERI

#ifndef _chemistry_qc_libint2_tbosar_h
#define _chemistry_qc_libint2_tbosar_h

namespace sc {

  namespace detail {
    template <TwoBodyOper::type OperType>
    struct OSAR_CoreInts;

    template <>
    struct OSAR_CoreInts<TwoBodyOper::eri> {
	// Line below isn't legal
        //const static double small_T = 1E-15;       /*--- Use only one term in Taylor expansion of Fj(T) if T < small_T ---*/
        ::libint2::FmEval_Chebyshev3 Fm_Eval_;

        OSAR_CoreInts(unsigned int mmax, const Ref<IntParams>& params) :   Fm_Eval_(mmax) {}

        const double* eval(double* Fm_table, unsigned int mmax, double T, double rho = 0.0) const {
          static double oo2np1[] = {1.0,  1.0/3.0,  1.0/5.0,  1.0/7.0,  1.0/9.0,
            1.0/11.0, 1.0/13.0, 1.0/15.0, 1.0/17.0, 1.0/19.0,
            1.0/21.0, 1.0/23.0, 1.0/25.0, 1.0/27.0, 1.0/29.0,
            1.0/31.0, 1.0/33.0, 1.0/35.0, 1.0/37.0, 1.0/39.0,
            1.0/41.0, 1.0/43.0, 1.0/45.0, 1.0/47.0, 1.0/49.0,
            1.0/51.0, 1.0/53.0, 1.0/55.0, 1.0/57.0, 1.0/59.0,
            1.0/61.0, 1.0/63.0, 1.0/65.0, 1.0/67.0, 1.0/69.0,
            1.0/71.0, 1.0/73.0, 1.0/75.0, 1.0/77.0, 1.0/79.0};
	  
	  const static double small_T = 1E-15;
          if(T < small_T){
            return oo2np1;
          }
          else {
            Fm_Eval_.eval(Fm_table, T, mmax);
            return Fm_table;
          }
        }
    };

    namespace {
      // this does common work for geminal-depenent kernels
      struct OSAR_CoreInts_G12Base {
          IntParamsG12::ContractedGeminal g12_;

          OSAR_CoreInts_G12Base(const Ref<IntParams>& p) {

            Ref<IntParamsG12> params = require_dynamic_cast<IntParamsG12*>(p, "");

            const bool braonly = (params->ket() == IntParamsG12::null_geminal);
            g12_ = braonly ? params->bra() : IntParamsG12::product(params->bra(), params->ket());
          }

          struct softcomparer {
              softcomparer(double tol = DBL_EPSILON) : tol_(tol) {}
              bool operator()(double a, double b) const { return a + tol_ < b; }
              double tol_;
          };

          // reduce terms with common exponents
          static IntParamsG12::ContractedGeminal
          reduce(const IntParamsG12::ContractedGeminal& g12) {
            softcomparer comp;
            std::map<double, double, softcomparer> g12red(comp);
            typedef std::map<double, double, softcomparer>::iterator iter;

            const size_t np = g12.size();
            for(size_t p=0; p<np; ++p) {
              iter i = g12red.find(g12[p].first);
              if (i != g12red.end()) {
                i->second += g12[p].second;
              }
              else
                g12red[g12[p].first] = g12[p].second;
            }

            IntParamsG12::ContractedGeminal result;
            for(iter i=g12red.begin();
                i != g12red.end();
                ++i) {
              result.push_back(*i);
            }

            return result;
          }
      };
    } // namespace <anonymous>

    template <>
    struct OSAR_CoreInts<TwoBodyOper::r12_0_g12> : OSAR_CoreInts_G12Base {
        ::libint2::GaussianGmEval<double, 0> Gm_Eval_;
        OSAR_CoreInts(unsigned int mmax, const Ref<IntParams>& params) :
          Gm_Eval_(mmax, 1e-14), OSAR_CoreInts_G12Base(params)
          {
            g12_ = reduce(g12_);
          }
        double* eval(double* Fm_table, unsigned int mmax, double T, double rho = 0.0) {
          Gm_Eval_.eval(Fm_table, rho, T, mmax, g12_);
          return Fm_table;
        }
    };
    template <>
    struct OSAR_CoreInts<TwoBodyOper::r12_m1_g12> : OSAR_CoreInts_G12Base {
        ::libint2::GaussianGmEval<double, -1> Gm_Eval_;
        OSAR_CoreInts(unsigned int mmax, const Ref<IntParams>& params) :
          Gm_Eval_(mmax, 1e-14), OSAR_CoreInts_G12Base(params)
        {
          g12_ = reduce(g12_);
        }
        double* eval(double* Fm_table, unsigned int mmax, double T, double rho = 0.0) {
          Gm_Eval_.eval(Fm_table, rho, T, mmax, g12_);
          return Fm_table;
        }
    };
    template <>
    struct OSAR_CoreInts<TwoBodyOper::g12t1g12> : OSAR_CoreInts_G12Base {
        ::libint2::GaussianGmEval<double, 2> Gm_Eval_;
        OSAR_CoreInts(unsigned int mmax, const Ref<IntParams>& params) :
          Gm_Eval_(mmax, 1e-14), OSAR_CoreInts_G12Base(params)
        {
          // [exp(-a r_{12}^2),[T1,exp(-b r_{12}^2)]] = 4 a b (r_{12}^2 exp(- (a+b) r_{12}^2) )
          // i.e. need to scale each coefficient by 4 a b
          Ref<IntParamsG12> p = require_dynamic_cast<IntParamsG12*>(params, "");
          const IntParamsG12::ContractedGeminal& gbra = p->bra();
          const IntParamsG12::ContractedGeminal& gket = p->ket();
          for(size_t b=0, bk=0; b<gbra.size(); ++b)
            for(size_t k=0; k<gket.size(); ++k, ++bk)
              g12_[bk].second *= 4.0 * gbra[b].first * gket[k].first;
          g12_ = reduce(g12_);
        }
        double* eval(double* Fm_table, unsigned int mmax, double T, double rho = 0.0) {
          Gm_Eval_.eval(Fm_table, rho, T, mmax, g12_);
          return Fm_table;
        }
    };

    template <>
    struct OSAR_CoreInts<TwoBodyOper::delta> {
        OSAR_CoreInts(unsigned int mmax, const Ref<IntParams>& params) {}
        const double* eval(double* Fm_table, unsigned int mmax, double T, double rho = 0.0) const {
          const static double one_over_two_pi = 1.0 / (2.0 * M_PI);
          const double G0 = exp(-T) * rho * one_over_two_pi;
          //const double G0 = 0.5 * sqrt(M_PI / rho);   // <- 1 (=> product of 1-e overlaps), instead of delta function
          std::fill(Fm_table, Fm_table+mmax+1, G0);
          return Fm_table;
        }
    };
  }

  namespace {

    inline void
    swtch(GaussianBasisSet* &i,GaussianBasisSet* &j)
    {
      GaussianBasisSet *tmp;
      tmp = i;
      i = j;
      j = tmp;
    }

    inline void
    pswtch(void**i,void**j)
    {
      void*tmp;
      tmp = *i;
      *i = *j;
      *j = tmp;
    }

    inline void
    iswtch(int *i,int *j)
    {
      int tmp;
      tmp = *i;
      *i = *j;
      *j = tmp;
    }

  }; // end of namespace <anonymous>

/** TwoBodyOSARLibint2 is a specialization of Int2eLibint2 that computes integrals using general Obara-Saika-Ahlrichs recurrence */
template <TwoBodyOper::type OperType>
class TwoBodyOSARLibint2: public Int2eLibint2 {
  private:

    static bool store_pair_data() { return true; } // hardwired for now

    // Storage for target integrals
    double *target_ints_buffer_;

    /*--- Intermediate scratch arrays (may be used in new[] and delete[]) ---*/
    double *cart_ints_;       // cartesian integrals, in by-contraction-quartet order
    double *sphharm_ints_;    // transformed integrals, in by-contraction-quartet order
    double *perm_ints_;       // redundant target integrals in shell quartet order, shells permuted

    /*--- Pointers to scratch arrays (never used in new[] and delete[]) ---*/
    double *prim_ints_;       // this points to the appropriate location for raw integrals
    double *contr_quartets_;
    double *shell_quartet_;

    /*--- Precomputed data ---*/
    Ref<ShellPairsLibint2> shell_pairs12_;
    Ref<ShellPairsLibint2> shell_pairs34_;

    /*--- Internally used "interfaces" ---*/
    struct {
      int p12, p34, p13p24;           // flags indicating if functions were permuted
      ShellPairLibint2 *shell_pair12, *shell_pair34;   // Shell pairs corresponding to the original
                                                     // (before permutation) order of shell
      int *op1, *op2, *op3, *op4;     // pointers to the primitive indices in the original order
      /////////// The rest of data has been permuted according to p12, p34, p13p24
      double A[3], B[3], C[3], D[3];
      double AB2, CD2;
      int gc1, gc2, gc3, gc4;
      int p1, p2, p3, p4;
      int am;
    } quartet_info_;
    typedef Libint_t prim_data;
    void quartet_data_(prim_data *Data, double scale);
    /*--- Compute engines ---*/
    std::vector<Libint_t> Libint_;
    double* Fm_table_;
    detail::OSAR_CoreInts<OperType> coreints_;
  
  public:
    TwoBodyOSARLibint2(Integral *,
	     const Ref<GaussianBasisSet>&,
	     const Ref<GaussianBasisSet>&,
	     const Ref<GaussianBasisSet>&,
	     const Ref<GaussianBasisSet>&,
	     size_t storage,
	     const Ref<IntParams>& oper_params);
    ~TwoBodyOSARLibint2();

    double *buffer(unsigned int t = 0) const {
      return target_ints_buffer_;
    }

    static size_t storage_required(const Ref<GaussianBasisSet>& b1,
				   const Ref<GaussianBasisSet>& b2 = 0,
				   const Ref<GaussianBasisSet>& b3 = 0,
				   const Ref<GaussianBasisSet>& b4 = 0);
    
    // evaluate integrals
    void compute_quartet(int*, int*, int*, int*);
    
};

template <TwoBodyOper::type OperType>
TwoBodyOSARLibint2<OperType>::TwoBodyOSARLibint2(Integral *integral,
           const Ref<GaussianBasisSet>& b1,
           const Ref<GaussianBasisSet>& b2,
           const Ref<GaussianBasisSet>& b3,
           const Ref<GaussianBasisSet>& b4,
           size_t storage,
           const Ref<IntParams>& oper_params) :
  Int2eLibint2(integral,b1,b2,b3,b4,storage),
  coreints_(b1->max_angular_momentum() +
          b2->max_angular_momentum() +
          b3->max_angular_momentum() +
          b4->max_angular_momentum(),
          oper_params)
{
  // The static part of Libint's interface is automatically initialized in libint.cc
  int l1 = bs1_->max_angular_momentum();
  int l2 = bs2_->max_angular_momentum();
  int l3 = bs3_->max_angular_momentum();
  int l4 = bs4_->max_angular_momentum();
  int lmax = std::max(std::max(l1,l2),std::max(l3,l4));
  if (lmax > LIBINT2_MAX_AM_ERI) {
    throw LimitExceeded<int>("TwoBodyOSARLibint2::TwoBodyOSARLibint2() -- maxam of the basis is too high,\
 not supported by this libint2 library. Recompile libint2.",__FILE__,__LINE__,LIBINT2_MAX_AM_ERI,lmax);
  }

  /*--- Initialize storage ---*/
  const int max_num_prim_comb = bs1_->max_nprimitive_in_shell()*
    bs2_->max_nprimitive_in_shell()*
    bs3_->max_nprimitive_in_shell()*
    bs4_->max_nprimitive_in_shell();
  // need one Libint_t object for each primitive combination
  // if Libint2 does not support contractions, just allocate 1
#if LIBINT2_CONTRACTED_INTS
  Libint_.resize(max_num_prim_comb);
#else
  Libint_.resize(1);
#endif
  ConsumableResources::get_default_instance()->consume_memory(Libint_.size() * sizeof(Libint_[0]));

  const int max_cart_target_size = bs1_->max_ncartesian_in_shell()*bs2_->max_ncartesian_in_shell()*
    bs3_->max_ncartesian_in_shell()*bs4_->max_ncartesian_in_shell();
  const int max_target_size = bs1_->max_nfunction_in_shell()*bs2_->max_nfunction_in_shell()*
    bs3_->max_nfunction_in_shell()*bs4_->max_nfunction_in_shell();

  size_t storage_needed = LIBINT2_PREFIXED_NAME(libint2_need_memory_eri)(lmax) * sizeof(LIBINT2_REALTYPE);
  LIBINT2_PREFIXED_NAME(libint2_init_eri)(&Libint_[0],lmax,0);  // only need to initialize stack of the first Libint_t object
  manage_array(Libint_[0].stack, storage_needed/sizeof(LIBINT2_REALTYPE));

  target_ints_buffer_ = allocate<double>(max_target_size);
  cart_ints_ = allocate<double>(max_cart_target_size);
  if (bs1_->has_pure() || bs2_->has_pure() || bs3_->has_pure() || bs4_->has_pure() ||
      bs1_->max_ncontraction() != 1 || bs2_->max_ncontraction() != 1 ||
      bs3_->max_ncontraction() != 1 || bs4_->max_ncontraction() != 1) {
    sphharm_ints_ = allocate<double>(max_target_size);
    storage_needed += max_target_size*sizeof(double);
  }
  else {
    sphharm_ints_ = 0;
  }
  if (l1 || l2 || l3 || l4) {
    perm_ints_ = allocate<double>(max_target_size);
    storage_needed += max_target_size*sizeof(double);
  }
  else
    perm_ints_ = 0;

  // See if can store primitive-pair data
  size_t primitive_pair_storage_estimate = (bs1_->nprimitive()*bs2_->nprimitive() +
    bs3_->nprimitive()*bs4_->nprimitive())*sizeof(prim_pair_t);
  //  ExEnv::errn() << scprintf("need %d bytes to store primitive pair data\n",primitive_pair_storage_estimate);

  if (store_pair_data()) {
    shell_pairs12_ = new ShellPairsLibint2(bs1_,bs2_);
    if ( (bs1_ == bs3_ && bs2_ == bs4_) /*||
             // if this is (ab|ba) case -- should i try to save storage?
         (bs1_ == bs4_ && bs2_ == bs3_)*/ )
      shell_pairs34_ = new ShellPairsLibint2(shell_pairs12_);
    else
      shell_pairs34_ = new ShellPairsLibint2(bs3_,bs4_);
    storage_needed += primitive_pair_storage_estimate;
  }

  storage_used_ = storage_needed;
  // Check if storage_ > storage_needed
  check_storage_();

  int mmax = bs1_->max_angular_momentum() +
    bs2_->max_angular_momentum() +
    bs3_->max_angular_momentum() +
    bs4_->max_angular_momentum();
  Fm_table_ = new double[mmax+1];
}

template <TwoBodyOper::type OperType>
TwoBodyOSARLibint2<OperType>::~TwoBodyOSARLibint2()
{
  unmanage_array(Libint_[0].stack);
  LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&Libint_[0]);
  Libint_[0].stack = 0;
  ConsumableResources::get_default_instance()->release_memory(Libint_.size() * sizeof(Libint_[0]));
  deallocate(target_ints_buffer_);
  deallocate(cart_ints_);
  if (sphharm_ints_)
    deallocate(sphharm_ints_);
  if (perm_ints_)
    deallocate(perm_ints_);
#ifdef DMALLOC
  dmalloc_shutdown();
#endif
  delete[] Fm_table_;
}

template <TwoBodyOper::type OperType>
size_t
TwoBodyOSARLibint2<OperType>::storage_required(const Ref<GaussianBasisSet>& b1,
               const Ref<GaussianBasisSet>& b2,
               const Ref<GaussianBasisSet>& b3,
               const Ref<GaussianBasisSet>& b4)
{
  Ref<GaussianBasisSet> bs1 = b1;
  Ref<GaussianBasisSet> bs2 = b2;
  Ref<GaussianBasisSet> bs3 = b3;
  Ref<GaussianBasisSet> bs4 = b4;

  if (bs2 == 0)
    bs2 = bs1;
  if (bs3 == 0)
    bs3 = bs1;
  if (bs4 == 0)
    bs4 = bs1;

  int l1 = bs1->max_angular_momentum();
  int l2 = bs2->max_angular_momentum();
  int l3 = bs3->max_angular_momentum();
  int l4 = bs4->max_angular_momentum();
  int lmax = std::max(std::max(l1,l2),std::max(l3,l4));

  size_t storage_required = storage_required_(bs1,bs2,bs3,bs4);

  const int max_num_prim_comb = bs1->max_nprimitive_in_shell()*
    bs2->max_nprimitive_in_shell()*
    bs3->max_nprimitive_in_shell()*
    bs4->max_nprimitive_in_shell();
#if LIBINT2_CONTRACTED_INTS
  storage_required += max_num_prim_comb * sizeof(Libint_t);
#else
  storage_required += sizeof(Libint_t);
#endif

  const int max_cart_target_size = bs1->max_ncartesian_in_shell()*bs2->max_ncartesian_in_shell()*
    bs3->max_ncartesian_in_shell()*bs4->max_ncartesian_in_shell();
  const int max_target_size = bs1->max_nfunction_in_shell()*bs2->max_nfunction_in_shell()*
    bs3->max_nfunction_in_shell()*bs4->max_nfunction_in_shell();

  storage_required += LIBINT2_PREFIXED_NAME(libint2_need_memory_eri)(lmax) * sizeof(LIBINT2_REALTYPE);

  if (bs1->has_pure() || bs2->has_pure() || bs3->has_pure() || bs4->has_pure() ||
      bs1->max_ncontraction() != 1 || bs2->max_ncontraction() != 1 ||
      bs3->max_ncontraction() != 1 || bs4->max_ncontraction() != 1) {
    storage_required += max_target_size*sizeof(double);
  }

  if (l1 || l2 || l3 || l4) {
    storage_required += max_target_size*sizeof(double);
  }

  // See if can store primitive-pair data
  size_t primitive_pair_storage_estimate = (bs1->nprimitive()*bs2->nprimitive() +
    bs3->nprimitive()*bs4->nprimitive())*sizeof(prim_pair_t);
#if STORE_PAIR_DATA
  storage_required += primitive_pair_storage_estimate;
#endif

  return storage_required;
}

template <TwoBodyOper::type OperType>
void
TwoBodyOSARLibint2<OperType>::quartet_data_(prim_data *Data, double scale)
{

  /*----------------
    Local variables
   ----------------*/
  double P[3], Q[3], PQ[3], W[3];

  int p1 = quartet_info_.p1;
  int p2 = quartet_info_.p2;
  int p3 = quartet_info_.p3;
  int p4 = quartet_info_.p4;

  double a1 = int_shell1_->exponent(quartet_info_.p1);
  double a2 = int_shell2_->exponent(quartet_info_.p2);
  double a3 = int_shell3_->exponent(quartet_info_.p3);
  double a4 = int_shell4_->exponent(quartet_info_.p4);

  prim_pair_t* pair12;
  prim_pair_t* pair34;
  if (!quartet_info_.p13p24) {
    pair12 = quartet_info_.shell_pair12->prim_pair(*quartet_info_.op1,*quartet_info_.op2);
    pair34 = quartet_info_.shell_pair34->prim_pair(*quartet_info_.op3,*quartet_info_.op4);
  }
  else {
    pair12 = quartet_info_.shell_pair34->prim_pair(*quartet_info_.op3,*quartet_info_.op4);
    pair34 = quartet_info_.shell_pair12->prim_pair(*quartet_info_.op1,*quartet_info_.op2);
  }

  double zeta = pair12->gamma;
  double eta = pair34->gamma;
  double ooz = 1.0/zeta;
  double ooe = 1.0/eta;
  double ooze = 1.0/(zeta+eta);
#if LIBINT2_DEFINED(eri,roz)
  Data->roz[0] = eta*ooze;
  double rho = zeta*Data->roz[0];
#else
  double rho = zeta * eta * ooze;
#endif

  double pfac_norm = int_shell1_->coefficient_unnorm(quartet_info_.gc1,p1)*
  int_shell2_->coefficient_unnorm(quartet_info_.gc2,p2)*
  int_shell3_->coefficient_unnorm(quartet_info_.gc3,p3)*
  int_shell4_->coefficient_unnorm(quartet_info_.gc4,p4);
  double pfac = 2.0*sqrt(rho*M_1_PI)*scale*pair12->ovlp*pair34->ovlp*pfac_norm;

  P[0] = pair12->P[0];
  P[1] = pair12->P[1];
  P[2] = pair12->P[2];
  Q[0] = pair34->P[0];
  Q[1] = pair34->P[1];
  Q[2] = pair34->P[2];
  PQ[0] = P[0] - Q[0];
  PQ[1] = P[1] - Q[1];
  PQ[2] = P[2] - Q[2];
  double PQ2 = PQ[0]*PQ[0];
  PQ2 += PQ[1]*PQ[1];
  PQ2 += PQ[2]*PQ[2];
  double T = rho*PQ2;

  if (!quartet_info_.am) {
    const double* Fm = coreints_.eval(Fm_table_, 0, T, rho);
    Data->LIBINT_T_SS_EREP_SS(0)[0] = Fm[0]*pfac;
  }
  else {
#if LIBINT2_DEFINED(eri,oo2ze)
    Data->oo2ze[0] = 0.5*ooze;
#endif
#if LIBINT2_DEFINED(eri,roe)
    Data->roe[0] = zeta*ooze;
#endif
#if LIBINT2_DEFINED(eri,oo2z)
    Data->oo2z[0] = 0.5/zeta;
#endif
#if LIBINT2_DEFINED(eri,oo2e)
    Data->oo2e[0] = 0.5/eta;
#endif
    W[0] = (zeta*P[0] + eta*Q[0])*ooze;
    W[1] = (zeta*P[1] + eta*Q[1])*ooze;
    W[2] = (zeta*P[2] + eta*Q[2])*ooze;

    const double* Gm = coreints_.eval(Fm_table_, quartet_info_.am, T, rho);
    assign_FjT(Data,quartet_info_.am,Gm,pfac);

    /* PA */
#if LIBINT2_DEFINED(eri,PA_x)
    Data->PA_x[0] = P[0] - quartet_info_.A[0];
#endif
#if LIBINT2_DEFINED(eri,PA_y)
    Data->PA_y[0] = P[1] - quartet_info_.A[1];
#endif
#if LIBINT2_DEFINED(eri,PA_z)
    Data->PA_z[0] = P[2] - quartet_info_.A[2];
#endif
    /* QC */
#if LIBINT2_DEFINED(eri,QC_x)
    Data->QC_x[0] = Q[0] - quartet_info_.C[0];
#endif
#if LIBINT2_DEFINED(eri,QC_y)
    Data->QC_y[0] = Q[1] - quartet_info_.C[1];
#endif
#if LIBINT2_DEFINED(eri,QC_z)
    Data->QC_z[0] = Q[2] - quartet_info_.C[2];
#endif
    /* WP */
#if LIBINT2_DEFINED(eri,WP_x)
    Data->WP_x[0] = W[0] - P[0];
#endif
#if LIBINT2_DEFINED(eri,WP_y)
    Data->WP_y[0] = W[1] - P[1];
#endif
#if LIBINT2_DEFINED(eri,WP_z)
    Data->WP_z[0] = W[2] - P[2];
#endif
    /* WQ */
#if LIBINT2_DEFINED(eri,WQ_x)
    Data->WQ_x[0] = W[0] - Q[0];
#endif
#if LIBINT2_DEFINED(eri,WQ_y)
    Data->WQ_y[0] = W[1] - Q[1];
#endif
#if LIBINT2_DEFINED(eri,WQ_z)
    Data->WQ_z[0] = W[2] - Q[2];
#endif

    // using ITR?
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_x)
    Data->TwoPRepITR_pfac0_0_0_x[0] = - (a2*(quartet_info_.A[0]-quartet_info_.B[0]) + a4*(quartet_info_.C[0]-quartet_info_.D[0]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_y)
    Data->TwoPRepITR_pfac0_0_0_y[0] = - (a2*(quartet_info_.A[1]-quartet_info_.B[1]) + a4*(quartet_info_.C[1]-quartet_info_.D[1]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_z)
    Data->TwoPRepITR_pfac0_0_0_z[0] = - (a2*(quartet_info_.A[2]-quartet_info_.B[2]) + a4*(quartet_info_.C[2]-quartet_info_.D[2]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_x)
    Data->TwoPRepITR_pfac0_1_0_x[0] = - (a2*(quartet_info_.A[0]-quartet_info_.B[0]) + a4*(quartet_info_.C[0]-quartet_info_.D[0]))/eta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_y)
    Data->TwoPRepITR_pfac0_1_0_y[0] = - (a2*(quartet_info_.A[1]-quartet_info_.B[1]) + a4*(quartet_info_.C[1]-quartet_info_.D[1]))/eta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_z)
    Data->TwoPRepITR_pfac0_1_0_z[0] = - (a2*(quartet_info_.A[2]-quartet_info_.B[2]) + a4*(quartet_info_.C[2]-quartet_info_.D[2]))/eta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_x)
    Data->TwoPRepITR_pfac0_0_1_x[0] = (a1*(quartet_info_.A[0]-quartet_info_.B[0]) + a3*(quartet_info_.C[0]-quartet_info_.D[0]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_y)
    Data->TwoPRepITR_pfac0_0_1_y[0] = (a1*(quartet_info_.A[1]-quartet_info_.B[1]) + a3*(quartet_info_.C[1]-quartet_info_.D[1]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_z)
    Data->TwoPRepITR_pfac0_0_1_z[0] = (a1*(quartet_info_.A[2]-quartet_info_.B[2]) + a3*(quartet_info_.C[2]-quartet_info_.D[2]))/zeta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_x)
    Data->TwoPRepITR_pfac0_1_1_x[0] = (a1*(quartet_info_.A[0]-quartet_info_.B[0]) + a3*(quartet_info_.C[0]-quartet_info_.D[0]))/eta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_y)
    Data->TwoPRepITR_pfac0_1_1_y[0] = (a1*(quartet_info_.A[1]-quartet_info_.B[1]) + a3*(quartet_info_.C[1]-quartet_info_.D[1]))/eta;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_z)
    Data->TwoPRepITR_pfac0_1_1_z[0] = (a1*(quartet_info_.A[2]-quartet_info_.B[2]) + a3*(quartet_info_.C[2]-quartet_info_.D[2]))/eta;
#endif
#if LIBINT2_DEFINED(eri,eoz)
    Data->eoz[0] = eta * ooz;
#endif
#if LIBINT2_DEFINED(eri,zoe)
    Data->zoe[0] = zeta * ooe;
#endif
  }

  return;
}

template <TwoBodyOper::type OperType>
void
TwoBodyOSARLibint2<OperType>::compute_quartet(int *psh1, int *psh2, int *psh3, int *psh4)
{
#ifdef EREP_TIMING
  char section[30];
#endif
  GaussianBasisSet *pbs1=bs1_.pointer();
  GaussianBasisSet *pbs2=bs2_.pointer();
  GaussianBasisSet *pbs3=bs3_.pointer();
  GaussianBasisSet *pbs4=bs4_.pointer();
  int int_expweight1; // For exponent weighted contractions.
  int int_expweight2; // For exponent weighted contractions.
  int int_expweight3; // For exponent weighted contractions.
  int int_expweight4; // For exponent weighted contractions.
  int size;
  int ii;
  int size1, size2, size3, size4;
  int tam1,tam2,tam3,tam4;
  int i,j,k,l;
  int pi, pj, pk, pl;
  int gci, gcj, gck, gcl;
  int sh1,sh2,sh3,sh4;
  int osh1,osh2,osh3,osh4;
  int am1,am2,am3,am4,am12,am34;
  int minam1,minam2,minam3,minam4;
  int redundant_index;
  int e12,e13e24,e34;
  int p12,p34,p13p24;
  int eAB;

  osh1 = sh1 = *psh1;
  osh2 = sh2 = *psh2;
  osh3 = sh3 = *psh3;
  osh4 = sh4 = *psh4;

  /* Test the arguments to make sure that they are sensible. */
  if (   sh1 < 0 || sh1 >= bs1_->nbasis()
      || sh2 < 0 || sh2 >= bs2_->nbasis()
      || sh3 < 0 || sh3 >= bs3_->nbasis()
      || sh4 < 0 || sh4 >= bs4_->nbasis() ) {
    ExEnv::errn() << scprintf("compute_erep has been incorrectly used\n");
    ExEnv::errn() << scprintf("shells (bounds): %d (%d), %d (%d), %d (%d), %d (%d)\n",
                              sh1,bs1_->nbasis()-1,
                              sh2,bs2_->nbasis()-1,
                              sh3,bs3_->nbasis()-1,
                              sh4,bs4_->nbasis()-1);
    throw sc::ProgrammingError("", __FILE__, __LINE__);
  }

  /* Set up pointers to the current shells. */
  int_shell1_ = &bs1_->shell(sh1);
  int_shell2_ = &bs2_->shell(sh2);
  int_shell3_ = &bs3_->shell(sh3);
  int_shell4_ = &bs4_->shell(sh4);

  /* Compute the maximum angular momentum on each centers to
   * determine the most efficient way to invoke the building and shifting
   * routines.  The minimum angular momentum will be computed at the
   * same time. */
  minam1 = int_shell1_->min_am();
  minam2 = int_shell2_->min_am();
  minam3 = int_shell3_->min_am();
  minam4 = int_shell4_->min_am();
  am1 = int_shell1_->max_am();
  am2 = int_shell2_->max_am();
  am3 = int_shell3_->max_am();
  am4 = int_shell4_->max_am();
  am12 = am1 + am2;
  am34 = am3 + am4;

  // This condition being true is guaranteed by the constructor of IntegralLibint2
  //if (minam1 != am1 ||
  //    minam2 != am2 ||
  //    minam3 != am3 ||
  //    minam4 != am4 ) {
  //  ExEnv::errn() << scprintf("Int2eLibint2::comp_eri() cannot yet handle fully general contractions") << endl;
  //  fail();
  //}

  /* See if need to transform to spherical harmonics */
  bool need_cart2sph_transform = false;
  if (int_shell1_->has_pure() ||
      int_shell2_->has_pure() ||
      int_shell3_->has_pure() ||
      int_shell4_->has_pure())
    need_cart2sph_transform = true;


  /* See if contraction quartets need to be resorted into a shell quartet */
  bool need_sort_to_shell_quartet = false;
  int num_gen_shells = 0;
  if (int_shell1_->ncontraction() > 1)
    num_gen_shells++;
  if (int_shell2_->ncontraction() > 1)
    num_gen_shells++;
  if (int_shell3_->ncontraction() > 1)
    num_gen_shells++;
  if (int_shell4_->ncontraction() > 1)
    num_gen_shells++;
  if (am12+am34 && num_gen_shells >= 1)
    need_sort_to_shell_quartet = true;

  /* Unique integrals are needed only ?*/
  bool need_unique_ints_only = false;
  if (!redundant_) {
    e12 = 0;
    if (int_shell1_ == int_shell2_ && int_shell1_->nfunction()>1)
      e12 = 1;
    e34 = 0;
    if (int_shell3_ == int_shell4_ && int_shell3_->nfunction()>1)
      e34 = 1;
    e13e24 = 0;
    if (int_shell1_ == int_shell3_ && int_shell2_ == int_shell4_ && int_shell1_->nfunction()*int_shell2_->nfunction()>1)
      e13e24 = 1;

    if ( e12 || e34 || e13e24 )
      need_unique_ints_only = true;
  }


#ifdef EREP_TIMING
  sprintf(section,"erep am=%02d",am12+am34);
  tim_enter(section);
  tim_enter("setup");
#endif

  /* Convert the integral to the most efficient form. */
  p12 = 0;
  p34 = 0;
  p13p24 = 0;

  if (am2 > am1) {
    p12 = 1;
    iswtch(&am1,&am2);iswtch(&sh1,&sh2);iswtch(psh1,psh2);
    iswtch(&minam1,&minam2);
    pswtch((void**)&int_shell1_,(void**)&int_shell2_);
    swtch(pbs1,pbs2);
  }
  if (am4 > am3) {
    p34 = 1;
    iswtch(&am3,&am4);iswtch(&sh3,&sh4);iswtch(psh3,psh4);
    iswtch(&minam3,&minam4);
    pswtch((void**)&int_shell3_,(void**)&int_shell4_);
    swtch(pbs3,pbs4);
  }
  if (am12 > am34) {
    p13p24 = 1;
    iswtch(&am1,&am3);iswtch(&sh1,&sh3);iswtch(psh1,psh3);
    iswtch(&am2,&am4);iswtch(&sh2,&sh4);iswtch(psh2,psh4);
    iswtch(&am12,&am34);
    iswtch(&minam1,&minam3);
    iswtch(&minam2,&minam4);
    pswtch((void**)&int_shell1_,(void**)&int_shell3_);
    swtch(pbs1,pbs3);
    pswtch((void**)&int_shell2_,(void**)&int_shell4_);
    swtch(pbs2,pbs4);
  }
  bool shells_were_permuted = (p12||p34||p13p24);

  /* If the centers were permuted, then the int_expweighted variable may
   * need to be changed. */
  if (p12) {
    iswtch(&int_expweight1,&int_expweight2);
  }
  if (p34) {
    iswtch(&int_expweight3,&int_expweight4);
  }
  if (p13p24) {
    iswtch(&int_expweight1,&int_expweight3);
    iswtch(&int_expweight2,&int_expweight4);
  }

  /* Compute the shell sizes. */
  size1 = int_shell1_->ncartesian();
  size2 = int_shell2_->ncartesian();
  size3 = int_shell3_->ncartesian();
  size4 = int_shell4_->ncartesian();
  size = size1*size2*size3*size4;

  /* Compute center data for Libint */
  int ctr1 = pbs1->shell_to_center(sh1);
  int ctr2 = pbs2->shell_to_center(sh2);
  int ctr3 = pbs3->shell_to_center(sh3);
  int ctr4 = pbs4->shell_to_center(sh4);
  Libint_t& libint0 = Libint_[0];
  libint0.AB_x[0] = pbs1->r(ctr1,0) - pbs2->r(ctr2,0);
  libint0.AB_y[0] = pbs1->r(ctr1,1) - pbs2->r(ctr2,1);
  libint0.AB_z[0] = pbs1->r(ctr1,2) - pbs2->r(ctr2,2);
  libint0.CD_x[0] = pbs3->r(ctr3,0) - pbs4->r(ctr4,0);
  libint0.CD_y[0] = pbs3->r(ctr3,1) - pbs4->r(ctr4,1);
  libint0.CD_z[0] = pbs3->r(ctr3,2) - pbs4->r(ctr4,2);
  for(i=0;i<3;i++) {
    quartet_info_.A[i] = pbs1->r(ctr1,i);
    quartet_info_.B[i] = pbs2->r(ctr2,i);
    quartet_info_.C[i] = pbs3->r(ctr3,i);
    quartet_info_.D[i] = pbs4->r(ctr4,i);
  }
  quartet_info_.AB2  = libint0.AB_x[0]*libint0.AB_x[0];
  quartet_info_.AB2 += libint0.AB_y[0]*libint0.AB_y[0];
  quartet_info_.AB2 += libint0.AB_z[0]*libint0.AB_z[0];
  quartet_info_.CD2  = libint0.CD_x[0]*libint0.CD_x[0];
  quartet_info_.CD2 += libint0.CD_y[0]*libint0.CD_y[0];
  quartet_info_.CD2 += libint0.CD_z[0]*libint0.CD_z[0];

  /* Set up pointers to the current shell pairs. */
  quartet_info_.shell_pair12 = shell_pairs12_->shell_pair(osh1,osh2);
  quartet_info_.shell_pair34 = shell_pairs34_->shell_pair(osh3,osh4);

  /* Remember how permuted - will need to access shell pairs in grt_quartet_data_() using the original
     primitive indices */
  quartet_info_.p12 = p12;
  quartet_info_.p34 = p34;
  quartet_info_.p13p24 = p13p24;

  /* Remember the original primitive indices to access shell pair data
     Note the reverse order of switching, p13p24 first,
     then p12 and p34 - because we need the inverse mapping! */
  quartet_info_.op1 = &quartet_info_.p1;
  quartet_info_.op2 = &quartet_info_.p2;
  quartet_info_.op3 = &quartet_info_.p3;
  quartet_info_.op4 = &quartet_info_.p4;
  if (p13p24) {
    pswtch((void **)&quartet_info_.op1,(void **)&quartet_info_.op3);
    pswtch((void **)&quartet_info_.op2,(void **)&quartet_info_.op4);
  }
  if (p12)
    pswtch((void **)&quartet_info_.op1,(void **)&quartet_info_.op2);
  if (p34)
    pswtch((void **)&quartet_info_.op3,(void **)&quartet_info_.op4);

  /* Determine where integrals need to go at each stage */
  if (shells_were_permuted)
    if (need_sort_to_shell_quartet) {
      prim_ints_ = cart_ints_;
      if (need_cart2sph_transform)
        contr_quartets_ = sphharm_ints_;
      else
        contr_quartets_ = cart_ints_;
      shell_quartet_ = perm_ints_;
    }
    else {
      prim_ints_ = cart_ints_;
      if (need_cart2sph_transform) {
        contr_quartets_ = sphharm_ints_;
        shell_quartet_ = contr_quartets_;
      }
      else
        shell_quartet_ = cart_ints_;
    }
  else
    if (need_sort_to_shell_quartet) {
      prim_ints_ = cart_ints_;
      if (need_cart2sph_transform)
        contr_quartets_ = sphharm_ints_;
      else
        contr_quartets_ = cart_ints_;
      shell_quartet_ = target_ints_buffer_;
    }
    else {
      if (need_cart2sph_transform) {
        prim_ints_ = cart_ints_;
        contr_quartets_ = target_ints_buffer_;
        shell_quartet_ = target_ints_buffer_;
      }
      else {
        prim_ints_ = target_ints_buffer_;
        shell_quartet_ = target_ints_buffer_;
      }
    }

  /* Begin loops over generalized contractions. */
  int buffer_offset = 0;
  for (gci=0; gci<int_shell1_->ncontraction(); gci++) {
    tam1 = int_shell1_->am(gci);
    int tsize1 = INT_NCART_NN(tam1);
    quartet_info_.gc1 = gci;
    for (gcj=0; gcj<int_shell2_->ncontraction(); gcj++) {
      tam2 = int_shell2_->am(gcj);
      int tsize2 = INT_NCART_NN(tam2);
      quartet_info_.gc2 = gcj;
      for (gck=0; gck<int_shell3_->ncontraction(); gck++) {
        tam3 = int_shell3_->am(gck);
        int tsize3 = INT_NCART_NN(tam3);
        quartet_info_.gc3 = gck;
        for (gcl=0; gcl<int_shell4_->ncontraction(); gcl++) {
          tam4 = int_shell4_->am(gcl);
          int tsize4 = INT_NCART_NN(tam4);
          quartet_info_.gc4 = gcl;
          quartet_info_.am = tam1 + tam2 + tam3 + tam4;
          int size = tsize1*tsize2*tsize3*tsize4;

          /* Begin loop over primitives. */
          int num_prim_combinations = 0;
          for (pi=0; pi<int_shell1_->nprimitive(); pi++) {
            quartet_info_.p1 = pi;
            for (pj=0; pj<int_shell2_->nprimitive(); pj++) {
              quartet_info_.p2 = pj;
              for (pk=0; pk<int_shell3_->nprimitive(); pk++) {
                quartet_info_.p3 = pk;
                for (pl=0; pl<int_shell4_->nprimitive(); pl++) {
                  quartet_info_.p4 = pl;

                  // Compute primitive data for Libint
                  quartet_data_(&Libint_[num_prim_combinations], 1.0);

                  ++num_prim_combinations;
                }}}}

          if (quartet_info_.am) {
            // Compute the integrals
            Libint_[0].contrdepth = num_prim_combinations;
            LIBINT2_PREFIXED_NAME(libint2_build_eri)[tam1][tam2][tam3][tam4](&Libint_[0]);
            // Copy the contracted integrals over to prim_ints_
            const LIBINT2_REALTYPE* prim_ints = Libint_[0].targets[0];
            for(int ijkl=0; ijkl<size; ijkl++)
              prim_ints_[buffer_offset + ijkl] = (double) prim_ints[ijkl];

#if 0
            std::cout << *psh1 << " " << *psh2 << " " << *psh3 << " " << *psh4 << " " << std::endl;
            for(int ijkl=0; ijkl<size; ijkl++) {
              std::cout <<  "  " << prim_ints[ijkl] << std::endl;
            }
#endif
          }
          else {
            double ssss = 0.0;
            for(int p=0; p<num_prim_combinations; ++p)
              ssss += Libint_[p].LIBINT_T_SS_EREP_SS(0)[0];
            prim_ints_[buffer_offset] = ssss;
          }
          buffer_offset += size;
        }}}}

  /*-------------------------------------------
    Transform to spherical harmonics if needed
   -------------------------------------------*/
  if (need_cart2sph_transform)
    transform_contrquartets_(prim_ints_,contr_quartets_);

  //
  // If not CCA-compliant normalization -- re-normalize all integrals to 1
  //
#if INTEGRALLIBINT2_NORMCONV != INTEGRALLIBINT2_NORMCONV_CCA
  norm_contrcart1_(need_cart2sph_transform ? contr_quartets_ : prim_ints_);
#endif

  /*----------------------------------------------
    Resort integrals from by-contraction-quartets
    into shell-quartet order if needed
   ----------------------------------------------*/
  if (need_sort_to_shell_quartet)
    sort_contrquartets_to_shellquartet_(contr_quartets_,shell_quartet_);

  /*---------------------------------
    Permute integrals back if needed
   ---------------------------------*/
  if ((!permute_)&&shells_were_permuted) {
    // handle integrals first
    permute_target_(shell_quartet_,target_ints_buffer_,p13p24,p12,p34);
    // then indices
    if (p13p24) {
      iswtch(&sh1,&sh3);iswtch(psh1,psh3);
      iswtch(&sh2,&sh4);iswtch(psh2,psh4);
      iswtch(&am1,&am3);
      iswtch(&am2,&am4);
      iswtch(&am12,&am34);
      pswtch((void**)&int_shell1_,(void**)&int_shell3_);
      swtch(pbs1,pbs3);
      pswtch((void**)&int_shell2_,(void**)&int_shell4_);
      swtch(pbs2,pbs4);
      iswtch(&int_expweight1,&int_expweight3);
      iswtch(&int_expweight2,&int_expweight4);
    }
    if (p34) {
      iswtch(&sh3,&sh4);iswtch(psh3,psh4);
      iswtch(&am3,&am4);
      pswtch((void**)&int_shell3_,(void**)&int_shell4_);
      swtch(pbs3,pbs4);
      iswtch(&int_expweight3,&int_expweight4);
    }
    if (p12) {
      iswtch(&sh1,&sh2);iswtch(psh1,psh2);
      iswtch(&am1,&am2);
      pswtch((void**)&int_shell1_,(void**)&int_shell2_);
      swtch(pbs1,pbs2);
      iswtch(&int_expweight1,&int_expweight2);
    }
  }

  /*--- Extract unique integrals if needed ---*/
  if (need_unique_ints_only)
    get_nonredundant_ints_(target_ints_buffer_,target_ints_buffer_,e13e24,e12,e34);

}

}

#endif // header guard
#endif // if LIBINT2_SUPPORT_ERI

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
