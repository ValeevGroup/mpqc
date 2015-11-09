//
// Copyright (c) 2002--2010
// Toon Knapen, Karl Meerbergen, Kresimir Fresl,
// Thomas Klimpel and Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// THIS FILE IS AUTOMATICALLY GENERATED
// PLEASE DO NOT EDIT!
//

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_DRIVER_GEES_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_DRIVER_GEES_HPP

#include <boost/assert.hpp>
#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/detail/array.hpp>
#include <boost/numeric/bindings/detail/complex_utils.hpp>
#include <boost/numeric/bindings/is_column_major.hpp>
#include <boost/numeric/bindings/is_complex.hpp>
#include <boost/numeric/bindings/is_mutable.hpp>
#include <boost/numeric/bindings/is_real.hpp>
#include <boost/numeric/bindings/lapack/workspace.hpp>
#include <boost/numeric/bindings/remove_imaginary.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/traits/detail/utils.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/utility/enable_if.hpp>

//
// The LAPACK-backend for gees is the netlib-compatible backend.
//
#include <boost/numeric/bindings/lapack/detail/lapack.h>
#include <boost/numeric/bindings/lapack/detail/lapack_option.hpp>

namespace boost {
namespace numeric {
namespace bindings {
namespace lapack {

//
// The detail namespace contains value-type-overloaded functions that
// dispatch to the appropriate back-end LAPACK-routine.
//
namespace detail {

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * float value-type.
//
inline std::ptrdiff_t gees( const char jobvs, const char sort,
        external_fp select, const fortran_int_t n, float* a,
        const fortran_int_t lda, fortran_int_t& sdim, float* wr, float* wi,
        float* vs, const fortran_int_t ldvs, float* work,
        const fortran_int_t lwork, fortran_bool_t* bwork ) {
    fortran_int_t info(0);
    LAPACK_SGEES( &jobvs, &sort, select, &n, a, &lda, &sdim, wr, wi, vs,
            &ldvs, work, &lwork, bwork, &info );
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * double value-type.
//
inline std::ptrdiff_t gees( const char jobvs, const char sort,
        external_fp select, const fortran_int_t n, double* a,
        const fortran_int_t lda, fortran_int_t& sdim, double* wr, double* wi,
        double* vs, const fortran_int_t ldvs, double* work,
        const fortran_int_t lwork, fortran_bool_t* bwork ) {
    fortran_int_t info(0);
    LAPACK_DGEES( &jobvs, &sort, select, &n, a, &lda, &sdim, wr, wi, vs,
            &ldvs, work, &lwork, bwork, &info );
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<float> value-type.
//
inline std::ptrdiff_t gees( const char jobvs, const char sort,
        external_fp select, const fortran_int_t n, std::complex<float>* a,
        const fortran_int_t lda, fortran_int_t& sdim, std::complex<float>* w,
        std::complex<float>* vs, const fortran_int_t ldvs,
        std::complex<float>* work, const fortran_int_t lwork, float* rwork,
        fortran_bool_t* bwork ) {
    fortran_int_t info(0);
    LAPACK_CGEES( &jobvs, &sort, select, &n, a, &lda, &sdim, w, vs, &ldvs,
            work, &lwork, rwork, bwork, &info );
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<double> value-type.
//
inline std::ptrdiff_t gees( const char jobvs, const char sort,
        external_fp select, const fortran_int_t n, std::complex<double>* a,
        const fortran_int_t lda, fortran_int_t& sdim, std::complex<double>* w,
        std::complex<double>* vs, const fortran_int_t ldvs,
        std::complex<double>* work, const fortran_int_t lwork, double* rwork,
        fortran_bool_t* bwork ) {
    fortran_int_t info(0);
    LAPACK_ZGEES( &jobvs, &sort, select, &n, a, &lda, &sdim, w, vs, &ldvs,
            work, &lwork, rwork, bwork, &info );
    return info;
}

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to gees.
//
template< typename Value, typename Enable = void >
struct gees_impl {};

//
// This implementation is enabled if Value is a real type.
//
template< typename Value >
struct gees_impl< Value, typename boost::enable_if< is_real< Value > >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixA, typename VectorWR, typename VectorWI,
            typename MatrixVS, typename WORK, typename BWORK >
    static std::ptrdiff_t invoke( const char jobvs, const char sort,
            external_fp select, MatrixA& a, fortran_int_t& sdim,
            VectorWR& wr, VectorWI& wi, MatrixVS& vs, detail::workspace2<
            WORK, BWORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        BOOST_STATIC_ASSERT( (bindings::is_column_major< MatrixA >::value) );
        BOOST_STATIC_ASSERT( (bindings::is_column_major< MatrixVS >::value) );
        BOOST_STATIC_ASSERT( (boost::is_same< typename remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename remove_const< typename bindings::value_type<
                VectorWR >::type >::type >::value) );
        BOOST_STATIC_ASSERT( (boost::is_same< typename remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename remove_const< typename bindings::value_type<
                VectorWI >::type >::type >::value) );
        BOOST_STATIC_ASSERT( (boost::is_same< typename remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename remove_const< typename bindings::value_type<
                MatrixVS >::type >::type >::value) );
        BOOST_STATIC_ASSERT( (bindings::is_mutable< MatrixA >::value) );
        BOOST_STATIC_ASSERT( (bindings::is_mutable< VectorWR >::value) );
        BOOST_STATIC_ASSERT( (bindings::is_mutable< VectorWI >::value) );
        BOOST_STATIC_ASSERT( (bindings::is_mutable< MatrixVS >::value) );
        BOOST_ASSERT( bindings::size(wi) >= bindings::size_column(a) );
        BOOST_ASSERT( bindings::size(work.select(fortran_bool_t())) >=
                min_size_bwork( bindings::size_column(a), sort ));
        BOOST_ASSERT( bindings::size(work.select(real_type())) >=
                min_size_work( bindings::size_column(a) ));
        BOOST_ASSERT( bindings::size(wr) >= bindings::size_column(a) );
        BOOST_ASSERT( bindings::size_column(a) >= 0 );
        BOOST_ASSERT( bindings::size_minor(a) == 1 ||
                bindings::stride_minor(a) == 1 );
        BOOST_ASSERT( bindings::size_minor(vs) == 1 ||
                bindings::stride_minor(vs) == 1 );
        BOOST_ASSERT( bindings::stride_major(a) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column(a)) );
        BOOST_ASSERT( jobvs == 'N' || jobvs == 'V' );
        BOOST_ASSERT( sort == 'N' || sort == 'S' );
        return detail::gees( jobvs, sort, select, bindings::size_column(a),
                bindings::begin_value(a), bindings::stride_major(a), sdim,
                bindings::begin_value(wr), bindings::begin_value(wi),
                bindings::begin_value(vs), bindings::stride_major(vs),
                bindings::begin_value(work.select(real_type())),
                bindings::size(work.select(real_type())),
                bindings::begin_value(work.select(fortran_bool_t())) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename MatrixA, typename VectorWR, typename VectorWI,
            typename MatrixVS >
    static std::ptrdiff_t invoke( const char jobvs, const char sort,
            external_fp select, MatrixA& a, fortran_int_t& sdim,
            VectorWR& wr, VectorWI& wi, MatrixVS& vs, minimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        bindings::detail::array< real_type > tmp_work( min_size_work(
                bindings::size_column(a) ) );
        bindings::detail::array< fortran_bool_t > tmp_bwork( min_size_bwork(
                bindings::size_column(a), sort ) );
        return invoke( jobvs, sort, select, a, sdim, wr, wi, vs,
                workspace( tmp_work, tmp_bwork ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename MatrixA, typename VectorWR, typename VectorWI,
            typename MatrixVS >
    static std::ptrdiff_t invoke( const char jobvs, const char sort,
            external_fp select, MatrixA& a, fortran_int_t& sdim,
            VectorWR& wr, VectorWI& wi, MatrixVS& vs, optimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        real_type opt_size_work;
        bindings::detail::array< fortran_bool_t > tmp_bwork( min_size_bwork(
                bindings::size_column(a), sort ) );
        detail::gees( jobvs, sort, select, bindings::size_column(a),
                bindings::begin_value(a), bindings::stride_major(a), sdim,
                bindings::begin_value(wr), bindings::begin_value(wi),
                bindings::begin_value(vs), bindings::stride_major(vs),
                &opt_size_work, -1, bindings::begin_value(tmp_bwork) );
        bindings::detail::array< real_type > tmp_work(
                traits::detail::to_int( opt_size_work ) );
        return invoke( jobvs, sort, select, a, sdim, wr, wi, vs,
                workspace( tmp_work, tmp_bwork ) );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const std::ptrdiff_t n ) {
        return std::max< std::ptrdiff_t >( 1, 3*n );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array bwork.
    //
    static std::ptrdiff_t min_size_bwork( const std::ptrdiff_t n,
            const char sort ) {
        if ( sort == 'N' )
            return 0;
        else
            return n;
    }
};

//
// This implementation is enabled if Value is a complex type.
//
template< typename Value >
struct gees_impl< Value, typename boost::enable_if< is_complex< Value > >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixA, typename VectorW, typename MatrixVS,
            typename WORK, typename RWORK, typename BWORK >
    static std::ptrdiff_t invoke( const char jobvs, const char sort,
            external_fp select, MatrixA& a, fortran_int_t& sdim,
            VectorW& w, MatrixVS& vs, detail::workspace3< WORK, RWORK,
            BWORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        BOOST_STATIC_ASSERT( (bindings::is_column_major< MatrixA >::value) );
        BOOST_STATIC_ASSERT( (bindings::is_column_major< MatrixVS >::value) );
        BOOST_STATIC_ASSERT( (boost::is_same< typename remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename remove_const< typename bindings::value_type<
                VectorW >::type >::type >::value) );
        BOOST_STATIC_ASSERT( (boost::is_same< typename remove_const<
                typename bindings::value_type< MatrixA >::type >::type,
                typename remove_const< typename bindings::value_type<
                MatrixVS >::type >::type >::value) );
        BOOST_STATIC_ASSERT( (bindings::is_mutable< MatrixA >::value) );
        BOOST_STATIC_ASSERT( (bindings::is_mutable< VectorW >::value) );
        BOOST_STATIC_ASSERT( (bindings::is_mutable< MatrixVS >::value) );
        BOOST_ASSERT( bindings::size(w) >= bindings::size_column(a) );
        BOOST_ASSERT( bindings::size(work.select(fortran_bool_t())) >=
                min_size_bwork( bindings::size_column(a), sort ));
        BOOST_ASSERT( bindings::size(work.select(real_type())) >=
                min_size_rwork( bindings::size_column(a) ));
        BOOST_ASSERT( bindings::size(work.select(value_type())) >=
                min_size_work( bindings::size_column(a) ));
        BOOST_ASSERT( bindings::size_column(a) >= 0 );
        BOOST_ASSERT( bindings::size_minor(a) == 1 ||
                bindings::stride_minor(a) == 1 );
        BOOST_ASSERT( bindings::size_minor(vs) == 1 ||
                bindings::stride_minor(vs) == 1 );
        BOOST_ASSERT( bindings::stride_major(a) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column(a)) );
        BOOST_ASSERT( jobvs == 'N' || jobvs == 'V' );
        BOOST_ASSERT( sort == 'N' || sort == 'S' );
        return detail::gees( jobvs, sort, select, bindings::size_column(a),
                bindings::begin_value(a), bindings::stride_major(a), sdim,
                bindings::begin_value(w), bindings::begin_value(vs),
                bindings::stride_major(vs),
                bindings::begin_value(work.select(value_type())),
                bindings::size(work.select(value_type())),
                bindings::begin_value(work.select(real_type())),
                bindings::begin_value(work.select(fortran_bool_t())) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename MatrixA, typename VectorW, typename MatrixVS >
    static std::ptrdiff_t invoke( const char jobvs, const char sort,
            external_fp select, MatrixA& a, fortran_int_t& sdim,
            VectorW& w, MatrixVS& vs, minimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        bindings::detail::array< value_type > tmp_work( min_size_work(
                bindings::size_column(a) ) );
        bindings::detail::array< real_type > tmp_rwork( min_size_rwork(
                bindings::size_column(a) ) );
        bindings::detail::array< fortran_bool_t > tmp_bwork( min_size_bwork(
                bindings::size_column(a), sort ) );
        return invoke( jobvs, sort, select, a, sdim, w, vs,
                workspace( tmp_work, tmp_rwork, tmp_bwork ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename MatrixA, typename VectorW, typename MatrixVS >
    static std::ptrdiff_t invoke( const char jobvs, const char sort,
            external_fp select, MatrixA& a, fortran_int_t& sdim,
            VectorW& w, MatrixVS& vs, optimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        value_type opt_size_work;
        bindings::detail::array< real_type > tmp_rwork( min_size_rwork(
                bindings::size_column(a) ) );
        bindings::detail::array< fortran_bool_t > tmp_bwork( min_size_bwork(
                bindings::size_column(a), sort ) );
        detail::gees( jobvs, sort, select, bindings::size_column(a),
                bindings::begin_value(a), bindings::stride_major(a), sdim,
                bindings::begin_value(w), bindings::begin_value(vs),
                bindings::stride_major(vs), &opt_size_work, -1,
                bindings::begin_value(tmp_rwork),
                bindings::begin_value(tmp_bwork) );
        bindings::detail::array< value_type > tmp_work(
                traits::detail::to_int( opt_size_work ) );
        return invoke( jobvs, sort, select, a, sdim, w, vs,
                workspace( tmp_work, tmp_rwork, tmp_bwork ) );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const std::ptrdiff_t n ) {
        return std::max< std::ptrdiff_t >( 1, 2*n );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array rwork.
    //
    static std::ptrdiff_t min_size_rwork( const std::ptrdiff_t n ) {
        return n;
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array bwork.
    //
    static std::ptrdiff_t min_size_bwork( const std::ptrdiff_t n,
            const char sort ) {
        if ( sort == 'N' )
            return 0;
        else
            return n;
    }
};


//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. In
// addition, if applicable, they are overloaded for user-defined workspaces.
// Calls to these functions are passed to the gees_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for gees. Its overload differs for
// * User-defined workspace
//
template< typename MatrixA, typename VectorWR, typename VectorWI,
        typename MatrixVS, typename Workspace >
inline typename boost::enable_if< detail::is_workspace< Workspace >,
        std::ptrdiff_t >::type
gees( const char jobvs, const char sort, external_fp select, MatrixA& a,
        fortran_int_t& sdim, VectorWR& wr, VectorWI& wi, MatrixVS& vs,
        Workspace work ) {
    return gees_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( jobvs, sort, select, a, sdim, wr, wi,
            vs, work );
}

//
// Overloaded function for gees. Its overload differs for
// * Default workspace-type (optimal)
//
template< typename MatrixA, typename VectorWR, typename VectorWI,
        typename MatrixVS >
inline typename boost::disable_if< detail::is_workspace< MatrixVS >,
        std::ptrdiff_t >::type
gees( const char jobvs, const char sort, external_fp select, MatrixA& a,
        fortran_int_t& sdim, VectorWR& wr, VectorWI& wi, MatrixVS& vs ) {
    return gees_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( jobvs, sort, select, a, sdim, wr, wi,
            vs, optimal_workspace() );
}

//
// Overloaded function for gees. Its overload differs for
// * User-defined workspace
//
template< typename MatrixA, typename VectorW, typename MatrixVS,
        typename Workspace >
inline typename boost::enable_if< mpl::and_< is_complex<
        typename bindings::value_type< MatrixA >::type >,
        detail::is_workspace< Workspace > >,
        std::ptrdiff_t >::type
gees( const char jobvs, const char sort, external_fp select, MatrixA& a,
        fortran_int_t& sdim, VectorW& w, MatrixVS& vs, Workspace work ) {
    return gees_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( jobvs, sort, select, a, sdim, w, vs,
            work );
}

//
// Overloaded function for gees. Its overload differs for
// * Default workspace-type (optimal)
//
template< typename MatrixA, typename VectorW, typename MatrixVS >
inline typename boost::disable_if< mpl::or_< is_real<
        typename bindings::value_type< MatrixA >::type >,
        detail::is_workspace< MatrixVS > >,
        std::ptrdiff_t >::type
gees( const char jobvs, const char sort, external_fp select, MatrixA& a,
        fortran_int_t& sdim, VectorW& w, MatrixVS& vs ) {
    return gees_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( jobvs, sort, select, a, sdim, w, vs,
            optimal_workspace() );
}

//
// Overloaded function for gees. Its overload differs for
// * VectorW
// * User-defined workspace
//
template< typename MatrixA, typename VectorW, typename MatrixVS,
        typename Workspace >
inline typename boost::enable_if< mpl::and_< is_real<
        typename bindings::value_type< MatrixA >::type >,
        detail::is_workspace< Workspace > >,
        std::ptrdiff_t >::type
gees( const char jobvs, const char sort, external_fp select, MatrixA& a,
        fortran_int_t& sdim, VectorW& w, MatrixVS& vs, Workspace work ) {
    std::ptrdiff_t info = gees_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( jobvs, sort, select, a, sdim,
            bindings::detail::real_part_view(w), bindings::detail::imag_part_view(w),
            vs, work );
    bindings::detail::interlace(w);
    return info;
}

//
// Overloaded function for gees. Its overload differs for
// * VectorW
// * Default workspace-type (optimal)
//
template< typename MatrixA, typename VectorW, typename MatrixVS >
inline typename boost::disable_if< mpl::or_< is_complex<
        typename bindings::value_type< MatrixA >::type >,
        detail::is_workspace< MatrixVS > >,
        std::ptrdiff_t >::type
gees( const char jobvs, const char sort, external_fp select, MatrixA& a,
        fortran_int_t& sdim, VectorW& w, MatrixVS& vs ) {
    std::ptrdiff_t info = gees_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( jobvs, sort, select, a, sdim,
            bindings::detail::real_part_view(w), bindings::detail::imag_part_view(w),
            vs, optimal_workspace() );
    bindings::detail::interlace(w);
    return info;
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif