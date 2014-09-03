/* Fetching
 * Elaine Angelino, Eddie Kohler
 * Copyright (c) 2012-2014 President and Fellows of Harvard College
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, subject to the conditions
 * listed in the Fetching LICENSE file. These conditions include: you must
 * preserve this copyright notice, and you cannot mention the copyright
 * holders in advertising related to the Software without their permission.
 * The Software is provided WITHOUT ANY WARRANTY, EXPRESS OR IMPLIED. This
 * notice is a summary of the Fetching LICENSE file; the license in that file
 * is legally binding.
 */
#ifndef MW_MCMC_COMPILER_HH
#define MW_MCMC_COMPILER_HH

template <typename T, typename U = void> struct do_nothing;

/** @brief Binary function object that does nothing when called. */
template <typename T, typename U>
struct do_nothing {
    typedef T first_argument_type;
    typedef U second_argument_type;
    typedef void result_type;
    void operator()(const T&, const U&) const {
    }
};

/** @brief Unary function object that does nothing when called. */
template <typename T>
struct do_nothing<T, void> {
    typedef T argument_type;
    typedef void result_type;
    void operator()(const T&) const {
    }
};

#if defined(__GNUC__) && __GNUC__
# define LCDF_GCCVERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#else
# define LCDF_GCCVERSION 0
#endif
#if defined(__clang__) && __clang__
# define LCDF_CLANGVERSION (__clang_major__ * 10000 + __clang_minor__ * 100 + __clang_patchlevel__)
#else
# define LCDF_CLANGVERSION 0
#endif

#if defined(__cplusplus) && (__cplusplus >= 201103L || (defined(__GXX_EXPERIMENTAL_CXX0X__) && __GXX_EXPERIMENTAL_CXX0X__))
# define LCDF_HAVE_SOME_CXX11 1
#else
# define LCDF_HAVE_SOME_CXX11 0
#endif

#if !(LCDF_HAVE_SOME_CXX11 && (LCDF_CLANGVERSION >= 20900 || LCDF_GCCVERSION >= 40300))
# define static_assert(x, ...) switch ((int) (x)) case 0: case !!((int) (x)):
#endif

#if LCDF_HAVE_SOME_CXX11 && (LCDF_CLANGVERSION >= 20900 || LCDF_GCCVERSION >= 40300)
# define LCDF_HAVE_CXX_VARIADIC_TEMPLATES
#endif

#if defined(__ICC) && __ICC
# define ALLOW_FLOAT_EQUALITY _Pragma("warning(disable:1572)")
# define DISALLOW_FLOAT_EQUALITY _Pragma("warning(default:1572)")
#else
# define ALLOW_FLOAT_EQUALITY /* nothing */
# define DISALLOW_FLOAT_EQUALITY /* nothing */
#endif

#endif
