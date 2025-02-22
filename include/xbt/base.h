/* xbt.h - Public interface to the xbt (simgrid's toolbox)                  */

/* Copyright (c) 2004-2019. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#ifndef XBT_BASE_H
#define XBT_BASE_H

/* Define _GNU_SOURCE for getline, isfinite, etc. */
#ifndef _GNU_SOURCE
#  define _GNU_SOURCE
#endif

// Teach the compiler that some code path is unreacheable:
#if defined(__has_builtin)
#  if __has_builtin(__builtin_unreachable)
#    define XBT_UNREACHABLE() __builtin_unreachable()
#  else
#    include <stdlib.h>
#    define XBT_UNREACHABLE() abort()
#  endif
#elif (__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5)
#  define XBT_UNREACHABLE() __builtin_unreachable()
#else
#  include <stdlib.h>
#  define XBT_UNREACHABLE() abort()
#endif

/* On MinGW, stdio.h defines __MINGW_PRINTF_FORMAT and __MINGW_SCANF_FORMAT
   which are the suitable format style (either gnu_printf or ms_printf)
   depending on which version is available (__USE_MINGW_ANSI_STDIO): */
#ifdef __MINGW32__
#  include <stdio.h>

#define XBT_ATTRIB_PRINTF(format_idx, arg_idx)                                                                         \
  __attribute__((__format__(__MINGW_PRINTF_FORMAT, (format_idx), (arg_idx))))
#define XBT_ATTRIB_SCANF(format_idx, arg_idx) __attribute__((__MINGW_SCANF_FORMAT(__scanf__, (format_idx), (arg_idx))))
#else
#define XBT_ATTRIB_PRINTF(format_idx, arg_idx) __attribute__((__format__(__printf__, (format_idx), (arg_idx))))
#define XBT_ATTRIB_SCANF(format_idx, arg_idx) __attribute__((__format__(__scanf__, (format_idx), (arg_idx))))
#endif

#if defined(__cplusplus)
#if __cplusplus >= 201103L
#define XBT_ATTRIB_NORETURN [[noreturn]]
#endif
#if __cplusplus >= 201703L
#define XBT_ATTRIB_UNUSED [[maybe_unused]]
#endif
#if __cplusplus >= 201402L
#define XBT_ATTRIB_DEPRECATED(mesg) [[deprecated(mesg)]]
#endif
#elif defined(__STDC_VERSION__)
#if __STDC_VERSION__ >= 201112L
#define XBT_ATTRIB_NORETURN _Noreturn
#endif
#endif

#ifndef XBT_ATTRIB_NORETURN
#define XBT_ATTRIB_NORETURN __attribute__((noreturn))
#endif
#ifndef XBT_ATTRIB_UNUSED
#define XBT_ATTRIB_UNUSED  __attribute__((unused))
#endif
#ifndef XBT_ATTRIB_DEPRECATED
#define XBT_ATTRIB_DEPRECATED(mesg) __attribute__((deprecated(mesg)))
#endif

#define XBT_ATTRIB_DEPRECATED_v325(mesg)                                                                               \
  XBT_ATTRIB_DEPRECATED(mesg " (this compatibility wrapper will be dropped in v3.25)") /* Will be dropped in v3.25 */
#define XBT_ATTRIB_DEPRECATED_v327(mesg)                                                                               \
  XBT_ATTRIB_DEPRECATED(mesg " (this compatibility wrapper will be dropped in v3.27)") /* Will be dropped in v3.27 */

#if !defined(__APPLE__)
#  define XBT_ATTRIB_CONSTRUCTOR(prio) __attribute__((__constructor__(prio)))
#  define XBT_ATTRIB_DESTRUCTOR(prio) __attribute__((__destructor__(prio)))
#else
#  define XBT_ATTRIB_CONSTRUCTOR(prio) __attribute__((__constructor__))
#  define XBT_ATTRIB_DESTRUCTOR(prio) __attribute__((__destructor__))
#endif

#ifndef XBT_ALWAYS_INLINE /* defined also in libsosp */
#  if defined(__GNUC__)
#    define XBT_ALWAYS_INLINE inline __attribute__ ((always_inline))
#  else
#    define XBT_ALWAYS_INLINE inline
#  endif
#endif

/* Stringify argument. */
#define _XBT_STRINGIFY(a) #a

/* Concatenate arguments. _XBT_CONCAT2 adds a level of indirection over _XBT_CONCAT. */
#define _XBT_CONCAT(a, b) a##b
#define _XBT_CONCAT2(a, b) _XBT_CONCAT(a, b)
#define _XBT_CONCAT3(a, b, c) _XBT_CONCAT2(_XBT_CONCAT(a, b), c)
#define _XBT_CONCAT4(a, b, c, d) _XBT_CONCAT2(_XBT_CONCAT3(a, b, c), d)

/*
 * Expands to `one' if there is only one argument for the variadic part.
 * Otherwise, expands to `more'.
 * Works with up to 63 arguments, which is the maximum mandated by the C99 standard.
 */
#define _XBT_IF_ONE_ARG(one, more, ...)                                 \
    _XBT_IF_ONE_ARG_(__VA_ARGS__,                                       \
                     more, more, more, more, more, more, more, more,    \
                     more, more, more, more, more, more, more, more,    \
                     more, more, more, more, more, more, more, more,    \
                     more, more, more, more, more, more, more, more,    \
                     more, more, more, more, more, more, more, more,    \
                     more, more, more, more, more, more, more, more,    \
                     more, more, more, more, more, more, more, more,    \
                     more, more, more, more, more, more, more, one)
#define _XBT_IF_ONE_ARG_(a64, a63, a62, a61, a60, a59, a58, a57,        \
                         a56, a55, a54, a53, a52, a51, a50, a49,        \
                         a48, a47, a46, a45, a44, a43, a42, a41,        \
                         a40, a39, a38, a37, a36, a35, a34, a33,        \
                         a32, a31, a30, a29, a28, a27, a26, a25,        \
                         a24, a23, a22, a21, a20, a19, a18, a17,        \
                         a16, a15, a14, a13, a12, a11, a10, a9,         \
                         a8, a7, a6, a5, a4, a3, a2, a1, N, ...) N

/* Handle import/export stuff
 *
 * Rationale of XBT_PUBLIC:
 *   * This is for library symbols visible from the application-land.
 *     Basically, any symbols defined in the include/directory must be like this (plus some other globals).
 *
 *     UNIX coders should just think of it as a special way to say "extern".
 *
 *   * If you build the DLL, define the DLL_EXPORT symbol so that all symbols actually get exported by this file.
 *
 *   * If you link your application against the DLL or if you do a UNIX build, don't do anything special. This file
 *     will do the right thing for you by default.
 *
 * Rationale of XBT_EXPORT_NO_IMPORT: (windows-only)
 *   * Symbols which must be exported in the DLL, but not imported from it.
 *
 *   * This is obviously useful for initialized globals (which cannot be  extern or similar).
 *   * This is also used in the log mechanism where a macro creates the variable automatically. When the macro is
 *     called from within SimGrid, the symbol must be exported, but when called  from within the client code, it must
 *     not try to retrieve the symbol from the DLL since it's not in there.
 *
 * Rationale of XBT_IMPORT_NO_EXPORT: (windows-only)
 *   * Symbols which must be imported from the DLL, but not explicitly  exported from it.
 *
 *   * The root log category is already exported, but not imported explicitly when creating a subcategory since we
 *     cannot import the parent category to deal with the fact that the parent may be in application space, not DLL
 *     space.
 */

/* Build the DLL */
#if defined(DLL_EXPORT)
#  define XBT_PUBLIC            __declspec(dllexport)
#  define XBT_EXPORT_NO_IMPORT  __declspec(dllexport)
#  define XBT_IMPORT_NO_EXPORT
#  define XBT_PUBLIC_DATA       extern __declspec(dllexport)
#  define XBT_PRIVATE

/* Link against the DLL */
#elif (defined(_WIN32) && !defined(DLL_EXPORT))
#  define XBT_PUBLIC            __declspec(dllimport)
#  define XBT_EXPORT_NO_IMPORT
#  define XBT_IMPORT_NO_EXPORT  __declspec(dllimport)
#  define XBT_PUBLIC_DATA       extern __declspec(dllimport)
#  define XBT_PRIVATE

#elif defined(__ELF__)
#  define XBT_PUBLIC            __attribute__((visibility("default")))
#  define XBT_EXPORT_NO_IMPORT  __attribute__((visibility("default")))
#  define XBT_IMPORT_NO_EXPORT  __attribute__((visibility("default")))
#  define XBT_PUBLIC_DATA       extern __attribute__((visibility("default")))
#  define XBT_PRIVATE           __attribute__((visibility("hidden")))

#else
#  define XBT_PUBLIC            /* public */
#  define XBT_EXPORT_NO_IMPORT
#  define XBT_IMPORT_NO_EXPORT
#  define XBT_PUBLIC_DATA       extern
#  define XBT_PRIVATE           /** @private */

#endif

/* C++ users need love */
#ifndef SG_BEGIN_DECL
# ifdef __cplusplus
#  define SG_BEGIN_DECL() extern "C" {
# else
#  define SG_BEGIN_DECL()
# endif
#endif

#ifndef SG_END_DECL
# ifdef __cplusplus
#  define SG_END_DECL() }
# else
#  define SG_END_DECL()
# endif
#endif
/* End of cruft for C++ */

#endif
