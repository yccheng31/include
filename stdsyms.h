/* $Header: stdsyms.h,v 1.10.98.6 96/02/20 20:17:37 venkatp Exp $ */

#ifndef _SYS_STDSYMS_INCLUDED /* allows multiple inclusion */
#define _SYS_STDSYMS_INCLUDED

/* Namespace control symbols for open systems standards.  See stdsyms(5). */

#ifdef _HPUX_SOURCE
#  ifndef _INCLUDE__STDC__
#    define _INCLUDE__STDC__
#  endif /* _INCLUDE__STDC__ */
#  define _INCLUDE_POSIX_SOURCE
#  define _INCLUDE_POSIX2_SOURCE
#  define _INCLUDE_POSIX4_SOURCE
#  define _INCLUDE_XOPEN_SOURCE
#  define _INCLUDE_XOPEN_SOURCE_EXTENDED
#  define _INCLUDE_AES_SOURCE
#  define _INCLUDE_HPUX_SOURCE
#else
#  ifdef _AES_SOURCE
#    ifndef _INCLUDE__STDC__
#      define _INCLUDE__STDC__
#    endif /* _INCLUDE__STDC__ */
#    define _INCLUDE_POSIX_SOURCE
#    define _INCLUDE_XOPEN_SOURCE
#    define _INCLUDE_AES_SOURCE
#  else
#    if defined(_XOPEN_SOURCE) || defined(_XOPEN_SOURCE_EXTENDED)
#      ifndef _INCLUDE__STDC__
#        define _INCLUDE__STDC__
#      endif /* _INCLUDE__STDC__ */
#      define _INCLUDE_POSIX_SOURCE
#      define _INCLUDE_XOPEN_SOURCE
#    else
#      if defined(_POSIX_SOURCE) || defined(_POSIX2_SOURCE) || defined(_POSIX_C_SOURCE)
#        ifndef _INCLUDE__STDC__
#          define _INCLUDE__STDC__
#        endif /* _INCLUDE__STDC__ */
#        define _INCLUDE_POSIX_SOURCE
#      else
#        ifndef __STDC__
#          ifndef _INCLUDE__STDC__
#            define _HPUX_SOURCE
#            define _INCLUDE__STDC__
#            define _INCLUDE_POSIX_SOURCE
#            define _INCLUDE_POSIX2_SOURCE
#            define _INCLUDE_POSIX4_SOURCE
#            define _INCLUDE_XOPEN_SOURCE
#            define _INCLUDE_XOPEN_SOURCE_EXTENDED
#            define _INCLUDE_AES_SOURCE
#            define _INCLUDE_HPUX_SOURCE
#          endif /* not _INCLUDE__STDC__ */
#        else /* __STDC__ */
#          ifndef _INCLUDE__STDC__
#            define _INCLUDE__STDC__
#          endif /* _INCLUDE__STDC__ */
#        endif /* else __STDC__ */
#      endif /* _POSIX_SOURCE || _POSIX2_SOURCE || _POSIX_C_SOURCE */
#    endif /* _XOPEN_SOURCE || _XOPEN_SOURCE_EXTENDED */
#  endif /* _AES_SOURCE */
#endif /* _HPUX_SOURCE */


/* _CLASSIC_TYPES enables old semantics (e.g., char * instead of void *),
   as appropriate to the requested namespace */

#if defined(_CLASSIC_TYPES) || defined(_XPG2) || defined(_SVID2)
#  ifndef __STDC__
#    ifdef _HPUX_SOURCE
#      define _CLASSIC_ANSI_TYPES
#      define _CLASSIC_POSIX_TYPES
#      define _CLASSIC_XOPEN_TYPES
#      define _CLASSIC_ID_TYPES
#      if defined(_SVID2) && !defined(_XPG2)
#         define _XPG2
#      endif /* _SVID2 && not _XPG2 */
#    else /* not _HPUX_SOURCE */
#      if defined(_POSIX_SOURCE) && ! defined(_XOPEN_SOURCE)
#        define _CLASSIC_ANSI_TYPES
#      endif /* _POSIX_SOURCE && not _XOPEN_SOURCE */
#    endif /* else not _HPUX_SOURCE */
#  endif /* not __STDC__ */
#endif /* _CLASSIC_TYPES */


/* Enable extended X/Open namespace (_INCLUDE_XOPEN_SOURCE_EXTENDED) whenever
   the user defines _XOPEN_SOURCE_EXTENDED, or defines both _XOPEN_SOURCE and
    _XPG4_EXTENDED.  */

#if defined(_XOPEN_SOURCE_EXTENDED) || (defined(_XOPEN_SOURCE) && defined(_XPG4_EXTENDED))
#  ifndef _INCLUDE_XOPEN_SOURCE_EXTENDED
#    define _INCLUDE_XOPEN_SOURCE_EXTENDED
#  endif
#endif /* _XOPEN_SOURCE_EXTENDED || (_XOPEN_SOURCE && _XPG4_EXTENDED) */


/* Enable extended X/Open semantics (_XPG4_EXTENDED) whenever the user
   defines _XOPEN_SOURCE_EXTENDED.  */

#ifdef _INCLUDE_XOPEN_SOURCE_EXTENDED
#  ifndef _XPG4_EXTENDED
#    define _XPG4_EXTENDED
#  endif /* not _XPG4_EXTENDED */
#  define _INCLUDE_AES_SOURCE
#endif /* _XOPEN_SOURCE_EXTENDED */


/* Enable base X/Open semantics (_XPG4) whenever _XPG4_EXTENDED is on,
 * since some header files still use "#ifdef _XPG4". */
/* Sometime in the future, those #ifdefs should be removed, so _XPG4
 * and _XPG4_EXTENDED can be mutually exclusive. */

#ifdef _XPG4_EXTENDED
#  ifndef _XPG4
#    define _XPG4
#  endif /* not _XPG4 */
#endif /* _XPG4_EXTENDED */


/* X/Open namespace defaults to XPG4 semantics for X/Open functions */

#ifdef _XOPEN_SOURCE
#  if !defined(_XPG2) && !defined(_XPG3) && !defined(_XPG4)
#    define _XPG4
#  endif /* not _XPG2 && not _XPG3 && not _XPG4 */
#endif /* _XOPEN_SOURCE */


/* AES namespace defaults to XPG3 semantics for X/Open functions */

#ifdef _AES_SOURCE
#  if !defined(_XPG2) && !defined(_XPG3) && !defined(_XPG4)
#    define _XPG3
#  endif /* not _XPG2 && not _XPG3 && not _XPG4 */
#endif /* _AES_SOURCE */


/* HP-UX namespace defaults to XPG4 semantics for X/Open functions */
/* After 10.0, this should probably be changed to _XPG4_EXTENDED. */

#ifdef _INCLUDE_HPUX_SOURCE
#  if !defined(_XPG2) && !defined(_XPG3) && !defined(_XPG4)
#    define _XPG4
#  endif /* not _XPG2 && not _XPG3 && not _XPG4 */
#endif /* _INCLUDE_HPUX_SOURCE */


/* XPG3 semantics gets POSIX.1-1988 semantics by default */

#if defined(_XPG3) && !defined(_XPG4) && !defined(_POSIX_C_SOURCE)
#  ifndef _POSIX1_1988
#    define _POSIX1_1988
#  endif /* _POSIX1_1988 */
#endif /* _XPG3 && not _XPG4 && not _POSIX_C_SOURCE */


/* X/Open namespace with XPG4 semantics gets POSIX.2 namespace by default */

#ifdef _XPG4
#  ifndef _INCLUDE_POSIX2_SOURCE
#    define _INCLUDE_POSIX2_SOURCE
#  endif /* not _INCLUDE_POSIX2_SOURCE */
#endif /* _XPG4 */


/* _POSIX2_SOURCE is an obsolescent synonym for _POSIX_C_SOURCE being 2 */

#ifdef _POSIX2_SOURCE
#  ifndef _POSIX_C_SOURCE
#    define _POSIX_C_SOURCE 2
#  endif /* not _POSIX_C_SOURCE */
#endif /* _POSIX2_SOURCE */


/* _POSIX_C_SOURCE selects POSIX namespace depending on its value */

#if defined(_POSIX_C_SOURCE)
#  if (_POSIX_C_SOURCE >= 199309) && !defined(_INCLUDE_POSIX4_SOURCE)
#      define _INCLUDE_POSIX4_SOURCE
#  endif /* _POSIX_C_SOURCE >= 199309 && !_INCLUDE_POSIX4_SOURCE */
#  if (_POSIX_C_SOURCE > 1) && !defined(_INCLUDE_POSIX2_SOURCE)
#      define _INCLUDE_POSIX2_SOURCE
#  endif /* (_POSIX_C_SOURCE > 1) && _INCLUDE_POSIX2_SOURCE */
#endif /* defined _POSIX_C_SOURCE */


/* Define ANSI C prototypes for either Standard C or C++ compilations */

#if defined(__STDC__) || defined(__cplusplus)
#  ifndef _PROTOTYPES
#     define _PROTOTYPES
#  endif /* _PROTOTYPES */
#endif /* __STDC__ || __cplusplus */


/* Obsolescent symbols for workstation vs. Server I/O types */

#if !defined(_WSIO) && !defined(_SIO)
#if defined(__hp9000s700) || defined(__hp9000s300)
#  define _WSIO
#else /* not (__hp9000s700 || __hp9000s300) */
#  define _SIO
#endif /* not (__hp9000s700 || __hp9000s300) */
#endif /* not _WSIO && not _SIO */


/* Large (64-bit) Files symbols */

#if defined(__STDC_EXT__) || !defined(__STDC__) || defined(__LP64__)

# ifdef __LP64__
#   ifndef _FILE_OFFSET_BITS
#     define _FILE_OFFSET_BITS 64 
#   else
#    if _FILE_OFFSET_BITS != 64
      #error "_FILE_OFFSET_BITS definition incompatible with __LP64__."
#    endif
#   endif /* _FILE_OFFSET_BITS */
# endif /* __LP64__ */

# ifdef _FILE_OFFSET_BITS
#  if _FILE_OFFSET_BITS == 64
#    define _FILE64
#  else
#    if _FILE_OFFSET_BITS != 32
	#error "_FILE_OFFSET_BITS defined to invalid number!!"
#    endif
#  endif/* _FILE_OFFSET_BITS == 64 */
# endif /* _FILE_OFFSET_BITS */

# define _LFS_LARGEFILE   1
# define _LFS64_LARGEFILE 1

# ifdef _LARGEFILE64_SOURCE
#  ifndef _LARGEFILE_SOURCE
#    define _LARGEFILE_SOURCE
#  endif
# endif

#else  /* strict ANSI */

# ifdef  _FILE_OFFSET_BITS
#  if !defined(__LP64__) && _FILE_OFFSET_BITS == 64
    #error "Large Files (ILP32) not supported in strict ANSI mode."
#  endif
# endif

# ifdef _LARGEFILE64_SOURCE
    #error "Large File interfaces not supported in strict ANSI mode."
# endif

#endif /* __STDC_EXT__ || ! __STDC__  || __LP64__ */

#endif /* _SYS_STDSYMS_INCLUDED */
