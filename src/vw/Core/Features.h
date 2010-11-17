// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// This file contains some useful macros and definitions so they don't get
// scattered everywhere.

#ifndef __VW_CORE_FEATURES_H__
#define __VW_CORE_FEATURES_H__

#include <vw/config.h>

#if defined(VW_COMPILER_HAS_ATTRIBUTE_DEPRECATED) && (VW_COMPILER_HAS_ATTRIBUTE_DEPRECATED==1)
#define VW_DEPRECATED __attribute__((deprecated))
#else
#define VW_DEPRECATED
#endif

#if defined(VW_COMPILER_HAS_ATTRIBUTE_NORETURN) && (VW_COMPILER_HAS_ATTRIBUTE_NORETURN==1)
#define VW_NORETURN __attribute__((noreturn))
#else
#define VW_NORETURN
#endif

#if defined(VW_COMPILER_HAS_ATTRIBUTE_WARN_UNUSED_RESULT) && (VW_COMPILER_HAS_ATTRIBUTE_WARN_UNUSED_RESULT==1)
#define VW_WARN_UNUSED __attribute__((warn_unused_result))
#else
#define VW_WARN_UNUSED
#endif

#define VW_NOTHROW VW_IF_EXCEPTIONS(throw())

/// The master compile-time debugging level flag.  The default value
/// for VW_DEBUG_LEVEL is guessed based on whether or not NDEBUG
/// is defined if the user has not specified it explicitly.
#ifndef VW_DEBUG_LEVEL
#ifdef NDEBUG
#define VW_DEBUG_LEVEL 0
#else
#define VW_DEBUG_LEVEL 1
#endif
#endif

/// A quick macro for selectively disabling code in non-debug builds.
#if VW_DEBUG_LEVEL == 0
#define VW_DEBUG(x)
#else
#define VW_DEBUG(x) x
#endif

#endif
