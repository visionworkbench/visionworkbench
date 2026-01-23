// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


/// \file Features.h 
///       This file contains some useful macros and definitions so they don't get
///       scattered everywhere.

#ifndef __VW_CORE_FEATURES_H__
#define __VW_CORE_FEATURES_H__

#include <vw/vw_config.h>

#if defined(VW_COMPILER_HAS_ATTRIBUTE_DEPRECATED)
#define VW_DEPRECATED __attribute__((deprecated))
#else
#define VW_DEPRECATED
#endif

#if defined(VW_COMPILER_HAS_ATTRIBUTE_NORETURN)
#define VW_NORETURN __attribute__((noreturn))
#else
#define VW_NORETURN
#endif

#if defined(VW_COMPILER_HAS_ATTRIBUTE_WARN_UNUSED_RESULT)
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
