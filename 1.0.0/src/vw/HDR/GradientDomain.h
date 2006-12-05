//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//

/// \file GradientDomain.h
/// 
/// Implements the tone mapping method described in the paper,
/// "Gradient Domain High Dynamic Range Compression," by Fattal,
/// Lischinski, and Werman (Siggraph 2002).
///
/// Note: This algorithm does not yet work to our full satisfaction.
/// Use with caution.

#ifndef __VW_HDR_GRADIENTDOMAIN_H__
#define __VW_HDR_GRADIENTDOMAIN_H__


#include <vw/Image/ImageViewBase.h>

namespace vw { 
namespace HDR {

  template <class ImageT, class MethodT>
  ImageT tone_map(ImageViewBase<ImageT> &image, MethodT &method);
			
}}  // namespace vw::HDR

#endif  // __VW_HDR_GRADIENTDOMAIN_H__
