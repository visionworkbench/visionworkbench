// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
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
// __END_LICENSE__

/// \file Image.h
/// 
/// A convenience header that includes the header files in vw/Image.
/// 
#ifndef __VW_IMAGE_H__
#define __VW_IMAGE_H__

#include <vw/Image/PixelIterator.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelAccessors.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageIO.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/PixelTypeInfo.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/PerPixelAccessorViews.h>
#include <vw/Image/UtilityViews.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/ImageResource.h>
#include <vw/Image/ImageResourceView.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/Interpolation.h>
#include <vw/Image/Convolution.h>
#include <vw/Image/Filter.h>
#include <vw/Image/Transform.h>
#include <vw/Image/BlockProcessor.h>
#include <vw/Image/BlockRasterize.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/BlockCacheView.h>
#include <vw/Image/Palette.h>

#endif // __VW_IMAGE_H__
